;;-------------------------------------------------------------------------------------------------
;; file:   planar_new.ctl
;; brief:  Scheme configuration input file for the FDTD solver Meep simulating the scattering of a 
;;         Gaussian beam at a plane dielectric interface
;; author: Daniel Kotik
;; date:   2013-2018
;;
;; example invocation: meep s-pol\?=false planar_new.ctl
;;
;; coordinate system in meep (defines center of computational cell):  --|-----> x
;;                                                                      |
;;                                                                      |
;;                                                                      v y
;;------------------------------------------------------------------------------------------------ 

;;------------------------------------------------------------------------------------------------
;; physical parameters characterizing light source and interface characteristics 
;; (must be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-input-var interface "concave" 'string (lambda type (or (string=? type "planar" ) 
                                                              (string=? type "concave") 
                                                              (string=? type "convex" ))))
(define-param s-pol? true )                 ; true for s-spol, false for p-pol
(define-param ref_medium 0)                 ; reference medium whose wavenumber is used as inverse scaling length
                                            ; (0 - free space, 1 - incident medium, 2 - refracted medium)
                                            ; k is then equivalent to k_ref_medium: k_1 = k_0*n_1 or k_2 = k_0*n_2
(define-param n1  1.54)                     ; index of refraction of the incident medium
(define-param n2  1.00)                     ; index of refraction of the refracted medium
(define-param kw_0  10)                     ; beam width (10 is good)
(define-param kr_w  50)                     ; beam waist distance to interface (30 to 50 is good if
                                            ; source position coincides with beam waist)
(define-param kr_c 470)                     ; radius of curvature (if interface is either concave of convex)

(define Critical                            ; calculates the critical angle in degrees
    (cond
      ((> n1 n2) (* (/ (asin (/ n2 n1)) (* 2.0 pi)) 360.0))
      (else      (display "\nWarning: Critical angle is not defined, since n1 < n2!\n\n"))
    ))  

(define Brewster                            ; calculates the Brewster angle in degrees
        (* (/ (atan (/ n2 n1)) (* 2.0 pi)) 360.0))

(define-param chi_deg  (* 0.99 Critical))   ; define incidence angle relative to the Brewster or critical angle,
;(define-param chi_deg  45.0)               ; or set it explicitly (in degrees)

;;------------------------------------------------------------------------------------------------ 
;; specific Meep paramters (may need to be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param sx 5)                         ; size of cell including PML in x-direction
(define-param sy 5)                         ; size of cell including PML in y-direction
(define-param pml_thickness 0.25)           ; thickness of PML layer
(define-param freq    12)                   ; vacuum frequency of source (5 to 12 is good)
(define-param runtime 10)                   ; runs simulation for 10 times freq periods
(define-param pixel   10)                   ; number of pixels per wavelength in the denser
                                            ; medium (at least >10; 20 to 30 is a good choice)
(define-param source_shift -2.15)           ; source position with respect to the center (point of impact) in Meep
;(define-param source_shift (* -1.0 rw))    ; units (-2.15 good); if equal -rw, then source position coincides with
                                            ; waist position
(define-param relerr 0.0001)                ; relative error for integration routine (0.0001 or smaller)

;;------------------------------------------------------------------------------------------------
;; derived Meep parameters (do not change)
;;------------------------------------------------------------------------------------------------
(define k_vac (* 2.0 pi freq))
(define n_ref (cond ((= ref_medium 0) 1.0)  ; index of refraction of the reference medium
                    ((= ref_medium 1)  n1)
                    ((= red_medium 2)  n2)))
                    
(define rw  (/ kr_w (* n_ref k_vac)))
(define w_0 (/ kw_0 (* n_ref k_vac)))
(define r_c (/ kr_c (* n_ref k_vac)))
(define shift (+ source_shift rw))          ; distance from source position to beam waist (along y-axis)

;;------------------------------------------------------------------------------------------------
;; placement of the planar dielectric interface within the computational cell
;;------------------------------------------------------------------------------------------------
;; helper functions
(define (alpha _chi_deg)                    ; angle of inclined plane with y-axis
        (- (/ pi 2.0) (* (/ _chi_deg 360) 2 pi)))
(define (Delta_x _alpha)                    ; inclined plane offset to the center of the cell
        (* (/ sx 2.0) (/ (-(- (sqrt 2.0) (cos _alpha)) (sin _alpha)) (sin _alpha))))
(define (chi_rad _chi_deg)                  ; conversion degrees to radians
        (* (/ _chi_deg 360.0) (* 2.0 pi)))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(cond
    ((string=? interface "planar")
        (set! default-material (make dielectric (index n1)))
        (set! geometry (list
                        (make block         ; located at lower right edge for 45 degree tilt
                        (center (+ (/ sx 2.0) (Delta_x (alpha chi_deg))) (/ sy -2.0))
                        (size infinity (* (sqrt 2.0) sx) infinity)
                            (e1 (/ 1.0 (tan (alpha chi_deg)))  1 0)
                            (e2 -1 (/ 1.0 (tan (alpha chi_deg))) 0)
                            (e3 0 0 1)
                        (material (make dielectric (index n2)))))))
    ((string=? interface "concave")
        (set! default-material (make dielectric (index n2)))
        (set! geometry (list
                    (make cylinder
                    (center (* -1 (* r_c (cos (chi_rad chi_deg)))) (* r_c (sin (chi_rad chi_deg)))) 
                                          ; Mittelpunkt wird nach rechts verschoben,
                                          ; so dass der Auftreffpkt immer mittig liegt
                    (height infinity)
                    (radius r_c)
                    (material (make dielectric (index n1)))))))
    ((string=? interface "convex" ) )
)

;;------------------------------------------------------------------------------------------------
;; add absorbing boundary conditions and discretize structure
;;------------------------------------------------------------------------------------------------
(set! pml-layers 
    (list (make pml (thickness pml_thickness))))
(set! resolution                            ; set resolution in pixels per Meep distance unit
      (* pixel (* (if (> n1 n2) n1 n2) freq))) 

;;------------------------------------------------------------------------------------------------
;; beam profile distribution(s) (field amplitude) at the waist of the beam
;;------------------------------------------------------------------------------------------------
(define (Gauss W_y)
        (lambda (r) (exp (* -1.0 (expt (/ (vector3-y r) W_y) 2.0)))
        ))

;(define (Asymmetric W_y)
;        (lambda (r) ...
;        ))

;;------------------------------------------------------------------------------------------------
;; spectrum amplitude distribution(s)
;;------------------------------------------------------------------------------------------------
(define (f_Gauss W_y)
        (lambda (k_y) (* (/ W_y (* 2.0 (sqrt pi)))
                         (exp (* -1.0 (expt (* 0.5 k_y W_y) 2.0))))
        ))

;(define (f_asymmetric W_y)
;        (lambda (k_y) ...
;        ))

;;------------------------------------------------------------------------------------------------
;; plane wave decomposition 
;; (purpose: calculate field amplitude at light source position if not coinciding with beam waist)
;;------------------------------------------------------------------------------------------------
(define (integrand f y x k)
        (lambda (k_y) (* (f k_y)
                        (exp (* 0+1i x (sqrt (- (* k k) (* k_y k_y)))))
                        (exp (* 0+1i k_y y)))
        ))

;; complex field amplitude at position (x, y) with spectrum amplitude distribution f
;; (one may have to adjust the 'relerr' parameter value in the integrate function)
(define (psi f x k)
        (lambda (r) (car (integrate (integrand f (vector3-y r) x k)
                          (* -1.0 k) (* 1.0 k) relerr))
        ))

;;------------------------------------------------------------------------------------------------
;; display values of physical variables
;;------------------------------------------------------------------------------------------------
(print "\n")
(print "Values of specified variables:    \n")
(print "chi:   " chi_deg        " [degree]\n") ; angle of incidence
(print "incl.: " (- 90 chi_deg) " [degree]\n") ; interface inclination with respect to the x-axis
(print "kw_0:  " kw_0  "\n")
(print "kr_w:  " kr_w  "\n")
(print "k_vac: " k_vac "\n")
(print "polarisation: " (if s-pol? "s" "p") "\n")
(print "interface: " interface "\n")
(print "\n")
;(print "The value of our Gaussian spectrum amplitude is: " ((f_Gauss w_0) 20.0) "\n")
;(print "integrand " ((integrand 0.8 2.0 k_vac w_0) 20.0) "\n")
;(print "Field amplitude: " ((psi 1.0 k_vac w_0) 0.5) "\n")

(set! sources (list
                (make source
                    (src (make continuous-src (frequency freq) (width 0.5)))
                    (if s-pol? (component Ez) (component Hz))
                    (amplitude 3.0)
                    (size 0 2.0 0)
                    (center source_shift 0 0)
                    ;(amp-func (Gauss w_0)))
                    ;(amp-func (Asymmetric (/ w_0 (sqrt 3.0)))))
                    (amp-func (psi (f_Gauss w_0) shift (* n1 k_vac))))
                ))


;;------------------------------------------------------------------------------------------------
;; specify output functions and run simulation
;;------------------------------------------------------------------------------------------------
(use-output-directory interface)            ; put output files in a separate folder
(set! force-complex-fields? false)          ; default: false
(set! eps-averaging? true)                  ; default: true

(define (eSquared r ex ey ez)
        (+ (* (magnitude ex) (magnitude ex)) (* (magnitude ey) (magnitude ey))
           (* (magnitude ez) (magnitude ez))))

(define (output-efield2) (output-field-function (if s-pol? "e2_s" "e2_p")
                                                (list Ex Ey Ez) eSquared))

(run-until runtime
     (at-beginning output-epsilon)          ; output of dielectric function
     (if s-pol?
         (at-end output-efield-z)           ; output of E_z component (for s-polarisation)
         (at-end output-hfield-z))          ; output of H_z component (for p-polarisation)
     (at-end output-efield2))               ; output of electric field intensity
