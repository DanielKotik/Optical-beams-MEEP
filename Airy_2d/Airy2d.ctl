;;-------------------------------------------------------------------------------------------------
;; file:    Airy2d.ctl
;; brief:   Scheme configuration input file for the FDTD solver Meep simulating the scattering of an 
;;          incomplete Airy beam at a planar dielectric interface
;; author:  Daniel Kotik
;; version: 1.2.0
;; date:    28.11.2019
;;
;; example invocations: a) launch the serial version of meep with specified polarisation (p)
;;
;;                              meep s-pol\?=false Airy2d.ctl
;;
;;                      b) launch the parallel version of meep using 8 cores
;;
;;                              mpirun -quiet -np 8 meep-mpi Airy2d.ctl
;;
;; coordinate system in meep (defines center of computational cell):  --|-----> x
;;                                                                      |
;;                                                                      |
;;                                                                      v y
;;
;; example visualisation (square brackets contain optional arguments for overlaying the dielectric function):
;;
;;          h5topng -S2 -X scalex -c hot [-a yarg -A eps-000000000.h5] e2_s-000003696.h5 
;;
;;          (if necessary, scale the x dimension of the image by scalex)
;;------------------------------------------------------------------------------------------------

(print "\nstart time: "(strftime "%c" (localtime (current-time))) "\n")

;;------------------------------------------------------------------------------------------------
;; physical parameters characterizing light source and interface characteristics 
;; (must be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param s-pol? true )                 ; true for s-spol, false for p-pol
(define-param ref_medium 0)                 ; reference medium whose wavenumber is used as inverse scaling length
                                            ; (0 - free space, 1 - incident medium, 2 - refracted medium)
                                            ; k is then equivalent to k_ref_medium: k_1 = k_0*n_1 or k_2 = k_0*n_2
(define-param n1  1.00)                     ; index of refraction of the incident medium
(define-param n2  0.65)                     ; index of refraction of the refracted medium
(define-param kw_0   7)                     ; beam width (>5 is good)
(define-param kr_w   0)                     ; beam waist distance to interface (30 to 50 is good if
                                            ; source position coincides with beam waist)
(define-param M  0)                         ; center of integration window
(set-param!   W  4)                         ; width of integration window

(define (Critical n1 n2)                    ; calculates the critical angle in degrees
    (cond
      ((> n1 n2) (* (/ (asin (/ n2 n1)) (* 2.0 pi)) 360.0))
      ((< n1 n2) (print "\nWarning: Critical angle is not defined, since n1 < n2!\n\n") (exit))
      ((= n1 n2) (print "\nWarning: Critical angle is not defined, since n1 = n2!\n\n") (exit))
    ))

(define (Brewster n1 n2)                    ; calculates the Brewster angle in degrees
        (* (/ (atan (/ n2 n1)) (* 2.0 pi)) 360.0))

;; define incidence angle relative to the Brewster or critical angle, or set it explicitly (in degrees)
;(define-param chi_deg  (* 0.85 (Brewster n1 n2)))
(define-param chi_deg  (* 1.0 (Critical n1 n2)))
;(define-param chi_deg  45.0)

;;------------------------------------------------------------------------------------------------ 
;; specific Meep paramters (may need to be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param sx 10)                        ; size of cell including PML in x-direction
(define-param sy 10)                        ; size of cell including PML in y-direction
(define-param pml_thickness 0.25)            ; thickness of PML layer
(define-param freq    12)                   ; vacuum frequency of source (4 to 12 is good)
(define-param runtime 90)                   ; runs simulation for X times freq periods
(define-param pixel   15)                   ; number of pixels per wavelength in the denser medium
                                            ; (at least >10; 20 to 30 is a good choice)
;(define-param source_shift 0)              ; source position with respect to the center (point of impact) in Meep
;(define-param source_shift (* -1.0 rw))    ; units (-2.15 good); if equal -rw, then source position coincides with
                                            ; waist position
(define-param source_shift (* -0.4 (- sx (* 2 pml_thickness))))
(define-param relerr 1.0e-5)                ; relative error for integration routine (1.0e-4 or smaller)

;;------------------------------------------------------------------------------------------------
;; derived Meep parameters (do not change)
;;------------------------------------------------------------------------------------------------
(define k_vac (* 2.0 pi freq))
(define k1    (* n1  k_vac  ))              ; wave number inside the incident medium
(define n_ref (cond ((= ref_medium 0) 1.0)  ; index of refraction of the reference medium
                    ((= ref_medium 1)  n1)
                    ((= red_medium 2)  n2)))
                    
(define rw  (/ kr_w (* n_ref k_vac)))
(define w_0 (/ kw_0 (* n_ref k_vac)))
(define shift (+ source_shift rw))          ; distance from source position to beam waist (along y-axis)

;;------------------------------------------------------------------------------------------------
;; placement of the dielectric interface within the computational cell
;;------------------------------------------------------------------------------------------------
;; helper functions
(define (alpha _chi_deg)                    ; angle of inclined plane with y-axis
        (- (/ pi 2.0) (* (/ _chi_deg 360) 2 pi)))
(define (Delta_x _alpha)                    ; inclined plane offset to the center of the cell
        (* (/ sx 2.0) (/ (-(- (sqrt 2.0) (cos _alpha)) (sin _alpha)) (sin _alpha))))
(define (chi_rad _chi_deg)                  ; conversion degrees to radians
        (* (/ _chi_deg 360.0) (* 2.0 pi)))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! default-material (make dielectric (index n1)))
(set! geometry (list
                (make block                 ; located at lower right edge for 45 degree tilt
                (center (+ (/ sx 2.0) (Delta_x (alpha chi_deg))) (/ sy -2.0))
                (size infinity (* (sqrt 2.0) sx) infinity)
                (e1 (/ 1.0 (tan (alpha chi_deg)))  1 0)
                (e2 -1 (/ 1.0 (tan (alpha chi_deg))) 0)
                (e3 0 0 1)
                (material (make dielectric (index n2))))
                ))

;;------------------------------------------------------------------------------------------------
;; add absorbing boundary conditions and discretize structure
;;------------------------------------------------------------------------------------------------
(set! pml-layers 
    (list (make pml (thickness pml_thickness))))
(set! resolution                            ; set resolution in pixels per Meep distance unit
      (* pixel (* (if (> n1 n2) n1 n2) freq)))
(set! Courant                               ; set Courant factor (mandatory if either n1 or n2 is smaller than 1)
      (/ (if (< n1 n2) n1 n2) 2))

;;------------------------------------------------------------------------------------------------
;; beam profile distribution (field amplitude) at the waist of the beam
;;------------------------------------------------------------------------------------------------
(define (Gauss W_y)
        (lambda (r) (exp (* -1.0 (expt (/ (vector3-y r) W_y) 2.0)))
        ))

;; incomplete Airy function
(define (Ai_inc W_y M W)
        (lambda (r) (car
        (integrate (lambda (xi) (exp (* 0+1i (+ (* (/ -1 3) (expt xi 3)) (* xi (/ (vector3-y r) W_y))))))
                   (- M W) (+ M W) relerr))
        ))

;; simple test outputs
;(print "w_0: " w_0 "\n")
;(print "Airy function 1: " ((Ai_inc w_0 0 4) (vector3 1 -0.3 1)) "\n")
;(exit)

;;------------------------------------------------------------------------------------------------
;; spectrum amplitude distribution
;;------------------------------------------------------------------------------------------------
(define Heaviside
        (lambda (x) (cond ((<  x 0) 0) ((>= x 0) 1))
        ))

(define (f_Gauss W_y)
        (lambda (k_y) (exp (* -1.0 (expt (* 0.5 k_y W_y) 2)))
        ))

(define (f_Airy W_y M W)
        (lambda (k_y) (* W_y (exp (* 0+1i (* (/ -1 3) (expt (* k_y W_y) 3))))
                             (Heaviside (- (* W_y k_y) (- M W))) (Heaviside (- (+ M W) (* W_y k_y))))
        ))

;; simple test outputs
;(print "Airy spectrum: " ((f_Airy w_0 0 4) 0.2) "\n")
;(exit)

;;------------------------------------------------------------------------------------------------
;; plane wave decomposition 
;; (purpose: calculate field amplitude at light source position if not coinciding with beam waist)
;;------------------------------------------------------------------------------------------------
(define (integrand f x y)
        (lambda (k_y) (* (f k_y)
                        (exp (* 0+1i x (sqrt (- (* k1 k1) (* k_y k_y)))))
                        (exp (* 0+1i k_y y)))
        ))

;; complex field amplitude at position (x, y) with spectrum amplitude distribution f
;; (one may have to adjust the 'relerr' parameter value in the integrate function)
(define (psi f x)
        (lambda (r) (car (integrate (integrand f x (vector3-y r))
                          (* -1.0 k1) (* 1.0 k1) relerr))
        ))

;(print "Airy function 2: " ((psi (f_Airy w_0 0 4) 0) (vector3 1 -0.3 1)) "\n")
;(exit)

;;------------------------------------------------------------------------------------------------
;; display values of physical variables
;;------------------------------------------------------------------------------------------------
(print "\n")
(print "Specified variables and derived values: \n")
(print "chi:   " chi_deg        " [degree]\n") ; angle of incidence
(print "incl.: " (- 90 chi_deg) " [degree]\n") ; interface inclination with respect to the x-axis
(print "kw_0:  " kw_0  "\n")
(print "kr_w:  " kr_w  "\n")
(print "k_vac: " k_vac "\n")
(print "polarisation: " (if s-pol? "s" "p") "\n")
(print "\n")

;;------------------------------------------------------------------------------------------------
;; specify current source, output functions and run simulation
;;------------------------------------------------------------------------------------------------
(use-output-directory)                      ; put output files in a separate folder
(set! force-complex-fields? false)          ; default: false
(set! eps-averaging? true)                  ; default: true

(set! sources (list
                  (make source
                      (src (make continuous-src (frequency freq) (width 0.5)))
                      (if s-pol? (component Ez) (component Ey))
                      (size 0 9 0)
                      (center source_shift 0 0)
                      ;(amp-func (Gauss w_0)))
                      ;(amp-func (Ai_inc w_0 M W)))
                      (amp-func (psi (f_Airy w_0 M W) shift)))
                  ))

;; calculates |E|^2 with |.| denoting the complex modulus if 'force-complex-fields?' is set to true, otherwise |.|
;; gives the Euclidean norm
(define (eSquared r ex ey ez)
        (+ (expt (magnitude ex) 2) (expt (magnitude ey) 2) (expt (magnitude ez) 2)))

(define (output-efield2) (output-real-field-function (if s-pol? "e2_s" "e2_p")
                                                     (list Ex Ey Ez) eSquared))

(run-until runtime
     (at-beginning (lambda () (print "\nCalculating inital field configuration. This will take some time...\n\n")))
     (at-beginning output-epsilon)          ; output of dielectric function
     (if s-pol?
         (at-end output-efield-z)           ; output of E_z component (for s-polarisation)
         (at-end output-efield-y))          ; output of E_y component (for p-polarisation)
     (at-end output-efield2))               ; output of electric field intensity

(print "\nend time: "(strftime "%c" (localtime (current-time))) "\n")
