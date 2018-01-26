;;-------------------------------------------------------------------------------------------------
;; file:   LaguerreGauss3d.ctl
;; brief:  Scheme configuration input file for the FDTD solver Meep simulating the scattering of a polarised
;;         Laguerre-Gaussian beam at a planar dielectric interface (3d)
;; author: Daniel Kotik
;; date:   17.01.2018
;;
;; example invocations: a) launch the serial version of meep with specified polarisation (p)
;;
;;                              meep s-pol\?=false LaguerreGauss3d.ctl
;;
;;                      b) launch the parallel version of meep using 8 cores
;;
;;                              mpirun -np 8 meep-mpi LaguerreGauss3d.ctl
;;
;; coordinate system in meep (defines center of computational cell):  --|-----> x
;;                                                                      |
;;                                                                      |
;;                                                                      v y
;;
;; visualisation: h5topng -S2 -0 -z 0 -c  hot       -a yarg -A eps-000000000.h5 e2_s-000001232.h5
;;                h5topng -S2 -0 -z 0 -Zc dkbluered -a gray -A eps-000000000.h5   ez-000001232.h5
;;------------------------------------------------------------------------------------------------

;;------------------------------------------------------------------------------------------------
;; physical parameters characterizing light source and interface characteristics 
;; (must be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param s-pol? true )                 ; true for s-spol, false for p-pol
(define-param e_z        1)                 ; z-component of Jones vector (s-polarisation: e_z = 1, e_y = 0)
(define-param e_y        0)                 ; y-component of Jones vector (p-polarisation: e_z = 0, e_y = 1)
                                            ;                      (circular-polarisation: ...             )
(define-param m_charge   2)                 ; vortex charge (azimuthal quantum number, integer number)
(define-param ref_medium 0)                 ; reference medium whose wavenumber is used as inverse scaling length
                                            ; (0 - free space, 1 - incident medium, 2 - refracted medium)
                                            ; k is then equivalent to k_ref_medium: k_1 = k_0*n_1 or k_2 = k_0*n_2
(define-param n1  1.00)                     ; index of refraction of the incident medium
(define-param n2  1.00)                     ; index of refraction of the refracted medium
(define-param kw_0   8)                     ; beam width (>5 is good)
(define-param kr_w   0)                     ; beam waist distance to interface (30 to 50 is good if
                                            ; source position coincides with beam waist)

(define Critical                            ; calculates the critical angle in degrees
    (cond
      ((> n1 n2) (* (/ (asin (/ n2 n1)) (* 2.0 pi)) 360.0))
      (else      (display "\nWarning: Critical angle is not defined, since n1 < n2!\n\n"))
    ))  

(define Brewster                            ; calculates the Brewster angle in degrees
        (* (/ (atan (/ n2 n1)) (* 2.0 pi)) 360.0))

;(define-param chi_deg  (* 0.99 Critical))   ; define incidence angle relative to the Brewster or critical angle,
(define-param chi_deg  45.0)               ; or set it explicitly (in degrees)

;;------------------------------------------------------------------------------------------------ 
;; specific Meep paramters (may need to be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param sx 5)                         ; size of cell including PML in x-direction
(define-param sy 5)                         ; size of cell including PML in y-direction
(define-param sz 5)                         ; size of cell including PML in z-direction
(define-param pml_thickness 0.25)           ; thickness of PML layer
(define-param freq     4)                   ; vacuum frequency of source (default 4)
(define-param runtime 10)                   ; runs simulation for 10 times freq periods
(define-param pixel   10)                   ; number of pixels per wavelength in the denser
                                            ; medium (at least >10; 20 to 30 is a good choice)
(define-param source_shift 0)           ; source position with respect to the center (point of impact) in Meep
;(define-param source_shift (* -1.0 r_w))   ; units (-2.15 good); if equal -r_w, then source position coincides with
                                            ; waist position
(define-param relerr 0.0001)                ; relative error for integration routine (0.0001 or smaller)
(define-param maxeval 1000)                ; maximum evaluations for integration routine

;;------------------------------------------------------------------------------------------------
;; derived Meep parameters (do not change)
;;------------------------------------------------------------------------------------------------
(define k_vac (* 2.0 pi freq))
(define n_ref (cond ((= ref_medium 0) 1.0)  ; index of refraction of the reference medium
                    ((= ref_medium 1)  n1)
                    ((= red_medium 2)  n2)))
                    
(define r_w  (/ kr_w (* n_ref k_vac)))
(define w_0 (/ kw_0 (* n_ref k_vac)))
(define shift (+ source_shift r_w))         ; distance from source position to beam waist (along y-axis)

;;------------------------------------------------------------------------------------------------
;; placement of the planar dielectric interface within the computational cell
;;------------------------------------------------------------------------------------------------
;; helper functions
(define (alpha _chi_deg)                    ; angle of inclined plane with y-axis
        (- (/ pi 2.0) (* (/ _chi_deg 360) 2 pi)))
(define (Delta_x _alpha)                    ; inclined plane offset to the center of the cell
        (* (/ sx 2.0) (/ (-(- (sqrt 2.0) (cos _alpha)) (sin _alpha)) (sin _alpha))))

(set! geometry-lattice (make lattice (size sx sy sz)))
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

;;------------------------------------------------------------------------------------------------
;; 2d-beam profile distribution (field amplitude) at the waist of the beam
;;------------------------------------------------------------------------------------------------
(define (Gauss W_y)
        (lambda (r) (exp (* -1.0 (/ (+ (* (vector3-y r) (vector3-y r)) (* (vector3-z r) (vector3-z r))) 
                                    (* W_y W_y))))
        ))

;; some test outputs
;(print "Gauss 2d beam profile: " ((Gauss w_0) (vector3 0 0.5 0.2)) "\n")
;(exit)


;;------------------------------------------------------------------------------------------------
;; spectrum amplitude distribution(s)
;;------------------------------------------------------------------------------------------------
(define (f_Gauss W_y)
        (lambda (k_y k_z) (* (/ W_y (sqrt (* 2 pi))) ;TODO: remove unnecessary prefactor
                             (exp (* -1 (* (* W_y W_y) (/ (+ (* k_y k_y) (* k_z k_z)) 4)))))
        ))

;; spherical coordinate transformation in k-space
(define (phi k)
        (lambda (k_y k_z) (atan (/ k_y k) (/ (* -1 k_z) k))
        ))

(define (theta k)
        (lambda (k_z) (acos (/ (* -1 k_z) k))
        ))

(define (f_Laguerre_Gauss W_y k)
        (lambda (k_y k_z) (* ((f_Gauss W_y) k_y k_z) (exp (* 0+1i m_charge ((phi k) k_y k_z))) 
                             (expt ((theta k) k_z) (abs m_charge)))
        ))

;; some test outputs
;(print "Gauss 2d spectrum: "          ((f_Gauss 20) 0.1 0.1)                "\n")
;(print "Laguerre-Gauss 2d spectrum: " ((f_Laguerre_Gauss 20 k_vac) 0.1 0.1) "\n")
;(exit)

;;------------------------------------------------------------------------------------------------
;; plane wave decomposition 
;; (purpose: calculate field amplitude at light source position if not coinciding with beam waist)
;;------------------------------------------------------------------------------------------------
(define (integrand f x y z k)
        (lambda (k_y k_z) (* (f k_y k_z)
                             (exp (* 0+1i x (real-part (sqrt (- (* k k) (* k_y k_y) (* k_z k_z))))))
                             (exp (* 0+1i y k_y))
                             (exp (* 0+1i z k_z)))
        ))


;; complex field amplitude at position (x, y) with spectrum amplitude distribution f
;; (one may have to adjust the 'relerr' and 'maxeval' parameter values in the integrate function)
(define (psi f x k)
        (lambda (r) (car (integrate (integrand f x (vector3-y r) (vector3-z r) k)
                         (list (* -1.0 k) (* -1.0 k)) (list (* 1.0 k) (* 1.0 k)) relerr 0 maxeval))
        ))

;(print "Gauss 2d beam profile: " ((psi (f_Gauss w_0) 0.0 (* n1 k_vac)) (vector3 0 0.0 0.0)) "\n")
;(print "Gauss 2d beam profile: " ((Gauss w_0)                          (vector3 0 0.0 0.0)) "\n")
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
(print "vortex charge: " m_charge "\n")
(print "Jones vector components (e_z, e_y): ("e_z", "e_y")" "\n")
(print "degree of linear polarisation at pi/4: " (* 2 (imag-part (* (conj e_y) e_z)))   "\n")
(print "degree of circular polarisation: " (* 2 (real-part (* (conj e_y) e_z)))   "\n")
(print "polarisation: " (if s-pol? "s" "p") "\n")
(print "\n")
(exit)
;(print "The value of our Gaussian spectrum amplitude is: " ((f_Gauss w_0) 20.0) "\n")
;(print "integrand " ((integrand 0.8 2.0 k_vac w_0) 20.0) "\n")
;(print "Field amplitude: " ((psi 1.0 k_vac w_0) 0.5) "\n")

;;------------------------------------------------------------------------------------------------
;; specify current source, output functions and run simulation
;;------------------------------------------------------------------------------------------------
(use-output-directory)                      ; put output files in a separate folder
(set! force-complex-fields? false)          ; default: false
(set! eps-averaging? false)                  ; default: true

(set! sources (list
                  (make source
                      (src (make continuous-src (frequency freq) (width 0.5)))
                      (if s-pol? (component Ez) (component Hz))
                      (amplitude 1.0)
                      (size 0 2 2)
                      (center source_shift 0 0)
                      ;(amp-func (Gauss w_0)))
                      ;(amp-func (psi (f_Gauss w_0) shift (* n1 k_vac))))
                      (amp-func (psi (f_Laguerre_Gauss w_0 (* n1 k_vac)) shift (* n1 k_vac))))
                  ))

;; exploiting symmetries to reduce computational effort:
;; The plane of incidence (x-y-plane) is a mirror plane which is characterised to be orthogonal to the z-axis
;; (symmetry of the geometric structure). Symmetry of the sources must be ensured simultaneously, which is possible 
;; for certain cases by adding a phase. This works for example for pure s- or p-polarisation, where either only
;; the Ez or Hz component is specified.
(set! symmetries (list (make mirror-sym (direction Z) (phase -1))))

(define (eSquared r ex ey ez)
        (+ (* (magnitude ex) (magnitude ex)) (* (magnitude ey) (magnitude ey))
           (* (magnitude ez) (magnitude ez))))

(define (output-efield2) (output-real-field-function (if s-pol? "e2_s" "e2_p")
                                                     (list Ex Ey Ez) eSquared))

(run-until runtime
;     (at-beginning output-epsilon)          ; output of dielectric function
     (if s-pol?
         (at-end output-efield-z)           ; output of E_z component (for s-polarisation)
         (at-end output-hfield-z))          ; output of H_z component (for p-polarisation)
     (at-end output-efield2))               ; output of electric field intensity
