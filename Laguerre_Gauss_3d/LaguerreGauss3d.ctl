;;-------------------------------------------------------------------------------------------------
;; file:   LaguerreGauss3d.ctl
;; brief:  Scheme configuration input file for the FDTD solver Meep simulating the scattering of a polarised
;;         Laguerre-Gaussian beam at a planar dielectric interface (3d)
;; author: Daniel Kotik
;; date:   17.01.2018
;;
;; example invocations: a) launch the serial version of meep with specified polarisation (p)
;;
;;                              meep e_z=0 e_y=1 LaguerreGauss3d.ctl
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
;; example visualisations (square brackets contain optional arguments for overlaying the dielectric function)
;;      - slice within plane of incidence (x-y plane)
;;          h5topng -S2 -0 -z 0 -c  hot       [-a yarg -A eps-000000000.h5] e2_s-000001232.h5
;;          h5topng -S2 -0 -z 0 -Zc dkbluered [-a gray -A eps-000000000.h5]   ez-000001232.h5

;;      - slice transversal to the incident propagation axis (INDEX specifies slice index)
;;          h5topng -S2 -x IDNEX -c  hot       [-a yarg -A eps-000000000.h5] e2_s-000001232.h5
;;          h5topng -S2 -x INDEX -Zc dkbluered [-a gray -A eps-000000000.h5]   ez-000001232.h5

;;      - full 3D simulation (creating a VTK file to be opened e.g., with MayaVi)
;;          h5tovtk e2_s-000001232.h5 
;;------------------------------------------------------------------------------------------------

(print "\nstart time: "(strftime "%c" (localtime (current-time))) "\n")

;;------------------------------------------------------------------------------------------------
;; physical parameters characterizing light source and interface characteristics 
;; (must be adjusted - either here or via command line)
;;------------------------------------------------------------------------------------------------
(define-param e_z        1)                 ; z-component of Jones vector (s-polarisation: e_z = 1, e_y = 0)
(define-param e_y        0)                 ; y-component of Jones vector (p-polarisation: e_z = 0, e_y = 1)
                                            ;                      (circular-polarisation: ...             )
(define-param m_charge   2)                 ; vortex charge (azimuthal quantum number, integer number)
(define-param ref_medium 0)                 ; reference medium whose wavenumber is used as inverse scaling length
                                            ; (0 - free space, 1 - incident medium, 2 - refracted medium)
                                            ; k is then equivalent to k_ref_medium: k_1 = k_0*n_1 or k_2 = k_0*n_2
(define-param n1  1.00)                     ; index of refraction of the incident medium
(define-param n2  1.54)                     ; index of refraction of the refracted medium
(define-param kw_0   8)                     ; beam width (>5 is good)
(define-param kr_w   0)                     ; beam waist distance to interface (30 to 50 is good if
                                            ; source position coincides with beam waist)

(define Critical                            ; calculates the critical angle in degrees
    (cond
      ((> n1 n2) (* (/ (asin (/ n2 n1)) (* 2.0 pi)) 360.0))
      ((< n1 n2) (print "\nWarning: Critical angle is not defined, since n1 < n2!\n\n"))
      ((= n1 n2) (print "\nWarning: Critical angle is not defined, since n1 = n2!\n\n"))
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
(define-param source_shift -2.15)           ; source position with respect to the center (point of impact) in Meep
;(define-param source_shift (* -1.0 r_w))   ; units (-2.15 good); if equal -r_w, then source position coincides with
                                            ; waist position
(define-param relerr 0.0001)                ; relative error for integration routine (0.0001 or smaller)
(define-param maxeval 10000)                ; maximum evaluations for integration routine (we recommend 1000 for testing
                                            ; purposes and 10000 or higher for a final simulation run)

;;------------------------------------------------------------------------------------------------
;; derived Meep parameters (do not change)
;;------------------------------------------------------------------------------------------------
(define k_vac (* 2.0 pi freq))              ; vacuum wave number
(define k1    (* n1  k_vac  ))              ; wave number inside the incident medium
(define n_ref (cond ((= ref_medium 0) 1.0)  ; index of refraction of the reference medium
                    ((= ref_medium 1)  n1)
                    ((= red_medium 2)  n2)))
                    
(define r_w (/ kr_w (* n_ref k_vac)))
(define w_0 (/ kw_0 (* n_ref k_vac)))
(define shift (+ source_shift r_w))         ; distance from source position to beam waist (along y-axis)

(define s-pol?                              ; true if s-polarised
       (if (and (= e_z   1 ) (= e_y   0 )) true false))
(define p-pol?                              ; true if p-polarised
       (if (and (= e_z   0 ) (= e_y   1 )) true false))
(define a-pol?                              ; true if arbitrary (complex) polarised
       (if (and (not s-pol?) (not p-pol?)) true false))
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

;;cartesian coordinates (not recommended) ---------------------------------------------
(define (f_Gauss_cartesian W_y)
        (lambda (k_y k_z) (exp (* -1 (* (* W_y W_y) (/ (+ (* k_y k_y) (* k_z k_z)) 4))))
        ))

;; spherical coordinate transformation in k-space
(define (phi k)
        (lambda (k_y k_z) (atan (/ k_y k) (/ (* -1 k_z) k))
        ))

(define (theta k)
        (lambda (k_y k_z) (acos (/ (real-part (sqrt (- (* k k) (* k_y k_y) (* k_z k_z)))) k))
        ))

(define (f_Laguerre_Gauss_cartesian W_y m)
        (lambda (k_y k_z) (* ((f_Gauss_cartesian W_y) k_y k_z) (exp (* 0+1i m ((phi k1) k_y k_z))) 
                             (expt ((theta k1) k_y k_z) (abs m)))
        ))

;; spherical coordinates --------------------------------------------
(define (f_Gauss_spherical W_y)
        (lambda (theta . phi)               ; phi is an optional argument that can be passed (but is not used here)
                (exp (* -1 (expt (/ (* k1 W_y theta) 2) 2)))
        ))

(define (f_Laguerre_Gauss_spherical W_y m)
        (lambda (theta phi) (* ((f_Gauss_spherical W_y) theta) (expt theta (abs m)) (exp (* 0+1i m phi)))
        ))

;; some test outputs
(print "Gauss spectrum (cartesian): " ((f_Gauss_cartesian w_0) 1.0 5.2)          "\n")
(print "Gauss spectrum (spherical): " ((f_Gauss_spherical w_0) (/ pi 3))         "\n\n")

(print "L-G spectrum   (cartesian): " ((f_Laguerre_Gauss_cartesian w_0 m_charge) 1.0 5.2) "\n")
(print "L-G spectrum   (spherical): " ((f_Laguerre_Gauss_spherical w_0 m_charge) (/ pi 3) (/ pi 4)) "\n\n")

;;------------------------------------------------------------------------------------------------
;; plane wave decomposition 
;; (purpose: calculate field amplitude at light source position if not coinciding with beam waist)
;;------------------------------------------------------------------------------------------------
(define (integrand_cartesian f x y z)
        (lambda (k_y k_z) (* (f k_y k_z)
                             ;(exp (* 0+1i x (real-part (sqrt (- (* k1 k1) (* k_y k_y) (* k_z k_z))))))
                             ;(exp (* 0+1i y k_y))
                             ;(exp (* 0+1i z k_z)))
                             (exp (* 0+1i (+ (* x (real-part (sqrt (- (* k1 k1) (* k_y k_y) (* k_z k_z))))) 
                                             (* y k_y) (* z k_z)))))
        ))

(define (integrand_spherical f x y z)
        (lambda (theta phi) (* k1 k1 (sin theta) (cos theta) (f theta phi)
                               ;(exp (* 0-1i k1 z (sin theta) (cos phi)))
                               ;(exp (* 0+1i k1 y (sin theta) (sin phi)))
                               ;(exp (* 0+1i k1 x (cos theta))))
                               (exp (* 0+1i k1 (+ (* (sin theta) (- (* y (sin phi)) (* z (cos phi)))) 
                                                  (* (cos theta) x)))))
        ))

;; complex field amplitude at position (x, y) with spectrum amplitude distribution f
;; (one may have to adjust the 'relerr' and 'maxeval' parameter values in the integrate function)
(define (psi_cartesian f x)
        (lambda (r) (car (integrate (integrand_cartesian f x (vector3-y r) (vector3-z r))
                         (list (* -1.0 k1) (* -1.0 k1)) (list (* 1.0 k1) (* 1.0 k1)) relerr 0 maxeval))
        ))

(define (psi_spherical f x)
        (lambda (r) (car (integrate (integrand_spherical f x (vector3-y r) (vector3-z r))
                         (list 0 0) (list (/ pi 2) (* 2 pi)) relerr 0 maxeval))
        ))

(print "integrand      (cartesian): " ((integrand_cartesian (f_Laguerre_Gauss_cartesian w_0 m_charge)
                                                            -2.15 0.3 0.5)    4.0      0.0  )   "\n")
(print "integrand      (spherical): " ((integrand_spherical (f_Laguerre_Gauss_spherical w_0 m_charge)
                                                            -2.15 0.3 0.5) (/ pi 3) (/ pi 4)) "\n\n")

(print "psi            (cartesian): " ((psi_cartesian (f_Laguerre_Gauss_cartesian w_0 m_charge) -2.15)
                                       (vector3 0 0.3 0.5)) "\n")

(print "psi            (spherical): " ((psi_spherical (f_Laguerre_Gauss_spherical w_0 m_charge) -2.15) 
                                       (vector3 0 0.3 0.5)) "\n")

;(print "psi       (origin, simple): " ((Gauss w_0) (vector3 0 0.2 0.2)) "\n")
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
(print "Jones vector components: (e_z=" e_z ", e_y=" e_y ")")
(print " ---> " (cond (s-pol? "s-") (p-pol? "p-") (a-pol? "mixed-")) "polarisation" "\n")

(print "degree of linear   polarisation at pi/4: " (* 2 (imag-part (* (conj e_y) e_z))) "\n")
(print "degree of circular polarisation: "         (* 2 (real-part (* (conj e_y) e_z))) "\n")
(print "\n")

;;------------------------------------------------------------------------------------------------
;; exploiting symmetries to reduce computational effort
;; (only possible for beams without intrinsic orbital angular momentum, i.e. no vortex charge)
;;------------------------------------------------------------------------------------------------

;; The plane of incidence (x-y-plane) is a mirror plane which is characterised to be orthogonal to the z-axis
;; (symmetry of the geometric structure). Symmetry of the sources must be ensured simultaneously, which is only
;; possible for certain cases. If I am not mistaken this can only be achieved for vortex free beams with pure s- or
;; p-polarisation, i.e. where either the Ez or Ey component is specified.
(if (equal? m_charge 0)
        (cond (s-pol?                       ; s-polarisation
              (set! symmetries (list (make mirror-sym (direction Z) (phase -1)))))
              (p-pol?                       ; p-polarisation
              (set! symmetries (list (make mirror-sym (direction Z)           ))))
        )
)

;;------------------------------------------------------------------------------------------------
;; specify current source, output functions and run simulation
;;------------------------------------------------------------------------------------------------
(use-output-directory)                      ; put output files in a separate folder
(set! force-complex-fields? false)          ; default: false
(set! eps-averaging? false)                  ; default: true

(set! sources (filter (compose not unspecified?)
              (list
                  (if (not (equal? e_z 0))
                      (make source
                          (src (make continuous-src (frequency freq) (width 0.5)))
                          (component Ez)
                          (amplitude e_z)
                          (size 0 3 3)
                          (center source_shift 0 0)
                          ;;(amp-func (Gauss w_0))
                          ;;(amp-func (psi_cartesian (f_Laguerre_Gauss_cartesian w_0 m_charge) shift))
                          (if (equal? m_charge 0)
                              ; if vortex charge is zero use Gauss spectrum distribution (improves perfomance)
                              (amp-func (psi_spherical (f_Gauss_spherical w_0) shift))
                              (amp-func (psi_spherical (f_Laguerre_Gauss_spherical w_0 m_charge) shift))
                          )
                       )
                  )
                  (if (not (equal? e_y 0))
                      (make source
                          (src (make continuous-src (frequency freq) (width 0.5)))
                          (component Ey)
                          (amplitude e_y)
                          (size 0 3 3)
                          (center source_shift 0 0)
                          ;;(amp-func (Gauss w_0))
                          ;;(amp-func (psi_cartesian (f_Laguerre_Gauss_cartesian w_0 m_charge) shift))
                          (if (equal? m_charge 0)
                              ; if vortex charge is zero use Gauss spectrum distribution (improves perfomance)
                              (amp-func (psi_spherical (f_Gauss_spherical w_0) shift))
                              (amp-func (psi_spherical (f_Laguerre_Gauss_spherical w_0 m_charge) shift))
                          )
                      )
                  )
              ))
)

(define (eSquared r ex ey ez)
        (+ (* (magnitude ex) (magnitude ex)) (* (magnitude ey) (magnitude ey))
           (* (magnitude ez) (magnitude ez))))

(define (output-efield2) (output-real-field-function (cond (s-pol? "e2_s") (p-pol? "e2_p") (a-pol? "e2_mixed"))
                                                     (list Ex Ey Ez) eSquared))

(init-fields)

(run-until runtime
      (at-beginning (lambda () (print "Calculating inital field configuration. This will take some time...\n")))
;     (at-beginning output-epsilon)          ; output of dielectric function
;     (at-end output-efield-x)               ; output of E_x component 
;     (at-end output-efield-y)               ; output of E_y component 
;     (at-end output-efield-z)               ; output of E_z component
     (at-end output-efield2))               ; output of electric field intensity

(print "\nend time: "(strftime "%c" (localtime (current-time))) "\n")
