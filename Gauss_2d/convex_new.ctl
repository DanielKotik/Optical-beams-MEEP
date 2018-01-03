;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;; command line argument: meep s-pol\?=false ....ctl 
;; coordinate system in meep:
;; --|-----> x
;;   |
;;   |
;;   v y
;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

(set! force-complex-fields? true)

(define-param sx 5) ; size of cell in X direction
(define-param sy 5) ; size of cell in Y direction

(define-param runtime 10) ; time to run for

(define-param s-pol? true) ; true for s-spol, false for p-pol

(define-param n1 1.54) ; index of refraction of the denser medium
(define-param n2 1.00) ; index of refraction of the thinner medium

(define (Critical _n1 _n2)  ; calculates critical angle in degrees
        (* (/ (asin (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define (Brewster _n1 _n2)  ; calculates Brewster angle in degrees
        (* (/ (atan (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define (chi_rad _chi_deg) ; Umrechung Grad --> Radiant
        (* (/ _chi_deg 360.0) (* 2.0 pi)))

(define-param chi_deg  (* 1.0 (Brewster n1 n2))) ; incidence angle in degrees
;(define-param chi_deg  40.0)
(define-param krw 0.0)  ; beam waist distance to interface (30 is good)
(define-param kw_0 5)   ; beam waist (10 is good)
(define-param kr_c 470)  ; radius of curvature
(define-param freq 12)   ; vacuum frequency of source (10 is good)
(define-param pixel 20)  ; number of pixels per wavelength in the denser medium
                         ; (at least >10; 20 is a good choice)

; Calculation of resolution parameter assuming n1 > n2                         
(define-param resol (* pixel (* n1 freq)))   

(define-param k_vac (* 2.0 pi freq))
(define-param rw (/ krw (* 1.00 k_vac)))
(define-param w_0 (/ kw_0 k_vac))
(define-param r_c (/ kr_c k_vac))

(define-param source_shift -2.15)  ; source position with respect to the center
                                   ; (point of impact) in Meep units (-2.15 good)
                                   ; if equal -rw, then source position 
                                   ; coincides with waist position
                                   
; distance from source position to beam waist                                   
;(define-param shift (/ source_shift 2.0))
(define-param shift (+ source_shift rw))

;(define-param kdir (vector3 1 0 0)) ; muss Einheitsvektor sein

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! default-material (make dielectric (index n1)))

(set! geometry (list
				(make cylinder
				(center (* r_c (cos (chi_rad chi_deg))) (* -1.0 (* r_c (sin (chi_rad chi_deg))))) 
				                      ; Mittelpunkt wird nach rechts verschoben,
                                      ; so dass der Auftreffpkt immer mittig liegt
				(height infinity)
				(radius r_c)
				(material (make dielectric (index n2))))))

;; ***** Definition of beam profiles at origin of source *****

(define (Complexfactor k sigma shift)
		(+ 1.0 (/ (* 0+2i  shift ) (* k sigma sigma))))

(define (Gaussian sigma k)
		(lambda (r) (* (/ 1.0 (sqrt (Complexfactor k sigma shift))) 
                (exp (- (* 0+i k shift)
				(/(/ (* (vector3-y r) (vector3-y r)) (* sigma sigma)) 
	            (Complexfactor k sigma shift)))))
		))

(define (Uniform sigma)
        ;(lambda (r) (if (< (abs (r)) sigma) (1.0) (0.0))
        (lambda (r) (if (< (abs (vector3-y r)) sigma) 1.0  0.0)
		))
		
(define (Academic sigma)
        (lambda (r) (if (< (abs (vector3-y r)) sigma) 
					    (- 1.0 (abs (/ (vector3-y r) sigma))) 0.0) 
		))
		
(define (Asymmetric sigma)
        (lambda (r) (if (< (vector3-y r) (* -1.5 sigma)) 0.0
					    (* (+ (/ (vector3-y r) sigma) 1.5) 
		                   (exp (* -1.0 (+ (/ (vector3-y r) sigma) 1.5)))))
        ))

;; ***** definition of different spectrum amplitudes f *****

(define (f_Gauss w_0)
        (lambda (k_y) (* (/ w_0 (* 2.0 (sqrt pi)))
                        (exp (* -1.0 (expt (* 0.5 k_y w_0) 2.0))))
        ))
		
(define (f_asymmetric a b)
        (lambda (k_y) (* (/ w_0 (* 2.0 (sqrt pi)))
                        (exp (* -1.0 (expt (* 0.5 k_y w_0) 2.0))))
        ))

(define (f_uniform w_0)
        (lambda (k_y) (/ (sin (* k_y w_0)) (* pi k_y))
		))

(define (integrand f y x k)
        (lambda (k_y) (* (f k_y)
                        (exp (* 0+1i x (sqrt (- (* k k) (* k_y k_y))))) 
                        (exp (* 0+1i k_y y))) 
        ))

; complex field amplitude at position (x, y) with a spectrum aplitude f
; (one may have to adjust the 'relerr' value of the integrand function)
(define (psi f x k)
        (lambda (r) (car (integrate (integrand f (vector3-y r) x k) 
						            (* -1.0 k) (* 1.0 k) 0.0001))  ;1.49e-8
        )) 

(print "w_0 " w_0 "\n")
(print "k_vac " k_vac "\n")
(print "The value of our Gaussian spectrum amplitude is: "
       ((f_Gauss w_0) 20.0) "\n")


;(print "integrand " ((integrand 0.8 2.0 k_vac w_0) 20.0) "\n")
;(print "Field amplitude: " ((psi 1.0 k_vac w_0) 0.5) "\n")


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; one may have to adjust the y-component of the size vector 
(set! sources (list
				(make source
					(src (make continuous-src (frequency freq) (width 0.5)))
					(if s-pol? (component Ez) (component Ey))
					(amplitude 3.0)
					(size 0 3.0 0)
					(center source_shift 0.0)
					;(amp-func (Gaussian w_0 (* n1 k_vac))))
					(amp-func (psi (f_Gauss w_0) shift (* n1 k_vac))))
				(make source
					(src (make continuous-src (frequency freq)))
					(if s-pol? (component Hy) (component Hz))
					(if s-pol? (amplitude -3.0) (amplitude 3.0))
					(size 0 3.0 0) 
					(center source_shift 0.0) 
					;(amp-func (Gaussian w_0 (* n1 k_vac))))
					(amp-func (psi (f_Gauss w_0) shift (* n1 k_vac))))
				))

(set! pml-layers (list (make pml (thickness 0.25))))

(set! resolution resol)

(run-until runtime
	(if s-pol?
			(at-end output-efield-z)
			(at-end output-efield-y)))

;; (run-until runtime
;;  (at-beginning output-epsilon)
;;  (at-end output-efield-z))