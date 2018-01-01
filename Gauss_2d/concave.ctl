(set! force-complex-fields? true)

(define-param sx 5) ; size of cell in X direction
(define-param sy 5) ; size of cell in Y direction

(define-param runtime 20) ; time to run for

(define-param s-pol? true) ; true for s-spol, false for p-pol

(define-param n1 1.54) ; index of refraction of the denser medium
(define-param n2 1.00) ; index of refraction of the thinner medium

(define (Critical _n1 _n2)  ; calculates critical angle in degrees
        (* (/ (asin (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define (Brewster _n1 _n2)  ; calculates Brewster angle in degrees
        (* (/ (atan (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define-param chi_deg (* 0.5 (Critical n1 n2))) ; incidence angle in degrees
(define-param krw 0.5) ; beam waist distance to interface (dimensionless)
(define-param kr_c 3.0) ; radius of curvature (dimensionless)
(define-param kw_0 0.212) ; beam waist (dimensionless)
(define-param shift 0.5)

(define-param kdir (vector3 1 0 0)) ; muss Einheitsvektor sein

(define (chi_rad _chi_deg) ; Umrechung Grad --> Radiant
        (* (/ _chi_deg 360.0) (* 2.0 pi)))

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! default-material (make dielectric (index n2)))

(set! geometry (list
				(make cylinder
				(center (* -1 (* kr_c (cos (chi_rad chi_deg)))) (* kr_c (sin (chi_rad chi_deg)))) 
				                      ; Mittelpunkt wird nach rechts verschoben,
                                      ; so dass der Auftreffpkt immer mittig liegt
				(height infinity)
				(radius kr_c)
				(material (make dielectric (index n1))))))

(define (Complexfactor k r sigma shift)
		(+ 1.0 (/ (* 0+1i (- shift (vector3-x r))) (* pi (vector3-norm k) sigma sigma))))

(define ((Gaussian sigma k) r)
		(* (/ 1.0 (sqrt (Complexfactor k r sigma shift))) (exp (- (* 0+2i pi (vector3-dot k r))
				(/(/ (vector3-dot r r) (* sigma sigma)) (Complexfactor k r sigma shift))))))

(set! sources (list
				(make source
					(src (make continuous-src (frequency 15)))
					(if s-pol? (component Ez) (component Ey))
					(size 0 1.0 0)
					(center (* -1.0 krw) 0.0)
					(amp-func (Gaussian kw_0 kdir)))))

(set! pml-layers (list (make pml (thickness 0.25))))

(set! resolution 100)

(run-until runtime
	(if s-pol?
			(at-end output-efield-z)
			(at-end output-efield-y)))
