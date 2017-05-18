(set! force-complex-fields? true)

(define-param sx 5) ; size of cell in X direction
(define-param sy 5) ; size of cell in Y direction

(define-param runtime 10) ; time to run for

(define-param s-pol? false) ; true for s-spol, false for p-pol

(define-param n1 1.54) ; index of refraction of the denser medium
(define-param n2 1.00) ; index of refraction of the thinner medium

(define (Critical _n1 _n2)  ; calculates critical angle in degrees
        (* (/ (asin (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define (Brewster _n1 _n2)  ; calculates Brewster angle in degrees
        (* (/ (atan (/ _n2 _n1)) (* 2.0 pi)) 360.0))

(define-param chi_deg  (* 1.0 (Brewster n1 n2))) ; incidence angle in degrees
(define-param krw 0.5)   ; beam waist distance to interface (dimensionless)
(define-param kw_0 0.1)  ; beam waist (dimensionless)
(define-param shift 0.0)
(define-param freq 10)   ; frequency of source
(define-param pixel 20)  ; number of pixels per wavelength in the denser medium
                         ; (at least >10; 20 is a good choice)
; Calculation of resolution parameter assuming n1 > n2                         
(define-param resol (* pixel (* n1 freq)))                        

(define (alpha _chi_deg) ; Winkel der schiefen Ebene
        (- (/ pi 2.0) (* (/ _chi_deg 360) 2 pi)))
(define (Delta_x _alpha) ; Versatz der schiefen Ebene zum Mittelpunkt
        (* (/ sx 2.0) (/ (-(- (sqrt 2.0) (cos _alpha)) (sin _alpha)) (sin _alpha))))

(define-param kdir (vector3 1 0 0))  ; muss Einheitsvektor sein

(set! geometry-lattice (make lattice (size sx sy no-size)))

(set! default-material (make dielectric (index n1)))

(set! geometry (list
				(make block
				(center (+ (/ sx 2.0) (Delta_x (alpha chi_deg))) (/ sy -2.0)) 
				                                    ; located at lower right
                                                    ; edge for 45 degree tilt
				(size infinity (* (sqrt 2.0) sx) infinity)
                    (e1 (/ 1.0 (tan (alpha chi_deg))) 1 0)
                    (e2 -1 (/ 1.0 (tan (alpha chi_deg))) 0)
                    (e3 0 0 1)
				(material (make dielectric (index n2))))))

(define (Complexfactor k r sigma shift)
		(+ 1.0 (/ (* 0+1i (- shift (vector3-x r))) (* pi (vector3-norm k) sigma sigma))))

(define (Gaussian sigma k)
		(lambda (r) * (/ 1.0 (sqrt (Complexfactor k r sigma shift))) (exp (- (* 0+2i pi (vector3-dot k r))
				(/(/ (* (vector3-y r) (vector3-y r)) (* sigma sigma)) (Complexfactor k r sigma shift))))))

(set! sources (list
				(make source
					(src (make continuous-src (frequency freq)))
					(if s-pol? (component Ez) (component Ey))
					(amplitude 3.0)
					(size 0 3.0 0) 
					(center (* -1.0 krw) 0 0) 
					(amp-func (Gaussian kw_0 kdir)))
				;(make source
				;	(src (make continuous-src (frequency freq)))
				;	(if s-pol? (component Hy) (component Hz))
				;	(if s-pol? (amplitude -3.0) (amplitude 3.0))
				;	(size 0 3.0 0) 
				;	(center (* -1.0 krw) 0.0) 
				;	(amp-func (Gaussian kw_0 kdir))))
				))


(set! pml-layers (list (make pml (thickness 0.25))))

(set-param! resolution resol)

(run-until runtime
	(if s-pol?
			(at-end output-efield-z)
			(at-end output-efield-y)))

;; (run-until runtime
;;  (at-beginning output-epsilon)
;;  (at-end output-efield-z))