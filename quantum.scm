;;
;; This file includes the quantum algorithms in the paper
;; "A lambda calculus for quantum computation". 
;; However, in this version linearity assumptions are left 
;; implicit, and it is up to the programmer to correctly use 
;; linear arguments.  Failing to do so may result in nonsensical 
;; results.
;;
;; Copyright 2003 - Andre van Tonder

(include "quantum-definitions.scm")

;==============================================================
; First let's apply a few Hadamard gates in sequence
; to a single qubit:

(qeval (H (H (H 0))))

; => (superposition (0.7071067811865 1) 
;                   (0.7071067811865 0))

(qeval (H (H (H (H 0)))))

; => (superposition (1.0 0))


;==============================================================
; Let's define a function that makes an EPR state:

(define (make-epr) (cnot (H 0) 0))

(qeval (make-epr))

; => (superposition (0.7071067811865 (1 1)) 
;                   (0.7071067811865 (0 0)))


;================================================================
; The Deutsch algorithm is easily defined:
; Here Uf is the oracle corresponding to an unknown 
; function f : Bit -> Bit.
; Note that all variables appear linearly.

(define (deutsch Uf)
  (bind* ([x (H 0)]
          [y (H 1)]
          [(x* y*) (Uf x y)])
    (list (H x*) (H y*))))

; For example, 
; Uf = cnot corresponds to f(0) = 0, f(1) = 1: 

(qeval (deutsch cnot))

; => (superposition (1.0 (1 1)))
;    The first qubit is indeed 1 = f(0) + f(1) mod 2.

; Another example,
; Uf = (lambda (x y) (list x y)) corresponds to f(0) = f(1) = 0:

(qeval (deutsch 
        (lambda (x y) (list x y))))

; => (superposition (1.0 (0 1)))
;    The first qubit is indeed 0 = f(0) + f(1) mod 2.


;==================================================================
; Quantum teleportation:  Once again all variables appear linearly.

(define (teleport x)
  (bind* ([(e1 e2)  (make-epr)]
          [(x* e1*) (alice x e1)])
    (bob x* e1* e2)))

(define (alice x e)
  (bind ([(x* e*) (cnot x e)])
    (list (H x*) e*)))

(define (bob x e1 e2)
  (bind* ([(e1* e2*) (cX e1 e2)]
          [(x* e2**) (cZ x e2*)])
    (list x* e1* e2**)))


; Let's teleport a qubit in the state (|0> - |1>)/sqrt(2):

(qeval (teleport (H 1)))

; => (superposition (-0.3535533905933 (1 0 1))
;                   (-0.3535533905933 (0 0 1))
;                   (-0.3535533905933 (1 1 1))
;                   (-0.3535533905933 (0 1 1))
;                   (0.3535533905933  (1 1 0))
;                   (0.3535533905933  (0 1 0))
;                   (0.3535533905933  (1 0 0))
;                   (0.3535533905933  (0 0 0)))
; Ignoring the first two bits, the third bit, belonging to Bob, is 
; now in the state (|0> - |1>)/sqrt(2)


;==================================================================
; Creating a uniform superposition:

(qeval (map (!suspend H) '(0 0 0)))

; => (superposition (0.3535533905933 (1 1 1))
;                   (0.3535533905933 (1 1 0))
;                   (0.3535533905933 (1 0 1))
;                   (0.3535533905933 (1 0 0))
;                   (0.3535533905933 (0 1 1))
;                   (0.3535533905933 (0 1 0))
;                   (0.3535533905933 (0 0 1))


;====================================================================
; Quantum Fourier transform:

(define (fourier lst)
  (reverse (fourier* lst)))

(define (fourier* lst)
  (list-match lst
    [()        '()]
    [(hd . tl) (bind ([(hd* . tl*) (phases (H hd) tl 2)])
                 (cons hd* (fourier* tl*)))]))

(define (phases target controls n)
  (list-match controls
    [() (list target)]
    [(control . tl) 
     (bind* ([(control* target*) ((cR n) control target)]
             [(target** . tl*)   (phases target* tl (add1 n))])
       (cons target** (cons control* tl*)))]))


; Testing the Fourier transform:

(qeval (fourier '(1 0)))

; => (superposition (-0.5 (1 1)) 
;                   (-0.5 (0 1)) 
;                   (0.5  (1 0)) 
;                   (0.5  (0 0)))

(qeval (fourier '(1 1 1)))

; => (superposition (0.25+0.25i            (1 1 1))
;                   (-0.25-0.25i           (0 1 1))
;                   (-0.25+0.25i           (1 0 1))
;                   (0.25-0.25i            (0 0 1))
;                   (0.0+0.3535533905933i  (1 1 0))
;                   (-0.0-0.3535533905933i (0 1 0))
;                   (-0.3535533905933      (1 0 0))
;                   (0.3535533905933       (0 0 0)))



















































