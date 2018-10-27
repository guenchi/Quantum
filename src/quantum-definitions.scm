;;
;; This file defines a simulator for testing the quantum algorithms 
;; in the paper "A lambda calculus for quantum computation". 
;; It relies on a quantum version of McCarthy's nondeterministic
;; operator amb, which ambiguously returns one of its arguments.
;; Our version, called quamb, keeps track of a complex amplitude
;; for each history in the quamb-tree.
;; The quantum-eval macro performs the quantum mechanical sum over 
;; histories.
;;
;; We then define some quantum primitives in terms of quamb.  
;; The user may similarly define more primitives and oracles
;; using the same techniques.  See discussion below for details.
;;
;; Copyright 2003 - Andre van Tonder


;=========================================================================
; Deconstructors for writing linear algorithms:

(define-syntax list-match
  (syntax-rules ()
    [(list-match exp 
       [()      exp1]
       [(h . t) exp2])
     (let ([lst exp])
       (cond [(null? lst) exp1]
             [(pair? lst) (let ([h (car lst)]
                                [t (cdr lst)])
                            exp2)]
             [else 'list-match-error]))]))

; A list-deconstructing version of let

(define-syntax bind
  (syntax-rules ()
    [(bind () body) 
     body]
    [(bind ((pat1 exp1) (pat2 exp2) . rest) body)
     (bind ((pat1 exp1)) (bind ((pat2 exp2) . rest) body))]
    [(bind ((() exp)) body)
     (list-match exp
       (()      body)
       ((h . t) 'bind-error))]
    [(bind (((pat1 . pat2) exp)) body)
     (list-match exp 
       (()      'bind-error)
       ((h . t) (bind ((pat1 h) 
                       (pat2 t)) 
                    body)))]
    [(bind ((identifier exp)) body)
     (let ((identifier exp)) body)]))

(define-syntax bind*
  (syntax-rules ()
    [(bind* . rest) (bind . rest)]))

;==============================================================

; Simulator infrastructure:
; First we define a quantum version of McCarthy's amb-operator:

(define current-history-amplitude 1)

(define (quamb-fail) 
  (error "History tree exhausted")(newline))

; Included for portability:

(define (error msg)
  (display msg) (car 1))

; Macro that provides programmer interface of quamb from 
; implementation

(define-syntax quamb
  (syntax-rules ()
    [(quamb (coeff label) ...)
     (quamb-proc (list coeff label) ...)]))
    
(define (quamb-proc . branches)
  (let ([prev-quamb-fail quamb-fail]
        [prev-amplitude current-history-amplitude])
    (call-with-current-continuation
     (lambda (cont)
       (for-each (lambda (branch)
                   (bind (((coeff label) branch))
                     (call-with-current-continuation
                      (lambda (cont*) 
                        (set! current-history-amplitude (* coeff prev-amplitude))
                        (set! quamb-fail (lambda () (cont* '())))
                        (cont label)))))
                 branches) 
       (set! current-history-amplitude prev-amplitude)
       (prev-quamb-fail)))))

; The following macro collects all branches of an evaluation 
; to obtain a quantum superposition: 
;
; quantum-eval : Expr -> Superposition

(define-syntax qeval
  (syntax-rules ()
    [(_ expr) 
     (begin
       (let ([prev-quamb-fail quamb-fail]
             [result '()])
         (call-with-current-continuation 
          (lambda (cont)
            (set! quamb-fail
                  (lambda () (cont '())))  ; will fall through after last history
            (let ([value expr]) 
              (set! result (add (list current-history-amplitude value)
                                result)))
            (quamb-fail)))                 ; evaluate next history or fall through
         (set! quamb-fail prev-quamb-fail)
         (cons 'superposition (superposition-round result))))]))
                      
; Infrastructure for manipulating superpositions:
     
; Branch = (cons Number Obj)
; Superposition = [Branch]

; add : (Branch, Superposition) -> Superposition 
     
(define (add branch superposition)
  (list-match superposition
    [() (list branch)]
    [(branch* . rest)
     (bind ([(coeff label)   branch]
            [(coeff* label*) branch*])
       (if (equal? label label*) 
           (cons (list (+ coeff coeff*) label*) rest)
           (cons branch* (add branch rest))))]))

; superposition-round : Superposition -> Superposition

(define (superposition-round superpos)
  (foldl (lambda (branch accum)
           (bind* (((coeff label) branch)
                   (rounded-coeff (complex-round coeff *epsilon*)))
             (if (zero? rounded-coeff)
                 accum
                 (cons (list rounded-coeff label)
                       accum))))
         '()
         superpos))

; This may not be included in your favorite Scheme:

(define (foldl f base lst)
  (list-match lst 
    [()      base]
    [(h . t) (foldl f (f h base) t)]))

(define (real-round x epsilon) 
  (* epsilon (round (/ x epsilon))))

(define (complex-round z epsilon)
  (make-rectangular (real-round (real-part z) epsilon)
                    (real-round (imag-part z) epsilon)))

; Accuracy 

(define *epsilon* 0.0000000000001)


;====================================================================
;; This completes the skeleton of the simulator.  
;; Now we can define some quantum primitives.  
;; These primitives should be regarded as part of the 
;; simulator.  In particular, they are not valid functions in the 
;; quantum lambda calculus, since their definitions make use of the 
;; simulator primitive quamb, may branch according to the value of
;; a qubit, and may use qubit arguments nonlinearly, all of which
;; are not part of the quantum lambda calculus.
;; We may call such functions "Oracles" to distinguish them from proper 
;; functions in the quantum lambda calculus.  The user may of course
;; define additional oracles to suit his or her purposes, but should 
;; ensure that they define unitary operations.  

;; Some useful definitions
 
(define 1/sqrt2 (/ 1 (sqrt 2)))
(define -1/sqrt2 (- 1/sqrt2))

(define (H a)
  (case a
    [(0) (quamb (1/sqrt2 0) (1/sqrt2 1))]
    [(1) (quamb (1/sqrt2 0) (-1/sqrt2 1))]))
   
(define (X a)
  (case a
    [(0) 1]
    [(1) 0]))

(define (Z a)
  (case a
    [(0) 0]
    [(1) (quamb (-1 1))]))

(define (controlled op)
  (lambda (a b)
    (list a 
          (case a
            [(0) b]
            [(1) (op b)]))))
  
(define cX (controlled X))
(define cZ (controlled Z))

(define cnot cX)

(define pi 3.141592653589793)

(define (cR n)
  (controlled 
    (lambda (a)
      (case a
        [(0) 0]
        [(1) (quamb ((exp (/ (* 2 pi 0+i) (expt 2 n))) 1))]))))

;=================================================================
; Constructor for creating !-suspensions (see paper)
; and extracting expression from under !. 

(define-syntax !suspend
  (syntax-rules ()
    [(!suspend exp) (lambda () exp)]))

(define (!resume susp)
  (susp))

;================================================================

;;
;; Here we provide some additional useful definitions:
;;

; Override primitive map, which does not
; interact well with call/cc in MzScheme.
; Parameter f occurs nonlinearly and must therefore be a
; !-suspension.  

(define (map* f lst)
  (list-match lst 
    [()        '()] 
    [(hd . tl) (cons ((!resume f) hd) (map* f tl))]))
       
(define map map*)


; Replace the primitive append with the following linear version:

(define (append* lst1 lst2)
  ((list-match lst1
     [()        (lambda (u) u)]
     [(hd . tl) (lambda (u) (cons hd (append* tl u)))]) 
   lst2))
   
(define append append*)

; Linear version of reverse:

(define (reverse* lst)
  (list-match lst
    [()        '()]
    [(hd . tl) (append (reverse* tl) (list hd))]))

(define reverse reverse*)


                
 
