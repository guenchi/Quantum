# Quantum

This project is forked form http://www.het.brown.edu/people/andre/qlambda/index.html


***Scheme simulator***

This page provides a simulator for a functional language based on Scheme for expressing and simulating quantum algorithms. Example code implementing some simple quantum algorithms is provided. For the theory, see my papers

A lambda calculus for quantum computation - André van Tonder.

Quantum computation, categorical semantics and linear logic - André van Tonder.

which develop a linear lambda calculus for expressing quantum algorithms. The code below implements the examples in the first one of these two papers.

***Instructions***

The code should work in any Scheme adhering to the R5RS standard. It has been tested in PLT's DrScheme and Petite Chez Scheme, both of which are freely available.

Instructions are as follows for DrScheme:

Install DrScheme . When asked during installation, choose the "Pretty Big" language level.
Save the two files
quantum.scm
quantum-definitions.scm
wherever you like, but they should be in the same directory.
Launch DrScheme, load the file quantum.scm, and click on the "execute" button to run the algorithms.
Examples:

Simple gate combinations

First let's apply a few Hadamard gates in sequence to a single qubit. The QEVAL performs the sum over histories:
  (qeval (H (H (H 0))))

      =>  (superposition (0.7071067811865 1) 
                         (0.7071067811865 0))
We see that the answer is a superposition of 0 and 1 with the indicated coefficients. Another simple check:
  (qeval (H (H (H (H 0)))))

    => (superposition (1.0 0))
Let's define a function that makes an EPR state:
  (define (make-epr) (cnot (H 0) 0))
Evaluating this function we get
  (qeval (make-epr))

    => (superposition (0.7071067811865 (1 1)) 
                      (0.7071067811865 (0 0)))
The Deutsch algorithm

The Deutsch algorithm is easily defined: Here Uf is the oracle (a unitary transformation of two qubits) corresponding to an unknown function f : Bit -> Bit. Note the use of the pattern matching in the third line of the BIND* form (which corresponds to LET* in the article) to linearly extract the constituent qubits.
  (define (deutsch Uf)
    (bind* ([x       (H 0)]
            [y       (H 1)]
            [(x* y*) (Uf x y)])     
      (list (H x*) 
            (H y*))))
For example, Uf = cnot corresponds to f(0) = 0, f(1) = 1:
  (qeval (deutsch cnot))

    => (superposition (1.0 (1 1)))
and the first qubit in the answer is indeed 1 = f(0) + f(1) mod 2.
Another example, Uf = (lambda (x y) (list x y)) corresponds to f(0) = f(1) = 0:

  (qeval (deutsch 
           (lambda (x y) (list x y))))

    => (superposition (1.0 (0 1)))
 
Now the first qubit in the answer is indeed 0 = f(0) + f(1) mod 2.
Quantum teleportation

Quantum teleportation (with deferred measurement) is easily defined by transcribing the corresponding quantum gate network. Once again, we use the pattern matching BIND* form to linearly extract the qubits.
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
Let's teleport a qubit in the state (0 - 1)/sqrt(2):
  (qeval (teleport (H 1)))

    => (superposition (-0.3535533905933 (1 0 1))
                      (-0.3535533905933 (0 0 1))
                      (-0.3535533905933 (1 1 1))
                      (-0.3535533905933 (0 1 1))
                      (0.3535533905933  (1 1 0))
                      (0.3535533905933  (0 1 0))
                      (0.3535533905933  (1 0 0))
                      (0.3535533905933  (0 0 0)))
       
Notice how linearity requires us to keep all the qubits in the answer (see the article for more discussion of this). Ignoring the first two bits, the third bit, belonging to Bob, is now in the state (0 - 1)/sqrt(2).
Higher order functions

Creating a uniform superposition is easy using the higher-order function MAP, which applies a given transformation to each element in a list.

  (define (map f lst)
    (list-match lst 
      [()        '()] 
      [(hd . tl) (cons ((!resume f) hd) 
                       (map f tl))]))
Note the use of pattern matching to linearly extract the constituents of the list. Also note that, since the parameter f appears nonlinearly in the MAP procedure, it is assumed to have been wrapped in a !-suspension. The !RESUME form extracts the wrapped function. To apply the H gate to each element of the superposition, we proceed as follows:
  (qeval (map (!suspend H) '(0 0 0)))

    => (superposition (0.3535533905933 (1 1 1))
                      (0.3535533905933 (1 1 0))
                      (0.3535533905933 (1 0 1))
                      (0.3535533905933 (1 0 0))
                      (0.3535533905933 (0 1 1))
                      (0.3535533905933 (0 1 0))
                      (0.3535533905933 (0 0 1))
Quantum Fourier transform

The Quantum Fourier transform may be defined as follows:
   
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
Testing the Fourier transform:
  (qeval (fourier '(1 0)))

    => (superposition (-0.5 (1 1)) 
                      (-0.5 (0 1)) 
                      (0.5  (1 0)) 
                      (0.5  (0 0)))

  (qeval (fourier '(1 1 1)))

    => (superposition (0.25+0.25i            (1 1 1))
                      (-0.25-0.25i           (0 1 1))
                      (-0.25+0.25i           (1 0 1))
                      (0.25-0.25i            (0 0 1))
                      (0.0+0.3535533905933i  (1 1 0))
                      (-0.0-0.3535533905933i (0 1 0))
                      (-0.3535533905933      (1 0 0))
                      (0.3535533905933       (0 0 0)))
    
Simulator primitives

The simulator relies on a quantum version of McCarthy's nondeterministic operator AMB, which ambiguously returns one of its arguments. Our version, called QUAMB, keeps track of a complex amplitude for each history in the quamb-tree. The quantum-eval macro then performs the quantum mechanical sum over histories.
QUAMB is easy to use. For example, the Hadamard and the Pauli X and Z operations are defined as (see the file quantum-definitions.scm above)

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
The user may define additional primitives and oracles using the quamb operator.

However, the user should be aware that functions that use QUAMB, that branch according to the value of a qubit, or that use arguments referring to values containing qubits nonlinearly, are not syntactically valid functions of the quantum lambda calculus and are not guaranteed to be unitary.

However, there is nothing wrong with regarding such functions as extensions of the simulator. We may call such functions "Oracles" to distinguish them from proper functions in the quantum lambda calculus. The user may of course define additional oracles to suit his or her purposes, but should ensure that they define unitary operations.
