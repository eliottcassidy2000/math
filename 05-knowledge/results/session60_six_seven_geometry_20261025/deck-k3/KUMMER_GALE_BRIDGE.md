# The central six-root K3 itself lies in the genus-two Kummer model

**Status: PROVED over C, using the CITED classical diagonal Kummer model.**
This is stronger than showing its symplectic quotient is Kummer. It makes
no untwisted assertion over Q.

## 1. General six-point conic criterion via an explicit matrix transform

Suppose six distinct branch-line normals in the dual P^2 lie on a
nonsingular conic. Choose a projective coordinate on that conic avoiding
infinity at all six points. After a base linear transformation and
individual nonzero column scalings, the six normal columns become

    v_i=(1,s_i,s_i^2)^t,    s_i distinct, i=1,...,6.

The corresponding root rescalings require square roots of those scalars,
which is why the assertion is made over C. Let

    P(T)=product_i(T-s_i),    w_i=1/P'(s_i).

Lagrange interpolation, by comparing the T^5 coefficient for the
polynomials 1,T,...,T^4, gives

    sum_i w_i s_i^k = 0,    k=0,...,4.                    (G1)

Consequently the three rows (w_i s_i^a)_i, a=0,1,2, form the entire
relation space among the six columns v_i: they annihilate the three
rows (s_i^b)_i, b=0,1,2, and have rank three. Thus the maximal six-root
surface has equations

    sum_i w_i s_i^a r_i^2=0,    a=0,1,2.

Put z_i=sqrt(w_i) r_i. Its equations become exactly

    sum_i z_i^2 = sum_i s_i z_i^2 = sum_i s_i^2 z_i^2 = 0.  (G2)

Garbagnati--Sarti, Section 5.1, equation (2), printed p.22, and
Section 10.2, printed p.55, identify (G2) as the
classical Kummer model associated with the genus-two curve
C: y^2=product_i(x-s_i). This primary model is the only imported
identification; (G1) is the explicit map carrying our surface into it.
Primary source:
https://sarti.pages.math.cnrs.fr/home/GarbagnatiSarti_RMJ.pdf
https://arxiv.org/abs/1305.3514

This does NOT reason that a surface is Kummer merely because it has a
Kummer quotient. It identifies its actual three defining quadrics with
an established Kummer model by invertible changes of coordinates.

## 2. Apply it to the hourglass central deletion

Deleting E leaves the six dual normals

    (1,1,0),(1,-1,0),(1,0,1),(1,0,-1),(1,1,1),(1,-1,-1).

They lie on the nonsingular conic

    A^2-B^2-C^2+BC=0.                                    (G3)

Its determinant is 3/4, hence nonzero. Therefore the central six-root
surface S_0 is a Kummer K3 over C, not just related to one by a finite
symplectic quotient.

An exact convenient parameter on (G3) in A=1 is the slope t of the line
through (B,C)=(1,0): C=t(B-1). The resulting parametrization is

    B=(t^2-1)/(t^2-t+1),
    C=t(t-2)/(t^2-t+1).

The six parameters, in the normal order displayed above, are respectively

    2, 0, -1, 1, infinity, 1/2.

The value t=2 is the tangent limit at (1,0). Apply the Mobius map
s=1/(t-3), obtaining six finite distinct parameters

    -1, -1/3, -1/4, -1/2, 0, -2/5.                       (G4)

The script verifies (G3), the point-parameter correspondence, and every
identity (G1) with the exact rational weights computed from (G4).
The separate symmetry/CM packet identifies this genus-two curve over C
with v^2=u^6-1 and computes its Jacobian's Picard contribution. The
square-root scalings above and that complex curve isomorphism need not
be defined over Q; rational twists and integral lattices are not discarded.

## 3. The six is structural

The two rank-three matrices in the proof are Gale dual, and the
self-duality here depends on six=2*3 and the degree-four vanishing range
in (G1). This is a genuine bridge from six branch normals on a conic to
three diagonal quadrics and a genus-two six-branch-point curve, not a
claim that every six-root K3 is Kummer. A general six-point configuration
in P^2 does not lie on a conic, so the required matrix transform fails.
