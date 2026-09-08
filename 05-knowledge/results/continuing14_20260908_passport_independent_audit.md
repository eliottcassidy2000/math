# Independent audit: equal full pointed torsion, different generic source curves

**Verdict: ACCEPTED analytically, with independent exact controls.** No
mathematical repair is requested. The all-degree generic-fibre
reconstruction theorem and the explicit degree-fifteen pair both pass.

Audited primary: `continuing14_20260908_same_torsion_different_fibres.md`,
SHA256 `601a55857454668fd85d902e3ee49345dd698d20212ef7ca88a454ada6e90fa7`;
its frozen source SHA256 is
`780fad38dd9865860f1415a11558e1f11120715af479d5329f27ead6a82e4033`.
Only the proof and formulas are inherited as the audit target; the
independent engine imports no producer or supplier engine. External
terminology and motivation do not enter the mathematical dependency graph.

## 1. The reconstruction iff is valid over the stated field

Let K=C(T), Q_i'=f_i^2, with f_i squarefree of degree q_i>=2. The inherited
mutation supplier identifies the original generic affine fibre exactly as

    U_i=Spec K[z,1/f_i(z),1/(T-Q_i(z))].

Its unique smooth projective completion is P1_K. The complement contains
q_i+1 distinct degree-one closed points, all constant (the q_i roots and
infinity), and one closed point of degree 2q_i+1. Indeed Q_i-T is prime
in C[z,T], since its quotient is C[z], and primitive over C[T]; Gauss's
lemma proves its irreducibility in K[z]. It is separable and disjoint
from the constant points. The argument counts closed points over K, not
geometric points after algebraic closure.

A K-isomorphism extends to a K-isomorphism of the smooth completions and
therefore preserves residue degrees of missing points. It permutes the
constant missing points and carries the unique larger-degree point to
its counterpart. In particular the two q_i agree. At least three
constant points exist. Their three distinct constant images determine
the extended Mobius transformation uniquely in PGL2(C), even though it
was initially allowed to have coefficients in K.

At the moving closed point the first residue field is C(z), with the
embedding of K specified by T=Q_1(z). Its image satisfies
Q_2(phi(z))=T, hence Q_2 o phi=Q_1 identically in C(z). No denominator
problem occurs: the pole of a constant Mobius map is a constant point,
and the moving closed point is not constant. The unique polynomial pole
then forces phi to fix infinity, so phi(z)=az+b with a,b constant and
a!=0. Conversely an affine identity Q_2(az+b)=Q_1(z) maps the derivative
root sets by differentiation and also maps the moving puncture; it
restricts to an isomorphism of the displayed rings.

Thus the exact iff is established, for the entire specified supplier
family. The absence of h and lambda from the generic ring does not
identify total fibrations or their source exceptional components. The
restriction to the same target coordinate T is essential. The proof
does not persist unchanged over an algebraic closure of K.

## 2. All hypotheses of the two actual W2 mutations are paid

Write S(w)=6w^5-15w^4+10w^3. For either root r of
6r^2-15r+10, put d=-1-r, W=dz^3+r, Q=S(W), and
f=s zW(W-1), where s^2=90d. Then Q'=f^2, deg Q=15 and
deg f=7. The zeros 0, W=0 and W=1 are disjoint and simple: r is
neither 0 nor 1, d!=0, and the relevant cubics have nonzero constant
terms. Hence f is squarefree and its zero at 0 is simple.

Q(0)=0. At the common actual source point z=1 one has W=-1,
Q(1)=-31 and f(1)=2s. For both signs of s, this is nonzero and not
-1, because 360d!=1. Thus h=-1/2 and lambda=1 simultaneously satisfy
all mutation hypotheses for both parameter roots and both square-root
choices. No source parameter is chosen separately to avoid an infinite
exceptional set.

The independent engine checks the actual W2 substitution

    z=-R[h^2(3-hR)+B(1-hR)^3], h=-1/2,

and the rational-mate chain identity J(T,1/(2A^2))=1. The inherited
all-point submersion proof therefore applies: on X=0 the normal
derivative is 2s(1+2s)!=0, and on the added divisor the tangential
derivative is -f'(0)!=0. All response statements remain in the original
C[x,t] quotient, with no implicit global Hamiltonian-ring assertion.

## 3. Full torsion really agrees, with the distinguished unit

The factorizations show four critical points of Q at value 0 and three
at value 1, each of local degree three. The remaining root counts in
Q=0 and Q=1 are three and six, respectively. The source point z=1
contributes the separate value -31. Under the mutation supplier this
gives respectively four, three and one full source torsion towers.
Each special source fibre also has its one regular component; the full
component counts are five, four and two. This does not identify the
ramification of Q with critical points of the globally smooth T.

The complete primitive parts have only double poles: coefficient 1/2
on every derivative-root component, coefficient (1+2s)^2/2 on X=0,
and zero on each regular component. At 0 and 1 the coefficient vectors
are the same for the two examples. At -31 the coefficient is nonzero,
so a constant scalar map on that one-dimensional coefficient space
carries one vector to the other. These maps commute with both the
target coordinate and its canonical derivative, and extend to every
pole order. They therefore give an isomorphism of the FULL torsion
modules which carries the distinguished unit and its embedded generated
submodule. This is stronger than merely matching scalar annihilators.

The full ambient rank is eight. The unit generates three towers, one
at each support, and has scalar annihilator

    [T(T-1)(T+31)]^2.

The local pole order is two; the unit-generated arm count is three;
the displayed annihilator has polynomial degree six. These three
different integers must not be called the same scalar order. Intrinsic
rational-pair degree is thirty and ordinary original-source puncture
count is twenty-three, as supplied and independently checked by the
degree formulas. Completed fibres are rational; this does not make
their punctured affine fibres constant in moduli.

## 4. The generic affine source curves differ

Every Q is centered of degree fifteen. An affine equivalence between
the two Q's must have zero translation by comparison of degree
fourteen, whose only possible contribution is the translated leading
term. After this, the equivalence would yield S o L=S for an affine
map L with L(r_+)=r_-. Differentiation makes L permute the critical
points 0 and 1 of S. Their distinct values force each point to be
fixed, so L is the identity, contradicting r_+!=r_-.

Independently, if a_15 and a_12 are the corresponding coefficients,

    a_12^5/a_15^4 = (15^5/6^4)(2r-1)^5

is invariant under coordinate scaling. The independent quotient-ring
engine finds the two values unequal, without approximating complex
roots. The critical-point proof explains this distinction for all
candidate coordinate scalings; the finite control is not its proof.
The tempting reflection obeys S(1-w)=1-S(w), so changes the target and
is not a fixed-target equivalence. A nonaffine constant Mobius map
would introduce a finite pole and also cannot evade the argument.

Combining this with Section 1 proves nonisomorphism over C(T). No
assertion after arbitrary algebraic base extension, no equality of
embedded pair fields across examples, no classification of arbitrary
degree-fifteen functions, and no Jacobian-conjecture conclusion is
introduced.

## 5. Independent reproduction

Run `python continuing14_20260908_passport_independent_audit.py` and
`python -O continuing14_20260908_passport_independent_audit.py`.
The engine uses polynomial remainders over Q[r]/(r^2-5r/2+5/3), complete
factorizations and discriminants, actual chart and bracket identities,
and constant module maps. It imports no primary engine and does not
use finite sampled fibres to prove the all-degree reconstruction iff.
Normal and optimized output and certificate are frozen byte-identically
with LF newlines; 92 always-active exact gates pass. The source locates its certificate beside itself or
under repository 05-knowledge/results after relocation.
