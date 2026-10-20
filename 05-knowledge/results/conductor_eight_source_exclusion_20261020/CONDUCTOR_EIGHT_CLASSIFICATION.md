# The conductor-eight source frontier

20 October 2026. **PROVED; independently audited PASS in the stated scope.**

## Outcomes and exact scope

Let k have characteristic zero and let A=k[p(t),q(t)] be a subalgebra of k[t]
with fraction field k(t) and delta dim_k(k[t]/A)=4. Then:

1. Its degree semigroup is exactly one of <2,9>, <3,5>, <4,5,6>.
2. All three semigroups are realized by immersed polynomial plane curves.
3. After ordinary triangular degree reductions over k, the only nonprincipal
   terminal coordinate degrees, up to swapping, are respectively
   (2,9) or (6,9); (3,5); and (4,6).
4. Each listed terminal degree pattern is genuinely polynomially realizable.
   Characteristic admissibility by itself is not used to infer this.
5. For R=C[x], no P,Q in R[y] with one positive-degree coordinate monic,
   J(P,Q) in C*, and generic degree semigroup <4,5,6> exist.
6. Any actual Keller pair with a monic quadratic coordinate is a polynomial
   automorphism, and in particular cannot be a nonprincipal terminal pair.

Statements 1--4 are about primitive polynomial plane algebras. Statement 5
uses the actual constant Jacobian and the inherited constant-leading reduction,
not an unsupported globalization of semigroup generators. Statement 6 does
not rule out the different terminal (6,9) embedding of a <2,9> algebra.

Consequently the already audited actual-source conductor-eight frontier is
reduced to the two terminal patterns (3,5) and (6,9). This note does not exclude
either of them, does not raise the full lower bound above eight, and does not
resolve the planar Jacobian conjecture.

## Inheritance and live concepts

Closest proved mechanism: the trace/residue symmetry proof and degree-gap
codimension identity in the frozen low-delta packet. The inherited
constant-leading gate permits target degree reduction over C for actual
Keller pairs. The fixed cubic/quartic packet supplies the degree-separated
free-module method.

Canonical hostile: a characteristic sequence can be planar in the abstract
without a polynomial parametrization of that particular sequence. The
Sathaye--Stenerson introduction specifically distinguishes these questions.
Here explicit parametrizations, not the converse planar-semigroup theorem,
prove every claimed polynomial realization.

Corrected near miss: a quadratic element over C(x) is not a quadratic member
of a terminal Keller pair; terminal (6,9) cannot be discarded by a
coordinate-quadratic argument. A concrete immersed realization is below.

Least-used sidecar: a quartic/sextic presentation has an exact four-element
free module; the first odd remainder degree reads its complete codimension.
The constant term of the *whole* Jacobian then makes its linear quartic
coefficient a polynomial unit.

Concept board: symmetric genus-four gaps; characteristic gcd drops;
polynomial realizability; odd remainder module bases; global coefficient
recovery from the actual terminal pair; constant-Jacobian unit constraints.

The anchor is exact conductor-eight classification. The niche is the
quartic/sextic terminal embedding. The wildcard is the nonic alternate
embedding of the quadratic semigroup. All three remain within one frontier.

## 1. Symmetric numerical semigroups of genus four

The frozen low-delta proof gives symmetry and largest gap 2 delta-1=7.
Thus exactly one of n and 7-n belongs to the semigroup for 0<=n<=7.
The integer 1 is a gap, since otherwise A=k[t], so 6 belongs.

If 2 belongs, closure forces 4 to belong, symmetry excludes 3 and 5,
and the gap set is {1,3,5,7}. This is <2,9>.
If 2 does not belong, 5 belongs. Choosing 3 gives gap set {1,2,4,7}
and semigroup <3,5>. Choosing 4 instead gives gap set {1,2,3,7}
and semigroup <4,5,6>.

These cases are exhaustive, not a finite search conjecture. The checker
independently lists all symmetric closed subsets below conductor eight.

## 2. The characteristic input and all nonprincipal patterns

The following classical one-place plane-curve facts are **CITED** from
Sathaye--Stenerson, Section 3.1, printed pages 8--9. For a characteristic
sequence (delta_0,...,delta_h), put

    d_i=gcd(delta_0,...,delta_(i-1)), n_i=d_i/d_(i+1).

The final gcd is one; n_i>1 for i>=2; n_i delta_i belongs to the semigroup
of earlier entries; and delta_i<n_(i-1) delta_(i-1) for i>=2.
Every primitive semigroup generator must occur in the sequence.
For nonprincipal initial orders, both are at most conductor+1=9,
and later entries are at most conductor-1=7.

The first two orders are the coordinate degrees of the given nonprincipal
plane presentation. This identification and triangular reduction are the
same classical input already checked in the frozen fixed cubic/quartic
packet. Constant field extension does not change a degree filtration:
filtered coefficient kernels commute with field extension. Therefore use
over k=C(x), after passage to its algebraic closure, is legitimate.
The normalization infinity point is residually rational and unique.

Here is a hand classification, independent of the enumerator:

- For <2,9>, the bounded positive entries are {2,4,6,8,9}.
  The primitive 9 must be one of the first two entries because later
  entries are at most seven. Pairing 9 with 2 gives the coprime sequence.
  Pairing it with 4 or 8 reaches gcd one without including primitive 2.
  Pairing it with 6 leaves gcd three, and the required next entry 2
  drops that gcd to one. Thus the patterns are (2,9) and (6,9,2),
  with either orientation of their first two entries.
- For <3,5>, the bounded entries are {3,5,6,8,9}. Both 3 and 5
  must occur. A coprime first pair therefore must be exactly 3,5.
  The only noncoprime nondividing alternatives are 6,8 and 6,9.
  The former has gcd two: adding either primitive 3 or 5 reaches
  gcd one before the other can occur. The latter has gcd three:
  3 cannot occur later since it makes no strict gcd drop, and
  adding 5 terminates the sequence without 3. No alternative remains.
- For <4,5,6>, the bounded entries are {4,5,6,8,9}. All of 4,5,6
  are primitive. A coprime first pair cannot include all three.
  The noncoprime nondividing possibilities are 4,6; 6,8; and 6,9.
  The first leaves gcd two and permits the required final 5.
  In 6,8, the missing primitive 4 makes no gcd drop; in 6,9,
  either missing primitive 4 or 5 would immediately reach gcd one
  before the other occurs. Only (4,6,5) remains.

The three-axiom tests verify both orientations of each listed pattern.
These conclusions classify characteristic sequences in this precise genus,
not polynomial curves of arbitrary conductor.

## 3. Explicit polynomial realizations, including the exceptional branches

The ordinary <2,9> and <3,5> examples are

    (t^2, t^9+t),    (t^3, t^5+t).

In each case the coprime degrees and power basis give exactly the displayed
semigroup. The derivatives have no common zero, since the first derivative
vanishes only at zero and the second derivative is one there.

For <4,5,6>, put

    p=t^4+t,  q=t^6+(3/2)t^3+t^2,
    r=q^2-p^3-2p^2+(3/4)q-p
      =-t^5+(1/8)t^3-(5/4)t^2-t.

The degrees of p and r are coprime, so k(p,q)=k(t). The monogenic
free k[p]-basis 1,q,q^2,q^3 transforms unimodularly to 1,q,r,qr.
Their degrees are 0,6,5,11, distinct modulo four. Thus their degree set
is precisely the union of those four translates of 4N, namely <4,5,6>.
Its delta is (0+6+5+11-0-1-2-3)/4=4. The derivative gcd is one.

The sextic/nonic alternative is also real. Put

    f=t^6+4t^2,  g=t^9+6t^5+6t.

The exact identity

    g^2=f^3+8f+4t^2

recovers t^2. Hence

    k[f,g]=k[t^2,t(t^8+6t^4+6)].

Its degree semigroup is <2,9>; the odd generator has degree nine, and
the even/odd direct sum over k[t^2] proves delta four. The derivative
gcd is one, or equivalently the quadratic normal form has p(0)=6!=0.
The quartic z^4+6z^2+6 is squarefree, so over C this example has four
ordinary nodes. In particular an attempted universal polynomial
nonrealizability obstruction for (9,6,2) is false.

See sextic_nonic/ for a separately derived exact certificate and family.

## 4. Universal degree-(4,6) coefficient classification

Assume k(p,q)=k(t) and dim_k(k[t]/k[p,q])=4, with p,q monic of
degrees four and six. Center the quartic and initially subtract constants
from p and q, over k. Write

    p0=t^4+a t^2+b t,
    q0=t^6+A t^5+c t^4+B t^3+e t^2+L t.

If A!=0, q0^2-p0^3 has degree eleven. The unimodular basis
1,q0,r,q0 r then has degrees 0,6,11,17, giving delta seven.
Therefore A=0.

Subtract (2c-3a)p0 q0 from q0^2-p0^3. Its degree-nine coefficient
is 2B-3b. If nonzero the same basis has degrees 0,6,9,15 and
delta six. Thus B=3b/2.

Next subtract sigma p0^2, where sigma=ac-c^2+2e. Its degree-seven
coefficient is 3ab/2-2bc+2L. If nonzero the basis degrees
0,6,7,13 give delta five. Thus

    L=b(c-3a/4).

Restoring scalar terms d,f gives the exact necessary normal form

    p=t^4+a t^2+b t+d,
    q=t^6+c t^4+(3b/2)t^3+e t^2+b(c-3a/4)t+f.       (N)

Every degree computation is made after subtracting d,f only for the
k-algebra comparison. Such x-dependent scalar subtraction is NOT claimed
to preserve a Jacobian. The original d,f remain in (N) and in Section 5.

For completeness this normal form has delta four exactly when

    b D != 0,   D=e-ac+3a^2/4.                       (D)

Indeed, after the preceding cancellations subtract kappa q0, with
kappa=-a^3+a^2c-ae-3b^2/4. The degree-five coefficient is -bD.
If it is nonzero, the basis degrees 0,6,5,11 prove delta four and
primitive fraction field (also gcd(4,5)=1).
If b=0 both coordinates are even and the fraction field is not k(t).
If b!=0 and D=0, remove the degree-four multiple of p0 too; the exact
remaining polynomial is

    (b^3/16)t(3a+2t^2).

It has degree three. The basis degrees 0,6,3,9 give delta three,
not four. This is the adjoining cubic/quartic boundary, not a
contradiction to it.

The mechanism is a *full module basis*, not an inference from the
existence of one isolated low-degree remainder.

## 5. No actual monic Keller source with semigroup <4,5,6>

Let R=C[x], K=C(x). Suppose one coordinate of P,Q is monic in y and
J(P,Q)=j in C*. The inherited top-coefficient gate makes the other
positive leading coefficient constant. Every ordinary degree reduction
then uses a constant coefficient and preserves the nonzero constant
Jacobian. The generic semigroup is unchanged.

A y-constant terminal coordinate would force a degree-one other coordinate
by its Jacobian, contradicting generic delta four. Section 2 therefore
forces terminal monic degrees (4,6) in the <4,5,6> sector.

Center the actual monic quartic by its cubic coefficient. This is a
polynomial source shear t=y+s(x), s in R, with determinant one.
All remaining quartic and sextic coefficients are in R. Apply the
necessary generic relations of Section 4 over K: they are literal
coefficient equalities, so the actual pair is exactly (N) with
a,b,c,d,e,f in R. No generic selected generator is globalized.

At t=0 the whole Jacobian is

    j=b((c-3a/4)d'-f').

Thus b is a polynomial unit, b in C*, and b'=0.
The coefficients of t^7,t^5,t^3,t^2 of the whole Jacobian form a linear
system forcing

    c'=3a'/2,
    d'=aa'/2,
    e'=ca',
    f'=ea'/2.                                      (E)

The simultaneous solution is legitimate because the t^2 equation may
be divided by b in C*. Explicit substitution of (E) into the ENTIRE
Jacobian gives the universal identity

    J(p,q)=-a'b(3a^2-4ac+6bt+4e)/8.                (J)

Its t coefficient is -3a'b^2/4. Since J is constant and b is a unit,
a'=0. Formula (J) now gives J=0, contradicting j!=0.

This excludes this full generic semigroup with a monic original coordinate.
It does not need constant fiber delta, a globally monic conductor,
saturation, finite embedding classification, or the moving-period theorem.
It does not exclude the other two genus-four semigroups.

## 6. The terminal-quadratic actual-Keller lemma

Suppose p is monic quadratic in y and J(p,q)=j in C*. A polynomial
source shear completes the square, so p=t^2+a(x).
At t=0,

    a'(x) q_t(x,0)=j.

Both factors are nonzero constants. Put a'=alpha in C*, so
a=alpha x+beta and (p,t) are polynomial coordinates. In these coordinates

    J(p,q)=alpha (partial q/partial t at fixed p).

Therefore q=(j/alpha)t+F(p), F in C[p]. This is a polynomial automorphism.
A constant-coefficient target shear removes F and leaves degree one,
contradicting a nonprincipal terminal (2,9) pair.

This theorem applies to an actual monic quadratic coordinate of the pair,
not simply to a quadratic element of K[p,q]. The explicit (6,9)
realization in Section 3 is the hostile control for confusing those two.

## 7. Verification and source boundary

Run from this packet's directory:

    python3 verify_frontier.py > results.json
    python3 -O verify_frontier.py > results-optimized.json
    cmp results.json results-optimized.json

The exact checker tests the exhaustive symmetric semigroup universe, all
bounded characteristic sequences using all three axioms, every symbolic
coefficient and remainder identity, all positive examples, the delta-three,
delta-five, delta-six, delta-seven hostile branches, and genuine
quadratic Keller positive controls. Explicit exceptions, not bare Python
assertions, enforce the checks.

The proofs in Sections 1,3--6 are universal algebraic arguments. The
bounded sequence enumeration is a corroboration of the exhaustive proof
in Section 2, whose finite scope follows from the cited conductor bound.
No polynomial finite search proves a general Keller exclusion.

Primary source freshly opened on 20 October 2026:
Sathaye--Stenerson, *On Plane Polynomial Curves*, author text dated
10 September 2014,
https://www.ms.uky.edu/~sohum/ma561/affgeom/papnew.pdf .
Section 3.1 supplies the characteristic axioms, bounds, and primitive-entry
rule; Section 2.2 supplies the one-place and field conventions. Its
introduction explicitly distinguishes plane and polynomial realization.
The frozen fixed cubic/quartic packet pins the first-entry identification
to Abhyankar--Sathaye (1996), Section 2, printed p.1062,
https://doi.org/10.1090/S0002-9939-96-03254-6 .
No new fresh retrieval of that second paper is claimed here.

