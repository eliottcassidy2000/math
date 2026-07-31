---
id: THM-2871
title: "Quartic leading-face survivor and integral-depression flag dichotomy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  At a simple
  quartic leading face with B*delta3 a unit, one chosen coordinate root
  has valuation -1 and the other three form the separable leading cubic.
  THM-2867's matching blow-up is exactly an affine copy of that x-finite
  cubic, but only in the associated graded.  At the next leading-cubic
  face, a unit quadratic coefficient gives one escaping sheet and no
  integral depression, whereas C=B*c gives an integral trace-zero
  coordinate and a paired slope-1/2 escape with
  delta3=-B(4*P3^3+27*B*Q3^2).  The sporadic Keller square law is one
  further non-universal sidecar.  A connected etale quasi-finite
  generic-S4 family realizes an arbitrary depressed cubic on its
  nonproper leading face, so none of the Keller/Jelonek ownership follows
  from the local graph data.
source: root/quartic-leading-face-flag-2026-07-28
audit: >
  thm2871-leading-face-flag-audit-2026-07-28 (independent Newton/root
  normalization, projective-versus-affine finiteness, cubic/orientation
  special-fibre, integral-depression, square-law, generic-S4,
  exact-nonproper-locus, graph-owner, and exact-evidence audit: ACCEPT)
depends_on:
  - THM-2867-homogeneous-quartic-orientation-sextic-and-cubic-leading-residual
related:
  - THM-2867-homogeneous-quartic-orientation-sextic-and-cubic-leading-residual
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2758-quartic-pair-sum-sextic-resolvent-pullback-and-discriminant-square
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - HYP-9027-twojet-disc-jelonek-odd-exponent-law
script: 04-computation/quartic_leading_face_integral_depression_thm2871.py
output: 05-knowledge/results/quartic_leading_face_integral_depression_thm2871.out
script_sha256: 40e8c43b044f328df66e23cd27e205155892ca8e5fe82fb6c4a3885d95321edf
output_sha256: 5d0925e6146faaba623bd0229060fd16999cb4b2c160d19b2f21ebfbde103dcc
hash_basis: LF-normalized bytes
---

# THM-2871 -- the cubic survivor still needs an integral depression flag

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2867 proves the striking identity

```text
U_0(t)=-64B^2 h(-(t+4C)/(4B)),               (1)
```

where `U_0` is the first matching blow-up of

```text
f=A X^4+B X^3+C X^2+D X+E                    (2)
```

at `A=0`, and

```text
h=B X^3+C X^2+D X+E.                          (3)
```

Equation `(1)` really does recover the leading cubic.  The point of the
present theorem is to type exactly what that recovery means.

There are two successive operations:

```text
quartic graph coordinate
   --simple A-leading face--> x-finite cubic survivor
   --B-leading flag--------> one escape or a paired escape. (4)
```

The first arrow is an associated-graded or normal-cone identification.  It
does not turn a generic quartic resolvent field into a subfield of a quartic
root field and does not construct a degree-three Keller map.  The second
arrow has a sharp divisibility gate: the cubic trace-zero translation
extends integrally exactly when its quadratic coefficient is divisible by
its leading coefficient.  Even after that gate, the special square law of
the fixed sporadic Keller map remains an additional identity.

The theorem is local algebra and invariant theory.  Its final family is a
connected etale generic-`S4` hostile showing why no stronger Keller
conclusion is formal.

## 1. First stage: one quartic coordinate pole and three finite roots

Let `R` be a henselian discrete valuation ring with fraction field `K`,
uniformizer `pi`, and residue characteristic different from `2` and `3`.
After strict henselization, assume the residue field is separably closed.
Let

```text
f(X)=A X^4+B X^3+C X^2+D X+E in R[X],         (5)
```

and put

```text
delta_3=Disc_X(BX^3+CX^2+DX+E)
       =C^2D^2-4BD^3-4C^3E-27B^2E^2+18BCDE.  (6)
```

Assume

```text
v(A)=1,                     v(B)=v(delta_3)=0. (7)
```

The Newton point `(4,1)` is joined to `(3,0)` by a terminal segment
of slope `1`.  Every coefficient point of degree at most three has
nonnegative height.  Hence precisely one root has valuation `-1`, while
the other three roots are integral.  More precisely, if `x_infty` is the
pole root, then

```text
v(x_infty)=-1,              A x_infty/B = -1 mod pi,       (8)
```

and the reductions of the three integral roots are the three distinct roots
of

```text
h_bar(X)=bar B X^3+bar C X^2+bar D X+bar E.   (9)
```

This is also a Hensel factorization statement.  Over the strict
henselization,

```text
one rank-one pole branch + three rank-one x-finite branches. (10)
```

No ramification is hidden in the pole.  The homogeneous quartic
discriminant satisfies

```text
Disc(f)|_(A=0)=B^2 delta_3,                  (11)
```

which is a unit under `(7)`.  Thus the normalization of the projective root
cover -- equivalently, the separated root branches after adjoining the
point at infinity for the chosen coordinate -- is etale at this leading
face.  This does **not** say that `Spec R[X]/(f)` is finite over `Spec R`:
`f` is nonmonic because `A` is not a unit.  The missing affine-coordinate
sheet is a pole of the chosen graph coordinate, not a quartic branch.

The same proof works without a separably closed residue field by replacing
the three individual finite roots with their rank-three finite-etale
algebra.  Strict henselization is used only to name the branches.

## 2. The matching blow-up is the cubic survivor, in the special fibre

Use the integral quartic covariants

```text
P=8AC-3B^2,
Q=B^3-4ABC+8A^2D,
R=-3B^4+16AB^2C-64A^2BD+256A^3E,

K=(P^2-R)/4,
T(v)=v^3+Pv^2+Kv-Q^2.                         (12)
```

Direct expansion gives

```text
T(B^2+At)=A^3 U(t),                           (13)

U(t)=t^3+8Ct^2+16(C^2+BD-4AE)t
    +64(BCD-B^2E-AD^2).                       (14)
```

At `A=0`,

```text
U_0(t)=t^3+8Ct^2+16(C^2+BD)t+64(BCD-B^2E),   (15)
```

and the denominator-free form of `(1)` is

```text
U_0(-4(C+Bx))=-64B^2h(x).                     (16)
```

Consequently, if `alpha_1,alpha_2,alpha_3` are the finite root residues,
then

```text
t_i=-4(C+B alpha_i)                           (17)
```

are exactly the roots of `U_0`.  The first matching blow-up is therefore
the characteristic polynomial of an affine coordinate on the rank-three
`x`-finite survivor algebra.

There is a parallel ordered-difference statement.  Put

```text
I_0=C^2-3BD,
tau_3=2C^3-9BCD+27B^2E.                       (18)
```

The natural orientation boundary is

```text
O_0(W)=W^6-2B^4I_0W^4+B^8I_0^2W^2-B^14delta_3
       =product_(i<j)(W^2-B^6(alpha_i-alpha_j)^2).        (19)
```

Thus the orientation blow-up is the ordered-pair or directed-difference
Galois-closure carrier of the same cubic.  It is squarefree exactly when

```text
B delta_3 tau_3 !=0.                         (20)
```

Equations `(16)` and `(19)` identify the **special fibres** of two Rees
or normal-cone constructions.  They do not give a generic field inclusion.
THM-2598's generic quartic root field and matching cubic field remain
incomparable in the live `A4/S4` cases.  Specialization can identify their
associated graded algebras without producing a rational map between their
generic covers.

## 3. Second stage: the integral-depression dichotomy

Now let `S` be a henselian DVR with uniformizer `b`, residue characteristic
different from `2` and `3`, and consider a cubic

```text
g(X)=bX^3+c_2X^2+dX+e.                         (21)
```

There are two sharply different leading faces.

### 3.1 A unit quadratic coefficient gives one escaping sheet

Assume

```text
v(c_2)=0,                 v(d^2-4c_2e)=0.      (22)
```

The terminal Newton segment joins `(2,0)` to `(3,1)`.  Hence one root has
valuation `-1`, while two roots are integral and reduce to the two distinct
roots of

```text
c_2X^2+dX+e.                                  (23)
```

The cubic discriminant is a unit:

```text
Disc(g)|_(b=0)=c_2^2(d^2-4c_2e).              (24)
```

The usual depressed coordinate would translate by

```text
c_2/(3b),                                      (25)
```

which has a pole.  Thus the trace-zero coordinate does not extend over the
leading face, and there is no paired-escape law.

### 3.2 Divisibility by the leading coefficient gives a paired escape

The trace-zero translation extends over `b=0` exactly when

```text
c_2=bc                                         (26)
```

for some `c in S`.  With the integral change

```text
X=Y-c/3,                                      (27)
```

one obtains

```text
g(X)=bY^3+P_3Y+Q_3,                           (28)

P_3=d-bc^2/3,
Q_3=e-cd/3+2bc^3/27.                          (29)
```

The exact discriminant factorization is

```text
delta_3
 =-4bP_3^3-27b^2Q_3^2
 =-b(4P_3^3+27bQ_3^2).                        (30)
```

If `P_3` is a unit, the Newton segment from `(1,0)` to `(3,1)` has
slope `1/2`.  Therefore

```text
one integral root + two roots of valuation -1/2.          (31)
```

The two pole roots form the tame quadratic pair, and

```text
v_b(delta_3)=1.                               (32)
```

This is the universal integral-depression and even-escape mechanism behind
a depressed cubic.  It is exactly where the primes `2` and `3` enter:
three roots admit a trace-zero translation, and the missing quadratic term
forces the two nonfinite roots into one conjugate pair.

The two cases exhaust a transverse DVR flag.  Either `c_2` is a unit, or
`c_2` lies in the principal ideal `(b)` and `(26)` holds.  Higher
valuations or failure of the displayed residual unit conditions belong to
intersections and require another Newton polygon.

## 4. The sporadic square law is an additional sidecar

Equation `(30)` gives an odd leading factor.  It does **not** make the
remaining cofactor a square.  The stronger square law is the separate
condition

```text
P_3^3+(27/4)bQ_3^2=S^2,                      (33)
```

equivalently

```text
delta_3=-4bS^2.                               (34)
```

The fixed sporadic Keller map in THM-2473 satisfies this condition for a
special reason.  Its cubic is

```text
L X^3+(4-3 beta gamma)X-2gamma,               (35)
```

where

```text
L=27alpha^2gamma^2-18alpha beta gamma+16alpha
  +beta^3gamma-beta^2,

S=27alpha gamma^2-9beta gamma+8.              (36)
```

Direct expansion gives the trisection identity

```text
(4-3 beta gamma)^3+27L gamma^2=S^2.           (37)
```

Taking

```text
b=L,             P_3=4-3 beta gamma,
Q_3=-2gamma                                      (38)
```

turns `(33)` into `(37)` and `(30)` into

```text
Disc=-4LS^2.                                   (39)
```

Thus the fixed map has three logically separate properties:

```text
integral depression;
paired leading escape with odd discriminant;
trisection cofactor square.                    (40)
```

Only the first two follow from `(26)` and the residual unit gate.  The third
is not universal cubic invariant theory.

## 5. Sharp local hostiles

The one-escape boundary is already visible in

```text
g_1=bX^3+X^2+1.                                (41)
```

Its quadratic residual is separable and

```text
Disc(g_1)=-4-27b^2.                            (42)
```

There are two finite roots and one valuation-`-1` root.  The discriminant is
a unit and the depression is not integral.

The paired-escape boundary is

```text
g_2=bX^3+X+1,                                  (43)
Disc(g_2)=-b(4+27b).                           (44)
```

It has an integral depressed coordinate, one finite root, two
valuation-`-1/2` roots, and odd discriminant.  Nevertheless `4+27b` is not
a square in `Q(b)`.  Therefore paired escape still does not imply `(33)`.

These examples isolate the first failed implication in both directions:
trace-zero depression is neither automatic at a cubic leading face nor
sufficient for the sporadic square factor.

## 6. A connected etale generic-`S4` hostile

There is a universal graph family satisfying all of the local leading-face
geometry while leaving the cubic arbitrary.  Put

```text
Phi=aX^4+X^3+pX+q,                            (45)

Delta=256a^3q^3-27a^2p^4-192a^2pq^2
      -6ap^2q-4p^3-27q^2.                    (46)
```

Let

```text
Y=Spec Q[a,p,q,Delta^(-1)],
X=Spec Q[a,p,q,Delta^(-1),X]/(Phi).           (47)
```

Then

```text
X -> Y is connected, etale, and quasi-finite
of generic degree four.                       (48)
```

Indeed:

1. `Phi` is irreducible.  Viewed as a polynomial of degree one in `a`, a
   nontrivial factor independent of `a` would divide both `X^4` and
   `X^3+pX+q`, whose gcd is one.  Gauss's lemma finishes the argument.
2. If `a!=0`, a common affine root of `Phi` and `dPhi/dX` makes `Delta`
   zero.  If `a=0`, such a root makes
   `-4p^3-27q^2=Delta|_(a=0)` zero.  Hence the relative derivative is a
   unit on `(47)`, proving etaleness.
3. Every fibre is the zero set of a nonzero polynomial of degree four or
   three, so the map is quasi-finite.

The generic Galois group is exactly `S4`.  At

```text
(a,p,q)=(2,-1,1)                             (49)
```

the polynomial

```text
2X^4+X^3-X+1                                 (50)
```

has discriminant `2673`, and exact finite-field factor degrees

```text
mod 5:   4,
mod 13:  1+1+2,
mod 17:  1+3.                                 (51)
```

Thus its Galois group is transitive and contains a four-cycle, a
transposition, and a three-cycle.  The unique such subgroup of `S4` is
`S4`.  A good specialization group embeds in the generic group, so the
generic group is `S4` as well.

The nonproper set of `(48)` is exactly

```text
{a=0}.                                       (52)
```

Over `a!=0`, division by `a` makes `(45)` monic, hence finite.  Over the
generic point of `a=0`, the fibre has the three roots of

```text
X^3+pX+q,                                    (53)
```

whereas the generic degree is four.  A proper quasi-finite morphism would be
finite and an etale finite rank could not drop from four to three.  Hence
the entire divisor `(52)` is nonproper.

Equivalently, at every point of `(52)` the reciprocal coordinate has the
simple omitted pole branch from Section 1; `Delta^(-1)` keeps the remaining
cubic separable.  This gives the same failure of the valuative criterion
pointwise along the divisor.

On that divisor, the survivor is the arbitrary separable depressed cubic
`(53)`.  The two blow-ups are

```text
U_0(u)=u^3+16pu-64q,                          (54)

O_0(w)=w^6+6pw^4+9p^2w^2+4p^3+27q^2.         (55)
```

At `(p,q)=(-1,1)` these become

```text
U_0=u^3-16u-64,
O_0=w^6-6w^4+9w^2+23.                        (56)
```

The family `(47)` is a sharp hostile, not a Keller counterexample.  Its
target is an open subset, its source is the spectral hypersurface rather
than affine space, and no constant-Jacobian polynomial coordinate system is
asserted.  It proves that the following package is still insufficient:

```text
connected + etale + quasi-finite + generic S4
+ an actual nonproper leading divisor
+ three surviving sheets
+ the exact cubic and ordered-difference blow-ups.        (57)
```

The missing content is global affine-source/Keller ownership, not another
quartic invariant identity.

## 7. What an actual Keller graph chart would transfer

Let

```text
F:A^n -> A^n                                  (58)
```

be a hypothetical degree-four Keller map, let `x` be a polynomial source
coordinate which is primitive over the target function field, and let
`(5)` be a cleared graph relation for `x`.  Write `Xbar` for the finite
normalization of the target in the source function field.

At a divisor satisfying `(7)`, Sections 1--2 identify four normalization
branches:

```text
one x-pole branch;                 three x-finite branches. (59)
```

Every branch belonging to the original affine source is `x`-finite, because
`x` is polynomial there.  The converse is false: another source coordinate
may have a pole on a branch where `x` remains finite.  Therefore

```text
affine-source branches subseteq x-finite branches,        (60)
```

and `h`, `U_0`, and `O_0` initially describe only the right side.

Equality in `(60)` is an additional reconstruction condition.  In a
primitive-element reconstruction of the remaining source coordinates as
rational functions of `x`, every denominator must be a unit at all three
finite branches.  Equivalently, the `x`-finite open of `Xbar` must agree
with the original affine source over the generic point of the divisor.

If that condition holds, the actual transfer is exactly:

```text
a rank-three finite-etale cover of a dense open in {A=0},
with characteristic polynomial h and ordered-difference closure O_0.    (61)
```

It is still not a degree-three Keller map between affine spaces.  Neither
the source divisor nor the target divisor in `(61)` is proved to be affine
space, and the restriction need not have a constant Jacobian in polynomial
coordinates.

There is also a divisor-type mismatch.  Under `(7)`, `Disc(f)` is a unit,
so `A=0` is an unramified coordinate-pole component.  THM-2473's
load-bearing divisor is instead the odd cubic-discriminant/Jelonek component
where the cubic leading coefficient vanishes and two roots escape.  The
first place that can reproduce that mechanism is the deeper flag

```text
A=B=0,                                        (62)
```

not the generic point of `A=0`.  On `(62)`, the exact tests are:

```text
(i)  C lies in (B): integral depression and paired escape;
(ii) 4P_3^3+27BQ_3^2 is four times a square;
(iii) B=0 is an actual Jelonek owner for the rank-three cover;
(iv) all source-coordinate reconstruction denominators are units.       (63)
```

Failure of `(i)` gives the one-escape boundary.  Failure of `(ii)` leaves
the universal odd cubic discriminant but loses the sporadic trisection law.
Failure of `(iii)` or `(iv)` means the polynomial carrier still has no
physical Keller owner.

This is the cheapest exact graph-chart decision procedure exposed by the
theorem.  It replaces a global comparison of two resolvent discriminants by
ideal membership, one square test, and a reconstruction-denominator
valuation audit.

## 8. Exact companion

Run

```bash
python 04-computation/quartic_leading_face_integral_depression_thm2871.py
python -O 04-computation/quartic_leading_face_integral_depression_thm2871.py
```

The dependency-free companion uses only the Python standard library and
contains no truth-bearing `assert`.  It implements:

- an exact sparse multivariate polynomial ring;
- a Sylvester determinant deriving `(46)` independently;
- the affine cubic identity `(16)`;
- the integral depression and discriminant factorization `(28)--(30)`;
- the sporadic square identity `(37)`;
- both hostile discriminants `(42)` and `(44)`;
- exact Newton polygons for the three displayed slope patterns;
- finite-field gcd/Frobenius factor degrees in `(51)`; and
- an independent enumeration of all `30` subgroups of `S4`, finding `S4`
  as the unique transitive subgroup with the three cycle types in `(51)`.

Normal and optimized runs byte-match the stored transcript.

## 9. Scope ledger

```text
PROVED:
  exact two-stage henselian Newton polygons;
  one quartic x-pole and the rank-three finite survivor;
  unramified nature of the simple A-leading face;
  exact U_0 affine identification with the leading cubic;
  ordered-difference orientation boundary;
  C-unit one-escape versus C=Bc paired-escape dichotomy;
  exact integral depression and delta_3 factorization;
  separation of odd discriminant from the sporadic square law;
  sharp local hostiles;
  connected etale quasi-finite generic-S4 hostile family;
  exact graph-owner/reconstruction boundary.

NOT PROVED:
  that a live Keller graph quartic enters the unit chart (7);
  that x-finite branches are affine-source branches;
  an affine-space or Keller realization of the cubic survivor;
  a Jelonek owner for B=0;
  the square law (33) for an arbitrary graph cubic;
  exclusion of A4 or S4 Keller monodromy;
  JC(2), DC(2), general JC, or any LRC statement.           (64)
```

The independent hostile audit rederived all three Newton polygons, the
projective/affine normalization split, both special-fibre identities, and
the square-law boundary.  It separately checked connectedness, etaleness,
the exact nonproper locus, and the specialization-to-generic `S4`
argument, as well as every Keller-owner caveat.  Normal, optimized, and
stored transcripts agree exactly; both declared LF-normalized hashes
match; the companion compiles without optional CAS, truth-bearing
`assert`, or floating-point decisions; and the documentation gate passes.

**QED.**
