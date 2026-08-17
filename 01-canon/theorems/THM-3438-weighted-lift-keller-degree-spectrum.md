---
id: THM-3438
title: "Weighted-lift Keller degree spectrum and the quartic G1 witness"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  A one-variable weighted lift produces a three-dimensional Keller map of
  every generic field degree n>=3.  Together with the classical degree-two
  exclusion, this gives KDeg(m)={1} union {3,4,5,...} for every m>=3.
  Its n=4 member is an explicit z-quadratic degree-four G1 witness.  The
  geometric monodromy of every weighted member is S_n, so each is a
  composition atom; every numerically factorable grade also contains an
  explicit composite.  The theorem does not classify all maps within a
  degree or settle the z-affine order-{1,3} conjecture or JC(2).
source: codex2 independent derivation from the Gallagher weighted-lift construction, 2026-08-15
audit: two independent algebra/monodromy referees plus normal/-O/stored exact-companion replay; coefficient-free ordered-two-root proof of S_n; primary-source formula audit
depends_on:
  - THM-1330-keller-monoid-exact-picture-inverse-jelonek-cusp-rule
related:
  - THM-1300-jacobian-counterexample-dixmier-A3-explicit
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-3440-weighted-lift-cyclic-infinity-torsor-and-7x13-character-grid
  - THM-3441-weighted-quartic-jelonek-components-and-inertia-parity
  - THM-3541-keller-grade-monoid-factorization-lengths-and-elasticity
  - HYP-9030-keller-degree-semigroup
script: 04-computation/jc_weighted_lift_degree_spectrum_thm3438.py
output: 05-knowledge/results/jc_weighted_lift_degree_spectrum_thm3438.out
script_sha256: f31dc1acc5925586802dc83ba41db0b38a43afdd5dedbd84c0a0a2a7cd643982
output_sha256: 110a3ed0526bd84e47924ed3f38b2367bedf4c688121f974bf8e850ed1a1bc57
semantic_sha256: fb0f262a5a3d98d2f334eadb53380729b388656d8aa3035f32446ab2cebba721
hash_basis: LF-normalized bytes
---

# THM-3438 -- weighted-lift Keller degree spectrum and the quartic G1 witness

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The construction was published with exact source calculations in Alexis
Gallagher's
[weighted-lift derivation](https://github.com/algal/jacobianfun/blob/main/RESEARCH.md).
The proof below is an independent function-field derivation.  This records a
mathematical verification, not a historical-priority claim.

## 1. General weighted lift

Let `k` be a characteristic-zero field.  Choose `b,c in k*` and a polynomial
`p in k[w]` of degree `d>=2` such that

```text
p(0)=0,              p(1)=-c,
R(1)=0,              R(w)=integral_0^w p(s) ds.          (1)
```

Put

```text
kappa=p'(1)/c,       kappa!=-2,
a=-(1+kappa)/(2+kappa),                                  (2)
q(0)=0,              q'(w)=w p'(w)/c.                    (3)
```

For source coordinates `(x,y,z)`, define

```text
v=xy,              t=x^2 z,             u=1+v,
gamma=1+a v+b t,   w=u gamma,                            (4)

beta=c+p(w)/gamma,
alpha=u+q(w)/gamma^2,                                   (5)

F_p(x,y,z)=(alpha/x^2, beta/x, x gamma).                (6)
```

Then the apparent quotients in `(5)--(6)` cancel, `F_p` is polynomial, and

```text
det J(F_p)=bc.                                           (7)
```

Moreover

```text
[k(x,y,z):k(F_p)]=d+1.                                  (8)
```

Thus every such seed gives a Keller nonunit of generic fibre degree `d+1`.

## 2. Polynomiality and the determinant identity

Equation `(3)` and `p(0)=0` show that `q` has a double zero at zero.  Hence
`p(w)/gamma` and `q(w)/gamma^2` are polynomials in `(u,gamma)`.  At
`(v,t)=(0,0)` one has `(u,gamma,w)=(1,1,1)`.  Integration by parts gives

```text
q(1)=(p(1)-R(1))/c=-1.                                  (9)
```

Consequently

```text
beta(0,0)=0,        alpha(0,0)=0.                       (10)
```

Differentiating `(5)` with respect to `v` at the origin gives

```text
alpha_v(0,0)=1+kappa(1+a)+2a
            =1+kappa+a(kappa+2)=0                       (11)
```

by `(2)`.  Every monomial of `beta(v,t)` therefore contains `v` or `t`, so it
is divisible by `x` after `(v,t)=(xy,x^2z)`.  Every monomial of `alpha(v,t)`
contains `t` or at least `v^2`, so it is divisible by `x^2`.  This proves
polynomiality of `(6)`.

For target invariants put

```text
P=beta gamma=c gamma+p(w),
Q=alpha gamma^2=w gamma+q(w).                            (12)
```

In coordinates `(w,gamma)`, direct differentiation gives

```text
Jac_(w,gamma)(P,Q)=-c gamma.                             (13)
```

The changes `(v,t)->(u,gamma)->(w,gamma)` have determinants `b` and
`gamma`.  Therefore

```text
Jac_(v,t)(P,Q)=-bc gamma^2.                              (14)
```

For every weighted lift of the form `(6)`, direct row factoring gives

```text
Jac_(v,t)(P,Q)=-det(JF_p) gamma^2.                       (15)
```

Equations `(14)--(15)` prove `(7)`.

## 3. The inverse equation proves the generic degree

Integrating `(3)` gives

```text
cq(w)=w p(w)-R(w).                                      (16)
```

Equations `(12)` therefore imply the inverse equation

```text
R(w)=wP-cQ.                                             (17)
```

It is not merely a degree bound.  On the open set `C!=0`, where `(A,B,C)`
denotes the three target coordinates, the change

```text
P=BC,              Q=AC^2                               (18)
```

is birational.  For a root of `(17)` with `gamma!=0`, reconstruction is

```text
gamma=(P-p(w))/c,       u=w/gamma,
x=C/gamma,              v=u-1,
y=v/x,                  t=(gamma-1-a v)/b,
z=t/x^2.                                                   (19)
```

Thus one root gives exactly one input on the generic open set.

The polynomial

```text
T(w)=R(w)-Pw+cQ in k[P,Q,w]                              (20)
```

has degree `d+1` in `w` and is irreducible: it is primitive and linear in
the independent variable `Q`, while

```text
k[P,Q,w]/(T) ~= k[P,w]
```

by eliminating `Q`.  Hence `(T)` is prime, so `(20)` is irreducible over
`k(P,Q)`.  Formula `(19)` identifies the source function field with the
simple extension generated by `w`.  This proves `(8)`, including separability
in characteristic zero and the exact generic sheet count.

## 4. One explicit seed for every degree

For every `n>=3`, put `d=n-1` and

```text
p_d(w)=2w-3w^2
       +w(1-w)(w^(d-2)-6/[d(d+1)]).                     (21)
```

For `d=2`, the final term is identically zero.  For `d>=3`, its leading
coefficient is `-1`, so `deg p_d=d`.  In every case,

```text
p_d(0)=0,              p_d(1)=-1,
integral_0^1 p_d(w)dw=0,
p_d'(1)=-5+6/[d(d+1)]!=-2.                              (22)
```

The integral in `(22)` follows from

```text
integral_0^1 w(1-w)w^(d-2)dw=1/[d(d+1)],
integral_0^1 w(1-w)dw=1/6.                              (23)
```

Apply Sections 1--3 with `b=c=1`.  The resulting rational polynomial map

```text
F_n:A^3 -> A^3
```

has determinant one and generic field degree exactly `n`.

Let `KDeg(m)` be the set of generic field degrees of Keller endomorphisms of
`A^m` over characteristic zero.  Stable extension by identity coordinates
preserves degree, while THM-1330 records that field degree two is impossible.
Therefore, for every `m>=3`,

```text
KDeg(m)={1} union {n:n>=3}=N_{>=1}\{2}.                 (24)
```

This is an exact classification of **degree values**, not of Keller maps up
to source/target automorphisms.

## 5. The quartic 2-jet witness

The cubic seed `p(w)=w-2w^3` gives a particularly small integer model.  Put

```text
u=1+3xy,                gamma=1-4xy-x^2z,                (25)

G(x,y,z)=((2u+u^2-3u^4 gamma^2)/x^2,
          (1+u-2u^3 gamma^2)/x,
          x gamma).                                      (26)
```

Both quotients cancel.  The coordinate degrees are `(12,11,4)`, the degrees
in `z` are `(2,2,1)`, and

```text
det JG=-6.                                               (27)
```

The collision certificate is

```text
G(1,0,0)=G(-1,0,2)=(0,0,1).                             (28)
```

For target `(A,B,C)`, put

```text
P=BC,             Q=AC^2,             w=u gamma.         (29)
```

Then

```text
w^4-w^2+2Pw-Q=0.                                        (30)
```

Conversely, on the generic open set a root reconstructs uniquely by

```text
gamma=P-w+2w^3,       u=w/gamma,       x=C/gamma,
y=(u-1)/(3x),
z=(1-(4/3)(u-1)-gamma)/x^2.                             (31)
```

The quartic in `(30)` is irreducible over `k(P,Q)` by the same linear-in-`Q`
argument as `(20)`.  At `(A,B,C)=(1,0,1)` it is

```text
w^4-w^2-1,                                               (32)
```

whose discriminant is `-400`; all four roots satisfy `gamma!=0` and reconstruct
four distinct inputs.  Thus `G` has field degree four independently by both
the function-field and explicit-fibre routes.

Because `(26)` is quadratic in `z`, it is a positive degree-four 2-jet G1
witness.  It supersedes the global-open verdict in THM-2465 while leaving all
of that theorem's stated local stratum exclusions intact.  It does **not**
refute the z-affine order-`{1,3}` conjecture, since `G` is not z-affine.

## 6. Full symmetric monodromy and the atom/decomposition grades

The inverse equation supplies more than its degree.  Its ramification locus
in the incidence cover over the `(P,Q)`-plane is

```text
T(w)=0,             partial_w T=p(w)-P=0.
```

The reduced branch curve is therefore parametrized by

```text
(P,Q)=(p(w),q(w)).                                    (33)
```

It is irreducible.  At a generic point `p'(w)!=0`, and

```text
dQ/dP=q'(w)/p'(w)=w/c.                                (34)
```

Thus the tangent slope recovers `w`: the parametrization `(33)` is
generically one-to-one, exactly one pair of sheets meets at a generic branch
point, and generic inertia is a transposition.

There is also an exact two-root argument.  For distinct ordered roots `r,s`
of `(20)`, subtraction and either root equation give

```text
P=(R(r)-R(s))/(r-s),
cQ=Pr-R(r).                                             (35)
```

The divided difference is polynomial in `(r,s)`.  Hence the ordered-two-root
incidence is isomorphic to `A^2_(r,s)\{r=s}` and is irreducible.  Therefore
geometric monodromy is 2-transitive.  A 2-transitive subgroup containing one
transposition contains every transposition, so

```text
Mon_geom(F_n)=S_n                  for every n>=3.       (36)
```

The point stabilizer `S_(n-1)` is maximal in `S_n`.  Therefore the source
function field of every weighted member has no proper intermediate field,
and every `F_n` is composition-irreducible.

For a fixed ambient dimension `m>=3`, define `AtomDeg(m)` to be the grades
containing at least one composition atom and `DecompDeg(m)` the grades
containing at least one composition of two nonunits.  Equations `(24)` and
`(36)`, degree multiplicativity, and composition of weighted members give

```text
AtomDeg(m)={n:n>=3},
DecompDeg(m)={ab:a,b>=3}.                               (37)
```

The second line is both necessary and sufficient.  Thus every numerically
nonfactorable grade contains only atoms, while every grade in `DecompDeg(m)`
is mixed: it contains the primitive `S_n` member `F_n` and a composite
`F_a o F_b`.  Degree nine is the first mixed grade.  The grades in which
**every** map is forced atomic are exactly the numerical nonproducts

```text
{odd primes} union {4,8} union {2p:p an odd prime}.      (38)
```

In particular, all maps of degree `3` through `8` are atoms, `(26)` is an
`S_4` atom, and the degree-five weighted member has monodromy `S_5`; its
generic inverse is not solvable by radicals.  This removes the old proposal
that an as-yet-unbuilt `A_5` quintic would be the first non-radical Keller
inverse.  It does not construct an `A_5` member: the explicit one is `S_5`.

Thus HYP-9030's proved monoid law and the subfamily `{3^k}` survive, but its
conjecture `KDeg(3)={3^k}` and atom prediction are refuted.  The tournament
strong-component analogy remains only a grammar for factorization; actual
blocks are controlled by intermediate fields and monodromy.

THM-3541 subsequently classifies the numerical monoid behind this paragraph.
Its atoms are exactly the odd primes, `4`, `8`, and twice an odd prime; all
factorization-length sets are explicit intervals and the numerical elasticity
is `3/2`.  This is a theorem about the grading.  It does not alter the fact
that the weighted `S_n` member is a map atom in every reducible grade.

## 7. Connection and loss ledger

| field | exact content |
|---|---|
| source | one-variable seed `p`, endpoint/integral conditions, and two nonzero lift scales |
| target | a polynomial Keller map in three variables and its inverse equation `(17)` |
| map | weighted invariants `(v,t)`, then `(P,Q,C,w)` |
| preserved | determinant, generic field degree, root multiplicity, weighted scaling, and stable degree under padding |
| destroyed by degree | formulas, which map in a mixed grade was chosen, Jelonek geometry, collisions, and composition factorization |
| restoring sidecars | inverse polynomial, reconstruction map, monodromy/block lattice, and boundary data |
| cheapest hostile to `{3^k}` | the quartic map `(26)` |
| G1 boundary | positive for z-quadratic 2-jets; z-affine order-`{1,3}` remains open |

The result gives no `JC(2)` or LRC consequence and no classification of all
counterexamples.  It closes the generic-degree **value** spectrum and the
existence spectra of atom/composite grades in every fixed dimension at least
three, but not equivalence or decomposition of an arbitrary map.

## 8. Exact companion

The standard-library companion independently checks polynomial cancellation,
both determinant routes, `(28)--(32)`, the seed identities, reconstruction,
generic degrees `3` through `100`, expanded lifts through degree `18`, and
hostile boundary controls.  Normal and optimized replays are byte-identical to
the stored output, with the hashes pinned in the frontmatter.  The monodromy
claim is proved structurally by `(33)--(36)` and was independently rederived
twice; the finite computation is evidence for the construction, not a
substitute for that proof.

Run

```bash
python -B 04-computation/jc_weighted_lift_degree_spectrum_thm3438.py
python -B -O 04-computation/jc_weighted_lift_degree_spectrum_thm3438.py
```
