---
id: THM-2748
title: "Single response-branch ramification and monomial-composition boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  A rational
  primitive on the source P1 that factors through a
  projective integral response component with one response pole forces the
  component normalization to be P1.  After placing the two pole points at
  infinity, both the response and the normalization map are pure powers.
  If the target pole has order P, the constant-U case forces P=degree(nu)=1;
  the monomial-U case forces m-1=P*degree(nu).  Applied to any one-branch
  escape from THM-2741, constant U is impossible and the response component
  is subject to this full monomial-composition normal form.  This does not
  prove that a reducible component contains a forced branch or close JC(2).
source: thm2694-full-lift-fibre-scout-2026-07-28
audit: jc-one-pole-audit-2026-07-28 (independent Luroth/function-field, pole-divisor, pure-power factorization, reduction/normalization-lift, fibre-pullback, THM-2741 table, and normal/optimized/hash audit: ACCEPT)
depends_on:
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2741-highest-odd-faber-response-pole-capacity-closure
related:
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
  - THM-2745-highest-odd-faber-componentwise-exact-prefix-closure
script: 04-computation/jc2_one_pole_response_composition_thm2748.py
output: 05-knowledge/results/jc2_one_pole_response_composition_thm2748.out
script_sha256: 0043fbbdd8101d14f9df7943b392606a045edcba6eb7562de34501343bf79a2a
output_sha256: b7220369b72b8b4ddeca97cedcb0e6699e5175fd98e123a64f1e0f95f465ba02
hash_basis: LF-normalized bytes
---

# THM-2748 -- one response pole forces a monomial composition

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2741 records the degree identity `m-1=eP` on a hypothetical reducible
component carrying only one of its response poles.  The same argument has a
rigid equality case: after the natural source and target coordinate changes,
the response and the normalization map must both be monomials.  This turns
the one-branch escape from a numerical congruence into an exact composition
normal form.

## 1. Abstract one-pole theorem

Work over `C`.  Let `Cbar` be a smooth projective integral curve, let

```text
nu:P1_x -> Cbar                                             (1)
```

be a nonconstant morphism, and let `R in C(Cbar)` be nonconstant.  Suppose
`R` has exactly one pole `P_infty`, of order `P>=1`.  On the affine source
line assume

```text
U(x) (R o nu)'(x)=kappa,
U in C[x]\{0},                       kappa in C*.          (2)
```

Assume the rational-primitive conclusion of THM-2723, equivalently one of

```text
U in C*,
R o nu=s_0+s_1 x;                                         (3a)

U=u_0(x-a)^m,
R o nu=s_0+s_1(x-a)^(1-m),          m>=2, s_1!=0.         (3b)
```

Then `Cbar` is isomorphic to `P1`.  More precisely:

1. In `(3a)`, necessarily

   ```text
   P=1,                  deg(nu)=1.                       (4)
   ```

2. In `(3b)`, put `X=(x-a)^(-1)` and `e=deg(nu)`.  There is a coordinate
   `Z` on `Cbar`, with its only pole at `P_infty`, and constants
   `z_0,rho,gamma in C`, `rho*gamma!=0`, such that

   ```text
   nu^*Z=z_0+gamma X^e,
   R=s_0+rho (Z-z_0)^P,
   rho gamma^P=s_1,
   m-1=eP.                                                (5)
   ```

In particular, the fibre over `P_infty` is the single source point `x=a`,
totally ramified with index `e`; `P` divides `m-1`; and the response map on
the target normalization has a unique finite critical point `z_0` after the
coordinate choice `(5)` when `P>=2` (and no finite critical point when
`P=1`).

The conclusion is sharp for every pair `P,e>=1`: choose `(5)`, put
`m=1+eP`, and choose `s_1=kappa/[u_0(1-m)]`.

## 2. Proof

The morphism `(1)` is nonconstant, so the induced inclusion of function
fields embeds `C(Cbar)` into `C(x)`.  Luroth's theorem makes this subfield
rational.  Since `Cbar` is smooth and projective,

```text
Cbar isomorphic to P1.                                    (6)
```

Choose a coordinate `Z` on `Cbar` whose unique pole is `P_infty`.  Then

```text
R=f(Z),                         f in C[Z], deg f=P.       (7)
```

The poles of `R o nu` are exactly the points of `nu^*(P_infty)`.  Either
formula in `(3)` has one pole on the projective source.  Therefore the fibre
of `P_infty` is a single point, with ramification index equal to
`e=deg(nu)`.  Choose a source coordinate `X` having its only pole at that
point.  Then

```text
nu^*Z=g(X),                         g in C[X], deg g=e.   (8)
```

In case `(3a)`, take `X=x`.  In case `(3b)`, take
`X=(x-a)^(-1)`.  Both cases become

```text
f(g(X))=s_0+s_1 X^N,                                    (9)
```

where `N=1` in `(3a)` and `N=m-1` in `(3b)`.

It remains to classify equality in `(9)`.  Every root `c` of `f-s_0` has a
preimage under the nonconstant complex polynomial `g`.  Equation `(9)` says
that every such preimage is `X=0`.  Hence `f-s_0` has the unique root
`z_0=g(0)`, and then `g-z_0` also has the unique root zero.  Degree and
factorization over `C` give

```text
f(Z)=s_0+rho (Z-z_0)^P,
g(X)=z_0+gamma X^e.                                    (10)
```

Substitution in `(9)` yields

```text
N=eP,                         rho gamma^P=s_1.          (11)
```

For `(3a)`, `N=1`, proving `(4)`.  For `(3b)`, `N=m-1`, proving `(5)`.

## 3. Highest-odd response application

Let `X_a` be a THM-2741 highest-odd response member, possibly reducible or
nonreduced, and let `D` be the closure in `(X_a)_red` of the generic image of
a physical source trajectory.  Because `P1_x` is reduced, the trajectory
kills target nilpotents and factors through `(X_a)_red`.  Moreover

```text
R_source'=kappa/U != 0,                                  (12)
```

so its generic image is not a point.  Since the response member is a
projective curve, `D` is a reduced irreducible curve component.  Let `Cbar`
be its smooth projective normalization.  The induced inclusion
`C(D) -> C(x)` lifts the generic trajectory to a rational map
`P1_x ---> Cbar`; properness extends it across every source point, and
nonconstancy makes the extension a finite surjective morphism

```text
nu:P1_x -> Cbar,             R_source=nu^*R,              (13)
```

where `R` is the affine response restricted to `Cbar`.

Suppose `Cbar` has a normalization point `P_infty` representing one of
THM-2741's forced formal infinity branches, with exact valuation

```text
ord_(P_infty)(R)=-P_j.                                   (14)
```

Write

```text
j=max{odd k:a_k!=0},       r=22-j,       g=gcd(r,6),
P_j=(150-6r)/g.                                          (15)
```

THM-2741 gives the exact table

```text
j       21  19  17  15  13  11   9   7   5   3   1
P_j    144  44 120 108  32  84  72  20  48  36   8.    (16)
```

Both rational-primitive alternatives in `(3)` have exactly one pole on
`P1_x`.  For every target pole `Q` of `R` and every `p in nu^(-1)(Q)`,

```text
ord_p(nu^*R)=e_p ord_Q(R)<0.                              (17)
```

Surjectivity makes each pole fibre nonempty, and fibres over distinct target
points are disjoint.  Hence a physical component can have at most one
response pole.  Since `P_infty` is already one, it is the unique pole and the
abstract theorem applies.  Every entry in `(16)` is at least eight, so the
constant-`U` conclusion `P_j=deg(nu)=1` is impossible.  Thus every surviving
component carrying such a branch satisfies `(5)` with

```text
m=1+eP_j,                       e=deg(nu).                (18)
```

This strictly strengthens the divisibility sentence in THM-2741: not only
the pole order, but the entire target response and normalization map are
pure powers after affine/Mobius gauge.

THM-2745 subsequently closed the complete degree-22 highest-odd physical
boundary by retaining the polynomial exact-prefix identities on every
component.  Thus the application above is no longer needed to close that
family.  The present theorem survives as a sharp abstract classification
after those lift equations are forgotten.

## 4. Boundaries

The one-pole hypothesis is load-bearing.  This theorem alone does not prove
that a reducible or nonreduced degree-22 member has a component containing
one of the forced THM-2741 branches.  THM-2745 supplies that separate global
`h=0` divisor analysis and closes the physical boundary.  A component
containing more than one pole is excluded here, not classified by `(5)`.

For a nonreduced response member, the theorem is applied to the normalization
of the reduced generic-image component: nilpotent and embedded structure is
killed by the reduced source and cannot change `R_source`.  What remains
unproved is that this component contains one of the forced formal branches;
it may instead meet the `h=0` divisor elsewhere.

The coordinate changes are also load-bearing.  Before placing the unique
source and target poles at infinity, the two maps need not look like literal
monomials.  In positive characteristic, inseparable powers lie in the
derivative kernel and both the rational-primitive and root-factor arguments
need separate typing.

No factorization locus is classified here, no arbitrary Keller pair is
placed in the degree-22 response family, and neither `JC(2)` nor `DC(2)` is
proved.

## 5. Exact companion

Run

```bash
python 04-computation/jc2_one_pole_response_composition_thm2748.py
python -O 04-computation/jc2_one_pole_response_composition_thm2748.py
```

Both executions byte-match the stored transcript
`05-knowledge/results/jc2_one_pole_response_composition_thm2748.out`.  The
companion verifies every pole order in `(16)`, all first five ramification
degrees, the exact rational-primitive derivative identity, and the pure-power
composition.  It also checks the constant-`U` degree-one boundary and the
positive control with nonzero target shift `z_0=7` in the inverted source
coordinate `X`.  The proof for arbitrary degrees is Section 2, not a finite
computation.

QED.
