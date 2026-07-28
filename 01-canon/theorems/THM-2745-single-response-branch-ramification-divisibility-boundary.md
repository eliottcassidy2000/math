---
id: THM-2745
title: "Single response-branch ramification and monomial-composition boundary"
status: >
  PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
  AUDIT.  A rational primitive on the source P1 that factors through a
  projective integral response component with one response pole forces the
  component normalization to be P1.  After placing the two pole points at
  infinity, both the response and the normalization map are pure powers.
  If the target pole has order P, the constant-U case forces P=degree(nu)=1;
  the monomial-U case forces m-1=P*degree(nu).  Applied to any one-branch
  escape from THM-2741, constant U is impossible and the response component
  is subject to this full monomial-composition normal form.  This does not
  prove that a reducible component contains a forced branch or close JC(2).
source: thm2694-full-lift-fibre-scout-2026-07-28
depends_on:
  - THM-2723-split-exact-square-prefix-rational-primitive-pole-capacity
  - THM-2741-highest-odd-faber-response-pole-capacity-closure
related:
  - THM-2214-nonsplit-terminal-quartic-spectral-curve-closure-through-degree-ten
  - THM-2726-a21-transverse-integral-split-response-three-pole-closure
script: 04-computation/jc2_one_pole_response_composition_thm2745.py
output: 05-knowledge/results/jc2_one_pole_response_composition_thm2745.out
script_sha256: 7abb11dfed2ee12e4367489ebff8acc305c1812aca37daba548ab493bf36d0cb
output_sha256: e891aadabbe9af77f69062470ffe9e3626bca6fa661d60bac91b1ce7e8671ad6
hash_basis: LF-normalized bytes
---

# THM-2745 -- one response pole forces a monomial composition

**PROOF-COMPLETE CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE
AUDIT.**

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

Let a physical trajectory land on the normalization of an irreducible
component of a THM-2741 highest-odd response member.  Suppose that component
contains one of THM-2741's forced infinity branches and that the response
has no other pole on the component.  Write

```text
j=max{odd k:a_k!=0},       r=22-j,       g=gcd(r,6),
P_j=(150-6r)/g.                                          (12)
```

THM-2741 gives the exact table

```text
j       21  19  17  15  13  11   9   7   5   3   1
P_j    144  44 120 108  32  84  72  20  48  36   8.    (13)
```

The source equation is precisely `(2)` by THM-2723.  Since every entry in
`(13)` is at least eight, case `(3a)` is impossible.  Thus every surviving
one-pole component must satisfy the monomial normal form `(5)` with

```text
m=1+eP_j.                                                (14)
```

If the component has two or more response poles, their pairwise disjoint,
nonempty fibres already contradict `(3)`.  Therefore `(5)` is the complete
equality boundary for any physical component that contains a forced branch.

This strictly strengthens the divisibility sentence in THM-2741: not only
the pole order, but the entire target response and normalization map are
pure powers after affine/Mobius gauge.

## 4. Boundaries

The one-pole hypothesis is load-bearing.  The theorem does not prove that a
reducible or nonreduced degree-22 member has a component containing one of
the forced THM-2741 branches.  A component avoiding that singular branch
requires a separate global `h=0` divisor analysis.  A component containing
more than one pole is excluded, not classified by `(5)`.

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
python 04-computation/jc2_one_pole_response_composition_thm2745.py
python -O 04-computation/jc2_one_pole_response_composition_thm2745.py
```

Both executions byte-match the stored transcript
`05-knowledge/results/jc2_one_pole_response_composition_thm2745.out`.  The
companion verifies every pole order in `(13)`, all first five ramification
degrees, the exact rational-primitive derivative identity, and the pure-power
composition.  It also checks the constant-`U` degree-one boundary and the
necessity of the source/target translations.  The proof for arbitrary
degrees is Section 2, not a finite computation.

QED (candidate; independent hostile audit pending).
