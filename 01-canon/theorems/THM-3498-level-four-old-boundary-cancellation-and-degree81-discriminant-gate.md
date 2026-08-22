---
id: THM-3498
title: "Level-four old-boundary cancellation and degree-81 discriminant gate"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.  For the fixed
  sporadic Keller map of THM-2473, the fourth norm has old-boundary valuation
  v_L(N(J))=-43.  Thus G=L^43 N(J) lies in Q[a,b,c] and is coprime to L.
  A squarefree full-degree 81 good-reduction fibre makes the fourth
  discriminant recursion lawful, and its exact field square class is
  [Delta_4]=[2G].  No irreducibility, image-divisor, fourth nonproperness-set,
  all-level, Jacobian-conjecture, or LRC conclusion is asserted.
source: codex/turnpike-atlas/2026-08-16
audit: >
  The candidate at 95e0c4c7e was independently reconstructed in an isolated
  worktree.  A separate SymPy face factorization, nested-Horner finite-sheet
  evaluation, and literal sign/exponent ledger verify the old-boundary and
  square-class claims.  A genuinely different F_101 implementation uses
  FLINT regular-representation matrices and multiplicative Fourier inversion
  on all 100 nonzero field points, with X=0 held out as a control; it exactly
  reproduces the submitted degree-81 coefficient hash and squarefree gcd.
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2582-odd-block-discriminant-tower-and-composite-jelonek-square-class
  - THM-3495-level-three-sporadic-keller-norm-divisor-and-three-component-nonproperness
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - MISTAKE-413
  - THM-3504-level-four-sporadic-keller-image-prime-and-four-component-nonproperness
scripts:
  - 04-computation/keller_level_four_old_L_boundary_norm_probe_20260816.py
  - 04-computation/keller_level_four_degree81_finite_field_probe_20260816.py
  - 04-computation/jc_level4_boundary_squareclass_independent_audit_20260816.py
  - 04-computation/jc_level4_degree81_fourier_flint_independent_audit_20260816.py
outputs:
  - 05-knowledge/results/keller_level_four_old_L_boundary_norm_probe_20260816.out
  - 05-knowledge/results/keller_level_four_degree81_finite_field_probe_20260816.out
  - 05-knowledge/results/jc_level4_boundary_squareclass_independent_audit_20260816.out
  - 05-knowledge/results/jc_level4_degree81_fourier_flint_independent_audit_20260816.out
script_sha256:
  - a52a992a2455f9a6d8b5a2949b1a956969846df69f9df3492573c5ab864d837c
  - 4039b4081c9f0d95b197d2e3a7581c66433382e53dac3b95fa2526c3a4ba4f2e
  - 53f63e59b5de3cd07645f3031225179105f16cc341d9062a124612361171f817
  - 3d0bee9dd97993160fc7275cb4a96e77893013c421029eada4bbd5b46ac5d3e6
output_sha256:
  - a730eb715177e1be945259c01e048c974be624c67e8882a6d4c4c66293d7b85a
  - f5498e42510641227052f578cd269697876746a6847ba6bb8cd382e726c35169
  - c499294d374851a4fa953cd825206b7d0517a1f098b01bd99920cdd3d3fe40fb
  - aef44b43c5c26a97349129653fa22518688e964cb6a7d2deceb53027484c4567
coefficient_ledger_sha256: 1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7
hash_basis: raw LF bytes for files; ascending exact F_101 coefficient ledger as stated
---

# THM-3498 -- the fourth discriminant gate for the fixed sporadic Keller map

**PROVED + VERIFIED-EXACT + INDEPENDENTLY PROOF-AUDITED.**

Let `F:C^3->C^3` be the fixed sporadic Keller map of THM-2473.  Put

```text
L=27a^2c^2-18abc+16a+b^3c-b^2,
T=4-3bc,
S=27ac^2-9bc+8,
E(w)=Lw^3+Tw-2c.
```

Let `N` be the cubic function-field norm induced by `F`, and let `J` be
THM-3495's primitive absolutely irreducible polynomial normalized by

```text
N(H)=J/(2^35 L^7).
```

For the generically squarefree degree-`3^r` `x`-eliminant of `F^r`, write
`Delta_r` for its discriminant; brackets denote its class in
`Q(a,b,c)^*/Q(a,b,c)^{*2}`.

## 1. Theorem

At the generic divisor `(L)`,

```text
v_L(N(J))=-43.                                           (1)
```

Consequently

```text
G:=L^43 N(J) belongs to Q[a,b,c],    gcd(G,L)=1.          (2)
```

The generic fourth `x`-eliminant has full degree `81` and is squarefree.
Equivalently, after passing to a splitting field, its `27` cubic blocks are
separable and pairwise coprime.  THM-2582's odd-block identity is therefore
applicable and gives the exact field square class

```text
[Delta_4]=[2G].                                          (3)
```

In particular, the old divisor `L` has even valuation at the fourth rung.
Equation (3) retains the rational unit class `[2]`; it is a field square-class
statement, not merely a divisor-parity statement.

## 2. Newton face and the two divergent sheets

At the generic DVR of `(L)`, THM-3495 proves that `c,T,S,D` are units, where

```text
D=18ac-3b^2c+2b,
```

and that `E` has one root of valuation zero and two roots of valuation
`-1/2`.  On either divergent sheet put `u=1/w`.  Then

```text
x=w=u^-1,
y=D/S+O(u),
z=-3(D/S)u+O(u^2),
3xz-2y=-11D/S+O(u).                                     (4)
```

An exact extraction from THM-3495's frozen `66,146`-term coefficient ledger
finds the exposed weight and complete face

```text
max(i-k)=43,
in(J)=-2^58 3^51 13^8 79^4 313^2
      *x^43(3xz-2y)^15.                                 (5)
```

The face has exactly `16` terms.  Substitution of (4) into (5) leaves a
nonzero unit times `u^-43`.  Hence each divergent sheet has valuation

```text
v_L(J(q))=-43/2.                                         (6)
```

There is no cancellation between those two contributions: a field norm is
their product, so valuations add.  The face computation separately prevents
cancellation among equal-weight monomials inside either factor.

## 3. The finite-sheet hostile and exact norm valuation

The point

```text
(a,b,c)=(2/27,1,1) on L=0
```

has `(c,T,S,D)=(1,1,1,1/3)`.  Its finite inverse root is `w=2`, with inverse
point

```text
q=(2,5/6,-7/8).
```

Exact nested-Horner evaluation gives `J(q)!=0`; the reduced rational value
is frozen by SHA256

```text
f4ad5e424ac9abc719f21fb92625e7c25a9e37cc614815dd1b54e43fa5936554.
```

Thus the residual function on the finite sheet is generically a unit and
contributes valuation zero.  Together with (6),

```text
v_L(N(J))=(-43/2)+(-43/2)+0=-43,
```

which proves (1).

## 4. Why denominator clearing gives a polynomial coprime to L

Let `A=Q[a,b,c]` and `U=Spec(A[L^-1])`.  THM-2473 proves that the restriction

```text
F^-1(U) -> U
```

is finite etale of degree three.  Since `J` is regular on the source, its
finite-algebra norm belongs to `A[L^-1]`.  Because `A` is a UFD and `L` is
irreducible, write that Laurent element in reduced form as `P/L^m` with
`P in A` and `L` not dividing `P`.  Equation (1) says exactly `m=43`.
Therefore `G=P` belongs to `A` and has `v_L(G)=0`, proving both assertions in
(2).

This argument proves neither that `G` is primitive over `Z` nor that it is
squarefree or irreducible.  It also does not identify `V(G)` as an image
component.

## 5. Independent degree-81 good-reduction witness

At target `(a,b,c)=(1,1,1)` over `F_101`, form the first three nested inverse
cubic algebras.  Their dimensions are `3`, `9`, and `27`.  Every leading
coefficient, defining-cubic derivative, and inverse-chart denominator is a
unit, so the `27`-dimensional algebra is finite etale.

Take the determinant norm of the fourth cubic core

```text
L(q)X^3+T(q)X-2q_z.                                     (7)
```

The submitted implementation evaluates `82` determinants, interpolates a
degree-at-most-`81` polynomial, and checks an unused determinant.  The
independent audit instead represents every algebra element by its FLINT
regular multiplication matrix and evaluates (7) at all `100` nonzero field
elements.  Multiplicative Fourier inversion on `F_101^*` gives coefficients
`0` in degrees `82,...,99`, a nonzero degree-`81` coefficient, and the same
ascending coefficient-ledger hash

```text
1c05c0fd5ee48fc2dd030aebdb9ad6ddd8185fb933eb91e7e39ff553424ef5a7. (8)
```

The held-out value `X=0` agrees, all direct Horner evaluations agree, and
FLINT gives

```text
gcd(P_4,P_4')=1 in F_101[X].                             (9)
```

Because every denominator used in this fibre is nonzero modulo `101`, (8)--
(9) are a lawful good-reduction witness.  They prove that the corresponding
characteristic-zero leading coefficient and discriminant are not identically
zero.  Over an algebraic closure the finite etale index algebra splits into
`27` points, and its norm is the product of the `27` cubic blocks.  Full
degree and squarefreeness therefore prove that every block is cubic and
separable and that distinct blocks are coprime.

## 6. Sign and power-of-two audit for the fourth square class

THM-2582 gives, for odd block degree `3^r`,

```text
[Delta_(r+1)]=[N(Delta_r)] [Delta_1].                    (10)
```

The level-three normalization is itself worth expanding.  THM-3495 gives

```text
[Delta_3]
 =[-L N(H)]
 =[-J/(2^35 L^6)]
 =[-2J].                                                 (11)
```

Here the exponent `-6` of `L` is even, while the exponent `-35` of `2` is
odd; inversion does not change a square class.  Thus the factor `[2]` in
(11) is mandatory.

Now use (10), the cubic degree of `N`, and `[Delta_1]=[-L]`:

```text
[Delta_4]
 =[N(-2J)][-L]
 =[(-2)^3 N(J)][-L]
 =[8 L N(J)].                                            (12)
```

The two minus signs cancel.  Substituting `N(J)=G/L^43` gives

```text
[8 L N(J)]
 =[8G/L^42]
 =[2G],                                                  (13)
```

because `L^42=(L^21)^2` and `8/2=4` are squares.  This proves (3) without
discarding a constant unit and completes the MISTAKE-413 audit.

## 7. Exact companions, replay, and scope

The submitted companions and both independent audit companions replay
identically under ordinary and optimized Python.  The independent paths
share the already-proved THM-3495 polynomial `J` and the fixed inverse chart,
but they do not share either level-four extraction algorithm.

This theorem concerns one fixed polynomial map in dimension three.  It does
not construct the global polynomial `G`, factor it, determine the generic
degree of `V(J)->V(G)`, prove `gcd(G,HJ)=1`, or show
`S_(F^4)=V(LHJ G)`; those further claims are supplied separately by
THM-3504.  It gives no all-level newest-factor induction, arbitrary
Keller-map classification, Jacobian-conjecture, Dixmier-conjecture, or
Lonely Runner consequence.

**QED.**
