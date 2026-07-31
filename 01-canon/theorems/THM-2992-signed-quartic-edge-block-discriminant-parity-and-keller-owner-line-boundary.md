---
id: THM-2992
title: "Signed quartic edge, block-discriminant parity, and Keller owner-line boundary"
status: PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED
source: codex-quartic-v4-origin-probe-2026-07-30
depends_on:
  - THM-2633-derangement-character-obstruction-and-d4-keller-exclusion
  - THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
related:
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary
script: 04-computation/quartic_signed_edge_block_parity_thm2992.py
output: 05-knowledge/results/quartic_signed_edge_block_parity_thm2992.out
script_sha256: 54165f8323a92bba8e5dcbfe7d9c0e617028e7d61609956e5309d78e8387fb8f
output_sha256: 440dd4b921a97403be2419817e6d8be7087dd1779ed9de2d45bd81f46287f938
hash_basis: LF-normalized bytes
---

# THM-2992 -- signed quartic edge and the intrinsic fixed block

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

[THM-2864](THM-2864-quartic-edge-orientation-sextic-resolvents-and-d8-radicand-product.md)
constructs the signed-edge sextic `E(Y)=S(Y^2)` above the matching cubic of a
depressed quartic.  [THM-2968](THM-2968-quartic-edge-and-oriented-cycle-s4-complements.md)
shows that a quotient transposition in `S4/V4~=S3` has both transposition and
four-cycle lifts.  [THM-2974](THM-2974-discriminant-cover-integral-order-smith-and-owner-boundary.md)
then gives the sharp warning that an inertia-fixed sextic sheet need not be a
fixed quartic owner.

The missing local coordinate is smaller than a quartic origin.  On a fixed
signed-edge lane, the two quadratic block discriminants recover the intrinsic
**unordered two-sheet fixed block**.  They do not order that block and do not
choose one affine owner sheet inside it.

Let `R` be a strictly henselian discrete valuation ring, let `K=Frac(R)`, and
let

```text
v:K* -> Z
```

be its normalized valuation.  Assume that the residue field is separably
closed and has characteristic different from `2`.  Let

```text
f(T)=T^4+pT^2+qT+r,                    p,q,r in R          (1)
```

be separable over `K`, and put

```text
S(U)=U^3+2pU^2+(p^2-4r)U-q^2.                           (2)
```

Assume that there is a fixed signed edge

```text
s in R*,                 u=s^2 in R*,
S(u)=0,                  S'(u) in R*.                    (3)
```

Define

```text
a=(p+u+q/s)/2,                   b=(p+u-q/s)/2,
Q_s(T)=T^2-sT+a,                 Q_-s(T)=T^2+sT+b,        (4)

delta_s=Disc(Q_s)=-2p-u-2q/s,
delta_-s=Disc(Q_-s)=-2p-u+2q/s.                          (5)
```

Then

```text
f=Q_s Q_-s,                                                   (6)

delta_s+delta_-s=-4p-2u,
delta_s-delta_-s=-4q/s,
delta_s delta_-s=16r-4pu-3u^2,                              (7)

Res_T(Q_s,Q_-s)=S'(u)=3u^2+4pu+p^2-4r,                     (8)

Disc_T(f)=delta_s delta_-s S'(u)^2.                         (9)
```

Let `B_s` and `B_-s` be the two unordered root sets of the quadratic blocks.
Their inertia action is completely determined by

```text
epsilon_s=v(delta_s) mod 2,
epsilon_-s=v(delta_-s) mod 2:                               (10)
```

| parity | quartic inertia | fixed quartic sheets |
|---|---|---|
| `(0,0)` | identity | all four |
| `(1,0)` | transposition inside `B_s` | exactly `B_-s` |
| `(0,1)` | transposition inside `B_-s` | exactly `B_s` |
| `(1,1)` | double transposition | none |

In the two middle rows, the two-sheet fixed set is intrinsic.  Replacing
`s` by `-s` exchanges the names `B_s,B_-s` and the two entries of `(10)`,
but it does not change the underlying fixed subset of quartic sheets.
Consequently `(10)` recovers no canonical `+/-` label.

## 2. Exact factorization and discriminants

Write `h=q/s`.  Before imposing the matching-resolvent equation, direct
multiplication gives

```text
Q_sQ_-s-f=((p+u)^2-h^2-4r)/4=S(u)/(4u).                  (11)
```

Here `u` is a unit, so `(3)` makes the right side zero and proves `(6)`.
The two formulas in `(5)` follow immediately from the discriminants of the
monic quadratics.  Their sum and difference give the first two identities in
`(7)`.  Moreover `S(u)=0` is equivalent to

```text
h^2=(p+u)^2-4r.                                          (12)
```

Substitution into

```text
delta_s delta_-s=(2p+u)^2-4h^2
```

gives the product in `(7)`.

The `4 x 4` Sylvester determinant of the two quadratics is, before using
`(12)`,

```text
Res_T(Q_s,Q_-s)=h^2+2pu+2u^2.                            (13)
```

Equation `(12)` turns `(13)` into exactly `(8)`.  Finally, the discriminant
of a product of two monic coprime polynomials is

```text
Disc(Q_sQ_-s)=Disc(Q_s)Disc(Q_-s)Res(Q_s,Q_-s)^2,
```

which proves `(9)`.  The unit assumption on `S'(u)` isolates the two blocks:
the local degeneration, if any, occurs only inside them.

## 3. Square classes and the inertia table

For every nonzero `x in K`,

```text
x is a square in K  iff  v(x) is even.                    (14)
```

Indeed, the forward implication is immediate.  Conversely write
`x=pi^(2m)c` with `c in R*`.  The residue of `c` has a square root because
the residue field is separably closed and has characteristic not `2`; the
simple root lifts by Hensel's lemma.  Thus

```text
K*/K*2 ~= Z/2,                                           (15)
```

generated by the class of a uniformizer.

The roots of `(4)` are

```text
(s +- sqrt(delta_s))/2,
(-s +- sqrt(delta_-s))/2.                                (16)
```

If both valuations in `(10)` are even, both square roots lie in `K` and
inertia is trivial.  If exactly one is odd, the corresponding quadratic is
the unique ramified quadratic extension and inertia swaps precisely its two
roots.  The other two roots are individually fixed.

If both valuations are odd, their quotient has even valuation and hence is
a square by `(14)`.  Thus both square roots generate the same ramified
quadratic extension.  Its nontrivial inertia element swaps both pairs at
once, giving a double transposition rather than two independent quadratic
inertia generators.  This proves the table.

## 4. The `S4/V4` lift census

Freeze the matching

```text
{{0,1},{2,3}}.                                           (17)
```

A transposition of the three-matchings quotient which fixes `(17)` has four
lifts in `S4`:

```text
(0 1),        (2 3),        (0 2 1 3),        (0 3 1 2). (18)
```

The first two are sheet transpositions.  They fix the opposite blocks
`{2,3}` and `{0,1}` respectively, and differ by the matching-direction
translation `(0 1)(2 3) in V4`.  The last two are four-cycles: they exchange
the two blocks and fix no quartic sheet.  Conjugating `(18)` through the
other two matchings gives all twelve lifts of the three quotient
transpositions:

```text
6 sheet transpositions,             6 four-cycles.       (19)
```

Thus the matching cubic and the quotient transposition do not alone select
a quartic fixed block.  On the fixed signed-edge lane, the parity table does:
a transposition has one even-discriminant block, and that block is exactly
its intrinsic fixed set.  Changing the sign of the edge merely exchanges
the two coordinate labels used to describe the same result.

## 5. Conditional Keller owner-line consequence

Let `D` be a sign-support Jelonek component of a degree-four `S4` Keller
cover.  Assume that its transverse local quartic lies in the separable fixed
signed-edge unit lane `(1)--(3)`.  THM-2633 proves that the local inertia is a
sheet transposition and that the number `k_D` of finite affine branches obeys

```text
1<=k_D<=2.                                               (20)
```

The inertia table therefore forces exactly one odd block discriminant.  Let
`B_fix` be the even-discriminant block.  It is an intrinsic unordered
two-sheet set, and every finite affine branch is inertia-fixed.  Hence

```text
empty != S_D subseteq B_fix,
|S_D|=k_D in {1,2}.                                      (21)
```

If `k_D=2`, `(21)` gives `S_D=B_fix`.  If `k_D=1`, the theorem does not
choose which one of the two fixed roots is present in the affine source.
That choice is a Zariski-main present/omitted coordinate, not a resolvent or
valuation-parity coordinate.

The conclusion is local to one divisor and one fixed signed-edge unit chart.
It supplies neither compatibility between divisors nor a global degree-four
Keller realization.  It also does not repair the double-transposition hostile
of THM-2974, which lies outside the transposition row of the table.

## 6. Sharp controls and stopping boundaries

### 6.1 Signed-edge reversal is only a label gauge

The simultaneous replacement

```text
(s,q/s) -> (-s,-q/s)                                    (22)
```

fixes `f,u,S,Disc(f)` and exchanges

```text
Q_s <-> Q_-s,                  delta_s <-> delta_-s.     (23)
```

Therefore no intrinsic `+/-` name survives.  In the transposition row,
however, `(23)` also exchanges the parity coordinates, so the geometric
even block remains the same fixed two-sheet set.  Sign reversal is a gauge
hostile to labels, not a hostile to the intrinsic block decoder.

### 6.2 Exact transposition control

Over `R=C[[t]]`, take

```text
f_t(T)=T^4-(3+t)T^2+2(1-t)T
      =(T^2-2T+1-t)(T^2+2T).                            (24)
```

Here `s=2`, `u=4`, and

```text
(delta_s,delta_-s)=(4t,4),
S'(4)=(t-9)(t-1) in R*.                                 (25)
```

Inertia swaps `1+-sqrt(t)` and fixes `0,-2`, so the intrinsic fixed block is
literal.

### 6.3 Exact double-transposition control

The family

```text
(T^2-2T+1-t)(T^2+2T+1-2t)                              (26)
```

has `s=2,u=4` and block discriminants `(4t,8t)`.  Since `2` is a square in
`C((t))`, both square roots lie in the same ramified quadratic extension and
its inertia swaps both blocks.  The block resultant specializes to `16`, so
this is a genuine within-block double collision, not a cross-block collision.

### 6.4 Eisenstein four-cycle boundary

Consider

```text
g_t(T)=T^4+tT+t.                                        (27)
```

It is Eisenstein and tamely totally ramified over `C((t))`, so local inertia
acts as a four-cycle.  Its matching cubic is

```text
S_t(U)=U^3-4tU-t^2.                                     (28)
```

Writing `U=tz` gives

```text
t^(-2)S_t(tz)=tz^3-4z-1.                                (29)
```

The simple residue root `z=-1/4` lifts, so the fixed matching value `u` lies
in `K` with `v_K(u)=1`.  Consequently no signed edge `s` with `s^2=u` lies in
`K`.  After adjoining `s`, normalize the extended valuation by
`v_L(t)=2`; then

```text
v_L(u)=2,                 v_L(q/s)=1,
v_L(delta_s)=v_L(delta_-s)=1.                            (30)
```

The two odd parities now describe the square of the original four-cycle,
namely the double-transposition inertia over the signed-edge extension.
They do not classify the original four-cycle over `K`.  This is why the
fixed signed-edge hypothesis in `(3)` is load-bearing.

## 7. Preserved and destroyed information

The exact bridge is

```text
source:     one fixed signed edge s over a strict-henselian DVR;
map:        quadratic factorization followed by
            (v(delta_s),v(delta_-s)) mod 2;
target:     local inertia type and, for a transposition,
            its intrinsic unordered fixed two-sheet block;
preserved:  the matching channel, fixed-root block, and local cycle type;
destroyed:  +/- block labels, the present sheet inside a one-owner block,
            divisor gluing, and global affine/Keller realization;
sidecar:    Zariski-main present/omitted data when k_D=1.
```

The first failed stronger implication is therefore

```text
intrinsic fixed two-sheet block
    does not imply
one canonically selected affine owner sheet.             (31)
```

## 8. Exact evidence

The standalone companion checks `(6)--(13)` in an exact symbolic quotient
ring, including the off-shell error in `(11)`, the resultant normalization,
the product-discriminant identity, and the sign-gauge exchange.  It enumerates
all `24` elements of `S4`, all three quotient transpositions, and all twelve
lifts in `(19)`.  It also checks the exact identity, transposition,
double-transposition, and Eisenstein four-cycle controls above.

Normal and optimized runs are LF-identical to the stored transcript.  Reproduce
with

```text
python 04-computation/quartic_signed_edge_block_parity_thm2992.py
python -O 04-computation/quartic_signed_edge_block_parity_thm2992.py
```

Frozen LF hashes are

```text
script  54165f8323a92bba8e5dcbfe7d9c0e617028e7d61609956e5309d78e8387fb8f
output  440dd4b921a97403be2419817e6d8be7087dd1779ed9de2d45bd81f46287f938
```

The proof, local valuation typing, intrinsic-block scope, four-cycle boundary,
and exact companion have passed an independent hostile audit.

**QED.**
