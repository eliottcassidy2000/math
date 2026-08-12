---
id: THM-3346
title: "U-spine prime-toggle root atlas, two-channel content splitter, and conjugation monodromy"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the
  consecutive-parameter Pythagorean U-spine, shared Gaussian prime powers
  split exactly between the difference and reflected-sum channels.  Fixed
  admissible hypotenuse roots form an explicit CRT prime-toggle cube whose
  antipodal quotient is the fixed-hypotenuse parent torsor.  Filling commuting
  toggle squares leaves exactly one conjugation Z/2 monodromy class.  A
  square-triangular Pell row gives two exact Gaussian compositions into its
  adjacent square-hypotenuse U-spine depths, and these bridge grades have
  unbounded prime-toggle rank.  No Berggren, LRC, tournament, or JC transfer
  follows.
source: codex-2026-08-12-u-spine-prime-toggle-atlas
audit: >
  An independent arithmetic hostile audit rederived all displayed Gaussian
  products, gcd channels, primitive radius laws, CRT/root reconstruction,
  Hensel densities, Pell compositions, folded metric, and unbounded-rank
  indexing.  An independent topology audit proved the free antipodal cover,
  projective two-skeleton model, integral homology, small-rank boundary, and
  weighted combinatorial systole, with separate Smith computations through
  rank seven.  A third inheritance audit found no canonical duplicate and
  fixed the dependency and quotient-loss boundary.  A post-promotion hostile
  audit found and repaired MISTAKE-371: equality of selected-grade and full
  contents is controlled by gcd(C_r,C_s)=N, not only by N=C_r.  Normal,
  optimized, and stored transcripts match after LF normalization; both
  recorded hashes match.
depends_on:
  - THM-3333-gaussian-square-farey-pythagorean-triangular-light-cone
  - THM-3334-berggren-parabolic-spine-gaussian-collision-torsor
  - THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation
  - THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors
related:
  - THM-2606-affine-v4-parity-channels-partial-cubes-and-feuerbach-origin
  - THM-2622-affine-torsor-holonomy-fixed-section-spectrum-and-v4-c13-dictionary
  - THM-2753-six-edge-parity-erasure-and-three-matching-resolvent-restoration
  - THM-3173-six-state-free-factor-actions-and-pointed-frame-cube
  - THM-3335-square-triangular-pell-markov-pythagorean-selector
  - THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost
  - THM-584-complement-is-antipodal-map-level-parity-spectrum
script: 04-computation/u_spine_prime_toggle_atlas_thm3346.py
output: 05-knowledge/results/u_spine_prime_toggle_atlas_thm3346.out
script_sha256: c65c85bedc1e9ab74aab2b71c26bec7ca06f3e0dd0f3cb9f21a9f8e328846eab
output_sha256: e9652fbf1a0f1afbfbdea020354bd85ffdbbec79b0326e99349ba1ecc528dd27
hash_basis: working-tree bytes (LF)
---

# THM-3346 -- U-spine prime toggles are arithmetic clocks and one topological bit

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

No literature-priority claim is made.  Classical sums-of-two-squares,
Hensel, Pell, and cubical-cover facts are assembled here as a repository
proof interface.

## 1. Inheritance, board, and conventions

The closest proved mechanisms are:

- [THM-3334](THM-3334-berggren-parabolic-spine-gaussian-collision-torsor.md),
  which identifies the U-spine, its fixed-hypotenuse Gaussian factor-choice
  torsors, and unbounded CRT ranks;
- [THM-3336](THM-3336-primitive-gaussian-multiplication-content-curved-farey-triangulation.md),
  which proves the general primitive-content cocycle, split-prime charge law,
  and section-dependent Boolean groupoid;
- [THM-3341](THM-3341-u-spine-square-hypotenuse-transplant-and-triangular-plane-torsors.md),
  which classifies square hypotenuses and adjacent Markov/Pell coordinates on
  this spine; and
- [THM-3345](THM-3345-prime-xor-ancestry-path-groupoid-and-source-dependent-berggren-cost.md),
  which proves that prime-XOR motion has source-dependent Berggren cost and
  flat ancestry lift.

THM-2606 and THM-2753 own the three C4 parity channels and six-edge matching
dictionary at rank three; THM-2622 owns the general affine-section warning;
THM-584 owns antipodal hypercube graph spectra; and THM-3173 owns the pointed
rank-three cube/owner quotient.  None fills the prime-toggle squares or
computes the quotient two-complex in Section 6.

The new object is the cross-index arithmetic that makes those prime toggles
visible directly on U-spine positions.  Put

```text
z_t=(t+1)+it,                 C_t=N(z_t)=2t^2+2t+1,       (1)
S_t=Phi(z_t)=(2t+1,2t(t+1),C_t),                         (2)
```

for `t>=0`.  On the positive branch `t>=1`, THM-3333's triangular potential
gives

```text
inradius(S_t)=t,            C_t=4T_t+1.                  (3)
```

Thus the U-spine parameter is simultaneously Farey-fan position, Berggren
depth plus one, triangle inradius, and a triangular-number index.  These are
typed views of one labelled spinor, not interchangeable quotient invariants.

The canonical hostile is THM-3345's result: knowing a prime toggle and even
its folded arithmetic weight does not determine a fixed Berggren word or
cost.  The least-used sidecar is the commuting-square topology implicit in
THM-3336's raw allocation cube; Section 6 computes what survives its
conjugation quotient.

## 2. Shared primes split into exactly two content channels

For `r,s>=0`, direct multiplication gives

```text
z_r z_s
 =(r+s+1)+i(2rs+r+s),                                   (4)
z_r conjugate(z_s)
 =(2rs+r+s+1)+i(r-s).                                   (5)
```

Let `cont(a+ib)=gcd(a,b)` and define

```text
d_+(r,s)=gcd(C_r,r+s+1),
d_-(r,s)=gcd(C_r,s-r).                                  (6)
```

Then the exact primitive contents are

```text
cont(z_r z_s)=d_+(r,s),
cont(z_r conjugate(z_s))=d_-(r,s),                       (7)
```

and

```text
gcd(C_r,C_s)=d_+(r,s)d_-(r,s),
gcd(d_+(r,s),d_-(r,s))=1.                               (8)
```

To prove (7), apply THM-3336's content formula to the primitive spinors.  A
direct Bezout check makes the equality transparent:

```text
(2r+1)(r+s+1)-(2rs+r+s)=C_r,
(2r+1)(r-s)+(2rs+r+s+1)=C_r.                           (7a)
```

Together with the obvious divisibility in (4)--(5), these give the two gcds
in (7).  For (8), use

```text
C_s-C_r=2(s-r)(r+s+1).                                  (9)
```

All `C_r` are odd.  If an odd prime divided both factors in (6), it would
divide both `s-r` and `r+s+1`, hence `2r+1`; but

```text
2C_r=(2r+1)^2+1                                         (10)
```

then gives a contradiction.  Prime-power valuations do not split between
the two factors: modulo every `p^e|C_r`, the two simple roots of `C_s=0` are

```text
s=r                 or                 s=-r-1.           (11)
```

Therefore each common prime power lies wholly in the difference channel
`d_-` (same local Gaussian orientation) or reflected-sum channel `d_+`
(opposite orientation).

The minimal control using both channels is

```text
r=6, s=23:
gcd(C_6,C_23)=gcd(85,1105)=85=17*5,
d_-=17, d_+=5,
z_6 z_23=30+305i=5(6+61i),
z_6 conjugate(z_23)=306-17i=17(18-i).                   (12)
```

Squarefreeness is not assumed.  At `r=3,s=21`, the common factor `25` lies
entirely in `d_+=25`.

Consecutive U-spine norms are pairwise coprime: putting `s=r+1` in (8) gives
both factors one.  Accordingly, adjacent Farey/Berggren motion carries no
shared odd-prime toggle even though distant positions can share arbitrarily
many split primes.

## 3. Primitive compositions have exact radius invoices

After dividing (4)--(5) by their contents, apply the Euclid inradius formula
`n(m-n)` in the positive ordered chamber.  For `r,s>=1`, the primitive product
channel has inradius

```text
rho_+(r,s)=((r+s+1)(2rs-1))/d_+(r,s)^2.                 (13)
```

For `0<=r<s`, the primitive conjugate-product channel has inradius

```text
rho_-(r,s)=((s-r)(2rs+2r+1))/d_-(r,s)^2.                (14)
```

Thus the prime toggle is metric: opposite and same local orientations remove
different square content factors from two different bilinear radius laws.
This specializes THM-3336's radius-matrix covariance; it is not a new V4
circle action.

Self-multiplication is a useful boundary.  Here `d_+(r,r)=1` and
`d_-(r,r)=C_r`, so `z_r conjugate(z_r)` reduces to the unit while `z_r^2`
usually leaves the U-spine.  Among positive `r`, only

```text
(2+i)^2=3+4i                                             (15)
```

has consecutive output coordinates.  Gaussian self-squaring on this lane
must not be identified with THM-3341's state-dependent Berggren-ray
transplant.

## 4. The fixed-grade U-spine root atlas

Let

```text
N=product_(j=1)^k p_j^e_j,       p_j=1 mod 4,            (16)
```

with the `p_j` pairwise distinct, every `e_j>=1`, and `N>1`.  Define the
modular U-spine roots

```text
R_N={t mod N : N divides C_t}.                           (17)
```

The word *modular* is load-bearing: (17) is not the usually much smaller set
of integers satisfying `C_t=N`.  With `x=2t+1`, (10) gives the bijection

```text
R_N <--> {x mod N : x^2=-1 mod N}.                       (18)
```

Each prime power in (16) has two Hensel roots.  Hence `R_N` is a free
transitive affine `F_2^k` torsor.  More explicitly, let `E_j` be the CRT
idempotent satisfying

```text
E_j=1 mod p_j^e_j,            E_j=0 mod N/p_j^e_j.       (19)
```

The prime-power toggle is

```text
T_j(x)=(1-2E_j)x,
T_j(t)=t-E_j(2t+1) mod N.                                (20)
```

The `T_j` commute, are involutions, and act freely transitively.  Their total
product is

```text
x -> -x,                       t -> N-1-t.               (21)
```

This is global Gaussian conjugation up to a unit.  Indeed the Gaussian gcd

```text
G_N(t)=gcd_(Z[i])(N,z_t)                                (22)
```

has norm `N`, is defined up to a unit, and reconstructs the selected factor
allocation.
Prime-power by prime-power, a root in (18) selects exactly one of the two
Gaussian factors above `p_j`; consecutive coordinates make `z_t` primitive,
so it cannot select both.  Their CRT product has norm `N`, which proves the
norm assertion in (22).
Equivalently one can use the gcd with `(2t+1)+i`, after a compatible
unit/conjugation choice, because `(1+i)z_t=i((2t+1)-i)`.
Conversely, for `g=m+in` of norm `N`, its root can be taken as

```text
x=n m^(-1) mod N                                        (23)
```

after a compatible unit choice.  Quotienting (17) by (21) therefore recovers
exactly THM-3334's fixed-hypotenuse parent torsor

```text
R_N/(t~N-1-t)  isomorphic to F_2^k/<1>.                  (24)
```

The residue `t` is a clock/root address, not itself a parent; (22) followed
by `Phi` is the map.

For two roots `r,s in R_N`, define the `N`-primary channels

```text
delta_-(r,s)=gcd(N,s-r),
delta_+(r,s)=gcd(N,r+s+1).                              (25)
```

They are coprime and satisfy

```text
delta_-(r,s)delta_+(r,s)=N.                             (26)
```

Locally, `p_j^e_j|delta_-` iff their signs agree at `j`, while it divides `delta_+`
iff they differ.  Thus the unordered pair

```text
{delta_-(r,s),delta_+(r,s)}                             (27)
```

is exactly THM-3336's intrinsic folded content weight `{P_S,N/P_S}`.
It gives the `N`-primary parts of the full contents in (7).  More exactly,

```text
delta_-(r,s)=gcd(N,d_-(r,s)),
delta_+(r,s)=gcd(N,d_+(r,s)).                            (27a)
```

Indeed `N|C_r`, so taking `gcd(N,-)` in (6) gives (27a).  Since both pairs
are coprime and their products are respectively `N` and `gcd(C_r,C_s)`,

```text
(delta_-,delta_+)=(d_-,d_+)
iff gcd(C_r,C_s)=N.                                     (27b)
```

The fixed-grade specialization `N=C_r` with `N|C_s` is sufficient, but not
necessary.  The smallest counterexample to the old necessity wording is
`N=5,r=3,s=6`: here `(C_r,C_s)=(25,85)` and both channel pairs equal `(1,5)`
although `N!=C_r`.

The smallest modular/full-content mismatch is `N=5,r=3,s=21`.  Both indices
lie in `R_5`, but

```text
(delta_-,delta_+)=(1,5),       (d_-,d_+)=(1,25).         (27c)
```

The distinct-prime hostile `N=5,r=6,s=23` remains useful: its selected pair
is `(1,5)` while the full pair from (12) is `(17,5)`.  Thus the modular atlas
sees the selected grade and must not be advertised as the full gcd of the two
larger norms.

The first nontrivial controls are

```text
R_65={23,28,36,41},
R_85={6,23,61,78},
R_1105={23,231,418,431,673,686,873,1081}.               (28)
```

The last raw cube has eight roots; its antipodal quotient is the first
four-parent affine V4 fibre.

## 5. Split primes are exact Hensel clocks on inradius

For a prime `p=1 mod 4` and `e>=1`, let the two roots in (18) modulo `p^e`
be `rho_(p,e)` and `-1-rho_(p,e)`.  Then

```text
v_p(C_t)>=e
iff t=rho_(p,e) or -1-rho_(p,e) mod p^e.                (29)
```

Each branch lifts uniquely.  Consequently the natural densities in the
integer U-spine are

```text
density(v_p(C_t)>=e)=2/p^e,
density(v_p(C_t)=e)=2(p-1)/p^(e+1).                     (30)
```

For `p=3 mod 4` both densities are zero.  Across distinct primes these clocks
are exactly independent by CRT.  Formula (11) says two inradii sharing a
prime-power clock are on the same or reflected Hensel branch, with no third
case.

This is an exact arithmetic clock and a possible sidecar, not an LRC runner
clock.  No map from these residues to a lawful 13-speed owner/phase packet is
proved.

## 6. Filling toggle squares leaves one conjugation monodromy

Assume `k>=3`.  Let `Q_k^(2)` be the two-skeleton of the raw allocation cube:
vertices are roots, edges are single prime toggles, and square faces express
commuting pairs.  The antipode (21) is free on all cells of dimension at most
two exactly because `k>=3`.  Define

```text
K_N=Q_k^(2)/<antipode>.                                  (31)
```

Every edge loop in `Q_k^(2)` reduces by cancelling backtracks and commuting
adjacent toggles.  Equivalently, attaching the cube's cells of dimension at
least three does not change its fundamental group, so the two-skeleton is
simply connected.  It is therefore the universal double cover of `K_N`, and

```text
pi_1(K_N)=Z/2,
H_1(K_N;Z)=Z/2,
H^1(K_N;Z)=0,
H^1(K_N;F_2)=F_2.                                       (32)
```

Conceptually, `Q_k^(2)` is the two-skeleton of the cubical boundary
`partial I^k`, and central reflection is the antipodal action on
`partial I^k isomorphic to S^(k-1)`.  Hence `K_N` is exactly the two-skeleton
of the induced quotient CW structure on `RP^(k-1)`.  The mod-two class in
(32) is the deck class, equivalently the restriction of the standard first
Stiefel--Whitney class.  This is not the one-cell-per-dimension CW structure
and does not make `K_N` a surface when `k>=4`.

A quotient loop lifts with coordinate-toggle parity vector either all zero
or all one.  Any one coordinate parity evaluates the class.  Total edge-count
parity equals `k` times this class and therefore misses it when `k` is even.

The cell counts are

```text
f_0=2^(k-1),
f_1=k 2^(k-2),
f_2=binom(k,2)2^(k-3).                                  (33)
```

Thus

```text
chi(K_N)=2^(k-4)(k^2-5k+8),
H_2(K_N;Z)=Z^(chi(K_N)-1).                              (34)
```

Here `H_2` is free because it is a subgroup of the free cellular group
`C_2`; its rank follows from Euler characteristic and the rational vanishing
of `H_1` in (32).

At `k=3`, `(f_0,f_1,f_2)=(4,6,3)` and `K_N` is a cubical model of
`RP^2`.  Its one-skeleton is `K4`; the three square boundaries are precisely
THM-2606's three translation-invariant C4 parity channels.  Their union is
the quotient surface, not three unrelated matching decorations.  For
`k>=4`, `K_N` is a two-complex with many faces at each edge, not a surface.

Small ranks are different.  At `k=1`, interval reflection has a fixed
midpoint.  At `k=2`, square rotation has a fixed centre; the topological
quotient of the filled square is a disk, although the quotient one-skeleton
is two parallel edges forming a circle.  The free-cover theorem (31)--(33)
starts sharply at `k=3`.

THM-3345's source-dependent Berggren arrows fill every commuting square and
the Berggren tree has no graph `H_1`.  Hence the nonzero class in (32) maps to
zero under the ancestry path functor.  It records only the obstruction to a
global conjugation section, not an ancestry current, owner, or flux.

## 7. Folded content is the quotient geodesic invoice

Give coordinate `j` the edge weight

```text
w_j=log(p_j^e_j).                                       (35)
```

Let `d_N` be the quotient metric on the parent torsor induced from this
weighted raw cube: minimize over the two conjugate lifts of either endpoint.
For `k>=3`, this is exactly the weighted graph metric on `K_N^(1)`; the
definition also remains meaningful at the small-rank boundaries.  If roots
differ on a subset `S`, then (25)--(26) give

```text
d_N([r],[s])
 =min(sum_(j in S)w_j,sum_(j notin S)w_j)
 =log min(delta_+(r,s),delta_-(r,s)).                   (36)
```

Thus the folded content is literally the weighted geodesic distance invoice
on the quotient one-skeleton, not merely an edge colour.  The filled squares
do not change this graph metric.  For `k>=3`, define the weighted
combinatorial systole `sys_1^comb` to be the minimum length of an edge loop in
`K_N^(1)` representing the nontrivial fundamental-group class, equivalently
the nonzero class in `H_1(K_N;F_2)`.  Every such loop lifts from a root to its
antipode, toggles every coordinate oddly, and has weighted length at least

```text
sys_1^comb(K_N)=sum_j w_j=log N,                         (37)
```

attained by toggling each prime once.  This is torsion energy: the stable real
homology norm is zero because `H_1(K_N;Z)` is torsion.

## 8. Pell rows compile two compositions into adjacent square depths

Let

```text
T_n=n(n+1)/2=R^2                                        (38)
```

be a positive square-triangular row, and put

```text
M_-=2n+1-2R,              M_+=2n+1+2R,
t_-=2R-n-1,               t_+=2R+n,
H=4R^2+1=M_-M_+.                                        (39)
```

THM-3341 proves that `t_-` and `t_+` are consecutive positions on the
square-hypotenuse selector:

```text
C_(t_-)=M_-^2,             C_(t_+)=M_+^2.               (40)
```

Two spinors of norm `H` are

```text
u=(n+1)+in,
v=2R+i.                                                     (41)
```

They give two distinguished fixed-hypotenuse parents

```text
Phi(u)=(2n+1,4R^2,H),
Phi(v)=(4R^2-1,4R,H).                                   (42)
```

The complete operation-level bridge is

```text
u conjugate(v)=M_+ z_(t_-),
u v            =M_- (t_+ + i(t_++1)).                  (43)
```

The second output differs from the `z_(t_+)` convention only by the retained
leg/unit gauge.  Its contents are exactly `M_+` and `M_-`, so primitive
Gaussian composition lands on the two adjacent square U-spine spinors.
Moreover,

```text
H=C_n,                 gcd(M_-,M_+)=1,
W_H([u],[v])={M_-,M_+}.                                  (44)
```

The first equality uses `2R^2=n(n+1)`, the second follows from
`M_-M_+-(2R)^2=1`, and the folded weight follows from the two contents in
(43); `W_H` denotes THM-3336's folded weight, not the complex `K_H`.  Since
`M_-<=M_+`, Section 7 turns the Pell carry itself into the exact
parent-space geodesic invoice

```text
d_H([u],[v])=log M_-.                                   (45)
```

Thus (43) is a factor-to-height compiler: a weighted move inside the grade
`H` fibre removes one complementary factor, and the two primitive outputs
have the adjacent square grades `M_-^2` and `M_+^2`.  This need not be one
prime toggle; the factor can contain several prime-power coordinates.

For the first nondegenerate control `n=8,R=6`, one gets

```text
H=145,  M_-=5, M_+=29,  t_-=3, t_+=20,
u=9+8i, v=12+i,
u v=5(20+21i),       u conjugate(v)=29(4+3i).            (46)
```

The two distinguished parents are `(17,144,145)` and `(143,24,145)`.
At `n=1,R=1`, they collapse after leg exchange to the same `3-4-5` triangle;
this is the minimal hostile to distinctness.  No claim is made that (42)
exhausts every parent of `H` when `H` has additional split-prime choices.
It is also the sharp metric boundary: `M_-=1`, so (45) gives zero distance.

These bridge grades have unbounded prime-toggle rank.  To see this without a
new recurrence argument, choose distinct odd rational primes
`ell_1,...,ell_a`, put `L=product ell_i`, and use THM-3341's square-selector
root

```text
M_+=pell_L.                                               (47)
```

Since `L` is odd, this is the upper coordinate of one positive bridge row;
write its preceding selector root as `M_-`.  Strong divisibility makes the
nontrivial integers `pell_(ell_i)` pairwise coprime divisors of `M_+`, while
Cassini gives `gcd(M_-,M_+)=1`.  Therefore

```text
omega(H)=omega(M_-M_+)>=omega(M_+)>=a.                  (48)
```

All these primes split in `Z[i]` by THM-3341's square-selector argument.
Consequently, for arbitrarily large `k=omega(H)`, the very Pell grades in
(43) support the complex `K_H` with

```text
H_1(K_H;Z)=Z/2,
rank H_2(K_H;Z)=2^(k-4)(k^2-5k+8)-1.                   (49)
```

The first invariant is persistent conjugation torsion; the second is
unbounded free square-filling topology.  Neither survives as a nonzero
Berggren-tree current.

Section 8 unifies four proved carriers without identifying their actions:
the square-triangular Pell compiler supplies `u`, its universal companion
supplies `v`, their primitive contents are adjacent negative-Pell coordinates,
and their outputs are U-spine square-hypotenuse roots.  It supplies no fixed
Berggren branch word; THM-3341 and THM-3345 already forbid that simplification.

## 9. Exact consequence and stopping boundary

The theorem proves:

1. exact two-channel factorization of every shared U-spine prime power;
2. a lossless CRT root-coordinate atlas for fixed-hypotenuse factor choices;
3. exact Hensel clocks and densities on the inradius parameter;
4. one conjugation `Z/2` class after all commuting prime squares are filled;
5. a weighted geodesic and torsion-systole interpretation of folded content;
6. exact primitive-composition radius invoices; and
7. a Pell/Markov two-composition bridge to adjacent square U-spine depths,
   occurring at unbounded prime-toggle rank.

It does not prove a global prime-toggle ancestry transducer, an LRC(14) row or
safety certificate, a tournament orientation, a planar-JC flux, or a new
literature theorem.  The root atlas retains prime phase but loses height;
the parent quotient additionally loses conjugation; the Berggren functor
kills the remaining topological class.  Any frontier use must reattach
height, owner, phase, saturation, labelled columns, and actual consumer data.
The unweighted topology depends only on `k`; `N`, prime labels, exponents,
and the metric in (35) are arithmetic sidecars.

## 10. Exact evidence

The companion audits `90,601` exact integer pair identities, `135,150`
composition-radius rows, `92` modular roots, `188` explicit CRT toggles, `184`
Gaussian reconstructions, `336` Hensel-lift rows, `154` Pell bridge identities,
and integral cubical Smith forms through rank seven (with the rank-eight
closed-form row recorded separately).  Normal and optimized runs match the
stored transcript after LF normalization.  Reproduce with

```bash
python3 04-computation/u_spine_prime_toggle_atlas_thm3346.py
python3 -O 04-computation/u_spine_prime_toggle_atlas_thm3346.py
```

The universal claims rest on the displayed proofs; the finite sweeps are
hostile and reproducibility evidence rather than replacements for them.
