---
id: HYP-2680
title: LRC(14) signed multi-far boundary hierarchy
status: OPEN; exact Newton/Stirling boundary formula, analytic resonance bound missing
source: codex-2026-06-20-S51
depends_on:
  - THM-548
  - HYP-2679
  - HYP-2678
  - HYP-2677
  - HYP-2676
  - HYP-2675
  - THM-547
  - THM-546
  - THM-531
  - HYP-2648
  - HYP-2639
related:
  - HYP-2638
  - HYP-2637
  - HYP-2636
  - HYP-2633
  - HYP-2614
  - HYP-2608
  - OPEN-Q-108
---

# HYP-2680 - Signed Multi-Far Boundary Hierarchy

## Claim Being Tested

The signed two-far bound of THM-548 and HYP-2679 is the second Newton
difference of one general boundary hierarchy.  If `B` is a bounded core and
`F` is a far set, define

```text
Delta_S(B) = sum_{T subset S} (-1)^{|S|-|T|} p0(B union T).
```

The exact finite expansion is the ordinary Newton/Mobius identity

```text
p0(B union F) = sum_{S subset F} Delta_S(B).
```

The earlier shorthand "|S| <= 6" is unsafe without qualification: duplicate
or redundant far-sector hits still contribute to higher Newton differences.
The six-sector structure instead appears in the coefficient formula below,
because `p_t(B)=0` for `t>6`.

Let `p_t(B)` be the exact measure of points where `B*x` misses exactly `t` of
the six inner sectors.  In the fully decorrelated far-runner limit, the signed
`s`-far Newton term should be

```text
Phi_s(B) = 7^-s * sum_{t=1}^s (-1)^(s-t) * t! * S(s,t) * p_t(B),
```

where `S(s,t)` is a Stirling number of the second kind.  This recovers the
known cases

```text
Phi_1(B) = p1(B)/7
Phi_2(B) = (2*p2(B) - p1(B))/49
Phi_3(B) = (p1(B) - 6*p2(B) + 6*p3(B))/343.
```

The corresponding fully decorrelated boundary value for `r` far runners is the
Newton sum

```text
P_r(B) = p0(B) + sum_{s=1}^{r} binom(r,s) * Phi_s(B),
```

equivalent to THM-548's formula

```text
sum_t p_t(B) * sum_{i=0}^t (-1)^i binom(t,i) (1 - i/7)^r.
```

## Proof Obligation

The remaining analytic content is not the formula above; it is the signed
resonance discrepancy

```text
R_S(B; f_1,...,f_s) = Delta_S(B) - Phi_s(B).
```

For `s=2`, THM-548 identifies frequencies `m*f_1+n*f_2` and a sharp
two-far constant `C_2=13/(4*7^3)`.  The next target is the three-far analogue:
bound

```text
Delta_{u,v,w}(B) - (p1 - 6*p2 + 6*p3)/343
```

by a signed Abel packet whose denominators are controlled by all small linear
forms

```text
m*u+n*v+l*w.
```

This is stricter than pairwise resonance control.  A triple can have no zero
pair relation and still carry a three-body relation such as `u - 2v + w = 0`.
That is the exact place where the user's "two things with opposite bounded
signs" becomes a relation-lattice statement: rank one gives a finite
scale-invariant ledger; higher rank should pay a signed dimension penalty.

## Expected Split

```text
No small relation among far speeds
  -> signed Abel / Koksma decay around Phi_s(B).

Rank-one far relation lattice
  -> Freiman/scale reduction to a finite boundary atlas.

Rank >= 2 relation lattice
  -> multi-form discrepancy bound with dimension penalty.
```

The rank split preserves the LRC predicate `p0(B union F) <= cap_k` through the
Newton expansion and destroys only the order in which far speeds are added.
That loss is acceptable because the mixed differences are symmetric after
summing over subsets.

## Challenged Assumption

Do not assume tournament vertices must be runners or far speeds.  In this
route the useful vertices can be:

- Newton orders `s>=1`;
- relation-lattice ranks;
- proof obligations `Phi_s`, `R_s`, and direct cap margin;
- small linear forms among far speeds;
- state-word packets of the bounded core.

The quotient preserves signed boundary corrections.  It destroys individual
runner identity except through the retained relation forms.

## Next Computation

Completed first scout:

- `04-computation/lrc14_signed_multifar_boundary_hierarchy_codex_s51.py`
- `05-knowledge/results/lrc14_signed_multifar_boundary_hierarchy_codex_s51.out`

The script:

1. verifies the Stirling coefficients beyond `s=6`;
2. scans exact three-far triples against `Phi_3(B)` on the known true-wide
   cores;
3. records pair and triple resonance distances separately;
4. runs Tournament Analysis on proof-obligation vertices rather than runners.

## Exact S51 Findings

The coefficient table begins

```text
s=1: [1] / 7
s=2: [-1, 2] / 49
s=3: [1, -6, 6] / 343
s=4: [-1, 14, -36, 24] / 2401
s=5: [1, -30, 150, -240, 120] / 16807
s=6: [-1, 62, -540, 1560, -1800, 720] / 117649
```

The script checks the identity

```text
P_r(B) = p0(B) + sum_{s=1}^r binom(r,s) Phi_s(B)
```

against THM-548's direct `sum_t p_t(B)c_t(r)` formula for three cores and
`r=1..6`; all checks are exact equalities.

Three-far deviations behave like a relation-lattice object.  Named rows:

```text
B=(0,4,6,8,10,12,14), F=(15,16,17):
  Delta_3=269/99960
  Phi_3=-22/252105
  dev=31753/11428760
  p0=23869/66640, margin=213303/866320
  triple relation -15+2*16-17=0

B=(0,4,6,8,10,12,14), F=(15,18,21):
  dev=18443/3025260
  p0=779/2940, margin=12973/38220
  triple relation -15+2*18-21=0

B=(0,1,2,3,4,5,6,7), F=(17,23,31):
  dev=421/71900346
  p0=1548425/3563574, margin=13469887/46326462
```

The separated triple is already nearly decorrelated, while exact low-height
triple relations create much larger signed deviations.

The small all-core bank with far triple `(15,16,17)` checks `3003` primitive
bounded cores.  Every row has the exact triple relation `-u+2v-w=0`, no row has
an exact pair relation at coefficient height `4`, and deviation signs split

```text
positive: 1999
negative: 1004
zero:        0
```

The top abs deviation is

```text
row=(0,5,6,7,11,13,14,15,16,17)
Delta_3=1280449/14294280
Phi_3=-114263/72102030
dev=40633081/445721640
p0=1843399/7147140
margin=2476301/7147140.
```

The top direct `p0` row in that bank is

```text
row=(0,9,10,11,12,13,14,15,16,17)
p0=2290763/5717712
margin=1164997/5717712.
```

So three-far resonance can be much larger than the fully decorrelated
`Phi_3`, but these tested rows still have large direct cap margins.

## Tail-Rank S51 Addendum

The follow-up scout

- `04-computation/lrc14_signed_multifar_tail_rank_codex_s51.py`
- `05-knowledge/results/lrc14_signed_multifar_tail_rank_codex_s51.out`

uses the simultaneous-peel form directly.  For a far block `F`, it verifies the
exact decomposition

```text
p0(B union F) - P_|F|(B)
  = sum_{s=1}^{|F|} sum_{S subset F, |S|=s} (Delta_S(B)-Phi_s(B)).
```

This changes the proof target.  It is too crude to bound every individual
residual absolutely and add them.  The useful object is the signed order sum

```text
R_s(B;F) = sum_{S subset F, |S|=s} (Delta_S(B)-Phi_s(B)).
```

Structured far blocks can have large order terms, but consecutive orders often
carry opposite bounded signs.

Exact named rows:

```text
B=(0,4,6,8,10,12,14), F=(15,16,17):
  p0-P_3=9243793/68572560
  order signs R1,R2,R3 = +-+
  exact relation rank at height 3: 1

B=(0,4,6,8,10,12,14), F=(15,16,17,18):
  p0-P_4=222536009/1440023760
  order signs R1..R4 = +--+
  exact relation rank at height 3: 2

B=(0,4,6,8,10,12,14), F=(15,16,17,18,19,20):
  p0-P_6=249775162037/1340662120560
  order signs R1..R6 = ++-+-+
  exact relation rank at height 3: 5

B=(0,4,6,8,10,12,14), F=(17,23,31):
  p0-P_3=4331573/239667820
  order signs R1,R2,R3 = ++-
  exact relation rank at height 3: 0
```

The all-core four-far bank for `F=(15,16,17,18)` checks `3003` primitive bounded
cores.  It finds:

```text
R2/R3 opposite signs: 1644/3003
R3/R4 opposite signs: 2053/3003
total residual positive: 2966
total residual negative:   37
top direct row: (0,9,10,11,12,13,14,15,16,17,18)
  p0=4671421/9529520
  margin=2240099/9529520
  sign word R1..R4 = +++-
top |R4| row: (0,5,6,7,12,13,14,15,16,17,18)
  R4=10105997/157313520
  sign word R1..R4 = ++++
```

This makes the "two opposite bounded signs" route precise: the three-far proof
should first control `R_3` in the `r=3` peel, then the general proof should
control the signed order sums `R_s(B;F)` with a relation-rank budget.  Exact
low-height relations route to finite Freiman/scale atlases; high-rank blocks
should pay signed Abel/Koksma cancellation and the apex-prime suppression.

The relation-rank budget cannot be a raw rank-only scalar.  Older HYP-2637 and
HYP-2639, plus the S51 explorer audit, point to a typed relation ledger:
summand shell, multiplicand clearance, observer visibility, support size, and
packet sign must remain attached.  Equal energy or equal low-height rank can
still have different signed corrections, so the lemma should be support-
stratified before any absolute values are taken.  The incoming adversarial
THM-548 verification also warns against fixed coefficient boxes: a pair can
hide its shortest relation just beyond `|m|,|n|<=7`.  The multi-far analogue
must weight relation height rather than pretending a small box exhausts the
resonance lattice.

## Revised Next Target

The sharp target is no longer "prove all higher orders vanish."  They do not.
Incoming THM-548 adds the right proof shape: simultaneous peel, not iterative
peel.  For `r=3`, the exact target is

```text
p0(B union {u,v,w}) =
  P_3(B)
  + sum_i       (Delta_i - Phi_1)
  + sum_{i<j}   (Delta_{ij} - Phi_2)
  +             (Delta_{uvw} - Phi_3).
```

The one-far residuals stay over bounded `B` and route to THM-547; the two-far
residuals route to THM-548's `C_2`/simultaneous-peel bound; the new obligation
is the three-far residual

```text
Delta_{uvw}(B) - Phi_3(B).
```

For `r>3`, the sharper object is not only each individual residual but the
signed order ledger `R_s(B;F)`.  These residuals need a signed Abel
relation-lattice packet bound with separate treatment of:

- exact low-height multi-relations such as `u-2v+w=0`;
- near-relations with small `|m*u+n*v+l*w|`;
- genuinely high-rank/nonresonant far sets, where the signed packet should
  decay by an extra apex-prime power.

No proof of LRC(14) is claimed here.
