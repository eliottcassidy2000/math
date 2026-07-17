---
id: THM-933
title: THE SHARP LOCAL-DENSITY BLOCK-GLUING THEOREM — exact primitive discrepancy inside blocks, component debt across gaps; LRC(14) lacunary threshold sharpened from 15 to 7 and the uniform all-n threshold from 19 to 8
status: PROVED (elementary circle proof in-file; exact-rational referee passes; Lean algebraic core formalized, sorry-free, standard axioms only)
source: codex-2026-07-16-S21 (HYP-7152; owner asked for the local-density block-gluing theorem, frequent pulls, concrete steps, and formalization)
depends_on: [THM-928 (comparison and certified packet), THM-930 and HYP-7104 (within-block certificate suppliers; not needed for the abstract proof)]
script: 04-computation/lrc14_local_density_block_gluing_codex_S21.py
output: 05-knowledge/results/lrc14_local_density_block_gluing_codex_S21.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCLocalDensityBlockGluing.lean
---

# THM-933 — The sharp local-density block-gluing theorem

**Concurrent synthesis.**  Opus S333 first pushed the `THM-932` namespace
stub while this proof commit was being prepared.  Their target records local
density at a chosen interval scale `ell_i` and pays a coarse `K ell_i` loss.
This independently completed theorem uses the *optimal all-interval loss*
`q_i`, proves the closed product/debt formula, and supplies the exact constants
and referee.  Thus THM-933 is the sharp primitive-discrepancy refinement and a
direct certificate supplier for the coarser THM-932 interface; the two files
are kept distinct rather than overwriting the first-pushed claim.

## 1. Statement

Write `T = R/Z`, with normalized Lebesgue measure `mu`.  For a positive
integer speed `x` and a radius `0 < rho < 1/2`, put

```text
D_rho(x) = {t in T : ||xt|| < rho},
S_rho(B) = T \ union_{x in B} D_rho(x)
```

for every nonempty finite speed block `B`.  Define

```text
delta(B) = mu(S_rho(B)),
M(B)     = sum_{x in B} x.
```

Lift the centered indicator periodically to `R` and define

```text
H_B(u) = integral_0^u (1_{S_rho(B)}(t) - delta(B)) dt,
q(B)   = sup_R H_B - inf_R H_B.
```

The function `H_B` is continuous, piecewise linear, and one-periodic.  Thus
`q(B)` is finite and is obtained at endpoints of the danger teeth.

**Theorem (local-density block gluing).**  Let `B_1,...,B_m` be nonempty
blocks and choose certified density floors

```text
0 <= d_r <= delta(B_r).
```

Put `P_{r-1} = sum_{i<r} M(B_i)`.  Then

```text
mu(intersection_{r=1}^m S_rho(B_r))
 >= product_{r=1}^m d_r
    - sum_{r=2}^m q(B_r) P_{r-1} product_{s=r+1}^m d_s.       (BG)
```

In particular, a positive right-hand side proves a common safe time for all
speeds in the blocks.

The same statement holds when block `r` uses a stronger radius
`rho_r >= rho`: replace `S_rho(B_r)` by `S_{rho_r}(B_r)`.  Its intersection
is contained in the target safe set, so positivity still proves the target
radius.  This is how THM-928(C)'s radius-`1/13` certificate can serve the
radius-`1/14` LRC problem.

## 2. The exact one-interval certificate

**Lemma 1 (primitive discrepancy).**  For every circular interval `I`,

```text
mu(I intersect S_rho(B))
 >= delta(B) mu(I) - q(B)
 >= d mu(I) - q(B)                              (1)
```

whenever `0 <= d <= delta(B)`.  Moreover, `q(B)` is the smallest possible
constant in the first inequality.

**Proof.**  Choose a positively oriented lift `I=[a,b]` with
`0 <= b-a <= 1`.  Periodicity and the definition of `H_B` give

```text
mu(I intersect S_rho(B)) - delta(B)(b-a) = H_B(b)-H_B(a).
```

The right side is at least `inf H_B - sup H_B = -q(B)`.  Replacing
`delta(B)` by the smaller `d` only lowers the desired right side.  Conversely,
choose `a` at a maximum of `H_B` and, after a periodic lift if necessary,
choose `b` at a minimum.  Then the discrepancy is exactly `-q(B)`.  Endpoint
choices have measure zero and do not affect the identity.  ∎

**Lemma 2 (local-to-component sum).**  If `W` is a union of `kappa` circular
intervals, disjoint modulo endpoints, then

```text
mu(W intersect S_rho(B)) >= d mu(W) - kappa q(B).             (2)
```

**Proof.**  Apply (1) to the `kappa` components and sum.  ∎

This is the essential point.  A global density number alone cannot be
multiplied through an arbitrary earlier survivor.  The primitive sidecar `q`
is precisely the missing local information.

## 3. The component ledger

**Lemma 3 (tooth complexity).**  If `X` is a nonempty speed set, then

```text
kappa(S_rho(X)) <= sum_{x in X} x.                            (3)
```

**Proof.**  Since `rho < 1/2`, `D_rho(x)` consists of exactly `x` disjoint
open circular teeth.  List all `N = sum_x x` teeth.  The complement of the
first tooth has one component.  Adding one further open arc to the deleted
set can split at most one surviving component, so it increases the number of
complement components by at most one.  After `N` teeth the complement has at
most `N` components.  Overlaps can only reduce this number.  ∎

Now let

```text
W_r = intersection_{i=1}^r S_rho(B_i),
w_r = mu(W_r),
kappa_r = number of components of W_r.
```

The first block is global, so `w_1 = delta(B_1) >= d_1`.  For `r>=2`,
Lemmas 2 and 3 give

```text
w_r >= d_r w_{r-1} - q(B_r) kappa_{r-1}
    >= d_r w_{r-1} - q(B_r) P_{r-1}.                         (4)
```

All `d_r` are nonnegative.  Inductively substituting (4) therefore yields

```text
w_m >= d_m ... d_1
       - q(B_2)P_1 d_3...d_m
       - q(B_3)P_2 d_4...d_m
       - ...
       - q(B_m)P_{m-1},
```

which is exactly (BG).  This proves the theorem.  ∎

## 4. Scale covariance: the two-scale composition

For an integer `c>=1`, write `cB={cx:x in B}`.  Multiplication by `c` on
the circle is measure-preserving and

```text
S_rho(cB) = {t : ct in S_rho(B)}.
```

Consequently,

```text
delta(cB) = delta(B),
M(cB)     = c M(B).                                          (5)
```

If `h_B=1_{S_rho(B)}-delta(B)`, then `h_{cB}(t)=h_B(ct)` and

```text
H_{cB}(u) = (1/c) H_B(cu).
```

Hence

```text
q(cB) = q(B)/c.                                              (6)
```

Equations (5)-(6) are the promised two-scale interface.  A certificate inside
a template block retains its density under dilation, while the gap scale
divides its boundary discrepancy.  Explicitly, for `B_r=c_r A_r`, (BG) reads

```text
mu(intersection_r S_rho(c_r A_r))
 >= product_r d(A_r)
    - sum_{r=2}^m [q(A_r)/c_r]
        [sum_{i<r} c_i M(A_i)] product_{s>r} d(A_s).          (7)
```

This is “certificates within blocks, cascade across gaps” without an
independence assumption.

## 5. The sharpened lacunary theorem

Take the LRC(`n`) radius `rho=1/n` and a singleton block `{x}`.  Its safe
density is

```text
a = 1 - 2/n = (n-2)/n.
```

Across one danger tooth the centered safe primitive drops by

```text
(2/(nx)) a = 2(n-2)/(n^2 x),
```

and it rises by the same amount across the safe part of the cell.  Therefore

```text
delta({x}) = a,
q({x})     = 2(n-2)/(n^2 x).                                (8)
```

Let `x_1<...<x_{n-1}` and suppose `x_{r+1}/x_r >= R>1`.  Applying (BG) to
singleton blocks and using

```text
(x_1+...+x_{r-1})/x_r
 <= sum_{j=1}^{r-1} R^{-j} < 1/(R-1)
```

gives the strict bound

```text
mu(common safe set)
 > a^(n-1) - [a/(R-1)](1-a^(n-2))
 = [a/(R-1)](R a^(n-2)-1).                                 (9)
```

Thus `R a^(n-2) >= 1` suffices.

### LRC(14)

Here `a=6/7`.  Since

```text
7(6/7)^12 = 6^12/7^11 > 1
```

with `6^12-7^11=199455593`, every packet with consecutive ratios at least
`7` is lonely.  The explicit lower bound from (9) is

```text
mu(common safe set)
 > 199455593/13841287201
 = 0.014410191... .                                         (10)
```

This sharpens THM-928(A)'s elementary threshold `R>=15` to **`R>=7`**.

### Uniform in `n`

Put `m=n-2`.  The standard inequalities

```text
(1+2/m)^m < e^2 < 8
```

give `8((n-2)/n)^(n-2)>1`.  Hence **`R>=8` proves every R-lacunary
LRC(`n`) family simultaneously for every `n>=3`**, sharpening THM-928(A)'s
uniform elementary threshold `19` to `8`.

## 6. A genuinely block-separated LRC(14) certificate

The exact referee uses three templates and scales:

```text
A = {1,2,3,4,5},                       c_A=1,
B = 30 {1,2,4,5}   = {30,60,120,150}, c_B=30,
C = 900 {1,3,5,7}  = {900,2700,4500,6300}, c_C=900.
```

They contain thirteen distinct speeds, but they are not runner-by-runner
7-lacunary: the first block is consecutive.  Exact endpoint arithmetic gives

```text
              delta             q                 M
A             53/105            53/735            15
B             3/5               1/350             360
C             386/735           193/2315250       14400.
```

The two boundary debts are

```text
q(B) M(A) delta(C) = 193/8575,
q(C) (M(A)+M(B))   = 193/6174.
```

Therefore (BG) gives

```text
mu(common safe set) >= 81253/771750 = 0.1052840945... > 0.  (11)
```

Direct exact union gives

```text
mu(common safe set) = 55063/330750 = 0.1664792139...,
```

confirming the theorem with substantial margin.

## 7. Existing certificates become glueable

THM-928(C)'s packet

```text
{300,406,511,652,862,963,1074,1357,1459,1571,1776,1991,2208}
```

has the exact radius-`1/13` BONF5 lower density

```text
d = 624500321285438209432944191
    /15959187221015235325636692600
  = 0.0391310856... .
```

The new exact primitive sidecar is

```text
q = 495710884720685026577242468200179
    /263539708132900820325543980292980400
  = 0.0018809722... .
```

Thus the packet is no longer only a global positive-measure certificate: it
is a local-density block certificate that can enter (BG), with `q/c` after a
scale dilation.  HYP-7104 and THM-930 fit the same interface.  A lower density
certificate alone is not enough; the primitive-discrepancy sidecar is the
required local datum.

## 8. Tournament Analysis and assumption challenge

The referee uses blocks, not runners, as tournament vertices.  For two blocks
`A,B`, define

```text
debt(A then B) = q(B) M(A),
A -> B iff debt(A then B) <= debt(B then A),
```

breaking equality by increasing physical scale.  On the exact three-block
certificate the tournament is `A -> B -> C` and `A -> C`: score histogram
`{0:1,1:1,2:1}`, zero directed triangles, three singleton SCCs, and one
Hamiltonian path `(A,B,C)`.  It has two edge flips against the density-only
ordering, demonstrating that density does not encode boundary debt.  Exhausting
all six block orders shows `(A,B,C)` uniquely maximizes the theorem bound.

The quotient preserves exactly the density product and the pairwise boundary
debt.  It destroys the within-block endpoint order and actual wall chronology;
those are retained, rather than guessed away, in `q` and the exact endpoint
referee.  Alternate vertices considered were runners, gaps, endpoints, and
proof obligations.  Blocks are the smallest carrier on which both local
certificates and gap suppression coexist.

## 9. Scope

THM-933 proves the sharp local-density composition theorem, improves the
elementary lacunary constants, and closes broad block-separated 13-speed
families.  It is **not by itself a proof of arbitrary LRC(14)**: an arbitrary
packet still needs either one positive whole-block certificate or a partition
whose right side in (BG) is positive.  Its role in the endgame is exact and
modular: every future within-block density certificate becomes composable once
it carries `q`.

## 10. Lean formalization

`TournamentH7.LRCLocalDensityBlockGluing` formalizes the proof-bearing
algebraic core:

- `local_to_component_sum`: a local interval/component certificate sums with
  the exact `card * q` loss;
- `local_to_complexity_sum`: a component cap turns that loss into `M*q`;
- `lowerBound_le_actual`: every nonnegative-density recurrence propagates a
  certified lower bound;
- `weightedDebt_eq_suffix_sum` and `lowerBound_eq_closed`: the recursive ledger
  is exactly the product minus the suffix-density-weighted debt in (BG);
- `three_block_gluing`: the explicit `d1*d2*d3-(e2*d3+e3)` theorem;
- kernel arithmetic for `6^12>7^11`, the positive `R=7` margin, the exact
  `81253/771750` three-block ledger, and its domination by the direct measure.

The module contains no `sorry` and no `native_decide`; its axiom audit reports
only `propext`, `Classical.choice`, and `Quot.sound`.  The remaining geometric
formalization layer is to instantiate the abstract local hypothesis from the
circle primitive `H_B` and to formalize the tooth-union component cap.  Those
are isolated from, and feed directly into, the now-checked closed recurrence.
