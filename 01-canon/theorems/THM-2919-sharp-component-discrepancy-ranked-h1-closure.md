---
id: THM-2919
title: "Sharp component-discrepancy ranked-H1 closure"
status: >
  PROVED + FINITE-EXACT + VERIFIED.  The componentwise THM-1094
  discrepancy constant brings 80 additional ranked-H1 branches under the
  15,000 cutoff; 79 close by the ordered literal-child test and the last
  by one exact maximum-tree Hunter repair.  One route root results, already
  covered by THM-2916.  A route-saturation audit isolates the unique
  still-additive high-cutoff root, whose missing branch at N=18,869 closes
  by the same ordered test.  Relative to the live 1,041-root union through
  THM-2920 this adds exactly one root, giving union 1,042 and residual
  2,390.  This is not LRC(14).
source: codex-sharp-h1-tail-2026-07-29
depends_on:
  - THM-1094-exact-two-comb-component-theorem
  - THM-2907-paircap-exception-h4-heavy-link-child-census
  - THM-2911-all-root-finite-ranked-h1-hunter-closure
  - THM-2912-one-h3-row-exact-child-top-four-closure
  - THM-2913-one-h3-row-pair-hunter-toothpick-closure
  - THM-2916-two-h3-row-dynamic-tail-child-top-four-closure
  - THM-2920-two-h3-row-pair-hunter-recursive-toothpick-closure
related:
  - THM-2915-all-open-centre-exact-child-top-four-closure
verification:
  - 04-computation/lrc14_j6_sharp_h1_tail_recomposition_codex_20260729.py
  - 05-knowledge/results/lrc14_j6_sharp_h1_tail_recomposition_codex_20260729.out
---

# THM-2919 -- sharp component-discrepancy ranked-H1 closure

**PROVED + FINITE-EXACT + VERIFIED.**

## 1. Statement

In the exact `14,806`-row scalar-hard branch universe used by THM-2911,
exactly `6,180` rows satisfy the strict ranked-four gate

```text
R_4=q_1+q_2+q_3+q_4<6h/7.
```

THM-2911 certified all `5,999` such rows whose old strict discrepancy
cutoff was at most `15,000`.  Replacing only that tail by the sharper
componentwise consequence of THM-1094 puts `6,079` rows under the same
cutoff.  Every one of the `80` newly finite rows closes: `79` by the
ordered literal-child certificate and one by an exact relative-Hunter
maximum-tree repair.

Their branchwise composition with the proved `G_5`, hostile-centre pivot,
and H4-exception routes closes the additional route root

```text
(2,4,8,10,11,13,14).                                    (1)
```

Root `(1)` is independently contained in the proved THM-2916 bank, so it
does not enlarge the live union.  Exact saturation of the sharp route
then isolates one still-additive missing branch.  Its direct exact replay
closes the root

```text
(1,3,5,7,10,11,14).                                     (2)
```

The proved union through THM-2920 has `1,041` roots and does not contain
`(2)`.  Therefore this theorem gives

```text
proved union                                      1,041+1=1,042
finite residual                         3,432-1,042=2,390.            (3)
```

This is a finite exact improvement, not a proof of LRC(14).

## 2. The sharp component tail and its equality layer

Let `C` be a literal carrier of mass `h` with interval components
`I_1,...,I_r`, and put

```text
c_C(w)=mu(C intersect D_w).
```

Section 2 of THM-1094 proves for every interval `I`

```text
mu(I intersect D_w)<=mu(I)/7+6/(49w).                   (4)
```

Endpoints have measure zero.  Summing `(4)` over the disjoint components
of `C` gives

```text
c_C(w)<=h/7+6r/(49w)
       =h/7+(6/7)r/(7w).                                (5)
```

For a ranked-four row set

```text
epsilon=6h/7-R_4>0,
H_1={w not in P:c_C(w)>=h-R_4}= {w:c_C(w)>=h/7+epsilon}.
```

Combining this membership inequality with `(5)` shows

```text
w<=(6/7)r/(7epsilon).                                   (6)
```

Because `(5)` is non-strict, the safe integer cutoff is

```text
N_sharp=floor((6/7)r/(7epsilon)),                        (7)
```

not `ceil(...)-1` at an integral ratio.  The verifier finds exactly three
integral sharp ratios:

```text
3718, 154440, 1536.
```

The first and third rows were already in the old finite H1 bank, while the
second remains above `15,000`.  None of the `80` new rows is an equality
case, so retaining the correct boundary leaves the exact `6,079/80`
census unchanged.

## 3. The 80 newly finite branches

The new cutoffs range from `9,107` to `14,970`; their old cutoffs range
from `15,028` to `24,701`.  The literal H1 cores on the `79` direct
closures have sizes `315` through `755`.  For every ordered earliest
label `x`, the computation subtracts `D_x` from the literal carrier and
checks that the exact four largest later-label coverages sum strictly
below the child mass.  The smallest positive direct margin is

```text
167/1,891,890.                                           (8)
```

Only one branch has an unresolved earliest pivot:

```text
E=(3,6,7,9,11,12,13), K=13, rank=3,
a=19, P=(15,30,19), N=10842, |H_1|=525, pivot=25.        (9)
```

Star pruning leaves exactly one hostile five-set,

```text
Q=(16,25,31,32,64).
```

For any five-set `Q`, Hunter's tree inequality gives

```text
mu(union_(w in Q)(C intersect D_w))
 <=sum_(w in Q)c_C(w)-max_T sum_(xy in T)i_C(x,y),       (10)
```

where `T` ranges over spanning trees and
`i_C(x,y)=mu(C intersect D_x intersect D_y)`.  Exact Kruskal
maximization on the sole hostile set gives

```text
Psi_5(Q)=7,543,517/37,519,300,
h-Psi_5(Q)=1,932,131/28,139,475>0.                       (11)
```

Thus `(9)` also closes, and all `80` newly finite branch keys are proved.

## 4. The route-selected extension

The fixed-window branch union closes `139` route roots, one more than the
old `138`.  To decide whether a larger blind scan was useful, the verifier
performs a set-theoretic saturation audit: pretend, only for target
selection, that all `6,180` positive-epsilon sharp-H1 keys were available.
This gives `146` route roots.

Seven roots lie in the saturated route but not the fixed-window route.
Five were already in the historical `351`-root union through THM-2913.
Exactly two could have enlarged that historical union:

```text
(1,2,3,4,5,6,14),       missing sharp cutoff 34222,
(1,3,5,7,10,11,14),     missing sharp cutoff 18869.     (12)
```

The first root in `(12)` is independently proved by THM-2916.  The second
is not in the live union through THM-2920, so it is the unique
load-bearing high-cutoff target.  Its missing branch is

```text
K=9, rank=1, a=25, P=(25), N=18869, |H_1|=831.           (13)
```

The exact ordered child test on `(13)` has

```text
824 parent-bound pivots,
3 literal-local pivots,
4 short suffixes,
0 unresolved pivots,
```

with smallest strict margin

```text
89,371/111,411,300.                                      (14)
```

The other two hard rows on this body have old finite-H1 cutoffs `482` and
`821`, so `(13)` closes the whole root `(2)`.

The `101` sharp-positive keys above the fixed cutoff contain only the two
keys in `(12)` that could add a root beyond the old `351`-root union.
Thus the other `99` high-cutoff keys cannot contribute a new whole root
there.  This is a target-selection ceiling, not a claim that unscanned
branches themselves are closed.

## 5. Exact recomposition

Keep the full branch key

```text
(body, gate size, rank, apex, excluded prefix).
```

Writing `E,G,H,P,S` for the H4 exception, G5, old finite-H1,
hostile-centre pivot, and new sharp-H1 certificates, respectively, the
nonempty atomic branch classes after the extension are

```text
E 52, G 49, GH 2909, GS 6, H 2875,
HP 215, P 60, PS 4, S 71.                               (15)
```

The old branch union has `6,170` keys and `138` route roots.  The
fixed-window sharp union has `6,240` keys and `139` route roots.  Adding
`(13)` gives `6,241` keys and `140` route roots.  Relative to the
historical union through THM-2913, the two route additions are exactly
`(1)` and `(2)`.  THM-2916 already contains `(1)`, and neither THM-2916
nor THM-2920 contains `(2)`.  Exact set difference against the current
`1,041`-root tuple therefore leaves only `(2)`, proving `(3)`.

The sorted `80` fixed-window key tuple has digest

```text
a6fe73463c89a889137e0cf25c6d1badbb792654862676c3186c5c9adc12018f,
```

and the `81`-key tuple including `(13)` has digest

```text
3f6280598288370f5acd7ce51e7439faa730edc0af46980a65da6095d8c25b57.
```

## 6. Verification

The verifier hash-pins THM-1094, the THM-2911 scout and join, the
relative-Hunter engine, and the exact outputs through THM-2920.  It
reconstructs all branch keys, carriers, exclusions, equality cutoffs, and
root sets from those sources.  All interval arithmetic and decisions use
exact rational arithmetic.

Ordinary and optimized eight-worker replays produce byte-identical output
and pass every locked control.  The internal `82`-row certificate ledger
has SHA-256

```text
6de74387474e14bbfd622c5371d2d246588532e29764938f6210277478950831.
```

LF-normalized SHA-256 values are

```text
source  a628808c6a7a22f94c0b387d45c6e8066858f0cfdd77392cecb6db40744f59df
output  8eb5b8af65466539f070418b1d405b6b8b7c8d7f9595c2a7fa069196f80c0d38.
```

Reproduction:

```bash
python3 04-computation/lrc14_j6_sharp_h1_tail_recomposition_codex_20260729.py --workers 8
python3 -O 04-computation/lrc14_j6_sharp_h1_tail_recomposition_codex_20260729.py --workers 8
```
