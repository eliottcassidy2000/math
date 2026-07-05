# The l ≥ 7 Assembly Lemma: 13-adic tower induction (HYP-4119)

**opus-2026-07-05-S82.** The corrected architecture for the l ≥ 7 leg of hdich lift
rigidity (post MISTAKE-105), assembled from: the shadow dichotomy at level 169 (this
session's search), gap descent (HYP-4116, LRCGapDescent.lean), the class-level kernel
rows (LRCLiftRowsL7.lean), the S80 structural lemmas (LRCLiftPigeonhole.lean), kps's
pair walk (HYP-4117) and compression ladder, and the six-top ceiling caps.

## Setting

W = {v_1, ..., v_12} primitive (gcd 1), residues mod 13 pinned to {1..12} (residue
pinning, S75), position r holds v_r = r + 13 k_r, lift count l = #{r : k_r >= 1} >= 7.
Goal: M(W) > 1/13 strictly (the only primitive tight family is {1..12}, which has l = 0;
dilated standards λ·{1..12} — which DO reach l >= 7 for λ >= 3 — are non-primitive and
excluded upstream by the primitivity split).

## The three moves

**(GRID / witness).** At t = a/169 the margin of runner v depends only on v mod 169,
i.e. on the pattern κ_r = k_r mod 13. If the pattern's kill sets do not cover all 156
cells a ∈ [1,168] \ 13ℤ, an uncovered a gives min_r dist(v_r a/169) >= 14/169 > 1/13:
strictly loose, certified by ONE kernel row for the whole congruence class mod 169
(strictLonely13_of_row).

**(SHADOW / ascend) -- REVISED after the census.** The dichotomy CONJECTURE (full covers
= dilation shadows) is REFUTED: the exhaustive DFS finds tens of thousands of cover
patterns (all sampled cores irredundant on ALL 12 positions, none shadow, none affine or
quadratic in r); the cover variety is large and unstructured-so-far. Consequently
per-pattern kernel rows CANNOT enumerate the stratum. What survives: (i) every cover-class
member other than its small representative (values <= 168) has a runner >= 170 =
descent-eligible; (ii) sampled small reps are comfortably loose (M = 7/23, 4/17, 2/13 --
all >= 2/13), witnessed at merge denominators, NOT at 169 -- covers are "loose elsewhere"
classes, not near-tight ones. The tower ascent survives for the SHADOW covers specifically
(primitivity forces deviation); the non-shadow covers are handled by the strata below
WITHOUT pattern classification.

**(DESCENT / dodge).** Runners >= 170 > 134 = 2/L clear the entry bar of gap descent
(L = 7/468 from citing LRC(<= 6) on the <= 5 unlifted runners, margin 1/6, Lipschitz
B = 12). Spread tops (consecutive ratios >= 26/11) are dodged in ANY number
(spread_tower). Between descent stages, the six-top ceiling caps what fees can absorb
(<= 6 per window — provably necessary, HYP-4116); same-scale pairs are killed by kps's
pair walk (min <= 22B); same-scale clusters of >= 3 are the alignment stratum (below).

## The case tree (corrected strata, constants at band 1/13)

Split by v_max, not by pattern:

1. **BOUNDED stratum (all values <= 168, i.e. every k_r <= 12):** finite: <= C(12,7+)
   choose lifted positions x k-ranges; with the pinning/pigeonhole/primitivity prunes this
   is a mac-mini-scale C sweep (~1e10 cells against their census harness; their S54 rig
   did 3.6e9 nodes in 12 min). Verdict expected: all loose with witnesses at merge
   denominators (the sampled small reps all had M >= 2/13). Status: COMP target,
   HANDED to mac-mini (their DFS-prune offer, S54 letter).
2. **UNBOUNDED, spread top (v_max > 12 * v_2nd):** the S76 window peel removes the top
   (threshold 12B at citation margin 1/12); recurse on the 11-family (induction on count;
   strictness via interior points, pending pass).
3. **UNBOUNDED, compressed top (169 <= v_max <= 12 * v_2nd):** the top block lives at
   scale >= 15 with bounded ratios. Sub-split: (a) top block spread within itself at
   ratios >= 26/11: gap descent (any count); (b) top block a same-scale PAIR: kps pair
   walk (min <= 22B, band analogue pending at 1/13); (c) top block a same-scale >= 3
   cluster inside ratio < 26/11: THE ALIGNMENT SLIVER -- the one genuinely open case,
   sharpened to: >= 3 runners with distinct residues mod 13, values within ratio 26/11 at
   scale >= 169, plus a bounded core. This is mac-mini's alignment-band CRT geometry
   (HYP-4094) territory; fees absorb such a block only when its scale >= ~1605
   (boundary terms 4rho/W < slack L/13).

## The l >= 7 specifics

Pigeonhole (S80, formal): >= 7 lifts hit a unique-multiple residue r ∈ {7..12} whose
sieve survival forces k_r >= r, value >= 14r >= 98. Note 98 < 134: height forcing alone
does NOT clear the descent bar — the gap between 98 and 134 is why the grid stratum is
load-bearing (a fact S78 missed). With the dichotomy, the assembly needs no per-l
case split beyond visible >= 7: the grid/shadow/descent trichotomy is uniform in l.

## Status ledger

- exists_gap_subinterval / spread_tower / spread_tower_13: GREEN (S81).
- speedOK13_congr / strictLonely13_of_row + rows A–D: GREEN (this session).
- Shadow dichotomy at 169: REFUTED this session (covers >> shadows; census in
  shadow_dichotomy_169_opus_S82.out); the machinery (rows, congruence transport,
  descent) survives via the v_max strata instead.
- OPEN (i): the BOUNDED stratum sweep (mac-mini harness; COMP status).
- OPEN (ii): the alignment sliver -- compressed >= 3-blocks, scale [169, ~1605].
- OPEN (iii): strictness bookkeeping through descent (interior points; mechanical).
- OPEN (iv): the ascent/termination formalization (Int induction on levels; small).

## The l <= 6 caps are theorems, not choices (for the proof map)

The six-top ceiling (HYP-4116) upgrades every l <= 6 cap in the corpus from design
choice to necessity: any placement-uniform fee dominates its placement mean 2ρL, so
Σ over l tops >= 2ρlL >= L exactly when 2ρl >= 1 — at ρ = 1/13, 1/14, AND 2/25 the wall
is l = 6.5. kps's C_l ladder stopping at l = 6, mac-mini's T_l = 156(13−l)/(13−2l)
poles, and klein's 6-top multi-far are all the same theorem seen from three sides.
No sharpened mass bound crosses it; only adaptivity (descent) or arithmetic (grid) does.


## S84 addendum: the composition lemma and the l >= 9 leg (HYP-4136)

**THE COMPOSITION LEMMA (descent with fee-paying passengers; verified 6/6,
composition_descent_opus_S84.out):** in a window of length L, remove the <= 1 tooth each
of p SUB-GAP passengers (values with 13*L*w <= 11; teeth split the window into <= p+1
pieces): the largest passenger-free piece has length >= (L - Sum_p 2rho/w_p)/(p+1), and
the gap-descent chain runs inside it (passengers stay sub-gap as gaps shrink -- the
sub-gap condition only improves). One point ends up rho-far for base + passengers + all
spread tops.

**THE l >= 9 BUDGET TABLE** (base = <= 3 unlifted, margin 1/4, B <= 12, L0 = 3/104):
p = 0: entry bar 69 | p = 1: bar 224 | p = 2: bar 789 | p >= 3: DEAD (worst-case fees
14,15,16 exceed L0). So at l >= 9 the descent leg tolerates AT MOST TWO small passengers.

**THE l >= 9 PARTITION:** lifted values split (<= 2 passengers in [14,29]) + (spread
tops >= bar(p), ratios >= 26/11). RESIDUAL: >= 3 values <= 29, or any value in the mid
zone [30, bar(p)), or ratio-clusters above the bar. All residual pieces are BOUNDED
below ~800 except clusters, which go to the level-3 arithmetic (HYP-4126 probes) --
the bounded residual is a mac-mini-prune-stack sweep, NOT naive (k-ranges to ~60).

**ANY-CARRIER (MISTAKE-109):** the >= 98 height forcing was self-carried-only; the
carrier table (min value r*(s*r^{-1} mod 13), minima 14..24 for r = 7..12) replaces it.
The l >= 9 box constraints come from covering-duty matching, not per-position forcing.
