---
id: HYP-2661
status: VERIFIED exhaustive — AP-tail (1,2-tail) AND general 12-cores [1,19] (50388 cores, 0 violations); the multi-tail mouth-retention rigidity as a CARRY CONSERVATION law; comb proof in progress
source: kind-pasteur-2026-06-19-S16
depends_on:
  - THM-541   # drop-6 collar (mouth = 7/858)
  - THM-543   # one-tail AP layer (the {5:+2} exception)
  - THM-544   # two-tail AP layer (every two-tail >= 426/35035)
  - HYP-2659  # the Euler-parity / odd-shell carry carrier
  - HYP-2654  # multi-tail mouth-damage rigidity (this IS its sharp statement)
related:
  - HYP-2656  # Glaisher dyadic-halving / odd base
  - HYP-2660  # Glaisher-Witt-even-graph bridge
  - OPEN-Q-108
---

# HYP-2661: the dyadic-1 tower {1,2,4,8} OWNS the mouth — a carry conservation law

User hint (made precise): "Euler/Glaisher carry says track `2^a·m` as `2^a` units in odd shell `m`;
tournament/even-graph equinumerosity says free binary choices can be moved into a parity-completed
carrier. In the AP-tail rows, that separates the one-tail exception `{5:+2}` with full mouth survival
from the corrected two-tail minimum `{1:-4,5:+2,23:+2}`, which damages the old mouth but already pays
the second threshold."

## The law (the sharp form of HYP-2654's multi-tail mouth-retention)

For AP-tail 12-cores `C` (`= ({1..13}∖holes) ∪ tails`, `tails ⊆ [14,∞)`), let `meas(G_C)` be the
lonely measure at gap `1/14`. With the odd-shell carry profile `carry(C)[m] = Σ_{2^a m ∈ C} 2^a`:

> **`meas(G_C) < 426/35035` (the AP one-hole second value) ⟹ the full DYADIC-1 TOWER `{1,2,4,8} ⊆ C`,
> i.e. `carry(C)[1] = 15` is intact.** Equivalently, deleting ANY of `{1,2,4,8}` forces
> `meas(G_C) ≥ 426/35035`.

Equivalently: the only sub-`426/35035` AP-tail rows are **`{drop-6 = 7/858, one-tail exception =
7/858 + 1/980}`**, and both keep the complete shell-1 tower `1111₂`.

## Verification (exact, exhaustive)

| family | rows | below 426/35035 | of those, shell-1 damaged |
|--------|------|-----------------|---------------------------|
| 0-tail (drop-`e`) | 13 | 2 (drop-6, drop-12=`426/35035`) | 0 |
| 1-tail (2 holes, tail≤60) | 3666 | 1 (the `{5:+2}` exception) | **0** |
| 2-tail (3 holes, 2 tails≤40) | 100386 | 0 | 0 |
(3-tail in progress.) The carry profiles confirm the hint EXACTLY:
```
  drop-6:            carry={1:15, 3:5, 5:3, 7:1, 9:1, 11:1, 13:1}   meas=7/858        (< thr2)
  one-tail {5:+2}:   carry={1:15, 3:5, 5:5, …}                       meas=3859/420420  (< thr2, shell-1 FULL)
  two-tail min:      carry={1:11, 3:5, 5:5, …, 23:2}                 meas=50189/3223220 (≥ thr2, shell-1 DAMAGED by {1:-4})
```

## GENERALIZES beyond AP-tail to ALL 12-cores (OPEN-Q-108): exhaustive [1,19]

The rigidity is NOT special to AP-tail rows. Over ALL primitive positive 12-cores in `[1,B]` (exact):
```
   B=16: 1820 cores,  1 below thr2 (drop-6), 0 tower-violations
   B=18: 18564 cores, 1 below thr2 (drop-6), 0 tower-violations
   B=19: 50388 cores, 1 below thr2 (drop-6), 0 tower-violations   (= HYP-2651's full atlas box)
```
So **every sub-`426/35035` primitive 12-core in `[1,19]` contains the full dyadic tower `{1,2,4,8}`**, and
drop-6 is the unique sub-threshold core there (matching HYP-2651). This sharpens OPEN-Q-108: the
sub-threshold census is confined to `{1,2,4,8}`-containing cores — a strong universal mouth-owner.

**Structural reason (why exactly `{1,2,4,8}`):** odd part 1 has the DEEPEST dyadic tower among speeds
`≤13` (`{1,2,4,8}`, 4 levels, since `2^4=16>13`); the other odd parts have shorter towers
(`3→{3,6,12}`, `5→{5,10}`, `7,9,11,13→` singletons). By the Glaisher halving (HYP-2656), the deepest
tower gives the tightest NESTED constraints `||t||,||2t||,||4t||,||8t||>1/14`, carving the smallest
("ground-state") lonely set = the mouth. The tower depth `⌊log₂(n−1)⌋ = ⌊log₂13⌋ = 3` is set by `n=14`.

## The mechanism — CORRECTED (workflow owner-geometry, STRONG-PARTIAL): GLOBAL clamping, not local walling

The tower does NOT wall the mouth. The drop-6 mouth's 4 surviving components are walled EXCLUSIVELY by
the three TOP speeds `{11,12,13}` (their danger-arc edges; the THM-541 owners `R(13,2)→L(12,2)` etc.).
The tower `{1,2,4,8}` owns NO mouth wall. So the carry conservation is **GLOBAL, not local**:

- The four nested tower constraints `||t||,||2t||,||4t||,||8t||>1/14` (dyadic refinements, each halving
  the band — HYP-2656) **clamp the REST of `[0,1)` to near-zero safe mass**, leaving only the tiny mouth
  (`7/858`, owned by `{11,12,13}`). **Deleting `2^a` leaves the mouths intact (4/4 survive) but UNCOVERS
  far safe mass** elsewhere (where `||2^a t||>1/14` was the sole binding constraint), and that uncovered
  mass instantly exceeds `426/35035`. So the tower is a global suppressor; spending any bit re-opens a
  far region worth `≥ thr2`. (This is why the bound is robust and the `v=4` deletion — the middle band —
  is the cheapest: it uncovers the least far mass.)
- **Even-graph parity-completion:** the holes/tails are free binary choices; for sub-threshold they are
  NOT free — the shell-1 word must be parity-completed to `1111` to maintain the global clamp. The
  `{5:+2}` exception is "pure carry inside an existing shell" (a gauge move that keeps the clamp); the
  `{1:-4}` move spends a clamp bit and pays `≥ thr2` via the uncovered far mass (+ a new shell `{23:+2}`).
- **Even-graph parity-completion:** the holes/tails are free binary choices; for sub-threshold they are
  NOT free — the shell-1 binary word must be parity-completed to `1111`. The `{5:+2}` exception is "pure
  carry inside an existing shell" (a gauge move that keeps the tower); the `{1:-4}` move spends a tower
  bit and is forced to pay the second threshold via a new shell (`{23:+2}`).

## Tower-deletion binding minima (exact) — the quantitative carry conservation

The minimum `meas(G_C)` over AP-tail 12-cores MISSING each tower bit `v` (the "spending" cost):
```
   v=1 missing:  min = 1333/30940  = 0.04308  (margin to thr2: 0.0309)   at holes(1,6,10) tails(17,20)
   v=2 missing:  min = 27/1547     = 0.01745  (margin 0.0053)            at holes(2,6) tail(17)
   v=4 missing:  min = 335/23023   = 0.01455  (margin 0.0024, TIGHTEST)  at holes(4,6) tail(46)
   v=8 missing:  min = 6163/336336 = 0.01832  (margin 0.0062)            at holes(6,8) tail(16)
```
ALL ≥ `426/35035`; `v=4` is the binding deletion (the middle tower bit, closest to the threshold). The
shell-1-full cores `drop-12 = 426/35035` (the threshold itself) and `drop-6 = 7/858` (the floor) achieve
the sub-threshold values, confirming the tower is the ground state and any deletion pays `≥ thr2` with
the `v=4` deletion the cheapest. (1-tail & 2-tail exhaustive; the `v=4` finite minimum `335/23023` is the
target of the comb proof below.)

## Proof target (OPEN)

Prove the quantitative step **"`2^a ∉ C` (a∈{0,1,2,3}) ⟹ `meas(G_C) ≥ 426/35035`"** via the rational
periodic-comb cutoff (the THM-543/544 technique): the missing `||2^a t||` band gives a comb lower bound
`meas(G_C) ≥ (1 − 1/(7·2^a))·meas(G_{C+2^a}) + …`, and a finite exact residue check. This would close the
multi-tail mouth-retention rigidity (HYP-2654) as a single carry-conservation lemma, reducing the
AP-tail sub-threshold census to the two known templates. Then route to the general 12-core via the
Glaisher tower word + CRT cell + endpoint-owner ledger (HYP-2660).

## Honest status
LRC(14) NOT proved. This sharpens HYP-2654 to a clean carry-conservation law (the dyadic-1 tower owns the
mouth), VERIFIED exhaustive for ≤2 tails, and identifies the proof as a rational-comb cutoff on the
missing tower band. Files: `04-computation/lrc14_{carry_conservation,shell1_rigidity}_kps.py`.
