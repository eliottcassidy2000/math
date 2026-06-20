---
id: HYP-2661
status: VERIFIED exhaustive (1-tail, 2-tail); the multi-tail mouth-retention rigidity stated as a CARRY CONSERVATION law; proof OPEN (rational comb cutoff)
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

## The mechanism (Glaisher + even-graph reading)

- **Why the tower owns the mouth:** the four tower constraints `||t||, ||2t||, ||4t||, ||8t|| > 1/14`
  are the dyadic refinement of one constraint — each doubling halves the admissible interval (HYP-2656).
  They carve the small surviving "mouth" set (`7/858`). **Deleting `2^a` relaxes `||2^a t|| > 1/14`,
  re-opening a `1/(7·2^a)`-band and ENLARGING the lonely set past the threshold.** Carry conservation:
  the `1111₂` tower is the minimal-measure "ground state"; spending any tower bit pays `≥ 426/35035`.
- **Even-graph parity-completion:** the holes/tails are free binary choices; for sub-threshold they are
  NOT free — the shell-1 binary word must be parity-completed to `1111`. The `{5:+2}` exception is "pure
  carry inside an existing shell" (a gauge move that keeps the tower); the `{1:-4}` move spends a tower
  bit and is forced to pay the second threshold via a new shell (`{23:+2}`).

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
