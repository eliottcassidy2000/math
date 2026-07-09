# THM-672 — The descent-window torsion-occupancy theorem: on k ∈ [15,28], blocking requires torsion classes; composite all-unit descents are never blocked; wall primes block iff every ±-class is hit

**Status:** PROVED (master ledger + T1 + T2 + T3, elementary, proofs below) + VERIFIED
(exhaustive: T2 exact characterization 0 violations at k = 17/19/23 over all 13-subsets incl.
497,420 at k=23; T1/T3/ledger 0 violations at every composite k — subsets exhaustive k ≤ 26,
300k sampled at k = 27, 28, plus random multisets; the k=28 row was corrected once by this very
verification, then 0/12,709). **PLUS the decisive negative: the [15,28] window alone CANNOT
close the large-domain statement — full window-dodgers exist among covering sets (dodger search,
7/8 restarts) precisely because covering supplies the torsion occupants; every dodger found is
caught at k ≥ 29.** The large-domain statement is thereby localized to the k ≥ 29 descent
structure, not proved.
**Source:** mac-mini-2026-07-09-S65 (cont. 2). Targets HYP-5730 obligation (2): the a-priori
"covering + sums > Q₀ ⟹ C1/C2 fires".
**Depends on:** THM-668 (pair-sum ruler theorem; certificate C2 = divisor descent).

**Setup.** Descent modulus `k` with `15 ≤ k ≤ 28`, so `⌈k/14⌉ = 2`: the closed band is
`[2, k−2]` and the danger set is `{0, ±1}` mod `k`. Residues `R = {v_l mod k}` (multiset,
13 entries). If `0 ∈ R` the descent is dead. "Blocked" = no multiplier `s ∈ [1, k−1]` puts all
of `R·s` in the band. C2 fires on a pair sum `q` through `k | q` iff `R` is NOT blocked (then
`p = (q/k)s` is a banded multiplier mod `q` — THM-668 addendum).

## Lemma 1 (bad-set structure)

For `r ∈ R` with `g = gcd(r, k)`, let `A_r = {s : r·s mod k ∈ {0, ±1}}`.
- `g = 1`: `A_r = {0, r⁻¹, −r⁻¹}` — exactly 2 nonzero elements; `A_r = A_{r'} ⟺ r' ≡ ±r`.
- `g > 1`: `r·s ≡ ±1` is unsolvable (`g | rs`, `g ∤ 1`), and `r·s ≡ 0 ⟺ (k/g) | s`; so
  `A_r = (k/g)·Z/k` — exactly `g − 1` nonzero elements — **determined by `g` alone**, and nested:
  `g | g' ⟹ A_{(g)} ⊆ A_{(g')}`.

*Proof.* Direct computation; for units, inversion. ∎

## Theorem (master ledger — necessity)

> Blocked ⟹ `2·U + N ≥ k − 1`, where `U` = number of distinct occupied unit ±-classes and
> `N = |⋃_{occupied g>1} (k/g)·Z/k \ {0}|`.

*Proof.* Blocked means `[1, k−1] ⊆ ⋃_r A_r`; count the right side by Lemma 1: unit classes
contribute ≤ 2 each (distinct classes only), non-units contribute their union `N`. ∎

## T1 (unit-pigeonhole: composite all-unit descents are NEVER blocked)

> If `k` is composite and every residue is a unit, the band is solvable.

*Proof.* `N = 0`, `U ≤ φ(k)/2`, so the ledger is `≤ φ(k)`. For composite `k`,
`φ(k) ≤ k − √k < k − 1`. By the master ledger, not blocked. ∎
(For prime `k`, `φ(k) = k−1` — exactly the wall; hence T2.)

## T2 (wall primes — exact characterization)

> For prime `k` (here `k ∈ {17, 19, 23}`), with no zero residue:
> **blocked ⟺ `R` hits EVERY ±-class mod `k`.**

*Proof.* (⟹) Master ledger: `2U ≥ k−1 ⟹ U ≥ (k−1)/2` = all classes.
(⟸) Inversion permutes ±-classes bijectively, so `{±r⁻¹ : all classes}` = all of `[1, k−1]`:
every `s` is bad. ∎

## T3 (per-k torsion-occupancy necessary conditions, composite k)

From the master ledger, blocked ⟹ `N ≥ d(k) := k − 1 − φ(k)` (using `U ≤ φ(k)/2`). The nested
non-unit unions then force, mechanically (derivation per k; machine-verified against exhaustive
blocked enumerations):

| k | d(k) | blocked ⟹ (torsion residues occupied) |
|---|---|---|
| 15 | 6 | some `v ≡ ±3, ±6` AND some `v ≡ ±5` (mod 15) |
| 16 | 7 | some `v ≡ 8` (mod 16) |
| 18 | 11 | some `v ≡ 9` AND some `v ≡ ±6` (mod 18) |
| 20 | 11 | some `v ≡ 10` AND some `v ≡ ±4, ±8` (mod 20) |
| 21 | 8 | some `v ≡ ±7` AND some `v ≡ ±3, ±6, ±9` (mod 21) |
| 22 | 11 | some `v ≡ 11` AND some even `v ≢ 0 (mod 11)` (mod 22) |
| 24 | 15 | some `v ≡ 12` AND some `v ≡ ±8` (mod 24) |
| 25 | 4 | some `v ≡ ±5, ±10` (mod 25) |
| 26 | 13 | some `v ≡ 13` AND some even `v ≢ 0 (mod 13)` (mod 26) |
| 27 | 8 | some `v ≡ ±9` (mod 27) |
| 28 | 15 | some `v ≡ 14` AND some `v ≡ ±4, ±8, ±12` (mod 28) |

**The occupied classes are the small-torsion points of `Z/k`** (`k/2`; `±k/3`; `±k/5`): blocking
requires the residue system to sit on torsion — the arithmetic signature of AP-likeness.
*(Provenance note: the k=28 row initially conflated the g=7 residues `{±7}` with the g=4 class's
A-set `7Z`; the exhaustive verification flagged 5299 violations, the row was corrected to the g=4
RESIDUES `±4, ±8, ±12`, and re-verification gives 0 violations on 12,709 blocked configs. All
other rows verified with 0 violations on their exhaustive enumerations — the
derive-then-machine-verify loop working as intended.)*

## Assembled conditional theorem (what IS proved toward the large-domain statement)

> **Let S be any 13-set. If some pair sum `q = v_i + v_j` has a divisor `k ∈ [15, 28]` such that
> (i) no `v_l ≡ 0 (mod k)`, and (ii) the T2/T3 blocking-necessary condition FAILS at `k`
> (composite: a required torsion class is empty; prime: some ±-class is missed), then C2 fires
> on `q` and `S` is lonely: `M(S) ≥ 1/14`, witness `t = (q/k)s/q`.**

No covering, primitivity, or scale hypotheses. This converts the large-domain statement into a
pure residue-occupancy question.

## What this does and does not prove (honest scope)

**The covering condition naturally SUPPLIES the torsion occupants.** A covering set has multiples
of 2..14; an odd multiple of 13 is `≡ 13 (mod 26)`; a multiple of 7 prime to 3 is `≡ ±7
(mod 21)`; a multiple of 9 that is odd is `≡ 9 (mod 18)`; etc. So covering and window-blocking
are COMPATIBLE — the tight AP is the proof of concept (it occupies every torsion class ≤ its
range). **DODGER SEARCH RESULT (decisive, honest): the window alone does NOT close the large
domain.** Hill-climbs at cap 120 found full window-dodgers in 7 of 8 restarts — covering sets
with EVERY `k ∈ [15,28]` descent blocked or dead (e.g. `{7,22,27,28,31,46,55,60,61,91,100,115,
120}`: 39/39 window descents blocked; `91 ≡ 13 (26)`; `7,28,91 ≡ ±7 (21)`; all 8 ±-classes mod
17 hit). Exactly as the alignment predicts: covering hands the adversary the occupants. **But
every dodger found is immediately caught above the window** — full C2 fired at `k = 29` (`q =
58`), and the dodgers carry many live rulers (`29, 34, 35, 38, 44, 49`). So the counting window
ends at `k = 28` (for `k ∈ [29,42]` the danger is `{0,±1,±2}`, unit cost 4, and `4·φ(k)/2 <
k−1` never holds — no T1 analog), while the `k ≥ 29` descents remain empirically decisive
through finer structure. Example of what that structure is: for prime `k = 29`, `2` is a
primitive root and `−1 = 2^14`, so ±-classes ≅ `Z/14` with the doubling map acting as `+1`;
blocked ⟺ the occupied inverse-class set `C ⊆ Z/14` satisfies `C ∪ (C+1) = Z/14` (no two
consecutive missing) — a cycle-domination condition. Each higher `k` has such an exact orbit
characterization; there is no uniform counting theorem beyond 28.

**Net effect on the large-domain statement ("covering + sums > Q₀ ⟹ C1/C2 fires"): it is now
PRECISELY LOCALIZED, not proved.** Proved here: the exact anatomy of window-blocking (master
ledger, T1, T2, T3) + the conditional theorem below it. Refuted here: any proof through the
`[15,28]` window alone. Remaining: the `k ≥ 29` descent structure (orbit/domination
combinatorics per k, empirically unbeaten) — equivalently, klein's mid-band realization in
adaptive-split coordinates. The theorem's value:
1. It makes the dodge requirements EXACT — a full window-dodger must simultaneously (a) occupy
   the torsion classes mod every `k ∈ [15,28]` dividing any of its ~91 pair sums (or keep those
   `k` from dividing any sum), (b) hit all 8/9/11 ±-classes mod each of 17/19/23 that divides a
   sum (11 of 13 residues committed for 23 alone), and (c) still evade C0/C1 and the `k > 28`
   exact descents. Each requirement consumes degrees of freedom of only 13 speeds.
2. It gives the Lean-shaped certificate: conditions (i)+(ii) are finite residue checks; the
   consumer is kps-S114's `mreach_ge_of_pairsum_band`.
3. It explains the blocking census structurally (blocked configs are torsion-anchored /
   AP-like), and why `Q₀ ≈ 36`: the cost-2 window ends at 28; `[29, 42]` costs 4 and the
   unit-pigeonhole fails there (`φ(k)/2 ≥ 13` becomes possible).

**Verification:** `04-computation/lrc14_torsion_occupancy_macmini_S65cont2.py` — T1/T2/T3
verified against exhaustive blocked enumerations per k (subsets exhaustive `k ≤ 26`, sampled
300k at `k = 27, 28`, plus random multisets), and the dodger search (see `.out`).

**Related:** THM-668 (+ addendum certificates), THM-420 (Lemma B = the `k = 2n−1` unit case of
this window), klein's adaptive-split THM-667 (mid-band localization — the dodgers' habitat),
monad THM-667/668 (grid-port; detuned dispatch), opus LEM-015 / kps LRCSchurRigidity (E3
rigidity — the global counterpart of per-k torsion anchoring), HYP-5730.
