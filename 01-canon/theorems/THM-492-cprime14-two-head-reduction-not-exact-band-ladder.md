# THM-492 — C′(14) does NOT reduce exactly to C′(5)-on-the-3-core ∪ the mod-2/mod-7 fiber; the band ladder and the stranger dodge

**Status:** PROVED (degeneration lemma; band criterion; fibered-shell reduction; the family
theorem, part-finite-verification) + VERIFIED (exact integer/ℚ arithmetic throughout; 688 joint
failures; the 936-instance family exhaustive). **Source:** claudebox-2026-06-11-S7, executing
t-0122 / HYP-2436 Test 3 (user dispatch: "does C′(14) reduce exactly to C′(5)-on-the-3-core ∪ the
THM-421 mod-2/mod-7 fiber? attack LRC 14 creatively").
**Companions:** THM-491 (ramification tower), THM-421 (divisor-clock peel), THM-420 (k-clock /
transversal), THM-398 (C′ reduction + dominance dodge), S622-era twisted-shell dodge, S643 fiber
windows. **Normalization:** canon (n=14 runners, 13 speeds, M(S)=max_t min_v‖vt‖; C′(14):
multiple-of-14 ⟹ M > 1/14).

## 0. Primitivity hygiene (a statement-level correction to C′)

C′ as literally stated in THM-398 §0 (no gcd condition) is **false by dilation**:
`S = 2·{1,…,13} = {2,4,…,26}` contains 14 and has `M(S) = M({1,…,13}) = 1/14` exactly (tight) —
verified exact, witness structure at t = 1/28. Since M is dilation-invariant, C′ must be read for
**primitive** sets (`gcd(S) = 1`), which is what every verification (THM-420's "gcd-1 configs",
HYP-2102's boxes) already used. The reduction C′ ⟹ LRC is unaffected: apply Lemma A / C′ to the
primitive representative. All sets below are primitive with a multiple of 14.

## 1. The degeneration lemma — the answer is NO, and why

**Lemma (proved).** Let `3 | v`, `27 ∤ v`, and `gcd(a, 27) = 1`. Then
`‖v a/27‖ = ‖(v/3)·a/9‖ ≥ 1/9 > 1/14`.
*Proof.* Write v = 3w with 9 ∤ wa (a is a unit, w ≢ 0 mod 9 since 27 ∤ v); a nonzero residue mod 9
has clock distance ≥ 1/9. ∎ (Exhaustively confirmed over all residues: 0 violations.)

**Consequence.** At level 1/14, the 3-core of the shell-27 problem carries **only the divisor
conditions** {3|v, 9|v, 27|v}: the THM-491 descent 27→9 is arithmetically exact, but the
**forbidden band does not rescale**. The descended shell-9 problem has band {0}; C′(5) proper
(level 1/5 on its shell 9 = 2·5−1) has band {0, ±1}. The bands agree iff `1/n ∈ [1/9, 2/9)`,
i.e. n ∈ {5,…,9}. At n = 14 the 3-core head **degenerates to divisibility — C′(5) never appears.**

**Non-necessity witness (verified exact).** `S = {1,…,12, 28}`: its 3-core ÷3 = {1,2,3,4} is the
**tight** n=5 AP (M = 1/5, blocked at C′(5)'s own band), yet S is loose: M(S) = 1/13, with a
shell-27 unit twist t = 2/27 certifying (min ‖vt‖ = 2/27 > 1/14). So "C′(5) on the 3-core" is not
necessary for C′(14)-looseness.

## 2. The band criterion and the fibered-shell reduction (the corrected heads)

**Band criterion (proved; generalizes the S622 twisted-shell dodge beyond band ±1).**
`t = a/q`, `gcd(a,q) = 1`, is a strict witness (M > 1/14) **iff** every v ∈ S has
`(va mod q) ∉ ±{0, 1, …, ⌊q/14⌋}`. (3000-sample exact cross-check: 0 mismatches.)
Band-k shells are `q ∈ [14k, 14(k+1)−1]`; the band-1 ceiling is `2n−1 = 27`, the band-2 ceiling
`3n−1 = 41`, generally `(k+1)n−1`: **the shell horizon is a ladder, not a wall.**

**Fibered-shell reduction (proved).** At `q = d·m` with `d | 14`: a d-core runner `v = dw` sees
`‖va/q‖ = ‖wa/m‖` — the core/d enters its **own shell-m problem** with band `⌊m/14⌋`, which for
`m ≤ 13` is `{0}`, a **pure divisor condition** (no multiple of dm). Strangers see band
`⌊dm/14⌋`. The THM-421 clocks (m = 1), the S643 perturbation windows (m → ∞), and the 3-tower
(d = 1, m ∈ {3, 9, 27}) are all faces of the one lattice `Q = {d·m : d | 14, m ≤ 27}`
(76 moduli). Example: at q = 91 = 7·13, a 7-core needs only "no multiple of 91", i.e. the core/7
is handled by the **proven** LRC(13) clock; the stranger rides band ±6.

## 3. Non-sufficiency: the joint failures and the family theorem

**Heads as asked.** HEAD A = THM-421/S643 fiber (d ∈ {2,7,14}, all base clocks b/d, exact R-safe
window, exact core-arc coverage inside it). HEAD B = the 3-tower (shell-27 unit twists + shells
9, 3). Implemented exactly; **688 primitive multiple-of-14 configs found failing BOTH heads**
across three families — all 688 loose (C′ holds on them), with witnesses outside both heads:
minimal witness moduli {4,…,25} ∪ {40, 41}; the fibered lattice Q certifies **688/688**.
So the union [C′(5)-3-core head] ∪ [THM-421 fiber head] is **not sufficient** either. NO is the
answer to the exactness question in both directions.

**The family theorem (proved).** `S(r) = 7·{1,…,12} ∪ {r}` (r ∤≡ 0 mod 7, distinct): every valid
instance is loose.
- `r > 1092 = 13·84 = (n−1)·V′`: THM-398 Cor B2 — the dominance dodge applied to the **stranger**
  r (no divisibility hypothesis needed).
- `r ≤ 1092`: finite; exhaustively verified exact (936 instances; max minimal witness modulus 41).

**The five evaders (the sharp new fact).** r ∈ {611, 702, 793, 962, 1053} — all ≡ 0 (mod 13),
with `r mod 27 ∈ {0, ±10}` — additionally block **all twisted shells m ≤ 27** and the **B′ width
dodge of every multiple of 14**. They are the first explicit configs past the 2n−1 horizon: their
minimal witnesses sit at `q ∈ {40, 41}`, exactly one band rung up (band-2 range [28, 3n−1=41]).
Mechanism: 13 | r kills the 13-clock (the family's generic rescuer); r ≡ 0 mod 27 blocks shell 27
by the multiple, r ≡ ±10 completes the 9th unit ±-class cover (the core's units mod 27 cover
{±1,±2,±4,±5,±7,±8,±11,±13}, missing exactly ±10). The rescue: either band-2 shells, or B′ on the
stranger — the safe lattice of the core (the 1/7-scaled {1,…,12} safe components) dodges r's thin
arcs. **The dodge target must be allowed to be ANY runner — dodging the multiple of n fails here.**

## 4. Significance

- **Answers HYP-2436 Test 3 / t-0122: NO.** The hoped route `LRC(14) ⟸ LRC(5) + LRC(7)` is not
  available: the 3-adic head degenerates to divisor conditions (§1), and the two heads jointly
  leak (§3). The corrected route is the **fibered band ladder**: cores/d enter the proven
  LRC(≤13) regime as pure divisor conditions (m ≤ 13), strangers ride the widening bands, and the
  residual climbs exactly one rung per blocking layer — each layer consuming runner-residue
  resources (13 | r, r ≡ ±10 or 0 mod 27, …) from a budget of only 13 runners.
- **Corrects the S622-era toolkit:** [twisted shells m ≤ 2n−1] ∪ [B′ restricted to multiples of n]
  has a **nonempty residual** at n = 14 (the five evaders); the 46000-config sampling missed this
  measure-thin family. With B′ on any runner, no residual is known (HYP-2438).
- **The ladder quantifies "the cover always leaks" (S643):** no finite band horizon suffices in
  principle (a runner ≡ 0 mod lcm of all moduli blocks every clock — but is then B′-dominant);
  the ladder ∪ B′(any) is the candidate finite closure (HYP-2438).

## Artifacts

`04-computation/lrc14_two_head_exactness_cbx7.py`, `lrc14_fibered_lattice_cbx7.py`,
`lrc14_hyp2197_check_cbx7.py`, `lrc14_family_theorem_cbx7.py` (+ `.out` in
`05-knowledge/results/`). All arithmetic exact (integers / Fraction); every claimed witness
re-certified by direct evaluation. New: **HYP-2438** (the lattice-closure conjecture).
