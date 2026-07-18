# The bounded-spread classification is the apex-7 wall — an honest boundary

**death-star-2026-07-17-S56.** Following the census work (THM-996/997/999) and the far-element
covering bound (THM-1000), the natural next target is the **bounded-spread classification for n=14**:
prove that the only primitive tight families (`M = 1/14` exactly) are the AP `{1,…,13}` and the
Goddyn–Wong `{1,…,11,13,24}`. This note records, honestly, that this is **not** an independently
reachable lemma — it is the LRC(14) hard core itself — together with the exact point where the covering
machinery stalls, so no future session mistakes the wall for a gap.

**This is a boundary note, not a proof. The classification is NOT proved here.**

---

## 1. The landscape is consistent with bounded spread (computed)

- Every large candidate tested is **non-tight**: translated APs `{15,…,27}`, `{29,…,41}`, a GW-translate
  — all have `M > 1/14`. No "runaway" tight family appears.
- **Exhaustive:** the AP `{1,…,13}` is the *only* primitive tight family with `Vmax ≤ 21` (and, pending
  the running check, `{AP, GW}` are conjectured the only ones with `Vmax ≤ 24`). [scripts:
  `classify_probe`, `exhaustive24`.]
- THM-1000 (proved): a tight family has `Vmax ≤ (n−1)·v₂` — no dominant outlier; the deep-well far
  family is excluded from tightness.

So bounded spread *looks* true, and the far-element half is closed. What remains is the **clustered**
case: speeds split into a small cluster and a few outliers with no single dominant element.

## 2. Where the covering machinery stalls — the apex-7 wall

To bound the spread absolutely, split `V` into a small cluster `S` (`k` speeds) and large outliers `L`
(`13−k` speeds) separated by a gap. Tightness forces `L` to cover the good set `G_S` of `S` (measure
`μ_S > 0`, since `M(S) ≥ 1/(k+1) > 1/14`). The covering has a chance only if `L`'s total danger
measure exceeds `μ_S`:

| small `k` | large `13−k` | `L`-danger density `(13−k)/7` | union bound `μ_S ≥ 1−k/7` |
|---|---|---|---|
| ≤ 6 | ≥ 7 | ≥ 1 | positive — but `L` can tile, no obstruction |
| **7** | **6** | **6/7 < 1** | **collapses to 0** |
| ≥ 8 | ≤ 5 | `< 1` | 0 |

At `k = 7` the union bound on `μ_S` collapses to `0` (seven danger sets of measure `1/7` can, in
principle, cover everything), while the six outliers have density `6/7 < 1` — they cannot tile, yet
they *can* cover a good set of small measure. **Deciding whether six outliers cover a seven-speed
cluster's good set is exactly the Fraenkel/apex-7 tiling question** — the same wall that makes LRC(14)
hard in the first place (the factor-2 gap: `j = n−1` runners at threshold `1/n`; MISTAKE-122;
LEM-020's moment-order-3 necessity). Elementary covering, union bounds, and the single-element argument
of THM-1000 all live on the `≤ 6`-outlier side of this wall and cannot cross it.

## 3. Why the classification = the hard core (not an independent lemma)

The bounded-spread classification is not weaker than LRC(14): if `{AP, GW}` were proven to be the only
primitive tight families, then (both being non-covering — the AP omits a multiple of 14, GW is
sieve-covered at `q ≤ 14`) **no covering family would be tight**, so every covering family would have
`M > 1/14`, giving the strict margin that closes the covering case — i.e. LRC(14). This is precisely
THM-995's equality horn, which the repo labels "≥ LRC-hard, cannot be proved here," and OPEN-Q's inf-L
nucleus. So a genuine proof of the classification would *be* a proof of LRC(14).

## 4. What is honestly established, and what is not

**Established (rigorous):**
- No far element: `Vmax ≤ (n−1)·v₂` (THM-1000), excluding deep-well/ladder families.
- Finite verification: `{AP}` is the unique primitive tight family with `Vmax ≤ 21` (AP-only), and the
  classification holds up to the checked `Vmax` bound (a decidable, ever-extendable finite check via
  THM-999's bounded-denominator lemma, which makes tightness itself decidable per family).
- The residual is exactly the **clustered `≤ 6`-outlier** stratum — the apex-7 wall.

**Not established (and not honestly reachable here):**
- The absolute `Vmax` bound / the full classification `{AP, GW}`. This requires crossing the apex-7
  wall — the equidistribution/tiling nucleus that the whole program (both routes) has reduced to and
  that every elementary tool has failed on (klein-S283, opus-S266). It is the LRC(14) hard core.

**Verdict.** The far-element half of "primitive tight ⟹ bounded spread" is proved (THM-1000); the
clustered half is the apex-7 nucleus and is equivalent to LRC(14) itself. The session's chain
(THM-996 → 997 → 999 → 1000) is a clean *reduction* of the census residual down to this one object —
but it terminates *at* the wall, it does not pass through it. Any claim to have "proved the
classification" by elementary means should be treated as an error (cf. the overclaim entries in
`MISTAKES.md`).

→ THM-995 (VII), THM-996, THM-997, THM-999, THM-1000, LEM-020, MISTAKE-122, HYP-7300.
