# The Eisenstein/cusp dichotomy is the three-distance theorem: regularizable ⟺ single-rotation orbit (≤3 gaps)

*kind-pasteur-2026-07-01. The owner: "think hyperconcavity and extend HYP-3768's ι-even Eisenstein / ι-odd cusp into a checkable geometric form — regularizable/−1/12 ⟺ AP-interval (ordered); un-regularizable/cusp ⟺ generic (disordered)." opus HYP-3775 already made it a binary criterion ("witness residues = scaled interval"). This session identifies that criterion as the classical **three-distance (Steinhaus) theorem**, gives it a clean checkable count, ties it to last session's continued-fraction thread, and records an honest negative about the soft "hyperconcavity" route.*

## The finding: it is the three-distance theorem

For a covering speed-set `S`, take the binding witness rotation `a/D` (`M(S)=max_τ min_v‖vτ‖` attained there) and the **witness residues** `R = {v·a mod D}`. Computed exactly (`lrc_threedistance_eisenstein_cusp_kps.py`):

| object | witness | #distinct cyclic gaps of `R` | ⇒ |
|---|---|---|---|
| construction `{1..n−2, n(n−1)}`, n=7..14 | `n/Φ₆` (`D=Φ₆`) | **3** (every n) | three-distance / AP orbit ⇒ **regularizable** |
| beaters n=7,8,9,10 | `2/13,2/15,4/33,4/37` | **4,4,5,5** | not three-distance ⇒ **cusp** |

> **regularizable / Eisenstein / −1/12 ⟺ the witness residues are a SINGLE-ROTATION ORBIT** — `≤3` distinct cyclic gaps, the **Steinhaus / three-distance** signature of `{k·α mod 1}`, i.e. a scaled arithmetic progression.
> **un-regularizable / cusp / residual ⟺ the residues need `>3` gaps** — `≥2` rotations, a higher "gap dimension" — genuinely non-modular.

The construction's residues are the interval-core `{1,…,n−2}` plus the antipode `−1`: an arc with one hole, which has *exactly* three gap-lengths `{1, 2, Φ₆−n+1}` — the three-distance theorem's three gaps. That is opus HYP-3775's "scaled interval + antipode" seen as what it classically is: a one-rotation orbit.

## Why three gaps ⟺ regularizable (the mechanism)

A single-rotation orbit `{k·α}` is an AP, and **the sawtooth sum over an AP telescopes into a Dedekind sum** (Ostrowski / continued-fraction summation) — closed in `O(log)` by reciprocity (`Φ₆≡1 mod n`, one step) to `s(n,Φ₆) = −T/(12T+6) → −1/12 = ζ(−1)` (reproduced exactly). `>3` gaps means the residues are **not one orbit**, so there is no telescoping, no Dedekind form — the cusp.

This closes a loop with last session: the three-distance gaps of `{k·α}` are the **continued-fraction convergent denominators** of `α`. So **"one rotation" = "one continued fraction" = the `O(log)` reciprocity descent** (opus HYP-3770). The regularizable locus is exactly the single-continued-fraction locus; the cusp is where a single CF does not suffice — the "gap dimension" (number of rotations) is the geometric shadow of the ι-odd **genus / cusp index** (`X₀(14)`, curve 14a).

This also **partially rehabilitates opus HYP-3773's withdrawn "beaters = higher-dimensional" hope**: beaters *are* higher-dimensional — in the three-distance / rotation-count sense (`>3` gaps = `≥2` rotations) — even though (as opus correctly showed) they are *not* literal Zagier cotangent sums. The right notion of "dimension" here is the number of incommensurate rotations, not a higher-degree cotangent product.

## The honest negative: "hyperconcavity" does NOT threshold regularizability

The owner's "hyperconcavity" prompt pointed at a *soft/analytic* order parameter. I chased two and both **fail to cleanly separate** regularizable from cusp:
- **Cyclic additive energy** `E(R)=‖\hat 1_R‖₄⁴` (last session's best predictor): it is descaling-invariant and `→1` for the construction as `n→∞` (a perfect arc in the limit) — a real fact mirroring mac-mini HYP-3774's "finite `T` ↦ regularized `−1/12`". But at finite `n` it is noisy (the n=8 beater has a residue collision `16≡1 mod 15`; its `E_norm` even exceeds the construction's).
- **Log-concave (Fejér) autocorrelation** — the literal "hyperconcavity": the construction's autocorrelation *is* log-concave (`lc_viol=0`), but so is a dense arc's only sometimes (the arc `{1..6}` in `ℤ/7` gives `lc_viol=2` because it is dense), and ~4% of random sets score 0. Noisy.

The lesson is itself a result: **regularizability is EXACT / ARITHMETIC — the AP / three-distance property — not a soft concentration property.** Soft metrics confirm the construction is maximally ordered but cannot threshold the dichotomy, because the Eisenstein/cusp split is `≤3`-gaps-or-not (a discrete arithmetic fact), not a continuous "how ordered." This is consistent with last session's theme: the cusp/hard content is arithmetic/signed, out of reach of soft (concentration / positivity) methods.

## Honest scope and limits

- **Scope:** `≤3` gaps is the clean criterion **in the covering-min regime** (large witness denominator `D ~ Φ₆ ≫ #speeds`). Validated: among large-`D` random witnesses, `0/14` have `≤3` gaps — only the AP-orbit construction does. For small-`D` (non-covering) sets, `≤3` gaps is trivial (167/398 random) and not meaningful.
- **Binary, not depth:** the gap count marks cusp-or-not cleanly, but does **not** track the residual *depth* `n/Φ₆ − M` (gaps go 4,4,5,5 while depth decreases). So it quantifies "is there a residual," not "how big" — quantifying the depth remains the open OPEN-Q-108 core.
- **Not a new proof.** It re-expresses opus HYP-3775's criterion as the three-distance theorem (a cleaner, classical, checkable form), links it to the continued-fraction/reciprocity descent, and honestly bounds the soft-metric route.

## Net

The checkable geometric form of the ι-even Eisenstein / ι-odd cusp dichotomy is the **three-distance theorem**: `regularizable ⟺ witness residues are a single-rotation orbit (≤3 gaps = scaled AP = one continued fraction) ⟺ the sawtooth sum telescopes to the −1/12 Dedekind/E₂ anomaly`; `cusp ⟺ >3 gaps (≥2 rotations, higher gap-dimension = the genus)`. The soft "hyperconcavity" order parameter does not threshold it — the split is arithmetic, and the "rotation dimension" (not a continuous concentration) is the right geometric measure of the residual.

— Related: opus HYP-3775 (interval-vs-generic criterion, corrected Zagier hope), mac-mini HYP-3774 (ζ-regularization carrier), klein/mac-mini HYP-3768 (ι-even E₂ / ι-odd cusp, s(n,Φ₆)→−1/12), opus HYP-3770 (O(log) reciprocity descent), opus HYP-3773 (Φ₆=2T+1, "higher-dim" hope — rehabilitated here as gap-dimension), HYP-3586 (genus jump 0→1 at n=14), `the-barrier-residual-is-a-dedekind-sum-...-kps.md` (last session's continued-fraction / ζ(2)↔ζ(−1) thread), `the-cyclotomic-magic-function-is-the-fejer-kernel-kps.md` (Fejér = the arc autocorrelation), OPEN-Q-108. Scripts: `04-computation/lrc_threedistance_eisenstein_cusp_kps.py`, `lrc_hyperconcavity_cyclic_kps.py`, `lrc_hyperconcavity_eisenstein_cusp_kps.py` (line-vs-cyclic lesson). Not a HYP reservation — extends opus HYP-3775 (deferring to first-pusher).
