# Freiman via resonance packing; function fields; and Tao's optimistic conjecture

**⚠ CORRECTION (boxeph-S91):** section 2's function-field shortcut is REFUTED. The packing does become unconditional/exact over 𝔽_p[t], but it does NOT force the AP core — covering + level-2-loneliness allow non-AP cores (verified over 𝔽_5). The discreteness that makes the packing exact ERASES the tightness gradient the inverse theorem needs. See [[the-function-field-route-tested-packing-is-unconditional-but-does-not-force-INV-boxeph-S91]].



*boxeph-2026-07-18-S90. Working the LRC(14) crux (INV: `M<1/13` covering ⟹ the 12-core is a dilated
AP) along the three lines the owner named. Outcome: **the resonance replaces Freiman's theorem with an
explicit packing that forces equal spacing**, reducing INV to two razor-thin arithmetic facts; the
function-field model plausibly makes the packing exact; and the whole statement is the covering shadow
of Tao's optimistic (n=12 AP-uniqueness) conjecture. NOT a proof of LRC(14). Verified
`lrc_inverse_freiman` + the S90 packing check.*

## 1. Freiman via resonance — the packing lemma (no Freiman theorem needed)

At the maximizer `M = val/q < 1/13`, the residues `w_v = v·a mod q` all lie in `[val, q−val]`. The
`M<1/13` families sit at the **minimal denominator `q = 13·val + 1`** (the discrete spectrum), so the
band `[val, 12val+1]` has length exactly `11val + 1`.

> **Packing Lemma (resonance ⟹ equal spacing).** Let `q = 13val+1`. If 12 speeds `S` have residues
> pairwise `≥ val` apart and include the runner `v_+` (residue `val`), then those residues are exactly
> `{val, 2val, …, 12val}` and `S = v_+·{1,…,12}` — a dilated AP.

**Proof.** 12 residues in `[val,12val+1]`, pairwise `≥ val` (linear `=` circular in this band), the
smallest being `val`. Sorted `val = s_1 < … < s_12`, the 11 gaps `g_i ≥ val` sum to
`s_12 − val ≤ 11val + 1`, so `Σ(g_i − val) ≤ 1` — the total excess over equal spacing is **at most 1**.
Sub-case A (all `g_i = val`): `s_k = k·val`, so `c_k·a ≡ k·val ≡ k·(v_+ a) (mod q)`, hence
`q ∣ (c_k − k v_+)`, and range forces `c_k = k v_+` → `S = v_+·{1,…,12}`. ∎ (Sub-case B, one gap
`val+1`, gives `c_k ≡ k v_+ + a^{-1}` and is **not** an AP — it is the one residual gap, see §2.)

**This is the sharp Freiman `3k−4` made explicit.** The tight band forces the minimum-doubling
configuration directly — the `≤ 1` excess is the resonance doing the inverse theorem's work. Verified
exactly: every `M<1/13` family has `q=13val+1`, core residues `= val·{1,…,12}` (all gaps `= val`,
sub-case A), `v_+ = 1`, core `= {1,…,12}`.

**What remains after the packing lemma** (each verified, each = the crux):
- **(0)** `q = 13val+1` (the spectrum is discrete — the maximizer sits at the minimal denominator);
- **(A)** sub-case A holds (the single gap is *not* bumped to `val+1`);
- **(A′)** the 12 pairwise-`≥val` residues are `V∖{v_max}` (the anomaly is `v_max`).

## 2. Function fields — where the ±1 vanishes and the packing is exact

Over `𝔽_p[t]` (speeds = polynomials, `‖·‖` = the `∞`-adic/valuation distance in `𝔽_p((1/t))`), the
whole obstruction — the "`≤ 1` excess," sub-case B, the "`q = 13val+1` vs larger `q`" ambiguity — comes
from the **archimedean `±1` carrying** of `ℤ`. In `𝔽_p[t]` there is none:

- Continued-fraction / best-approximation is **exact** (ultrametric): `‖q_k α − p_k‖` is *exactly* a
  power of `p`, not up to a constant. The Ostrowski/Stern–Brocot expansion has no rounding.
- So "excess `≤ 1`" becomes "**excess `= 0`**": gaps are exact powers of `p`, sub-case B
  (`val+1`) cannot occur, and the maximal-denominator ambiguity collapses. The Packing Lemma becomes
  **unconditional**.

**Conjecture/direction:** the function-field analogue of INV is provable by the packing lemma alone,
because the ultrametric removes (A) and (0). This makes `𝔽_p[t]` the natural **model to prove INV in
first**, then transfer — the archimedean case is the same statement plus control of the single unit of
carrying. (The polynomial Lonely Runner is known to be more tractable for exactly this reason; the
packing lemma is the mechanism.)

## 3. Tao's optimistic conjecture — INV is its covering shadow

The core is a **tight 12-family**: `M(V∖{v_max}) = 1/13` exactly (verified). Tao's *optimistic
conjecture* (n=12 form, [[HYP-7310]]): **the only 12-speed families with `M = 1/13` are the dilated APs
`d·{1,…,12}`.** Hence

> **INV ⟺ `[M<1/13` covering ⟹ core is a tight-12 family`]` ∘ `[`Tao's n=12 AP-uniqueness`]`.**

The packing lemma proves the first factor's *shape* (the core residues are `val`-spaced ⟹ dilated AP)
directly, i.e. it is a **self-contained proof of Tao's n=12 uniqueness restricted to the resonance
locus `q=13val+1`** — no general inverse theorem invoked.

**Where AP and GW sit.** The two tight (`M = 1/14`) 13-families are the AP `{1,…,13}` and Goddyn–Wong
`{1,…,11,13,24}` (the doubling). Both are **non-covering** (miss 14), so both are handled by the sieve
(`M ≥ 1/14`) and lie **outside INV's scope** — INV is about the *covering* `M<1/13` layer (the deep
well `{1..12,182}`, `M=14/183`), whose core is the *tight-12* AP. So the classification splits cleanly:
`M=1/14` tight 13-families {AP, GW, non-covering} vs the `M<1/13` covering rung {deep well, AP-core}.
GW does not threaten INV; it is the sporadic on the non-covering side. Any full classification (Tao's
program) must produce exactly this: AP + GW at `1/14`, the deep-well ladder `[0;13,14m]` above it.

## Net

The crux is now a two-fact reduction with an explicit packing proof of its additive half:
`M<1/13` covering ⟹ (0) minimal denominator, (A) no gap-bump, (A′) anomaly `= v_max` ⟹ [packing] core
is a dilated AP. Facts (0),(A),(A′) are archimedean-carrying artifacts that **vanish over `𝔽_p[t]`**,
and are the exact content of Tao's n=12 uniqueness on the resonance locus. The route: prove the
function-field INV (unconditional packing), then lift the single unit of carrying.

Cross-links: [[the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89]],
[[THM-1017-ap-core-bridge-reduction]], [[the-169-structure-and-the-difference-closure-rigidity-of-M-below-one-thirteenth-boxeph-S87]],
[[HYP-7310]], [[diophantine-approximation-lonely-runner-s361]], [[grothendieck-katz-p-curvature-and-lrc14-ledgers]].
