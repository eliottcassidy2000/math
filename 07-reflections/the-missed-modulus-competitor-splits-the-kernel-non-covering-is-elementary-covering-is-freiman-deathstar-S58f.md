# The missed-modulus competitor forces q′∣(vᵢ+vⱼ) — and splits the AP-extraction kernel

*death-star-2026-07-19-S58f. Following the S58e handoff — "force the foreign-denominator
competitor `q′∣(vᵢ+vⱼ)` of the close pair." The competitor is explicit; it proves the kernel in the
non-covering case outright and reduces the rest to the covering (Freiman) case. This does NOT close
the kernel, but it isolates its elementary half from its hard half and reconciles the S58e covering
question. Script: `lrc14_missed_modulus_competitor_deathstar_S58f.py`.*

## The competitor (PROVED, elementary)

> **Missed-modulus lemma.** For a speed set `V`, let `k′ = min{ k ≥ 2 : k ∤ vᵢ for all i }` be the
> smallest modulus `V` misses. Then `t = 1/k′` is a competitor time with `min_i ‖vᵢ/k′‖ ≥ 1/k′`,
> so **`M(V) ≥ 1/k′`**. The two runners realizing the minimum sit at residues `j, k′−j (mod k′)`,
> hence the tied pair satisfies **`k′ ∣ (vᵢ+vⱼ)`** — exactly the `q′∣(vᵢ+vⱼ)` form.

*Proof.* If `k′ ∤ vᵢ` for every `i`, then `vᵢ mod k′ ∈ {1,…,k′−1}`, so
`‖vᵢ·(1/k′)‖ = min(vᵢ mod k′, k′−vᵢ mod k′)/k′ ≥ 1/k′`. Taking the min over `i` gives the bound; the
runner nearest `0` from each side pairs a residue `j` with `k′−j`, whose speeds sum to `0 (mod k′)`.
∎

This is the classical covering reduction, but stated as the *explicit* foreign-denominator
competitor the kernel needs, and it lands on the exact `q′∣(vᵢ+vⱼ)` shape. Verified: 400/400 random
13-sets missing a modulus have `M ≥ 1/k′`, zero violations.

## Consequence — the kernel splits cleanly

Contrapositive of the lemma at `k′ ≤ 13`:

> **`M(V) < 1/13 ⟹ V covers `{2,…,13}`** (a multiple of every `k` in that range).
> (And `M<1/14 ⟹ V` covers `{2,…,14}` = Cover14 — the classical reduction.)

So **every strict-interior family `(1/14 < M < 1/13)` automatically covers `2..13`.** The
AP-extraction kernel ("≤ 1 small gap") therefore splits:

- **Non-covering half — DONE.** If `V` misses any `k ≤ 13`, then `M ≥ 1/13`: `V` is not strict
  interior at all. The missed-modulus competitor `t=1/k` (with `k∣(vᵢ+vⱼ)`) is the explicit witness.
  The `S58d` fold-back `{1,…,11,13,24}` is *not* this case — it covers `2..13` and lives at the
  boundary `M=1/14`, `val=1`; it is excluded by strictness (`val≥2`), consistent with S58e.
- **Covering half — the Freiman wall.** For families covering `2..13`, the competitor is no longer a
  small-denominator time (those are all killed by the covering). It sits at a **pair-sum denominator**
  `q′ = vᵢ+vⱼ` (THM-724), and the value only *barely* clears `1/13`. Verified on covering non-AP
  cores: `{1,…,11,13,156}` (core `{1..13}∖{12}`, `156=12·13`) has `M = 13/161 ≈ 0.0808 > 1/13` at
  `q′ = 161 = 156+5`; `{1,…,10,12,13,110}` has `10/113`; `{1,…,11,13,1092}` has `91/1097`. Every
  covering non-AP core tested has `M > 1/13`, each with `≥ 2` small gaps — the deep well
  `{1,…,12,182}` (`14/183`) is the only strict-interior one, with an AP core and one small gap.

## What this settles, and what it doesn't

- **Reconciles S58e.** Covering `2..13` is *not* an extra hypothesis on the kernel — it is *forced*
  by `M < 1/13`. So the earlier "is covering needed?" question is dissolved: the kernel is exactly a
  statement about covering families, automatically.
- **Isolates the hard half.** The kernel's non-covering half is elementary (missed modulus); its
  covering half is precisely the sharp Freiman `3k−4` / E₃-Schur inverse (boxeph-S89, THM-730):
  *among 13-sets covering `2..13`, only the dilated-AP core is strict interior.* No elementary
  small-denominator competitor exists there — the pair-sum competitor's margin over `1/13` is the
  additive-energy deficit, and bounding it is the resummation wall (macmini-S76).
- **Does NOT prove the kernel.** The covering half is the genuine crux and remains open.

## Next

The covering half now has a crisp target: for a covering `2..13` family with non-AP core, show the
pair-sum competitor `q′ = vᵢ+vⱼ` clears `1/13`. The value there is `val′/q′` with `q′ = vᵢ+vⱼ`; the
task is to lower-bound `val′` by the core's distance-from-AP — i.e. THM-730's *quantitative* Schur
inverse, now localized to a single explicit competitor denominator rather than the whole spectrum.

— Related: `the-ap-extraction-kernel-is-global-not-local-...-deathstar-S58e.md`,
`the-ap-extraction-crux-is-a-residue-gap-rigidity-...-deathstar-S58d.md`,
`the-lrc14-crux-is-sharp-freiman-...-boxeph-S89.md`, `the-last-inch-is-third-order-...-macmini-S76.md`.
THMs: THM-724 (pair-sum), THM-730 (E₃/Schur), THM-523 (covering reduction), THM-1017 (bridge).
HYP-7310 (crux), HYP-7740/7742 (residue-gap reduction, maximizer lemma).
