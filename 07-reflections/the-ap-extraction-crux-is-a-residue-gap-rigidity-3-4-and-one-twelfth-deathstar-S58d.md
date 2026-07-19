# The AP-extraction crux is a residue-gap rigidity — how 3, 4, and 1/12 braid into it

*death-star-2026-07-19-S58d. Working the owner's directive: attack the n=12 AP-uniqueness
inverse theorem (the sole open half of LRC(14) after THM-1017) by reading the three
thoroughly-explored threads — the **third order** (3), the **Freiman `3k−4`** (4), and **`1/12`**
— against each other. Outcome: a rigorous **sharpening of boxeph-S87's difference-closure**
into a single local kernel, and a clean identification (with a witness family) of exactly where
the covering hypothesis becomes essential. **LRC(14) is not closed.** Scripts:
`lrc14_ap_residue_gap_reduction_deathstar_S58d.py`, `lrc14_val1_foldback_obstruction_deathstar_S58d.py`.*

## The crux, in residue coordinates

The open supplier (THM-1017; `HYP-7310` = Tao's optimistic conjecture): a primitive **Cover14**
family `V` with `M(V) < 1/13` has `V∖{v_max} = d·{1,…,12}`. Write `M = val/q` (THM-999), maximizer
`t = a/q`. Since `1/14 < M < 1/13`, THM-999 gives **`13·val < q < 14·val`**. Every residue
`r_i = v_i·a mod q` lies in the band `[val, q−val]` (because `‖v_i t‖ = |r_i|_q/q ≥ val/q`), a band
of length `q − 2val < 12·val`.

Sort the 13 residues `s_1 < ⋯ < s_13`. Their 12 consecutive gaps sum to
`span = s_13 − s_1 ≤ q − 2val < 12·val`.

## The rigorous kernel (sharpening boxeph-S87)

boxeph-S87's **difference-closure lemma** (proved): pigeonhole on 13 points in a band `< 12val`
forces **at least one gap `< val`**, and that gap is a resonance-aligned non-speed difference.
Its stall was "upgrade *one* aligned gap to the *global* AP." The residue picture collapses that
stall to one line of arithmetic:

> **Gap-excess pinning (elementary, proved).** Suppose **at most one** of the 12 gaps is `< val`.
> Then the other 11 gaps are each `≥ val`, so `span ≥ 11·val`; combined with `span < 12·val`, the
> 11 large gaps have **total excess `Σ(g−val) < val`**. Hence the residues are the perfect
> arithmetic progression `{val, 2val, …, 12val}` plus one edge residue, **up to a total
> perturbation of less than one step `val`** — exactly the Freiman `3k−4` near-AP regime, with the
> `< val` slack being boxeph-S90's single unit of "archimedean carrying."

Verified exactly on the deep wells `{1,…,12, 182m}`: residues `{m·14·j}_{j=1..12} ∪ {q−val}`,
gaps `[val×11, g₀]` with `g₀ = q − 13val` (`=1` at every deep well). Eleven core gaps **exactly**
`val`, one far gap `g₀ < val`.

So the entire global AP-extraction reduces to a **single local statement**:

> **KERNEL.** For a Cover14 family with `M < 1/13`, the twelve **core** residues have all eleven
> internal gaps `≥ val` (equivalently: no aligned non-speed difference lies *inside* the core; the
> one `< val` gap guaranteed by difference-closure is the far element's proximity to the core).

Difference-closure already delivers the one small gap. The kernel is "**no second small gap**,"
after which gap-excess pinning + integrality give the exact AP.

## Where covering becomes essential — the `val = 1` fold-back

`M < 1/13` **alone does not force an AP core.** Witness (this session):

> `V = {1,…,11, 13, 24}` has `M = 1/14 < 1/13`, core `{1,…,11,13}` is **not** an AP.

It escapes because it is **not Cover14** (no multiple of 14) and sits at the degenerate scale
**`val = 1`, `q = 14`**: the band is all of `[1,13]`, the residues fill it, and the "far" element
`24 ≡ 10 (mod 14)` **folds back** onto the residue of the speed `10` (a repeated residue, gap `0`).
The clean `{j·val} ∪ {q−val}` structure needs the far element at the band **edge** `q−val`, not
folded into the interior — and that is precisely what covering (a forced multiple of 14, pushed to
the `lcm(13,14)=182` scale) supplies. A directed search finds **no** `val ≥ 2`, non-Cover14,
non-AP-core family with `M < 1/13` (one-hole cores, `w < 120`): the rigidity switches on exactly at
`val ≥ 2`, and covering is what forces `val` up off the degenerate `val = 1` floor. **This locates
the covering hypothesis inside the proof:** it is not decoration — it removes the `val=1` fold-back
sheet on which the inverse theorem is false.

## The 3–4–1/12 braid, made precise

- **3 (third order / Schur, THM-730).** The core residues, once pinned to the AP `{j·val}`, carry
  the maximal `\binom{12}{2}` additive relations `r_i + r_j ≡ r_k`; that E₃/Schur maximality is the
  *content* of "these residues are an AP," and it is invisible at orders 1–2 (translation-blind,
  macmini-S76). The kernel above is the order-3 statement in gap language.
- **4 (Freiman `3k−4`).** The `< val` total excess in gap-excess pinning is exactly the one unit of
  carrying in the sharp Freiman `3k−4` inverse (`|C−C| = 2·12−1 = 23` iff AP; boxeph-S89/S90).
  Killing that last unit — near-AP residues `⟹` exact-AP integer speeds — is the `4`-flavored step.
- **1/12 (LRC(12) floor).** `M(core) ≥ 1/12` (settled LRC(12)) is the scale that makes `val ≥ 2`
  meaningful and underwrites the Hamming rungs THM-1004/1005/1006 (their good interval exists
  because `1/12 > 2/25`). It is the inverse-side `1/12`; the Bernoulli `1/12 = B₂/2` (disc_v,
  THM-732) is the *forward*-side engine that would resum the Schur deficit into `L>0`. The
  `two-twelves` split (kps) is why these must not be conflated.

## Honest status and the next attempt

- **Proved this session:** gap-excess pinning (elementary); the `val=1` fold-back witness; the
  reduction of AP-extraction to the local KERNEL "no second small gap" (+ integrality carry).
- **Not proved:** the KERNEL itself for Cover14 families. This is genuinely metric/covering — the
  `val≥2` search supports it but no argument yet forbids a second interior aligned gap.
- **Next attempt:** prove the KERNEL by forbidding *two* interior aligned non-speed differences
  `δ_a, δ_b` (both `|δ·a|_q < val`): their combinations `δ_a ± δ_b` land in `(−2val, 2val)`, and if
  covering forces any such combination (or a nearby multiple of `13` or `14`) to *be* a speed, then
  `|speed·a|_q < val` contradicts `val = min`. That is the same "covering supplies the missing
  difference" mechanism boxeph-S87 §4 flagged — now aimed at a single, crisp, second-gap target
  rather than the whole AP.

— Related: `the-169-structure-and-the-difference-closure-rigidity-of-M-below-one-thirteenth-boxeph-S87.md`,
`the-last-inch-is-third-order-...-macmini-S76.md`, `two-twelves-...-n14-kps.md`,
`the-lrc14-crux-is-sharp-freiman-additive-energy-and-a-discrete-markoff-spectrum-boxeph-S89.md`,
`freiman-via-resonance-packing-function-fields-and-taos-optimistic-conjecture-boxeph-S90.md`.
THMs: THM-1017 (bridge), THM-730 (E₃/Schur), THM-1004/1005/1006 (Hamming ≤3), THM-999/724 (val/q),
THM-732 (Bernoulli disc_v). Hypotheses: HYP-7310 (crux), HYP-4382 (12-speed equality).
