---
id: THM-758
title: The far-count split — covering families with at most 3 elements above 14 are closed by THM-738; the f>=4 complement remains open after per-core capped-envelope, coherent-pack, cluster, and exact-certificate dispatch. Absolute far count is not a hardness classifier: scale-coherent covering rays have f=13 and M=1/13
status: PARTIAL REDUCTION. Claim A (f<=3 implies at least 10 speeds in {1..14}, hence THM-738) is PROVED. Claim B is NOT globally finite-decided: THM-755 gives a finite interval only for each fixed core, the claimed terminal cutoff near 500 is sampled, and the S105 bank is capped and structurally restricted. HYP-6780 proves the cutoff is scale-covariant and supplies an unbounded primitive covering f=13 ray below it. The remaining target is a scale-normal coherent/cluster/incoherent classification, not raw bounded enumeration
source: klein-2026-07-14-S309
depends_on:
  - THM-738   # kps: every ≥10-in-{1..14} family is lonely (PROVED) — Claim A's engine
  - THM-746   # opus: exact density floor for two named tail shapes; not a universal large-diameter dispatch
  - THM-726   # multi-killer M≥1/13 (≥2 far outliers) — supports Claim B
related:
  - THM-755   # opus capped-envelope (H proved v>v*) — the analytic twin of this structural split
  - THM-753   # mac-mini safe-peel — the other structural reducer
  - HYP-6720  # klein-S309 (this reduction)
external: LRC(≤13) SETTLED.
---

# THM-758 — the far-count split (covering reduces to the loose ≥4-far case)

The structural reduction that dodges the equidistribution. The whole endgame's hard object — the disc /
k=7 / equidistribution wall — lived in the *tight* (low-`M`, near-covering-min) families. This theorem
shows those families are **all** in a half of the covering class that kps has already **proved**, and the
complementary half is uniformly **loose**.

## The split

Let `S` be a covering 13-set and `f = #{s ∈ S : s > 14}` (the number of far elements).

**Claim A — `f ≤ 3`: PROVED.** Then `S` has `≥ 13 − 3 = 10` elements `≤ 14`, i.e. `|S ∩ {1,…,14}| ≥ 10`,
so `M(S) ≥ 1/14` by **kps THM-738** (every 13-speed family with `≥10` speeds in `{1..14}` is lonely,
PROVED via the exact Bonferroni tree on all 1001 ten-element bodies). *Pure counting + a proved theorem.*

**Claim B — `f ≥ 4`: `M(S) ≥ 1/14` (OPEN in general).**  The original sample suggested
`M≥0.097`, but HYP-6780 refutes that stronger claim: the primitive covering near-dilates
`V_c={c,2c,…,12c,13c+1}` have `f=13` and exact `M=1/13<0.097` at arbitrarily large scale.  They still
satisfy LRC(14), by THM-757, and illustrate the missing coherent-pack dispatch.  Samples of genuinely
incoherent `f≥4` bodies are loose, but neither far count nor raw diameter proves that classification.

## Why the original wall-dodging conclusion was too strong

The deep well and the named undilated binding families do lie in Claim A.  That does not imply every
scale-coherent or low-margin family does: dilation preserves `M` while changing `f`.  The correct split
must first quotient coherent packs and translated/hierarchical clusters.  Only the normalized
incoherent residual is a plausible loose-margin class.

The S308 redirect "low-`M` covering implies near-AP or safe" may still have a scale-normal version, but
its absolute far-count form is false: `V_c` is low-`M` with `f=13`.  THM-738 proves only the forward
statement `f≤3 => lonely`; it does not classify every low-`M` family.

## Per-core decidability does not give a global finite band (codex-S1 correction)

For an `≥4`-far covering `S`, let `v=max(S)`, `P=S\{v}`, and
`v*=r_P/(π|G'_P|)` (THM-755).  This gives the following per-core dichotomy, not a global finite check:

- **`v > v*`:** the capped-envelope `disc_v ≤ 4r_P|G'_P|/(πv) + 2|G'_P|²` gives `disc_v < 6|G'_P|²`
  (⟺ `v > v*`), so THM-731's `L(S) = (6/7)|G'_P| − ε_v > 0`, i.e. `M(S) > 1/14`. **PROVED in one peel.**
- **`v ≤ v*`:** this bounds `v` only relative to that core.  Under `P→cP`, good-set measure is invariant,
  the component count and `v*` both multiply by `c`, so no uniform raw bound follows.  The `497` figure is
  the maximum in 120 sampled families.  S105 exact-verified 8260 generated interval-core cases but capped
  at 4000 per `k` and restricted the outlier pool; it did not execute the global residual.

Thus Claim B currently equals **[`v>v*`: THM-755, proved per peel] + [coherent/cluster routes, partly
proved] + [scale-normal residual, open]**.  HYP-6750's `q≤25` good-period bank is useful exact evidence
on 120 sampled residuals, not yet a uniform bounded-denominator theorem.

## The band residual closes via a bounded-`q` RATIONAL WITNESS, not a crude bound (klein-S312)

The `v ≤ v*` residual (band residual: `≥4`-far, bounded diameter, not capped-envelope-certifiable) has **no
crude analytic bound** — tested and refuted (HYP-6750). Both natural forms fail at `δ=1/14`:
- **Bonferroni** (`B1..B7`, the odd-truncation lower bounds on `G(1/14)=P(N=0)`): all **negative**
  (`B5≈−2`), the inclusion–exclusion oscillates despite `G_true≈0.12`.
- **Absolute relation-lattice**: `G=(6/7)¹³+Σ_{rel}∏ĝ`, but the absolute sum `rel_abs≈27000` (the
  absolute-symbol `b(θ)=Σ|ĥ_m|e(mθ)` has a `−log|2sin πθ|` singularity). The series converges only
  **conditionally** — `G≈0.12` is a signed cancellation of `~10⁴`-size terms. This is the known
  **"signed not absolute" cancellation wall**; no unsigned/truncated certificate can capture it.

**Evidence and candidate resolution.** A good period is a small-`q` rational lonely witness `a/q` with
all `(c·a mod q) ∈ [q/14,13q/14]`; it is a rigorous, cheap certificate for each family where it is found.
The S312 bank found one with `q∈[15,25]` for 120/120 sampled residuals (median 17).  This motivates the
uniform conjecture "every normalized incoherent residual has such a good period", but does not prove it:
the raw residual is not globally bounded, the sample is not exhaustive, and scale-coherent rays require
their separate pack/witness routes.  The proposed loose-good-period versus tight-Bonferroni dichotomy is
therefore a useful target, not yet the covering endgame theorem.

## What remains

The main obligation is a **scale-normal structure theorem** for `f≥4`: coherent dilation packs should
route to THM-668/737, translated or hierarchical clusters to THM-739/740, capped incoherent peels to
THM-755/731/732, and the remaining normalized shape-plus-residue atlas must be shown finite or
recursively decreasing.  Independently, running the still-CLAIMED THM-741 four-slot census would extend
the exact near-AP closure from `f≤3` to `f≤4`.

**Status ledger.** PROVED: Claim A (`f≤3`).  PROVED per fixed peel/core: THM-755.  PROVED on named
coherent and cluster families: the corresponding pack/cluster theorems.  VERIFIED only: broad `f≥4`
looseness and bounded-`q` witness banks.  OPEN: completeness of those routes on the scale-normal
`f≥4` residual.

*Files: `04-computation/lrc14_claimAB_klein_S309.py`, `lrc14_far_split_klein_S309.py` (+.out). HYP-6720.
The structural twin of opus's analytic capped-envelope (THM-755) and mac-mini's safe-peel (THM-753).*
