---
id: LEM-007
title: The triple-overlap mass law and the covering-variance structure — (I) the j-fold arc-overlap mass E[|A_{i1} ∩ … ∩ A_{ij}|] equals the IID value L^j (L = 1/7) plus THM-641-style Bernoulli corrections that VANISH for generic differences; in particular the triple iid mass is (1/7)³ = 1/343, and the consecutive triple {0,1,n} has mass exactly L²/n = 1/(49n) (crossing iid at n=7). (II) Via the exact measure inclusion-exclusion W = 1 − k/7 + S₂ − S₃ + ⋯ (Sⱼ = Σ j-fold overlaps), Var(W) has leading additive-energy coefficient Var(ov)/2 = 7.64·10⁻⁴ (ov = the pair-overlap trapezoid), but the NET coefficient is c = Var(W)/R2 ≈ 5.85·10⁻⁵ (± 5%), an ~92% reduction from the S₃⁺ cancellation — the covering resonance is dominated by the higher inclusion-exclusion layers, not a pair correction
status: (I) PROVED for the iid values (L^j; elementary) and derived via THM-641's offset-sweep for the correction structure (the {0,1,n} = L²/n law and the vanishing at generic differences machine-verified exactly); (II) the leading coefficient Var(ov)/2 = (2/1029 − 1/2401)/2 is EXACT; the net c ≈ 5.85·10⁻⁵ is EMPIRICAL (stable, 30-shape sample min/mean/max = 5.34/5.85/6.38·10⁻⁵) — a rigorous Var(W) ≤ c·R2 with the correct small c requires controlling the S₃⁺ cancellation (the barely-covers wall in variance form), NOT proved here. HONEST CORRECTION to klein-S181: the brick-(B) c is ≈ 6.4·10⁻⁵ (worst-case), not the optimistic 5.67·10⁻⁵, so brick (B) via PZ needs E[W] ≥ 0.1393 (margin +0.004 vs min E[W] ≈ 0.143), thin — the D3 route (margin +0.134) is the robust one.
source: klein-2026-07-08-S182 (HYP-5387)
depends_on:
  - THM-641   # the offset-sweep + Bernoulli-representation method for overlap masses (extended here to triples)
  - THM-638   # the pair-overlap mass (the j=2 case: L² + correction)
  - THM-657   # W = uncovered measure; the covering reformulation this expands
related:
  - THM-656   # the tent-variance resonance (this is the covering analog, but with the full I-E hierarchy)
  - THM-662   # brick (B) — this supplies the exact c the brick's PZ/D3 bound needs
  - THM-660   # PZ = E[W]²/E[W²]; Var(W) = E[W²] − E[W]² is what this bounds
external: Bernoulli polynomials / the three-distance theorem; Stevens iid circle covering (L^j overlaps).
---

# LEM-007 — the triple-overlap mass law and the covering variance

## (I) The j-fold overlap mass law

For arcs `A_i(x) = [frac(e_i x), frac(e_i x) + L)`, `L = 1/7`, the `j`-fold overlap
`|A_{i1} ∩ … ∩ A_{ij}|(x) = (L − span_j(x))_+`, `span_j =` circular diameter of the `j` phases.

> **IID value: `E[|∩_j arcs|] = L^j`** (each arc independently covers a fixed point with prob
> `L`, so `E[|∩|] = ∫ P(y ∈ all j) dy = L^j`). Triple: `(1/7)³ = 1/343`.

For actual phases, `E[|∩|] = L^j + (correction)`, and by THM-641's offset-sweep the correction is
a sum of periodized-Bernoulli terms in the pairwise differences that **vanishes for generic
(large, non-resonant) differences** (machine-verified: the triples `(3,7),(10,21),(50,97),(1,100)`
all hit `1/343` exactly). The corrections are `O(1)`-small and signed:

> **The consecutive triple `{0, 1, n}` (differences `1, n, n−1`) has overlap mass exactly
> `L²/n = 1/(49n)`** — above iid for `n < 7`, equal at `n = 7` (`= L³`), below for `n > 7`.

(Verified: `n=2,3,4 → 1/98, 1/147, 1/196`.) This is the triple analog of THM-638's signed
pair-mass law; the general closed form is the THM-641 Bernoulli sum over the difference lattice.

## (II) The covering variance = additive energy through the I-E hierarchy

The measure inclusion-exclusion is EXACT (unlike the divergent PROBABILITY I-E): with
`S_j(x) = Σ_{|T|=j} |∩_{i∈T} A_i|(x)`,

> `W(x) = 1 − meas(∪ A_i) = 1 − k/7 + S₂(x) − S₃(x) + S₄(x) − ⋯`.

`S₂ = Σ_{i<j} ov(frac(d_ij x))`, `ov(s) = (L − |s|)_+` the pair-overlap trapezoid, has EXACTLY
THM-656's structure: `∫ov = L² = 1/49`, `∫ov² = 2L³/3 = 2/1029`, so `Var(ov) = 2/1029 − 1/2401`
and the leading additive-energy term of `Var(S₂)` is `(R2/2)·Var(ov)`, i.e. per unit `R2`

> **leading `c_lead = Var(ov)/2 = (2/1029 − 1/2401)/2 = 7.64·10⁻⁴`** (EXACT).

But `Var(W) ≠ Var(S₂)`: the higher layers `S₃, S₄, …` (with their triple-overlap `L³ + …`
masses from part I) **cancel ~92%** of it. The net is empirically stable:

> **`c := Var(W)/R2 ≈ 5.85·10⁻⁵`** (30-shape k=11 sample: min/mean/max `5.34/5.85/6.38·10⁻⁵`,
> std `2.7·10⁻⁶`) — the covering "resonance lemma," but the small `c` lives in the I-E
> cancellation, not the pair layer. `c_lead / c ≈ 13`.

So the covering variance is additive-energy controlled (`Var(W) ≈ c·R2`, `c ≈ 6·10⁻⁵`), but unlike
the tent (THM-656, single pair layer, `Var = R2·V1 + small resonance`), here the resonance (the
`S₃⁺` layers) is DOMINANT. A rigorous `Var(W) ≤ c·R2` with the correct `c ≈ 6.4·10⁻⁵` requires
controlling that cancellation — the barely-covers wall in second-moment form. The exact leading
coefficient and the triple-mass law (part I, the first cancellation layer) are the tools for it.

## Consequence for brick (B) (honest correction to klein-S181)

Brick (B) via `PZ = 1/(1 + Var(W)/E[W]²)` needs `Var(W)/E[W]² ≤ (1−bar)/bar = 2.02`, i.e. with
`R2 ≤ 614` and worst-case `c = 6.4·10⁻⁵`: `E[W] ≥ sqrt(6.4·10⁻⁵·614/2.02) = 0.1393` (not the
`0.1313` from S181's optimistic `c = 5.67·10⁻⁵`). Margin `+0.004` vs `min E[W] ≈ 0.143` — thin;
the **D3 route (margin +0.134) is the robust brick-(B) closure**, tolerating the true `c`.

## Files
`lrc14_triple_overlap_varW_klein_S182.out` (iid `L^j`, triple `= 1/343`, `Var(S₂)` vs `Var(W)`),
`lrc14_triple_bernoulli_klein_S182.out` (the `{0,1,n} = L²/n` law, net-c stability 5.34–6.38·10⁻⁵).
