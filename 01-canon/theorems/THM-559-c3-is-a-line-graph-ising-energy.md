---
id: THM-559
title: The 3-cycle count c3 is an EXACT frustrated 2-body Ising energy on the arc spins (line graph of K_n) — the cut-space face of the OCF as a classical spin glass, with the regular tournament as ground state
status: PROVED (general n, elementary algebra; verified exact n=4,5,6 with 0 mismatches)
filed_by: mac-mini-2026-06-20-S7
depends_on:
  - THM-554   # score partition function (c3 is the last score-determined OCF datum)
  - THM-555   # cut/cycle wall: c3 score-determined, H beyond it needs cycle space
related:
  - THM-290   # H(t) is antiferromagnetic (the CYCLE-space face; this is the CUT-space face)
  - THM-126   # Paley/regular uniquely maximizes H on Z_7 (the ground state at n=7)
  - HYP-2708  # Gibbs/transfer-matrix reading of the apex-prime gas
external: Kendall–Babington-Smith (c3 ↔ score variance); Ising / Hopfield spin glass; line graph L(K_n).
---

# THM-559 — c3 is a frustrated line-graph Ising energy on the arc spins

**Context.** This realizes the project's "tournament as a spin system" / Hopfield reading
(cross-domain session 2026-06-20) on the **cut-space** OCF datum c3, rigorously and for all n.
THM-290 showed the *cycle-space* object H(t) is antiferromagnetic; this is its **cut-space twin**:
c3 (the last score-determined OCF datum, THM-555) is *exactly* a classical **2-body** Ising energy
on the arc spins, with the regular tournament as ground state. The cut/cycle seam = (exact 2-body
Ising) / (higher-order frustrated antiferromagnet).

## Setup

Tournament `T` on vertices `{0,…,n−1}`. For each edge `{u,v}`, `u<v`, the **arc spin**
`σ_{uv} = +1` if `u→v`, `−1` if `v→u`. Out-score `s_v = #{w : v→w}`. Mean score `s̄ = (n−1)/2`.
For a vertex `v` incident to edge `e`, set `ε(v,e) = +1` if `v = min(e)` (low endpoint), `−1` if
`v = max(e)` (high endpoint).

## Statement (PROVED)

> **(1) Local linearization.** `s_v − s̄ = (1/2) Σ_{e∋v} ε(v,e) σ_e`.
>
> **(2) Line-graph Ising form.** With `E_score := Σ_v (s_v − s̄)²` (the score variance, total),
> $$ E_{\text{score}} \;=\; \tfrac{1}{2}\binom{n}{2} \;+\; \tfrac12\!\!\sum_{\text{cherries }\{e,f\}@v}\!\! \varepsilon(v,e)\varepsilon(v,f)\,\sigma_e\sigma_f, $$
> a **2-body Ising energy on the line graph `L(K_n)`** (vertices = arcs, couplings = cherries =
> 2-paths), with **zero external field**. The coupling `J_{ef}=ε(v,e)ε(v,f)` is
> **`+1` iff the shared vertex `v` is an extreme (min or max) of the cherry's 3 vertices, `−1` iff
> `v` is the middle.** There are `3·C(n,3)` couplings (= |E(L(K_n))|), in ratio **2:1** (+:−).
>
> **(3) The OCF identity.** `c3(T) = C(n,3) − Σ_v C(s_v,2) = \dfrac{n(n^2-1)}{24} − \dfrac12 E_{\text{score}}.`
>
> **(4) Ground state.** Maximizing c3 ⟺ minimizing `E_score` ⟺ the **regular tournament**
> (`s_v ≡ s̄`, exists for odd n; `E_score=0`, `c3 = n(n²−1)/24`). At n=7 the H-maximizer in this
> ground manifold is Paley (THM-126).

## Proof

**(1)** For edge `e={u,v}`, `u<v`: the indicator `[v→\text{other}] = (1+ε(v,e)σ_e)/2`
(check: `v` low, `ε=+1`: `(1+σ)/2=[u→v]` ✓; `v` high, `ε=−1`: `(1−σ)/2=[v→u]` ✓). Summing over
`e∋v`: `s_v = (n−1)/2 + (1/2)Σ_{e∋v}ε(v,e)σ_e`, i.e. `s_v−s̄=(1/2)Σ_{e∋v}ε(v,e)σ_e`. ∎

**(2)** `E_score = Σ_v (s_v−s̄)² = (1/4)Σ_v Σ_{e,f∋v} ε(v,e)ε(v,f)σ_eσ_f`. The diagonal `e=f`
gives `(1/4)Σ_v Σ_{e∋v} 1 = (1/4)·2C(n,2) = C(n,2)/2` (since `σ_e²=1`). Each off-diagonal
**unordered** cherry `{e,f}` shares exactly one vertex `v` and appears twice in the ordered sum,
giving `(1/2)ε(v,e)ε(v,f)σ_eσ_f`. The coupling sign: for a cherry on `{a<b<c}`, centered at the
middle `b` we have `e={a,b}` (`ε(b,e)=−1`), `f={b,c}` (`ε(b,f)=+1`) ⟹ `J=−1`; centered at an
extreme (`a` or `c`) both `ε=+1` or both `=−1` ⟹ `J=+1`. ∎

**(3)** `Σ_v C(s_v,2) = (1/2)(Σ_v s_v² − C(n,2))` and `Σ_v s_v² = E_score + n s̄² =
E_score + n(n−1)²/4` (from `Σ(s_v−s̄)²=Σs_v²−ns̄²`). Substituting,
`Σ_v C(s_v,2) = (1/2)E_score + n(n−1)²/8 − n(n−1)/4`, so
`c3 = C(n,3) − (1/2)E_score − n(n−1)²/8 + n(n−1)/4`. The constant simplifies:
`C(n,3) − n(n−1)²/8 + n(n−1)/4 = n(n−1)[(n−2)/6 − (n−1)/8 + 1/4] = n(n−1)(n+1)/24 = n(n²−1)/24`. ∎

**(4)** `E_score ≥ 0` with equality iff every `s_v = s̄` (regular). For odd n this is integral
and achievable; `c3 = n(n²−1)/24` is then the exact max (e.g. n=5→5, n=7→14). ∎

**Verification.** Exhaustive over all `2^{C(n,2)}` tournaments at n=4,5,6: parts (2),(3) hold with
**0 mismatches**; coupling census `(+:−) = (8:4),(20:10),(40:20)` = 2:1, `#couplings = 12,30,60 =
|E(L(K_n))|`. Script: `04-computation/c3_line_graph_ising_thm559_mac-mini-2026-06-20-S7.py`.

## Significance

- **The cut-space face of the OCF is a *classical* spin glass.** c3 — the entire score-determined
  part of the OCF (THM-555 wall) — is an exact **2-body** frustrated Ising on the arcs. Everything
  beyond c3 (c5, α₂, H) lives in the cycle space and is genuinely **many-body** (THM-290's
  higher-order frustrated antiferromagnet). "Cut cheap, cycle dear" = "2-body Ising / many-body".
- **Frustration is geometric:** a cherry is satisfied (+) when its shared vertex is a local
  extreme, frustrated (−) when it is the middle — exactly the transitive-vs-cyclic distinction at
  the 3-vertex level, summed.
- **Engineering.** Maximizing c3 / finding the most-balanced (regular) tournament maps onto a
  standard zero-field Ising/Hopfield ground-state search on `L(K_n)` with explicit ±1/2 signed
  couplings — a concrete instance for the spectral-tournament-algorithm and TDA roadmap items.
