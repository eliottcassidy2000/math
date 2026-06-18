# B(k) via the Fourier identity μ(E)=F(k)+relation-correction (kind-pasteur-2026-06-18-S5)

The uniform-floor lemma `inf_E μ(E) > 0` (`μ(E)=meas{x: maxgap{frac(e_i x)}>2/7}`, `0∈E`, `k=|E|≤13`)
has a clean exact governing identity and a finite-pattern reduction. This is the analytic core of B(k).

## The exact identity (Parseval — rigorous)
Write `e=(e_1,…,e_{k-1})` (drop `e_0=0`) and `g:T^{k-1}→{0,1}`, `g(v)=1[maxgap{0,v_1,…,v_{k-1}}>2/7]`.
The orbit `x↦(e_i x mod 1)` is a closed curve, so
> **μ(E) = ∫₀¹ g(e x) dx = Σ_{m∈ℤ^{k-1}} ĝ(m)·1[m·e=0] = F(k) + Σ_{m≠0, m·e=0} ĝ(m),**
where `ĝ(0)=∫_{T^{k-1}}g = F(k)` is the iid ceiling (PROVED exact: `F(k)=Σ_{j≥1}(−1)^{j+1}C(k,j)(1−2j/7)_+^{k−1}`;
F(13)=3132376013/13841287201≈0.2263). The correction sums `ĝ` over the **nonzero relation lattice**
`Λ(E):={m∈ℤ^{k-1}: m·e=0}` (a rank-`k−2` sublattice). **CONSEQUENCE: μ(E) depends ONLY on Λ(E).**

## Validation (exact-ish grid)
- Relation-FREE large spread → μ≈F(k): k=4 `{0,1,97,9409}`→0.997=F(4); `{0,3,101,997}`→0.997;
  k=6 `{0,7,53,311,1009,4999}`→0.897=F(6). The correction →0.
- A SINGLE short relation → small dip: `{0,1,500,501}` (1+500=501) → 0.979 (one term `ĝ(1,1,−1)`).
- MANY short relations (low complexity) → large dip: consecutive `{0,1,2,3}`→0.905, `{0..12}`→0.179.
So the correction's size is governed by the **number of SHORT vectors of Λ(E)**, which is large only for
low-complexity (AP-like / small-spread) shapes — the precise form of mac-mini's "extremal = bounded-spread."

## The finite-pattern reduction (proof strategy for B(k))
1. **μ = F(k) + Σ_{m∈Λ\0} ĝ(m)** (above, exact).
2. **Fourier decay of the max-gap indicator:** g is a BV/polytope indicator (boundaries are the hyperplanes
   `v_i−v_j=±2/7`, `v_i=±2/7`, etc.), so `|ĝ(m)| ≤ C(k)/∏_i(1+|m_i|)` (CONSTANT C(k) TO COMPUTE — the crux
   technical estimate; or replace g by a Selberg–Beurling smooth minorant `h≤g` with `ĥ(0)>0` and rapid decay,
   giving μ ≥ ĥ(0) − Σ_{m∈Λ\0}|ĥ(m)|).
3. **Split** Λ into short (`|m|≤M`) and tail. **Tail** `Σ_{|m|>M, m∈Λ}|ĝ(m)| ≤ ε(M)` uniformly (Fourier decay
   + lattice-point count; rank k−2). **Short part:** `|m|≤M` ⟹ `e` satisfies a bounded-coefficient relation,
   so only FINITELY many short-relation PATTERNS occur; each pattern fixes the short-vector contribution
   exactly (a computable rational). Bounded-spread shapes realize every pattern.
4. **inf μ = min over the finitely many short-relation patterns of [F(k)+short-correction] − ε(M)**, a FINITE
   exact computation with an explicit tail. Positive ⟹ B(k) proved ⟹ (after ∩G_P) LRC(14).

## Status / what remains
- The identity (1) is PROVED (Parseval). The decay constant C(k) / the smooth-minorant construction (2) and the
  uniform tail (3) are the remaining ANALYTIC work (= the Erdős–Turán estimate, the genuine crux of B(k)).
- The finite-pattern enumeration (4) is the computational backbone (= mac-mini's bounded-spread floor, now with
  an explicit reason for finiteness: short-relation patterns of Λ).
- This REPLACES "spread bound B(k)" (a metric notion) by "short-relation patterns of the relation lattice"
  (the intrinsic notion) — `{0,1,500,501}` has large spread but ONE short relation (small dip), showing spread
  was a proxy; the relation lattice is the right object. → THM-527 (mac-mini), HYP-2586, OPEN-Q-108.
