---
id: LEM-008
title: The overlap-mass Fourier kernels (triple/quad and general) — for the covering
  reformulation (THM-657), the mean S-arc overlap is E[L_S] = sum_{m in Lambda_S} prod_a
  c_{m_a}, where c_m = int_0^{1/7} e(-2 pi i m t) dt and Lambda_S = {m in Z^S : sum m_a = 0,
  sum m_a e_a = 0} is the BALANCED ADDITIVE-RELATION LATTICE of S (rank |S|-2). Triples give
  a rank-1 sum over the PRIMITIVE TRIANGLE of differences; quads a rank-2 sum encoding the
  quad's additive relations (triangle + Sidon-violation vectors). Apex-7 (c_m = 0 for 7|m)
  and dilation-invariance are built in. This is the |S|>=3 extension of THM-641's pair-mass law.
status: PROVED (elementary Fourier expansion of the arc indicator; the lattice constraint is
  the y-integral balance sum m=0 plus the x-constant-frequency sum m_a e_a = 0). VERIFIED exact
  (direct x-integration vs kernel, 0 mismatch after the primitivity fix): triples {0,1,3},
  {0,1,7}, {0,3,12} (gcd 3), {0,7,20}; quads {0,1,2,3}, {0,2,4,6} (2·AP), {0,3,6,9} (3·AP),
  {0,2,3,5} (Sidon-violation 0+5=2+3) — see lrc_overlap_mass_kernels_kps_S83.out.
source: kind-pasteur-2026-07-08-S83 (HYP-5397), extending THM-641 (mac-mini pair-mass law)
depends_on:
  - THM-657   # the covering reformulation (W = uncovered measure, L_S = S-arc overlap)
related:
  - THM-641   # the |S|=2 weighted anchored pair-mass law (this is the |S|>=3 sibling)
  - THM-637   # apex-7 invisibility (here: c_m = 0 for 7|m)
  - HYP-5397  # kps-S82: the resonance sign is non-perturbative; these kernels are its building blocks
external: none (Fourier series of a box indicator on the circle).
---

# LEM-008 — the overlap-mass Fourier kernels

## Setting

On the covering reformulation (THM-657): each speed `e_i` places a `1/7`-arc
`A_i(x) = [frac(e_i x), frac(e_i x) + 1/7)`, and `W = 1 - meas(union A_i)` is the uncovered
measure. For a subset `S ⊆ {1..k}`, the **S-overlap** is `L_S(x) = meas(∩_{i∈S} A_i)`; the
inclusion–exclusion `W = Σ_S (-1)^|S| L_S` is the basis of every covering moment. This lemma
gives `E_x[L_S]` in closed Fourier form.

Let `chi = 1_{(0,1/7)}` with Fourier coefficients `c_m = ∫_0^{1/7} e(-2πi m t) dt`
(`c_0 = 1/7`, `|c_m|^2 = sin^2(πm/7)/(π^2 m^2)`, so **`c_m = 0` whenever `7 | m`** — the
apex-7 invisibility, THM-637/638).

## Statement

`L_S(x) = ∫_0^1 ∏_{i∈S} chi(frac(y - e_i x)) dy`. Expanding each `chi` in Fourier and
integrating in `y` (which forces `Σ_{i∈S} m_i = 0`) leaves, in `x`, only the constant
frequency (`Σ_i m_i e_i = 0`) after taking the mean:

> **`E_x[L_S] = Σ_{m ∈ Λ_S} ∏_{a∈S} c_{m_a}`,  `Λ_S := { m ∈ Z^S : Σ_a m_a = 0 and
> Σ_a m_a e_a = 0 }`** — the **balanced additive-relation lattice** of `S`, of rank `|S| - 2`.

**Specializations.**
- **`|S| = 2` (pair):** `Λ = {0}` (since `m_i = -m_j` and `m_i(e_i - e_j) = 0 ⟹ m_i = 0`), so
  `E[L_{ij}] = c_0^2 = 1/49` — the pair-overlap mean (the THM-641 pair-mass law is the
  *anchored/weighted* refinement of this).
- **`|S| = 3` (triple), rank 1:** `Λ = Z · (d_1,d_2,d_3)/g`, where `(d_1,d_2,d_3) =
  (e_j-e_k, e_k-e_i, e_i-e_j)` is the **triangle of differences** (`d_1+d_2+d_3 = 0`) and
  `g = gcd(d_1,d_2,d_3)` (the dilation). Hence
  > `E[L_{ijk}] = (1/7)^3 + Σ_{t≠0} c_{t d_1/g} c_{t d_2/g} c_{t d_3/g}`.
  The mass depends only on the **primitive** triangle — dilating the triple leaves it fixed.
- **`|S| = 4` (quad), rank 2:** `Λ` is spanned by the triangle vectors **plus** any
  Sidon-violation vector (e.g. `(1,-1,-1,1)` when `e_i + e_l = e_j + e_k`). So
  `E[L_{ijkl}] = Σ_{m ∈ Λ} ∏ c_{m_a}` encodes the quad's additive relations. (`2·AP`, `3·AP`,
  and the primitive `AP` all give the same mass — dilation-invariance.)
- **general `|S|`:** rank `|S| - 2` lattice of balanced additive relations.

## Proof

`chi(t) = Σ_m c_m e(mt)`, so `∏_{i∈S} chi(frac(y - e_i x)) = Σ_{(m_i)} (∏ c_{m_i})
e((Σ m_i) y) e(-(Σ m_i e_i) x)`. Integrating over `y ∈ [0,1]` kills all terms with
`Σ_i m_i ≠ 0`: `L_S(x) = Σ_{m: Σ m_i = 0} (∏ c_{m_i}) e(-(Σ m_i e_i) x)`. Averaging over
`x ∈ [0,1]` kills all terms with `Σ_i m_i e_i ≠ 0`. What remains is exactly the sum over
`Λ_S`. The rank is `|S|` minus the two independent constraints (`Σ m = 0`, `Σ m e = 0`),
i.e. `|S| - 2` (the constraints are independent whenever the `e_i` are not all equal). ∎

The `|c_m| = O(1/m)` decay makes the lattice sum absolutely convergent.

## Why it matters (and its limit)

These are the *mean-mass* building blocks of the covering variance
`Var(W) = Σ_{|S|,|T|≥2} (-1)^{|S|+|T|} Cov(L_S, L_T)` — the object of the brick-(B) resonance
lemma (`Var(W) ≤ c·R2`, the k=11 mile). The single-overlap masses `E[L_S]` are exact and
additive-energy/relation-structured (the triple = a triangle sum, the quad = an
additive-relation lattice). **The variance itself has its own kernels** (distinct from these
mean masses): opus-S152 (HYP-5417) derived the exact **j-fold variance kernels** `c_j` via a
tent-power recursion and — crucially — showed the ~96% cancellation is the **inactive-arc
damping** `(1−θ)^{2(k−|supp|)}` (the spectator coordinates `m_i = 0`, weight `1−θ` per arc not
in the overlap), reproducing the tight `Var(W)/R2 ≈ 6.1e-5` at k=11. (This corrects kps-S82's
"non-perturbative" reading: the Fourier series converges once the `m_i = 0` spectators are
kept — I had excluded them.) So LEM-008 supplies the exact **mean**-overlap kernels; opus's
`c_j` supply the **variance** kernels; the remaining brick-(B) piece is the additive-energy
multipliers `E_3, E_4` and the uniform bound `Var ≤ c·R2` over the tail.

## Update (kps-S85): the resummation, E_3, and the confirmed non-perturbative wall

**opus-S153 resummed these per-subset masses into ONE full-set relation lattice**
(HYP-5427): `E[W] = Σ_{n∈L} ∏_i ψ̂(n_i)`, `L = {n∈Z^k : Σn=0, Σn e=0}`, `ψ̂(0)=1−θ`,
`ψ̂(m)=−c_m`. My per-subset `Λ_S` are exactly the **support-slices** of this single `L`
(the `(−1)^{|S|}` inclusion–exclusion signs are absorbed by the `ψ̂(0)=1−θ` spectator
weights). The same resummation gives `Var(W) = Σ_{n∈L} K(n)`, `K(0)` the Poisson diagonal,
`K(n≠0)` the resonance of relation `n`. So `E_2 = R2` **exactly** (the support-4 relations
`n=(1,−1,−1,1)`, additive quadruples `e_i+e_l=e_j+e_k`, once the `d=0` diagonal `k²` is
placed in the Poisson term), and `E_3` = the support-6 relations (additive sextuples).

**My S85 Fourier-support cross-check of `E_3`.** Indexing the *variance* `Σ_{ν≠0}|Ŵ(ν)|²`
by the common support size `j` of the matched Fourier vectors, the dominant order-3 form is
the **balanced double-difference** `(1,1,−2)`: `E_i+E_j−2E_l`, whose `L²`-mass `Σ_ν Q(ν)²`
(`Q(ν)=#{i<j,l: E_i+E_j−2E_l=ν}`) is the triple analog of `R2=Σ_d r_d²`. This is the same
object mac-mini LEM-007 found as the leading resonance (their `(1,−2,1)` 3-AP form). `E_3` is
an order-3 additive energy — verified tracking `Q`-energy at `k≤8`.

**The non-perturbative wall — now confirmed from three sides** (kps-S82 instinct, opus-S153,
mac-mini LEM-007, klein-S183). The per-level resonances **GROW** and only cancel coherently:
at the clustered block `k=11` I get `E_3=22k, E_4=608k, E_5=8M` (opus: `|W_j|²=0.077, 0.226,
0.932`), while the true `Var=0.047` survives by massive inter-level cancellation. So
`Var ≤ c·R2` is **NOT** provable by bounding `E_j` term-by-term — the honest `S82` reading
(the resonance sign is non-perturbative) was right *for the per-level expansion*; opus-S152's
damping converges the *mean/coherent* object, not the per-level variance series.

**The brick-(B) census (exact, no truncation).** Over 960 spread primitive 11-sets
(diam ∈ [16,44], `R2≤614`): the **decisive** degree-3 floor `D3` stays `min D3 = 0.4917 ≫
bar 0.331` — brick (B) confirmed (and klein-S184 proved it exhaustively for diam≤24, min
`D3=0.4356`). BUT the constant `c=Var/R2` reaches **7.2e-5 > the block's 6.1e-5**, so
`6.1e-5` is *not* the sup and a uniform-`c` proof cannot work. The real lever is `E[W]`:
spread sets have larger uncovered mass, which lifts the floor faster than the larger `c`
lowers it — this is what klein's LEM-009 decorrelation route exploits (`L1=(6/7)E[W]`).

## Files
`04-computation/lrc_overlap_mass_kernels_kps_S83.py` (+ `.out`): derivation + exact
verification (direct x-integration vs lattice kernel) for triples and quads, including the
gcd/dilation and Sidon-violation cases; `lrc_resonance_fourier_kps_S82.py` (the
non-perturbativity of the full `Var(W)` series); `lrc14_E3_triple_resonance_kps_S85.py`
(+ `.out`): the `E_3` Fourier-support derivation, the per-level divergence, and the exact
`Var/R2` + `D3` brick-(B) census.
