# THM-918 — The Coning Tower: the sink is a black hole, the source is a mirror; spectrally invisible pairs exist at every n ≥ 8

**Status:** PROVED (transform laws + persistence induction for the panel (cpA, cpL, H, τ_in, τ_out)); VERIFIED (the cpK leg, exhaustively at n = 8, 9; census-scale counts at n ≤ 9)
**Instance:** klein-2026-07-16-S315 (cont.3)
**Files:** `04-computation/coning_tower_klein_S315.py`, `05-knowledge/results/coning_tower_klein_S315.out`
**Depends on:** the cospectral-tie census (klein-S315 cont.2, `cospectral_tie_census_klein_S315.py`, HYP-7026); matrix-tree; LEM-003 orbit-stabilizer context not needed here.
**Cocylinder context:** HYP-5097 / mac-mini-S51 (the cylinder tiling model).

---

## Definitions

For a tournament `T` on `n` vertices: the **cone** `C(T)` adds a sink (every old vertex beats the new one); the **cocone** `C*(T)` adds a source; the **suspension** `Σ(T) = C*(C(T))` adds both. `cpA, cpK, cpL` are the characteristic polynomials of the adjacency `A`, skew `K = Aᵀ − A`, and out-Laplacian `L = D_out − A`. `H` = number of Hamiltonian paths. `τ_in(T)` = sorted vector of spanning **in**-arborescence counts (matrix-tree minors of `L`); `τ_out(T) = τ_in(T^op)` (win-trees).

## The transform laws (PROVED, one line each; machine-verified exactly, 30 random tournaments n = 4..7)

1. **cpA(C(T))(x) = x · cpA(T)(x)** — the sink's row of `A` is zero: block triangular.
2. **cpL(C(T))(x) = x · cpL(T)(x − 1)** — `L(C(T)) = [[L + I, −𝟙], [0, 0]]`: block triangular after the diagonal shift (each old vertex gains one out-arc).
3. **H(C(T)) = H(T)** — the sink extends every Hamiltonian path uniquely at its end and can appear nowhere else.
4. **τ_in(C(T)) = (0, …, 0, det(L_T + I))** — no in-arborescence can be rooted off the sink (the sink has no out-arc), and the minor at the sink is `L_T + I`. **The sink crushes the in-tree vector to a single determinant.**
5. **τ_in(C*(T)) = (0, n · τ_in(T))** — the source is a universal leaf: block-triangular minors pick up the factor `n` (its out-degree). **The source scales; it deletes nothing.** Dually `τ_out(C(T)) = (0, n · τ_out(T))` via `C(T)^op = C*(T^op)`.

Laws 1–3 + 5 show the cone is **spectrally silent** (A- and L-spectra and H transform functorially, independent of everything but the base's) while law 4 shows it is **destructive for exactly one witness**: the in-tree resolution. This asymmetry is invisible to crossing geometry (HYP-5097: both cylinder ends price transitivity identically) — the arborescence panel sees the orientation of the cylinder that the crossing number cannot.

## Census-scale consequences (VERIFIED, exhaustive)

Let a **deep tie** be an (cpA, cpK, cpL)-cospectral pair of iso classes. At n = 7 there are 27 equal-H pairs inside deep-tie groups.

- **27/27** of their cones are invisible at n = 8 to the census panel (cpA, cpK, cpL, H, τ_in). The blind stratum multiplies: 4 pairs at n = 7 → ≥ 27 at n = 8. (No pair was saved by cpK: 0/27.)
- **The census panel was one-eyed:** τ_out splits **all 4** of the n = 7 "invisible" pairs. Joint-τ (τ_in, τ_out) blind pairs at n = 7: **0/27**.
- **The double-blind four:** exactly 4 of the 27 base pairs (H = 23, 29, 31, 43) are split *only* by τ_in (τ_out already tied). Coning launders precisely that witness: their cones tie at n = 8 on the **extended panel (cpA, cpK, cpL, H, τ_in, τ_out)**.
- **Suspensions:** 4/4 suspensions of the n = 6 equal-H deep ties are panel-invisible at n = 8.
- **Persistence:** all 4 double-blind pairs remain extended-panel-tied at n = 9 after another cone.

## The tower theorem

**For every n ≥ 8 there exist non-isomorphic tournaments on n vertices agreeing on cpA, cpL, H, τ_in, and τ_out** (PROVED: seed = the double-blind four at n = 8; induction step = laws 1–5: cpA and cpL shift functorially, H is fixed, τ_in stays collapsed at the shared value `det(L+I)` since cpL is tied, τ_out scales by the shared factor n). **Including cpK the statement is VERIFIED at n = 8, 9** and conjectured for all n (HYP-7096): the border determinant gives `cpK(C(T))(x) = x · cpK(T)(x) + 𝟙ᵀ adj(xI − K) 𝟙`, so the cpK leg persists iff the **d-moment sequence** `𝟙ᵀ K^j 𝟙` propagates its equality up the tower — it held 27/27 at the base and through both verified rungs (`𝟙ᵀK²𝟙 = −Σ d_v²` is the x-level, forced; higher j empirical).

## Reading

The cone = the cylinder (two-spine cocylinder model) with one circle collapsed to a point; the suspension = both collapsed. Coning is the **observer insertion** (Rédei = LRC + one baseline runner): the observer that everyone beats is spectrally silent yet erases the in-tree witness — a one-way mirror. To manufacture invisibility, find a pair distinguished by a single collapsible witness and apply the operation that collapses it. The **vertex-deleted deck still separates everything** (the sink card of `C(T)` is `T` itself), so these pairs quantify exactly the gap between spectral/tree invariants and reconstruction.
