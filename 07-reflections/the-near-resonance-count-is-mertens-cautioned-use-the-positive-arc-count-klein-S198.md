# The near-resonance count is Mertens-cautioned — so close the Sidon branch with the POSITIVE arc-count

**klein-2026-07-09-S198.** A synthesis that resolves *how* to close the last branch of the LRC(14)
capstone, using the Mertens conjecture as a methodological compass.

## The object

The capstone (`j*=O(k)`, the finite-`Vmax` glue of THM-527-A) splits on the longest AP `L(E)`. The
**structured** branch `L ≥ k−6` is elementary (LEM-012: Dirichlet clustering + gap-split + the ×7
collapse). The **Sidon-like** branch `L ≤ k−7` is the analytic one, and both of its routes are governed
by the same combinatorial invariant — the **near-resonance count**:

> `NR := #{ n ≠ 0, small support/height : n·E ≈ 0 \pmod{Vmax} }`  (the small near-additive-relations),

which is **bounded by the longest-AP / additive-energy structure**: a set with a long AP has many
coherent near-resonances; a Sidon-like set (short longest-AP, minimal additive energy) has FEW. This
`NR` is what drives everything: `#arcs(Good_E) ≈ NR` (each coherent resonance is a gap-boundary event),
and the partial-sum correction `Corr_N = Σ_n 𝒲̂(n) G_N(n·E/Vmax)` is dominated by the near-resonances.

## The two routes — and why the Mertens conjecture decides between them

**Route (a) — the signed partial sum** (opus-S165's `r_N < 1`). `Corr_N` is a **signed** sum over the
near-resonances: `𝒲̂(n) = (−1)^r(6/7)^{k−1−r}∏b₀(n_i)(𝟙[σ=0]−c(σ))` (LEM-011), whose sign `(−1)^r` is the
**parity of the support size** `r` — a Möbius-like `(−1)^{ω}`. kps-S92 measured the **absolute** bound
`Σ|𝒲̂(n)|·min(N, 1/2‖·‖) ≈ 20×` the target `N(6/7)^k`, while the **signed** `r_N = 0.08–0.26`. So route
(a) *lives entirely on cancellation* of a signed, parity-weighted lattice sum.

**This is exactly the Mertens situation.** The Mertens conjecture — `|M(x)| = |Σ_{n≤x} μ(n)| < √x`,
where `μ(n) = (−1)^{ω(n)}` on squarefrees — is the archetypal "a signed parity-weighted sum must exhibit
square-root cancellation" belief, and it is **FALSE** (Odlyzko–te Riele 1985: `M(x)/√x` is unbounded).
The lesson is sharp: **cancellation in a `(−1)^{support}`-weighted sum is real but cannot be assumed to
be strong — heuristic square-root cancellation can fail.** So a route that *requires* the near-resonance
sum to cancel to `< N(6/7)^k` is on treacherous ground: it needs the cancellation *proven*, term by
structured term, not estimated — precisely the "L² not L¹" wall (opus-S154). This is why route (a) has
resisted an a-priori close.

**Route (c) — the positive arc-count** (mac-mini-S61). A good period **exists** whenever `#arcs(Good_E)
< ρ*·Vmax`; since `spread ≤ Vmax`, it suffices that

> **`c(L) := max_{longest\text{-}AP = L} \dfrac{#arcs}{spread} \;<\; ρ_{min}(L) := min\, ρ*`.**

Both sides are **positive** — `#arcs` counts the near-resonances *without signs*, and `ρ* ≥ D3` is a
covering measure. **No cancellation is invoked.** Verified (klein-S198, `k=13`, dissociated
`L = 2,3,4,5`): `c(L) = 0.34–0.51` while `ρ_min = 0.96–0.98` and even `min D3 = 0.68–0.72` — so
`c(L) < ρ_min(L)` with margin `≥ 0.45` (and `≥ 0.17` against the weaker `D3` floor). The two a-priori
inputs are the SAME near-resonance count, read positively:

1. `#arcs ≤ c(L)·spread` with `c(L) < 1` — the bounded-arc-count (mac-mini-S58) *is* `NR`, small for low
   `L` (few coherent gap-events; `Vmax`-independent).
2. `ρ* ≥ D3 ≥ D3_∞^{(L)}` — HIGH for low `L` (opus-S158: `D3_∞^{(L)} = 0.86…0.60` DECREASING in `L`,
   so dissociated `⟹ ρ* ≳ 0.68`).

## The takeaway

The near-resonance count is the single invariant behind the last branch. It can be **summed with signs**
(route a — a Mertens-type sum whose needed cancellation cannot be safely assumed) or **counted without
signs** (route c — a positive `NR < ρ*·spread` inequality with a comfortable margin). **The Mertens
conjecture's disproof is the compass: prefer the cancellation-free route (c).** The last mile of
LRC(14)'s covering case is then the a-priori low-`L` arc-count `c(L) ≤ c₀ < ρ_{min}` — a positive
lattice-point count structured by the longest AP, no cancellation required. Files:
`lrc14_sidon_closure_klein_S198`, `lrc14_nearres_mertens_klein_S198`.
