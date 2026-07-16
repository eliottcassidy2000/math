# THM-924 — Tournament walk reciprocity: cpA determines cpK, all walk moments, and all skew d-moments; HYP-7096 proved; THM-918's tower fully unconditional

**Status:** PROVED (four one-line lemmas; machine-verified exactly on every iso class n ≤ 6, all 456 classes at n = 7 (L3), and random n = 8..10)
**Instance:** klein-2026-07-16-S316
**Files:** `04-computation/tournament_walk_reciprocity_klein_S316.py`, `05-knowledge/results/tournament_walk_reciprocity_klein_S316.out`
**Resolves:** HYP-7096 (subsumed, vastly stronger). **Upgrades:** THM-918 (cpK tower leg now PROVED). **Explains:** the census fact cpK splits 0/116 cpA-ties (klein-S315 cont.2).
**Literature note:** L1/L2 are classical-flavored (main-eigenvalue/complement tricks); derived independently here — bibliography check queued in the backlog.

---

## Setup

Tournament `T` on `n` vertices, adjacency `A` (so `A + Aᵀ = J − I`), skew `K = A − Aᵀ = 2A − (J − I)`. Write `cpM(x) = det(xI − M)`, and the **walk resolvent** `R_M(z) = 𝟙ᵀ(zI − M)⁻¹𝟙 = Σ_{j≥0} (𝟙ᵀM^j𝟙) z^{−j−1}` — for `M = A` the coefficients are the walk counts `w_j`, for `M = K` the skew d-moments `μ_j` of HYP-7096.

## The four lemmas (each one line)

**L1.** `det(zI − A + J) = (−1)ⁿ · cpA(−1−z)`.
*Proof:* `A − J = −(Aᵀ + I)` (tournament identity), so `det(zI − A + J) = det((z+1)I + Aᵀ)` and `cp_{Aᵀ} = cp_A`. ∎

**L2 (walk moments are cpA-determined).** `1 + R_A(z) = (−1)ⁿ cpA(−1−z) / cpA(z)`.
*Proof:* matrix determinant lemma: `det(zI − A + 𝟙𝟙ᵀ) = cpA(z)·(1 + R_A(z))`; apply L1. ∎

**L3 (cpA ⟹ cpK, closed form).** `cpK(y) = 2^{n−1} [ cpA((y−1)/2) + (−1)ⁿ cpA((−y−1)/2) ]`.
*Proof:* `yI − K = 2(zI − A) + J` at `z = (y−1)/2`; determinant lemma + L2: `cpK(y) = 2ⁿ cpA(z)(1 + ½R_A(z)) = 2^{n−1}(cpA(z) + (−1)ⁿcpA(−1−z))`. ∎
Sanity: `cpK(−y) = (−1)ⁿ cpK(y)` (correct skew parity). **The skew spectrum is exactly the symmetrization of the A-spectrum about the line Re z = −1/2** (the line where `Re λ ≥ −1/2` with equality iff the eigenvector is orthogonal to 𝟙).

**L4 (skew d-moments are cpA-determined — HYP-7096, strong form).** `R_K(y) = 1 − 2ⁿ cpA((y−1)/2)/cpK(y)`.
*Proof:* reverse direction of L3 via `A = (K + J − I)/2`: `cpA(z) = 2^{−n} cpK(y)(1 − R_K(y))`. ∎
So `μ_j = 𝟙ᵀK^j𝟙` is a function of cpA alone — no deep-tie hypothesis needed. HYP-7096's conjecture (deep ties have equal d-moments) is a trivial corollary and the d-moment bookkeeping is unnecessary.

**L5 (panel reversal-closure).** `cpL_in(x)·(x−n) = (−1)ⁿ x · cpL_out(n−x)` (from `L_in = nI − J − L_out` and `L_out𝟙 = 0`). The census panel (cpA, cpL, H, τ_in, τ_out) is closed under reversal up to determined data.

## Consequences

1. **THM-918's coning tower is now fully PROVED, cpK leg included:** `cpA(C(T)) = x·cpA(T)`, and by L3 cpK of the cone is a function of that — cpA-cospectral pairs have cpK-cospectral cones, with no side condition. Invisible pairs on (cpA, cpK, cpL, H, τ_in, τ_out) exist at every n ≥ 8, unconditionally.
2. **cpK adds zero information beyond cpA** — the census's 0/116 was a theorem, not luck. The converse fails badly: at n = 7 there are only **11 distinct skew charpolys** across all 456 classes (the symmetrization collapses; note `e₁ = Σα² = C(n,2)` is forced since `tr K² = −2C(n,2)`, so cpK at n = 7 is just the pair `(e₂, e₃)` — 11 attained lattice points).
3. **Walk counts are spectral for tournaments** (L2): any two cpA-cospectral tournaments have equal numbers of directed j-walks for every j — unlike general digraphs, where walk moments (main angles) are extra data.
4. Engineering: cpK, all `w_j`, all `μ_j` computable from cpA in O(n²) coefficient arithmetic — drop K-spectrum computations from tournament pipelines (atlas columns, fingerprint codes).
