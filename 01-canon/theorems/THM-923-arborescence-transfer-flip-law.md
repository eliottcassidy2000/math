---
id: THM-923
title: THE ARBORESCENCE TRANSFER LAW (exact one-flip laws for the determinant side, all n) — flipping arc u→v to v→u in any tournament changes the out-arborescence counts τ_r by EXACT cofactor formulas: (T) TRANSFER: Δτ_u = −F(u,v), Δτ_v = +F(u,v), where F(u,v) = det(L with rows/cols {u,v} deleted) = #{spanning 2-forests of out-trees rooted at u and v} — the flip TRANSFERS exactly F(u,v) arborescences from the old winner's root to the new winner's (bijective proof: arborescences rooted at u through arc u→v ↔ 2-forests rooted {u,v}, attach/detach the arc); the pair contribution to τ_tot = Σ_r τ_r is CONSERVED; (B) BYSTANDER: for r ∉ {u,v}, Δτ_r = (e_u−e_v)ᵀ adj(L̃_r) (e_u+e_v) = adj_uu + adj_uv − adj_vu − adj_vv — four adjugate entries of the reduced in-Laplacian, each a 2-forest count by the all-minors matrix-tree theorem. MECHANISM: the arc flip is a RANK-ONE update of the in-Laplacian, ΔL = (e_u+e_v)(e_u−e_v)ᵀ, so the matrix-determinant lemma gives every Δτ_r in closed form — the determinant side of the metagraph has an EXACT flip atom at every n, in sharp contrast to H (permanental: drift affine only at n ≤ 5, THM-914). CONTEXT: klein-S314c9 (HYP-7024) opened τ-vs-H (τ-vector separations, skew-spectrum product, BEST counts); mac-mini THM-895/902 placed τ on the curl side and made the τ-vector a complete separator at n=7. This theorem supplies their missing DYNAMICS: how τ moves under the flip walk, exactly
status: PROVED (rank-one update + matrix-determinant lemma; transfer part also bijectively) + REFEREED EXACTLY on all 64+1024 tournaments n=4,5 and 300 random n=6..9 — zero failures (arb_flip_law_referee_kps_S128c33.py). Census/drift interaction engine running (arb_census_engine_kps_S128c33.py): does τ join or linearize the THM-914 census? Paley closed form τ_r = p^{(p−3)/2}((p+1)/4)^{(p−1)/2} (from klein's skew product) pending referee
source: kind-pasteur-2026-07-16-S128 (cont.33; owner: arborescences vs Hamiltonian paths/cycles)
depends_on:
  - all-minors matrix-tree theorem (adjugate entries = rooted-forest counts)
related:
  - HYP-7024 (klein: the determinant shadow — statics), THM-895/902 (mac-mini: curl side, τ-separator)
  - THM-914 (the H-drift census law — the permanental contrast), THM-833 (Δc₃ atom)
  - HYP-7105 (the τ-census program, this session)
---

# THM-923 — the arborescence transfer law

**Setup.** L = D_in − A (in-Laplacian); τ_r = det(L̃_r) counts spanning out-arborescences
rooted at r (Tutte). Flip one arc u→v to v→u.

**Rank-one mechanism.** The flip changes L only in rows u,v: ΔL = (e_u+e_v)(e_u−e_v)ᵀ.
Matrix-determinant lemma on each reduced minor yields:

- **(T)** Δτ_u = −F(u,v), Δτ_v = +F(u,v) with F(u,v) = det(L_{−{u,v}}).
- **(B)** Δτ_r = adj(L̃_r)_{uu} + adj(L̃_r)_{uv} − adj(L̃_r)_{vu} − adj(L̃_r)_{vv}, r ∉ {u,v}.

**Bijective proof of (T).** An out-arborescence rooted at u uses arc u→v iff deleting that arc
leaves a 2-forest of out-trees rooted at {u, v}; conversely any such forest plus u→v is an
arborescence at u. So exactly F(u,v) arborescences at u die with the arc — and the same
forests plus the new arc v→u are exactly the new arborescences at v. ∎

**Consequences.**
1. The pair term cancels in τ_tot: only bystander roots move the total — the τ_tot drift is
   a sum of four-cofactor forms over bystanders, all polynomial-time.
2. Contrast with H: both the level (Tutte det vs #P permanental) and now the DYNAMICS
   (exact closed-form atom at all n vs census-nonlinear beyond n=5) split along the
   determinant/permanent divide — the metagraph's integrable vs chaotic coordinates.

## Evidence log
- [x] referee: 1,388 tournaments (all n=4,5; 75 random each n=6..9), zero failures
- [ ] census engine verdicts (Q1–Q4: fibers, affinity, linearization) — running
- [ ] Paley closed form referee; BEST/Eulerian cross-check vs klein's 50176/48512
