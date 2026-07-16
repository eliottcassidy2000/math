---
id: THM-923
title: THE ARBORESCENCE TRANSFER LAW (exact one-flip laws for the determinant side, all n) — flipping arc u→v to v→u in any tournament changes the out-arborescence counts τ_r by EXACT cofactor formulas: (T) TRANSFER: Δτ_u = −F(u,v), Δτ_v = +F(u,v), where F(u,v) = det(L with rows/cols {u,v} deleted) = #{spanning 2-forests of out-trees rooted at u and v} — the flip TRANSFERS exactly F(u,v) arborescences from the old winner's root to the new winner's (bijective proof: arborescences rooted at u through arc u→v ↔ 2-forests rooted {u,v}, attach/detach the arc); the pair contribution to τ_tot = Σ_r τ_r is CONSERVED; (B) BYSTANDER: for r ∉ {u,v}, Δτ_r = (e_u−e_v)ᵀ adj(L̃_r) (e_u+e_v) = adj_uu + adj_uv − adj_vu − adj_vv — four adjugate entries of the reduced in-Laplacian, each a 2-forest count by the all-minors matrix-tree theorem. MECHANISM: the arc flip is a RANK-ONE update of the in-Laplacian, ΔL = (e_u+e_v)(e_u−e_v)ᵀ, so the matrix-determinant lemma gives every Δτ_r in closed form — the determinant side of the metagraph has an EXACT flip atom at every n, in sharp contrast to H (permanental: drift affine only at n ≤ 5, THM-914). CONTEXT: klein-S314c9 (HYP-7024) opened τ-vs-H (τ-vector separations, skew-spectrum product, BEST counts); mac-mini THM-895/902 placed τ on the curl side and made the τ-vector a complete separator at n=7. This theorem supplies their missing DYNAMICS: how τ moves under the flip walk, exactly
status: PROVED (rank-one update + matrix-determinant lemma; transfer part also bijectively) + REFEREED EXACTLY on all 64+1024 tournaments n=4,5 and 300 random n=6..9 — zero failures. ENGINE VERDICTS (exhaustive n=3..7, all tournaments via tilings): (Q1) τ is a NEW coordinate — the odd-cycle census does NOT determine τ_tot (splits 1/3, 3/8, 16/29, 135/203 of census fibers at n=4..7); (Q2) τ_tot not affine in the census beyond n=3; (Q3) HONEST NEGATIVE: τ_tot does NOT linearize the H-drift (drift not affine in (H,c5,c7,τ) at n=6,7); (Q4/AUTONOMY) the τ-drift E[Δτ_tot|T] is a function of τ_tot ALONE at n=3..6 — NONTRIVIALLY (11 τ-values cover 12 states at n=5, 43 cover 52 at n=6, zero inconsistencies; census-only FAILS with 344 collisions at n=6) — and AUTONOMY FAILS at n=7 (33 of 298 τ-values split), the SAME n=7 threshold where c₇ entered the H-drift (THM-914); (τ,H,c5,c7) has zero split fibers at n=7 (427 states); minimal n=7 partner being computed. Affinity: τ-drift affine at n=3 (E=3−(4/3)τ, pure OU) and n=4 (E=7−τ+H/3), nonlinear from n=5. PALEY CLOSED FORM REFEREED: τ_r(Paley_p) = p^{(p−3)/2}((p+1)/4)^{(p−1)/2}, root-independent, exact at p=3,7,11 (from klein's skew product n·τ = Π(n²+λ²)/4); BEST cross-check reproduces klein's ec = 50176 (Paley₇) / 48512 (C₇{1,2,3})
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

## Results (the τ/H contrast table, exhaustive n=3..7)

| question | H (permanental) | τ (determinantal) |
|---|---|---|
| exact flip atom | none known beyond THM-250 special cases | (T)+(B) above, all n, closed form |
| drift coordinate | census (H,c5,c7); minimal at n=7 | τ ALONE at n≤6; fails n=7 (33/298) |
| drift affine until | n=5 (THM-914) | n=4 |
| n=7 threshold | c₇ enters | autonomy breaks (same threshold) |
| census determines it? | — | NO (τ splits 135/203 census fibers, n=7) |
| linearizes the other's drift? | — | NO (Q3 negative) |

The n=7 threshold is SHARED: both drifts lose their small-n coordinate exactly when
7-cycles exist. RESOLVED (arb_n7_tuples run, 427 (H,c5,c7,τ) states): the minimal n=7
partner is **H itself** — (τ,H) DETERMINES the τ-drift (417 fibers, 0 splits) while
(τ,c5) fails (6), (τ,c7) fails (10), and even (τ,c5,c7) fails (by exactly 1 collision:
H carries odd-cycle COLLECTION information beyond the counts — the OCF's higher 2-adic
digits, visible dynamically). Census-only fails (224).

**THE ONE-WAY COUPLING LAW (n ≤ 7, exhaustive).** The flip walk's drift structure is
hierarchic: E[ΔH|T] = f_n(H, c₅, c₇) — the permanental block is autonomous and never
listens to τ (Q3 negative) — while E[Δτ|T] = g_n(τ, H) — the determinantal coordinate
listens ONLY to the permanental one (τ alone suffices for n ≤ 6; H joins at n = 7).
H is the environment; τ is driven, with no back-reaction at drift level.

## Evidence log
- [x] referee: 1,388 tournaments (all n=4,5; 75 random each n=6..9), zero failures
- [x] census engine (arb_census_engine_kps_S128c33.out): Q1–Q4 as in status
- [x] autonomy n=7: FAILS 33/298 (arb_autonomy_n7_kps_S128c33.out); n≤6 HOLDS (arb_subcoord_paley .out)
- [x] Paley closed form p=3,7,11 + BEST cross-check (arb_subcoord_paley_kps_S128c33.out)
- [x] minimal n=7 partner = H: (τ,H) determines (417 fibers); (τ,c5,c7) fails by 1
      (arb_n7_tuples_kps_S128c33.py; tuples saved to n7_tau_census_tuples.json)
- [ ] the g_n profile (τ-drift vs τ curve): shape/convexity, n=5,6 — small named follow-up
- [ ] n=8 rung ((τ,H) still closed? needs bit-packed engine) — named future
- [ ] class-level transfer on circulants (the parallel-class reflection's item 2)
