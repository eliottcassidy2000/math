# Hypothesis Log — Index

Every hypothesis tested in this project, whether confirmed, refuted, or open.
Organized by topic. Each hypothesis has a detail file.

**Status codes:** CONFIRMED | REFUTED | OPEN | PARTIALLY-TRUE

---

## By Status

### CONFIRMED
| ID | Statement | Why it works | Source |
|----|-----------|-------------|--------|
| HYP-001 | M[a,b] = M[b,a] (transfer matrix symmetric) | Unknown mechanism; verified n<=8 | INV-001 |
| HYP-002 | H(T) = H(T^op) | Even-r structure of W; THM-L | opus-S35c11 |
| HYP-003 | M(r) symmetric for ALL r, not just r=1/2 | Unknown; verified n<=5 | opus-S35c11 |
| HYP-004 | W(T,r) = (r-1/2)^{n-1} F(T, (2r+1)/(2r-1)) | Mobius transform; THM-K | opus-S35c11 |
| HYP-005 | Deletion-contraction: H(T) = H(T\e) + H(T/e) | Standard DC on Redei-Berge function | opus-S35c11 |
| HYP-006 | Var(c_0)/Var(H) -> 0 exponentially | Perpendicularity; Fourier hierarchy | opus-S35c11 |
| HYP-007 | Singleton cancellation in M | Algebraic swap proof | THM-055 |
| HYP-008 | GS flip orthogonality Cov(H_sym, H_anti)=0 | Trivial for ANY involution on uniform measure | opus-S35c11 |
| HYP-009 | Blueself impossible at odd n | Grid symmetry + endpoint constraint | INV-039 |
| HYP-010 | Every SC tournament has involution anti-aut | Moon's theorem + Cauchy | THM-024 |
| HYP-011 | F_k(T) = F_{n-1-k}(T^op) complement duality | Path reversal maps ascents to descents | worpitzky_F_at_2.py |
| HYP-012 | Sum_T F(T,x) = A_n(x) * 2^{C(n,2)-(n-1)} | Each perm is HP of 2^{extra} tournaments | worpitzky_restricted_eulerian.py |

### REFUTED
| ID | Statement | Why it fails | First failure | Source |
|----|-----------|-------------|---------------|--------|
| HYP-101 | Per-path identity holds for all n | 3-cycle-only formula misses longer cycles | n=6 | MISTAKE |
| HYP-102 | M(T) = M(T^op) | Fails for generic tournaments | n=4 | reversal_proof_attempt.py |
| HYP-103 | M(T) = M(T^op)^T | Also fails | n=4 | reversal_proof_attempt.py |
| HYP-104 | M is linear in A, A^T, A^2, etc. | No such combination works | n=4 | simple_pattern_search.py |
| HYP-105 | M commutes with A | Fails | n=4 | M_functional_equation.py |
| HYP-106 | M(T\e) symmetric (for deletion-contraction proof) | Deletion breaks tournament structure | n=4 | dc_symmetry_path.py |
| HYP-107 | M(T/e) symmetric (for contraction) | Also fails | n=4 | dc_symmetry_path.py |
| HYP-108 | Even-odd split equivalent to OCF | Split is CONSEQUENCE, not equivalent | - | MISTAKE-008 |
| HYP-109 | Omega(T) always claw-free | Fails at n=9 (90%) | n=9 | INV-032 |
| HYP-110 | I(Omega(T),x) always real-rooted | Fails n=9 (maximally rare) | n=9 | THM-025 |
| HYP-111 | Omega(T) always perfect | Fails n=8 (53.8%) | n=8 | INV-032 |
| HYP-112 | Blueself = SC maximizer | Fails n=6 | n=6 | INV-040 |
| HYP-113 | Fourier degree-4 "two type" decomposition generalizes to n=9 | Dimension >> 200, non-proportional coefficients | n=9 | INV-050 |
| HYP-119 | W(T,r) = (H/2^{n-1})(2r+1)^{n-1} (type B universal) | ARTIFACT of wrong W definition | worpitzky_typeB_verify.py |
| HYP-120 | F(T,x) is palindromic | FALSE: F_k(T) = F_{n-1-k}(T^op), not F_{n-1-k}(T) | worpitzky_ehrhart.py |
| HYP-114 | GS flip orthogonality explains perpendicularity | Trivially true for any involution; no content | - | opus-S35c11 |
| HYP-115 | Contraction approach proves symmetry | Dead end | - | T017 |
| HYP-116 | Contiguous block decomposition | Dead end | - | T035 |
| HYP-117 | Per-vertex decomposition of unmatched counts | Dead end | - | T045 |
| HYP-118 | Cycle bijection under arc reversal | Dead end | - | MISTAKE-005 |
| HYP-213 | H_2(T,T\v) = 0 for all tournaments T, all v | BUG in verification (MISTAKE-016); ker(d_2^rel) underestimated | n=4 | beta2_relative_debug.py |

### OPEN
| ID | Statement | Current evidence | Source |
|----|-----------|-----------------|--------|
| HYP-207 | β₂(T) = 0 for ALL tournaments T | 0 counterexamples in ~47k tests (exhaustive n≤6, sampled n≤9) | THM-100, beta2_vanishing.py |
| HYP-208 | Odd-n maximizers have nontrivial path homology | True n=3,5,7,9; split at even n=6,8 (some contractible) | THM-099 |
| HYP-209 | ALL deletions β_k>0 implies β_{k+1}(parent)>0 | Verified: β₃→β₄ (n=7), β₄→β₅ (n=9), β₁→β₁ (n=5,6) | beta4_classes_n7.py, n9_max_betti_quick.py |
| HYP-210 | β_{n-4}>0 for odd-n H-maximizers | β₁=1 (n=3,5), β₄=6 (n=7), β₅=10 (n=9) | THM-099 (S41 update) |
| HYP-211 | Odd-n H-maximizer hereditary: all del's give n-1 maximizer | n=3,5,7,9 ✓ (regular maximizers only) | n9_max_betti_quick.py |
| HYP-212 | β_top = (n-1) + δ where δ=0 for prime n, δ>0 for composite | P₇: 6=6+0, Z₉: 10=8+2. Eigenspace decomposition: trivial gives δ, non-trivial each give 1 | n9_beta5_eigenspace.py |
| HYP-213 | H_2(T,T\v) = 0 for all tournaments T, all vertices v | **REFUTED** — BUG in verification (MISTAKE-016). n=4: 16/256, n=5: 840/5120, n=6: 35328/196608 | beta2_relative_debug.py |
| HYP-214 | Paley P_p has β_{p-3} = p-1 (d=p-3 pattern) | Verified p=7,11 | paley_path_homology_theory.md |
| HYP-215 | DT alone fills all Z_2 at n=5 | 1024/1024 (100%) | beta2_hamiltonian_homotopy.py |
| HYP-216 | Tight β_2 cases (surplus=1) = β_3=1 tournaments | 80/80 at n=6 (100%) | beta2_tight_cases.py |
| HYP-201 | Char poly determines H exactly | 0 ambiguous at n=5; 3 at n=6; 36 at n=7 | char_poly_H.py |
| HYP-202 | Spectral det(I-uA) correlates with H | Corr > 0.94 for optimal u | ihara_deep.py |
| HYP-203 | M has algebraic formula in terms of cycle invariants | Partial: diagonal formula known | THM-053 |
| HYP-204 | F(T,x) ascent poly is ALWAYS unimodal | 100% at n=3-5 (exhaustive), n=6-8 (sampling) | worpitzky_roots.py |
| HYP-205 | F(T,x) is almost always log-concave | 100% except 4/1024 at n=5; 100% at n=6-8 sampling | worpitzky_roots.py |
| HYP-206 | Roots of F(T,x) are real ~89% of the time (stable rate) | n=5: 95%, n=6: 90%, n=7: 89%, n=8: 89% | worpitzky_n6_test.py |
| HYP-217 | β₂=0 for |S|=2 circulants iff doubling-closed | β₂=0 iff 2s₁≡s₂ or 2s₂≡s₁ (mod n); PERFECT at n=5,7,9,11,13; one exception n=8 (s₂-s₁=n/2) | beta2_nonzero_analysis.py |
| HYP-218 | Tang-Yau Conj 4.8 is FALSE | C₈^{1,5} has β₃=β₄=1 with S∩(-S)=∅; P₇ has β₄=6 | tang_yau_counterexample.py |
| HYP-219 | β₂=0 density threshold: β₂(C_n^S)=0 for |S|≥⌈n/3⌉ | n=9:|S|≥4✓, n=11:|S|≥4✓, n=13:|S|≥5✓, n=15:|S|=5 has 15 exceptions | beta2_threshold_analysis.py |
| HYP-220 | Arc-flip preserves β₂=0: surplus ≥ |drop| always | Exhaustive n=5 (10240 flips), n=6 (491520 flips), 0 violations | beta2_arcflip_proof.py |
| HYP-221 | Surplus=0 stable: max_drop=0 from surplus=0 (n=5), surplus=1 (n=6) | Joint (δΩ₃,δZ₂) has δΩ₃≥δZ₂ always from tight cases | beta2_surplus_zero_stability.py |
| HYP-222 | DT+cancellation fills Z₂ for ALL tournaments | Exhaustive n≤6 (1024+32768), 0 failures | THM-101, beta2_dt_cancel_filling.py |
| HYP-227 | delta_\|A_3\| = (n-3)*delta_\|A_2\| under arc flip (THM-100) | PROVED algebraically; verified 0 violations n=4-9 | beta2_delta_ratio_proof.py |
| HYP-228 | delta_\|A_2\| = 2(d_u - d_v - 1) under arc flip u->v | PROVED algebraically; verified 0 violations n=5-9 | beta2_delta_ratio_proof.py |
| HYP-229 | Transitive tournament surplus = C(n-1,4) | Verified n=3-9. O2=C(n,3), O3=C(n,4), Z2=C(n-1,3) | beta2_min_surplus.py |
| HYP-230 | Min surplus grows super-linearly: 0,1,9,<=25 for n=5,6,7,8 | Exact n=5,6; sampled n=7 (10k), n=8 (20k) | beta2_min_surplus.py |
| HYP-223 | Completing oriented graph to tournament kills β₂ | Exhaustive n=4 (all 729 oriented graphs), always | beta2_twin_obstruction.py |
| HYP-224 | Edge removal from tournament creates β₂>0 only if cycle uses both endpoints | 80/80 events at n=5: unfillable 2-cycle always involves u AND v | beta2_edge_removal_anatomy.py |
| HYP-225 | DT-only deficit at n=6 is exactly 1 (rk=9 in dim-10 Z₂) | ALL 960 deficit cases: score (1,2,2,3,3,4) or (2,2,2,3,3,3), |DT|=9 | beta2_dt_deficit_analysis.py |
| HYP-226 | ∃ vertex v with β₁(T\v) ≤ β₁(T) for every tournament T | Exhaustive n=5 (1024/1024) | beta2_relative_correct.py |
| HYP-231 | δ: H_2(T,T\v) → H_1(T\v) injective for all T,v (equiv to β₂=0 + induction) | Exhaustive n=4,5 (256+5120 pairs, 0 mismatches) | beta2_connecting_map.py |
| HYP-232 | DT deficit gap ≤ 2 at n=7,8 | Gap distribution: n=7: 91.6% gap=0, 7.7% gap=1, 0.6% gap=2; n=8 similar | beta2_dt_deficit_n7.py |
| HYP-233 | beta_2 is ARC-FLIP INVARIANT: delta(dim Z_2) = delta(rk d_3) under any arc flip | Exhaustive n=5 (10240 flips, 0 mismatches), sampled n=6 (15000), n=7 (2500), n=8 (500) | beta2_arcflip_exactness.py |
| HYP-234 | beta_2 = 0 via arc-flip invariance: beta_2(T_trans)=0 + HYP-233 => beta_2=0 for all | PROOF STRATEGY. Equivalent to: delta(rk d_3) = delta(dim Omega_2) + delta(beta_1) | beta2_arcflip_mechanism.py |
| HYP-235 | dim(Omega_2) = \|A_2\| - J_2 where J_2 = #{(a,c): c→a and A²[a,c]>0} | CONFIRMED exhaustive n=4,5,6. Junk pairs contribute exactly one linear constraint each | beta2_omega_formula.py |
| HYP-236 | ALL Z_2 cycles use ALL n vertices | Exhaustive n=5 (3600/3600 cycles), sampled n=6 (5000 tournaments). Full vertex support | beta2_filling_structure.py |
| HYP-237 | β_p = 0 for ALL p ≥ 2 at ALL n | **REFUTED**: β₃=1 at n=6 (320/32768), β₄=1 at n=7. TRUE at n≤5. ONLY β₂=0 is universal (HYP-249) | beta3_analysis.py |
| HYP-238 | χ = 1 - β₁ for all tournaments (Euler char from simplex deformation) | n=4: χ∈{0,1}, n=5: χ∈{0,1}, matches 1-β₁ exactly. NOT constant! | beta2_simplex_deformation.py |
| HYP-239 | DT+cancellation fills ALL Z₂ for ALL tournaments | Exhaustive n=5 (DT alone), n=6 (960 need cancel). Sampled n=7 (1000: 920 DT+80 cancel), n=8 (1000: 896+104). 0 failures. DT rate: 100/97/92/90% | beta2_filling_algebraic.py, beta2_dt_cancel_n7.py |
| HYP-240 | DT deficit only for scores (1,2,2,3,3,4) and (2,2,2,3,3,3) at n=6 | 720+240=960 deficit tours. All deficit cases have beta_1=0. Max deficit=2 | beta2_deficit_anatomy.py |
| HYP-241 | Max DT deficit grows slowly: 0 (n=5), 2 (n=6,7), 3 (n=8) | Sampled. Cancellation ALWAYS covers deficit | beta2_dt_n7_deficit.py |
| HYP-242 | TT-span for Omega_2 is FALSE at n>=6 | TRUE at n=5 (opus-S45), FAILS 94.6% at n=6 | beta2_tt_span_proof.py |
| HYP-243 | H₂(flag complex) = 0 for ALL tournaments | **REFUTED**: 40/1024 at n=5 (all c₃=4), 6480/32768 at n=6, up to H₂=4 at n=7 | beta2_flag_complex.py |
| HYP-244 | DT_c flow conservation fills simplicial H₂ holes | n=5: 40/40 simplicial holes filled by exactly 3 DT_c paths; NT cancels via 3-cycle flow | beta2_simp_hole_filling.py |
| HYP-245 | Ω₃ ⊋ DT+cancel extra elements never in Z₃ | 144/1024 n=5 tours have extras; ALL extras have nonzero ∂₃ | beta2_extra_omega3.py |
| HYP-246 | ∂₃∂₃ᵀ\|_{Z₂} = nI for transitive tournament | CONFIRMED n=4,5,6. λ=n dominant eigenvalue for all tournaments (68.9% at n=5, 43.9% at n=6) | beta2_laplacian.py, beta2_eigenvalue_n.py |
| HYP-247 | ∂₃∂₃ᵀ\|_{Z₂} = nI for ALL n=4 tournaments | CONFIRMED exhaustive 64/64. Single eigenvalue 4 on Z₂ always | beta2_laplacian.py |
| HYP-248 | δ: H₂(T,T\\v)→H₁(T\\v) injective for ALL (T,v) at n≥4 | FALSE: 4 failures n=5, 1 at n=6. BUT injective for SOME v in ALL T (100% at n=4,5,6) | beta2_delta_injectivity.py |
| HYP-249 | δ injective for non-source/non-sink v (1≤d⁺≤n-2) | CONFIRMED n=5 (0/600 interior failures). All failures at d⁺=0 or d⁺=n-1 | beta2_delta_sourcesink.py |
| HYP-250 | β₂=0 via LES induction: base β₂=0 at n=3, step uses δ-injectivity at interior v | PROOF STRATEGY. Reduces to HYP-249 for all n. Verified n=4,5,6 | beta2_delta_injectivity.py |
| HYP-246 | dim(Ω₃) = |A₃| - rk(C) where C = constraint matrix with rows for invalid d₁/d₂ faces | CONFIRMED exhaustive n=5 (0/1024), n=6 (0/32768). C has cross-links from NN paths | beta2_omega3_correct.py |
| HYP-247 | β₁ ≠ c₃ for tournament path homology | TRUE: β₁∈{0,1} at n=5 while c₃∈{0,...,5}. rk(bd₂\|_Ω₂) usually = C(n-1,2), not C(n-1,2)-c₃ | beta2_beta1_check.py |
| HYP-248 | β₃ > 0 first at n=6: exactly 320/32768 tournaments (blow-up of C₃) | CONFIRMED exhaustive. Two iso classes: scores (1,1,1,4,4,4) [80, c₃=2] and (2,2,2,3,3,3) [240, c₃=8] | beta3_analysis.py |
| HYP-249 | β₂(T) = 0 for ALL tournaments at ALL n | Extended HYP-207. Exhaustive n≤6, sampled n≤8 (0 failures). β₃,β₄ CAN be nonzero. β₂ is special. | beta2_surjectivity.py |
| HYP-250 | DD paths alone fill Z₂ only at n=5; non-DD Ω₃ needed at n≥6 | n=5: 100% DD fills. n=6: 960/32768 need non-DD. n=7: ~8%. n=8: ~8.5% | beta2_surjectivity.py |
| HYP-251 | Transitive tournament path complex = (n-1)-simplex: dim(Ω_p) = C(n,p+1) | CONFIRMED n=3-8. All betti=0 except β₀=1. Path chain complex isomorphic to simplicial chain complex of Δ^{n-1} | beta2_simplex_connection.py |
| HYP-252 | dim(Ω₂) = \|A₂\| - n_bp where n_bp = # backward pairs with intermediaries | CONFIRMED exhaustive n=5 (1024/1024). Junk face constraints have independent rows | beta2_exactness.py |
| HYP-253 | Edge removal creating β₂>0: always δ(Ω₃)=+3, δ(Z₂)=+2, δ(rk₃)=+3 | Exhaustive n=5: all 80 cases. Only scores (0,1,3,3,3) and (1,1,1,3,4) affected | beta2_edge_restore.py |
| HYP-254 | χ_path(T) = 1 - β₁ (Euler char for n=5 tournaments) | CONFIRMED exhaustive n=5 (1024/1024). At n=6: χ = 1 - β₁ - β₃ + β₄ (200 samples, 0 mismatches) | beta2_exactness.py |
| HYP-255 | β₂=0 verified at n=7 (2000 samples), n=8 (500), n=9 (100). 0 counterexamples. | Extends HYP-249. β₁ rate decreases: 4.3% at n=7, 0.6% at n=8 | beta2_n7_verify.py |
| HYP-256 | Circulant tournaments have β_m=0 for m≥2 (proved by Fourier decomposition) | CONFIRMED via Tang-Yau-Hess arXiv:2602.04140, Corollary 3.15 for S={1,...,(n-1)/2} | Literature |
| HYP-257 | dim H₂(T,T\v) ∈ {0,1} for all interior v | CONFIRMED exhaustive n=5,6; sampled n=7 (500). Always (h₂_rel, β₁) = (1,1) when nonzero | beta2_delta_dimensions.py |
| HYP-258 | **Every tournament has interior v with H₂(T,T\v) = 0** | CONFIRMED exhaustive n=5 (1024/1024), n=6 (32768/32768). KEY PROOF STRATEGY: gives β₂=0 via LES without δ | beta2_h2rel_zero_vertex.py |
| HYP-259 | δ: H₂(T,T\v) → H₁(T\v) injective for all interior v (1≤d⁺≤n-2) | CONFIRMED: n=5 exhaustive (0/600), n=6 exhaustive (0 interior failures), n=7 sampled (0/37) | beta2_delta_sourcesink.py, beta2_delta_n7.py |
| HYP-260 | δ fails ONLY at source/sink vertices | CONFIRMED n=5 (4 boundary failures, 0 interior). All failures at d⁺∈{0,n-1} | beta2_delta_failures.py |
| HYP-261 | Source v: all v-paths are TT, v only at position 0; NT/Ω₂ coupling absent | PROVED algebraically. Source v→everyone, so (v,b,c) always TT. ∂₁∘∂₂=0 automatic for TT | beta2_delta_source_proof.py |
| HYP-262 | **Σ_v h₂(T,T\v) ≤ 3 for ALL tournaments** | CONFIRMED: n=5 (max=3, exh.), n=6 (max=3, exh.), n=7 (max=3, 500 samples). UNIVERSAL BOUND! Gives HYP-258 by pigeonhole for n≥7 | beta2_h2rel_sum.py |
| HYP-263 | β₁(T) > 0 ⟹ Σ h₂_rel = 0 (all vertices trivial) | CONFIRMED n=5: β₁=1 ⟹ all 304 tournaments have Σ=0. If T has a 1-cycle, no relative H₂ | beta2_h2rel_sum.py |
| HYP-264 | "Good deletion" (∃v with β₁(T\v)=0) works for ALL tournaments | **REFUTED**: 24/1024 at n=5 (all regular), 960/32768 at n=6 (all (2,2,2,3,3,3)). Fraction decreases: ~1% n=7, ~0.1% n=8 | beta2_good_deletion.py |
| HYP-265 | β₂(T)=0 for tournaments is NEW (not in GLMY literature) | Searched Burfitt-Cutler, Fu-Ivanov, Tang-Yau. No tournament-specific H₂ results. Tournaments are multisquare-free (Fu-Ivanov basis applies) | S42 literature search |
| HYP-266 | dim(Ω_p) for regular n=5: [5,10,10,10,5,0] (symmetric, larger than simplex at p=3,4) | CONFIRMED. Extra Ω₃ elements exactly fill extra ker(d₂). Transitive=[5,10,10,5,1,0] | beta2_full_chain.py |
| HYP-267 | χ(Ω_*) = 1 - β₁ for ALL tournaments (all higher β_p cancel in alternating sum) | CONFIRMED n=3,4,5 exhaustive. Implies β₂-β₃+β₄-...=0. Combined with HYP-249 gives consistency | beta2_full_chain.py |
| HYP-268 | dim(A_2) = C(n,3) + 2*c₃ for ALL tournaments | CONFIRMED n=4,5 exhaustive. Each transitive triple gives 1 path, each 3-cycle gives 3 | beta2_omega_poincare.py |
| HYP-269 | Ω₂ constraint rank = #{non-allowed (a,c) with mediators} | CONFIRMED n=5. Constraints from different (a,c) pairs use disjoint 2-path sets => automatically linearly independent | beta2_omega_poincare.py |
| HYP-270 | Position cone (all 4 insertion positions) fills ker(d₂) from A₃ | CONFIRMED n=4,5. sum_{v,pos} w_{v,pos} s_v^{pos} = id on ker(d₂\|Ω₂). But image NOT in Ω₃! | beta2_homotopy_proof.py |
| HYP-271 | β₂=0 extended: confirmed n=7 (2000 samples), n=8 (500), n=9 (100). Zero failures. | Extends HYP-255. β₃≠0 at n≥7 (8.2% at n=7, 19.2% at n=8). β₂ is dimension-specific! | beta2_large_n_sample.py |
| HYP-272 | Paley T₇ has β₄=6, Betti=[1,0,0,0,6,0]. Most extreme non-vanishing at n=7 | CONFIRMED. Random n=7 never has β₄≠0 (0/2000). Paley structure creates 4-dimensional holes | beta2_large_n_sample.py |
| HYP-273 | im(d₃\|Ω₃) ≠ im(d₃\|A₃) for most tournaments | CONFIRMED: 904/1024 at n=5, 32048/32768 at n=6. Ω₃ restriction loses most of A₃ image. Yet β₂=0 persists. | beta2_omega3_filling.py |

---

## By Topic

### Transfer matrix symmetry
HYP-001, HYP-003, HYP-102, HYP-103, HYP-104, HYP-105, HYP-106, HYP-107, HYP-203

### Complement invariance
HYP-002, HYP-102, HYP-103

### Fourier / W-polynomial
HYP-004, HYP-006, HYP-008, HYP-113, HYP-114

### OCF structure
HYP-005, HYP-101, HYP-108, HYP-115, HYP-116, HYP-117, HYP-118

### Omega(T) properties
HYP-109, HYP-110, HYP-111

### Path homology / topology
HYP-207, HYP-208, HYP-209, HYP-210, HYP-211, HYP-212, HYP-213, HYP-214, HYP-215, HYP-216, HYP-217, HYP-218, HYP-219, HYP-220, HYP-221, HYP-222, HYP-223, HYP-224, HYP-225, HYP-226, HYP-227, HYP-228, HYP-229, HYP-230, HYP-231, HYP-232, HYP-233, HYP-234, HYP-235, HYP-236, HYP-237, HYP-238, HYP-239, HYP-240, HYP-241, HYP-242, HYP-243, HYP-244, HYP-245, HYP-246, HYP-247, HYP-248, HYP-249, HYP-250, HYP-257, HYP-258, HYP-259, HYP-260, HYP-261, HYP-262, HYP-263

### Self-complementary / blueself
HYP-009, HYP-010, HYP-112

### Worpitzky / Forward-edge polynomial
HYP-011, HYP-012, HYP-119, HYP-120, HYP-204, HYP-205, HYP-206

### Spectral
HYP-201, HYP-202

---

## Template for new hypothesis files

```markdown
# HYP-NNN: [Short title]

**Status:** OPEN | CONFIRMED | REFUTED | PARTIALLY-TRUE
**Filed:** [session-id]
**Resolved:** [session-id] (if applicable)

## Statement
[Precise mathematical statement]

## Motivation
[Why we thought this might be true]

## Test method
[How it was tested — script name, parameters, sample size]

## Outcome
[What happened]

## WHY it works/fails
[The crucial insight — what structural/algebraic reason makes this true or false]

## Related
- Variables: [links]
- Hypotheses: [links]
- Theorems: [links]
- Scripts: [links]

## Tags
#tag1 #tag2
```
