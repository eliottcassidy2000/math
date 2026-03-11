# Hypothesis Log — Index

Every hypothesis tested in this project, whether confirmed, refuted, or open.
Organized by topic. Each hypothesis has a detail file.

**Status codes:** CONFIRMED | REFUTED | OPEN | PARTIALLY-TRUE

---

## By Status

### CONFIRMED
| ID | Statement | Why it works | Source |
|----|-----------|-------------|--------|
| HYP-001 | M[a,b] = M[b,a] (transfer matrix symmetric) | PROVED by induction (THM-030, kind-pasteur-S25) | INV-001 |
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
| HYP-316 | beta_1*beta_3=0 (seesaw) — PROVED for n<=7 only | im(d_2) mediator: drops by 1 for beta_1>0, saturated for beta_3>0. Open for n>=8. | kind-pasteur-S45 |
| HYP-318 | ker(d_1) = C(n,2)-n+1 constant for all tournaments | im(d_1)=n-1 since tournaments weakly connected (Redei) | kind-pasteur-S45 |
| HYP-320 | beta_1 in {0,1} only (never >1) | im(d_2) takes exactly 2 values, gap = 1 | kind-pasteur-S45 |
| HYP-322 | beta_2=0 for all tournaments | Exhaustive n<=6, 1000 samples n=7,8: 0 violations | kind-pasteur-S45 |
| HYP-323 | beta_5 first appears at n=8 | betti=[1,1,0,0,0,1,0,0] found at n=8 (~0.2%) | kind-pasteur-S45 |
| HYP-324 | beta_4 can be >1 (value 5 observed at n=8) | Near-regular c3=20 tournament | kind-pasteur-S45 |
| HYP-325 | THM-119: Disjoint support at Omega_2 (PROVED) | Each 2-path has at most 1 NA face => constraint matrix full rank | kind-pasteur-S45 |
| HYP-326 | beta_2=0 at n=9 (0/500) and n=10 (0/100) | Extends HYP-249/276. Now confirmed n<=10 | kind-pasteur-S45 |
| HYP-327 | Removing 1 edge from tournament can create beta_2>0 | 13/500 at n=6 (2.6%). Completeness is SHARP condition | kind-pasteur-S45 |
| HYP-328 | 3-path NA face distribution converges to 25/50/25 | NOT per-tournament (transitive=100/0/0). Average over ensemble only. Max dev: 0.75(n=5)->0.21(n=8) | kind-pasteur-S45 |
| HYP-336 | dim(Omega_2) = C(n,3) + 2*c3 - e_cyc | PROVED: THM-119 gives dim=|A2|-#NA; |A2|=C(n,3)+2*c3 (HYP-268); #NA = e_cyc (edges in 3-cycles). Exhaustive n=4,5,6 | kind-pasteur-S45 |
| HYP-337 | e_cyc NOT determined by c3 alone | Multiple e_cyc values per c3 at n=5,6 (depends on cycle arrangement/sharing) | kind-pasteur-S45 |
| HYP-338 | Defect rate wave: beta_1 decreasing, beta_3 increasing with n | n=5: (29.7%,0%), n=6: (14.6%,1%), n=7: (5.8%,7.2%), n=8: (1%,21%) | kind-pasteur-S45 |
| HYP-339 | Adjacent-odd seesaw: beta_{2k-1}*beta_{2k+1}=0 for ALL tournaments | PERFECT: 0 violations in 1500+ samples n=6-8. Includes β₁β₃, β₃β₅, β₅β₇ | kind-pasteur-S45 |
| HYP-340 | At most one nonzero beta_p (p>=1) in generic tournament | 99.8% at n=7, 100% at n=8 (500 samples). Rare exceptions: [1,1,0,0,1,0,0] at n=7 (0.1%) | kind-pasteur-S45 |
| HYP-341 | beta_4 onset at n=7 (not n=8) | beta_4=6 (Paley T_7) and beta_4=1 found at n=7 with 1000 samples (0.2%). Missed in earlier 200-sample runs | kind-pasteur-S45 |
| ~~HYP-342~~ | ~~beta_{2k-1} in {0,1} for ALL k~~ | **MOVED TO REFUTED**: beta_3=2 at n=8. Boolean only holds for k=1. See REFUTED section. | kind-pasteur-S48 |
| HYP-343 | rank(d_2) takes exactly 2 values: {n-1, n} | EXHAUSTIVE n=6: {9,10}. Sampled n=7: {14,15}. rank(d_2)=n-1 iff beta_1=1. Equivalent to HYP-320 | kind-pasteur-S45 |
| HYP-344 | beta_3+beta_4 coexistence extremely rare at n=8 | 0/1000 in large sample (SVD crash after), 1/300 in separate run. Rate < 0.3%. When present: beta_5=0 always | kind-pasteur-S45 |
| HYP-345 | beta_3>0 forces beta_4=0 in generic tournaments | 43/43 at n=8 in first test, but 1 counterexample in 300 (beta_3=beta_4=1). MOSTLY true but not universal | kind-pasteur-S45 |
| HYP-346 | beta_4>0 requires self-complementary score sequence | ALL 25 beta_4>0 cases at n=8 (2000 samples) have SC scores. Very strong constraint | kind-pasteur-S45 |
| HYP-347 | beta_1+beta_5 coexistence at n=8 | Found 1/2000: betti=[1,1,0,0,0,1,0,0], scores=(3,3,3,3,4,4,4,4), c3=20, chi=-1. Extremely rare | kind-pasteur-S45 |
| HYP-348 | beta_3+beta_4 coexistence rate ~0.15% at n=8 | 3/2000 with profile [1,0,0,1,1,0,0,0]. All have SC scores | kind-pasteur-S45 |
| HYP-349 | rank(d_4) gap = ker(d_3) - rank(d_4) is always 0 or 1 | Exhaustive n=6 (32448 gap=0, 320 gap=1). Sampled n=7 (925/75), n=8 (255/45). Equivalent to beta_3 in {0,1} | kind-pasteur-S46 |
| HYP-350 | Good vertex existence for beta_3: exists v with beta_3(T\v)=0 when beta_3(T)>0 | TRUE when beta_3>0. But FAILS unconditionally: Paley T_7 has beta_3=0 with NO good vertex (240 counterexamples at n=7). Good vertex holds for non-Paley n=7 (exhaustive) | kind-pasteur-S46, corrected opus-S53 |
| HYP-351 | dim H_3(T,T\v) <= 1 for ALL tournaments T and ALL vertices v | **REFUTED at n=8** (HYP-401): dim=2 for all-good beta_3=2 tournaments. Exhaustive n=6: ALL 1920 pairs dim=1. n=7: dim in {0,1}. FAILS n>=8 | kind-pasteur-S46, REFUTED kind-pasteur-S49 |
| HYP-352 | LES isomorphism: beta_3(T\v)=0 implies beta_3(T) = dim H_3(T,T\v) | n=6 exhaustive: 1920/1920 perfect match. Follows from H_2(T\v)=0 + LES exactness | kind-pasteur-S46 |
| HYP-353 | beta_3 completely fragile at n=6: ALL 6 deletions give beta_3=0 | 320/320 exhaustive. At n=7: 24/44 fragile, 16/44 have 1 surviving, 4/44 have 2 surviving | kind-pasteur-S46 |
| HYP-354 | Quotient proportionality: all ker(d_3) basis vectors project proportionally to H_3 | Exhaustive n=6 Type B (240/240). Cokernel vector in ker(d_3) has dim 7 at n=6, varies by tournament | kind-pasteur-S46 |
| HYP-355 | H_3 generator at n=6 Type B uses ALL C(6,4)=15 vertex 4-subsets (36 paths) | 240/240 Type B. Type A uses 9 paths on 9 subsets. Structure rigid within isomorphism class | kind-pasteur-S46 |
| HYP-356 | beta_1=1 forces rank(d_4) = ker(d_3) exactly (perfect saturation, gap=0) | n=6: 4800/4800 (100%), n=7: 27/27 (100%). ker(d_3) ranges 3-18 (n=6) and 12-49 (n=7) but gap ALWAYS 0 when beta_1=1. This IS the seesaw: extra Omega_4 content fills extra ker(d_3) | kind-pasteur-S46 |
| HYP-357 | ker(d_3) = dim(Omega_3) - dim(Omega_2) + rank(d_2) (seesaw formula) | Exact formula, verified 32768/32768 at n=6. Follows from beta_2=0 (THM-108): rank(d_3)=ker(d_2) | kind-pasteur-S46 |
| HYP-358 | Good vertex selection rule: max c3(v) gives beta_3(T\v)=0 at 97.7% | n=7: 43/44 success. Bad vertices have LOW c3(v) (1-4 vs good's 4-7). Intuition: removing the most cyclic vertex disrupts H_3 | kind-pasteur-S46 |
| HYP-359 | Bad-vertex TT always rank-critical (for ALL n, not just zero-redundancy) | Exhaustive n=5 (120/120), n=6 (4320/4320). At n=6 redundancy=2-3 but bad-TT still RC. Witness: global 1-form w vanishing on all non-bad columns. | beta2_badtt_rc_n6.py |
| HYP-360 | At n=6 #bad=3: exactly 3 TTs are RC (1 BAD + 2 MIXED), rest redundant | Exhaustive: BAD=1.00 RC, MIXED=1.33 RC avg, GOOD=0.67 RC avg. Total RC=3.00 = redundancy complement | beta2_rc_mechanism.py |
| HYP-362 | Witness cocycle ζ = H₁ generator of flipped tournament T' | CONFIRMED 100%: n=5 (120/120), n=6 (200/200 sampled). ζ normalized to ⟨ζ,∂_bad⟩=1 is EXACTLY the H₁ class emerging after flip. | beta2_witness_cocycle.py |
| HYP-363 | #bad = Σ_v β₁(T\v) when β₁(T)=0 | CONFIRMED exhaustive n=5 (720/720), n=6 (27968/27968). Each bad vertex contributes exactly 1 to the sum. | beta2_witness_cocycle.py |
| HYP-364 | Hidden cycle lifts span exactly #bad dimensions (never more) | CONFIRMED n=5 (480/480), n=6 (19488/19488). dim(lift span) = #bad always. | beta2_hyp282_test.py |
| HYP-365 | Hidden cycle z_v defects are purely local: only TTs involving v violated | CONFIRMED n=5: 100% of nonzero defects involve v. z_v is cocycle of T\v → trivially satisfies non-v TTs. | beta2_why_3_max.py |
| HYP-366 | coker(res_v) = 1 for ALL bad vertices (exactly 1 missing cocycle direction) | CONFIRMED exhaustive n=5 (840/840 instances), n=6 (405/405 sampled). | beta2_hyp282_test.py |
| HYP-367 | Alternating sum of 3 hidden cycle lifts ∈ im(∂₂) (1-boundary) | CONFIRMED: z₁-z₂+z₃ is in im(∂₂) but NOT a coboundary. Being exact doesn't force dependency. | beta2_redei_telescope.py |
| HYP-368 | Non-bad vertex subtournaments have zero redundancy when #bad=3 (n=5) | CONFIRMED exhaustive n=5: 240/240 non-bad T\\w have redundancy=0. At n=6: mixed (0,1,2). | beta2_redei_telescope.py |
| HYP-369 | Bad vertices biased toward HP endpoints (positions 0, n-1) | CONFIRMED n=5,6: endpoints ~30% overrepresented vs uniform. | beta2_redei_telescope.py |
| HYP-370 | All hidden cycle lifts z_v ∈ B₁ = im(∂₂) | CONFIRMED exhaustive n=5 (480/480), n=6 (1000/1000): rank_C₁/B₁ = 0 for ALL. Each z_v is a linear combination of TT boundaries. The quotient C₁/B₁ sees NOTHING — lifts are independent WITHIN B₁, not outside it. | beta2_quotient_rank.py |
| HYP-372 | z_v MUST use bad-TT boundary (can't be expressed without it) | CONFIRMED 100%: n=5 (360/360), n=6 (84/84). Removing bad-TT column makes z_v inexpressible. Consistent with RC. | beta2_badtt_coeff.py |
| HYP-373 | Bad-TT coefficient |α_bad| = 1/√2 universal for all z_v | CONFIRMED n=5 (120 tournaments, all α = ±0.707107). |α_bad|/max|α| = 1.0 always. Sign pattern not determined by TT position. | beta2_badtt_coeff.py |
| HYP-371a | EXACT: H_3(T,T\v) = beta_3(T) - rank(i_*) when beta_2=0 | LES + beta_2(T\v)=0 gives delta=0, j_* surjective. Exact, not inequality. Verified n=6 exhaustive | opus-S52 |
| HYP-371b | beta_3(T)=2 impossible: requires ALL deletions beta_3=1, never occurs | **REFUTED at n=8** (kind-pasteur-S48): beta_3=2 CONFIRMED, profile (1,0,0,2,0,0,0,0), 4/5000 (0.08%). Scores: (2,3,3,3,4,4,4,5) and (3,3,3,3,4,4,4,4). All have good vertices. | opus-S52, REFUTED kind-pasteur-S48 |
| HYP-371c | Relative dim(R_1)=n-1 always, ker(d_1^rel)=n-2 always | R_1 = arcs incident to v (fixed for tournaments). d_1^rel has rank 1 (projects to R_0={v}). | opus-S52 |
| HYP-375 | beta_3 <= 1 at n=9 (200 random samples, 0 violations) | **REFUTED at n=8 AND n=9** (kind-pasteur-S48): b3=2 at n=8 (4/5000=0.08%), b3=2 at n=9 (1/2000=0.05%). beta_3<=1 proved ONLY at n<=7 | kind-pasteur-S47, REFUTED kind-pasteur-S48 |
| HYP-376 | 3 * H_3 cokernel generator has INTEGER coefficients in {-3,...,3} at n=6 | ALL 240 Type B tours: 3*coker rounds to exact integers. Multiplied by 3 = integer coefficients {-3,-2,-1,0,1,2,3} | kind-pasteur-S47 |
| HYP-377 | Type B cokernel |coefficient| pattern is UNIVERSAL (single pattern for all 240 tours) | Exactly ONE sorted |coeff| pattern across all Type B. 36 nonzero / 63 total. Values 1/3, 2/3, 1 | kind-pasteur-S47 |
| HYP-378 | Subtournament c3 count determines # paths per vertex set in cokernel: c3=1 → 2 nonzero, c3=2 → 3 nonzero | 100%: 720 instances of c3=1 with 2nz, 2160 of c3=2 with 3nz, at n=6 | kind-pasteur-S47 |
| HYP-379 | Hybrid numpy-int64 + mod-p Gauss gives 3.8x speedup at n=8 for beta_3 | Benchmark: SVD=175ms, hybrid=46.5ms. beta_1 TT-span 6.4x faster. | kind-pasteur-S47 |
| HYP-380 | i_* INJECTIVE when b3(T)=1 AND b3(T\v)=1: rank(i_*)=1, H_3^rel=0 | **REFUTED at n=8** (kind-pasteur-S48): rank(i_*)=0 found even with b4=0. Trials 613,627,877/5000 (seed 12345). Holds n<=7 only. NOT needed for beta_3<=1 proof (good vertex + Claim II suffices) | kind-pasteur-S47, REFUTED kind-pasteur-S48 |
| HYP-381 | LES dichotomy: b3(T\v)=0 => (rank=0,H3rel=1); b3(T\v)=1 => (rank=1,H3rel=0) | **PARTIALLY REFUTED at n=8**: bad vertex half fails (rank=0 possible). Good vertex half STILL HOLDS. n<=7: perfect. n=8: good vertex part still works | kind-pasteur-S47, updated kind-pasteur-S48 |
| HYP-382 | Relative dims at n=6 fully determined by type: Type A=(1,5,9,6,0,0), Type B=(1,5,12,14,8,3) | Exhaustive: 480 Type A pairs, 1440 Type B pairs. Rigid structure | kind-pasteur-S47 |
| HYP-383 | Bad vertices have δ(β₃)=0: adding v adds EQUAL kernel and im(d_4) | n=7: 71/71 bad vertices. Perfect saturation — T\v cycle orthogonal to new v-paths content | kind-pasteur-S47 |
| HYP-384 | Restriction res: Z₁(T)→⊕_v H₁(T\v) is SURJECTIVE (coker=0) for all β₁=0 tournaments | n=5 exhaustive (480/480), n=6 (85/85). rank(res)=#bad exactly. Each hidden cycle is detected by restriction. | beta2_global_les_test.py |
| HYP-385 | #bad ≤ 3 is SPECIFIC to β₁(T)=0; β₁=1 allows #bad up to n | n=5 EXHAUSTIVE: β₁=1→#bad∈{3,4,5}; β₁=0→#bad∈{0,1,2,3}. n=6: β₁=1→#bad∈{4,5,6}. When β₁=1, ALL n vertices can be bad. β₁=0 is essential. | beta2_image_dim_test.py |
| HYP-385b | #bad ≤ 3+β₁(T) (Grok's "compat-subspace" formula) | REFUTED: n=5 β₁=1 gives #bad=5 (24 tours), but 3+1=4. n=6 β₁=1 gives #bad=6, but 3+1=4. Actual max is n when β₁=1, not 3+β₁. | beta2_image_dim_test.py |
| HYP-385c | H₁(T)=0 imposes "3 global relations" on ⊕H₁(T\v), giving 3D compat-subspace S | REFUTED EXHAUSTIVELY with β₁=0 ONLY: d(T) spans ALL of R^n. All n standard basis vectors e_i appear (each vertex can be uniquely bad). n=5: 720 β₁=0 tours, span=5. n=6: 27968, span=6. n=7: 5000 sampled, span=7. No compat-subspace. | beta2_compat_subspace_test.py |
| HYP-387 | Score obstruction fails for beta_3=2 at n=7: 4 score seqs compatible | (2,2,2,3,4,4,4), (2,2,3,3,3,4,4), (2,3,3,3,3,3,4), (3,3,3,3,3,3,3) all allow all-deletions-beta_3=1 | opus-S53 |
| HYP-388 | Only Paley T_7 has all 7 deletions with beta_3=1 (EXHAUSTIVE) | 240/2097152 = all labelings of Paley. Betti=[1,0,0,0,6,0,0], all deletions=[1,0,0,1,0,0]. beta_3(T_7)=0 (not 2!) | opus-S53 |
| HYP-389 | Exactly 2 iso classes of beta_3=1 at n=6 | Type A: score(1,1,1,4,4,4), 80 tours, 2 c3, not SC. Type B: score(2,2,2,3,3,3), 240 tours, 8 c3, SC | opus-S53 |
| HYP-390 | HEREDITARY SEESAW: beta_3(T)=1 => beta_1(T\v)=0 for ALL vertices v | n=6 exhaustive: 1920/1920, n=7: 700/700. STRONGER than seesaw (b1*b3=0 per tournament). Forces H_2(T,T\v)=0 for ALL v | kind-pasteur-S47 |
| HYP-391 | chi_rel DICHOTOMY: good vertex => chi_rel=-1, bad vertex => chi_rel=0 | n=6: 1920/1920 chi_rel=-1 (all good). n=7: 316/316 good=-1, 34/34 bad=0. Follows from chi(T\v)=1 for good, 0 for bad | kind-pasteur-S47 |
| HYP-392 | Relative complex concentrated in degree 3: H_p^rel=0 for p!=3 when beta_3(T)=1 | H_0^rel=0 (i_* iso on H_0), H_1^rel=0 (j_* surj from H_1=0), H_2^rel=0 (beta_1_Tv=0), higher=0 at tested n. Forces H_3^rel = -chi_rel | kind-pasteur-S47 |
| HYP-404 | Score alone cannot predict beta_3(T\v): every score 0-6 has P(beta_3=1)>0 at n=7 | Sampled: all scores have 1-2% rate of beta_3=1 in deletion. No "guaranteed good" score | opus-S53 |
| HYP-405 | 3-cycle count is constant per score sequence for beta_3=1 at n=7 | 12 score seqs observed, each with unique c3 count. Follows from Rédei formula c3 = C(n,3) - Σ C(d_i,2) | opus-S53 |
| HYP-406 | Paley T_7 is the ONLY good-vertex-free tournament at n=7 (exhaustive) | 240/2097152 = all labelings of Paley. |Aut(T_7)|=21, 7!/21=240. All have beta_3=0. | opus-S53 |
| HYP-393 | beta_3 ≤ 1 at n=7 EXHAUSTIVE (2,097,152 tournaments) | Case 1: 2,096,912 with good vertex → LES gives ≤1. Case 2: 240 Paley → beta_3=0 directly | opus-S53 |
| HYP-394 | CONSECUTIVE SEESAW: beta_k * beta_{k+1} = 0 for ALL k≥1, ALL tournaments | **REFUTED at n=8** (kind-pasteur-S48): 3/2000 have beta_3=1 AND beta_4=1, confirmed by BOTH SVD and mod-p. Holds exhaustively n=6, sampled n=7 (3000). FAILS at n=8 (~0.15% rate). The opus proof architecture (LES reduction to H_4^rel=0) still works at n≤7 but not n≥8. | opus-S54, REFUTED kind-pasteur-S48 |
| HYP-395 | BAD vertex ACYCLICITY: H_p(T,T\v) = 0 for ALL p when b3(T)=1 and b3(T\v)=1 | 80 b3=1 tournaments at n=7: 60/60 bad vertices have ALL relative homology = 0. FAILS at n=8: when rank(i_*)=0, H_4^rel=1 and H_3^rel=1. | opus-S54, updated kind-pasteur-S48 |
| HYP-396 | H_4(T,T\v) = 0 for ALL vertices of beta_3=1 tournaments | Holds n=7 (560/560). FAILS at n=8: when rank(i_*)=0 and b4=0, H_4^rel=1 from LES. Paley T_7: H_4^rel=1 (exception, b3(T)=0) | opus-S54, REFUTED kind-pasteur-S48 |
| HYP-397 | CORRECTED: beta_3<=1 ONLY at n<=7. b3=2 possible at n=8 | beta_3=2 at n=8: profile (1,0,0,2,0,0,0,0), 4/5000 (0.08%). Good vertices always exist even for b3=2 tours. The "good vertex + Claim II" argument is CIRCULAR (H_3^rel=b3 for good vertex). beta_3<=1 proved exhaustively only at n<=7 | kind-pasteur-S48 |
| HYP-398 | NEW BOUNDARIES TARGET ONLY NEW CYCLES: im(d_4^new) ∩ span(old ker_d3) ⊂ im(d_4^old) | **REFUTED at n=8** but with crucial nuance: (1) Fails for GOOD vertices (b3(T\v)=0) at high rate — irrelevant to proof. (2) Fails for ~0.65% of BAD vertices (5/767, seed 12345) — EXACTLY those with rank(i_*)=0 (HYP-404). (3) HOLDS for 99.35% of BAD vertices. kind-pasteur-S49 saw failures at GOOD vertices; opus-S56 targeting_n8 (seed 42, 61 BAD) saw no failures (sample too small for 0.65%). The equivalence HYP-398⟺rank(i_*)=1 is exact. HOLDS universally at n≤7 (34/34). | opus-S55, REFUTED kind-pasteur-S49, reconciled opus-S56 |
| HYP-399 | Embedded H_3(T\\v) generator = H_3(T) generator (mod boundaries) | n=7: 44/44 BAD vertices. The true H_3 generator of T\\v, when embedded in T, IS the H_3 generator of T (up to scalar, mod im(d_4)). H_3 is "inherited" from the deletion. | opus-S55 |
| HYP-400 | rank(i_*^p) = 0 for all p >= 1 except p=3 for BAD vertices (n=7,8) | n=7: 560 vertices, n=8: 120 vertices. ALL rank(i_*) at all degrees match exactly. But this holds only when beta_3(T)=1 (not beta_3=2). | opus-S55 |
| HYP-401 | Claim II FAILS at n=8: dim H_3(T,T\v) = 2 for all-good beta_3=2 tournaments | beta_3=2 with ALL 8 deletions b3=0 (score 3,3,3,3,4,4,4,4) => H_3(T,T\v)=b3(T)=2 for every v. Claim II (H_3^rel<=1) was backbone of n<=7 proof. Fails EXACTLY when beta_3=2. | kind-pasteur-S49 |
| HYP-402 | beta_3 <= 2 for ALL n (floor(n/4) bound) | n=8: max=2 (0.07%), n=9: max=2 (0.07%), consistent with floor(n/4). beta_3=2 has SC score at n=9 | kind-pasteur-S49 |
| HYP-403 | Two structural types of beta_3=2 at n=8: all-good (SC score) and mixed (non-SC score) | All-good: score (3,3,3,3,4,4,4,4), 8/8 good vertices. Mixed: score (2,3,3,3,4,4,4,5), 6/8 good + 2/8 bad. Both types confirmed in 5000-sample run | kind-pasteur-S49 |
| HYP-404 | HYP-398 ⟺ i_*-injectivity: extra_kills=0 iff rank(i_*)=1 at n=8 | 762 injective BAD vertices: ALL extra_kills=0. 5 non-injective: ALL extra_kills=1. The equivalence is exact. extra_kills=1 means new d_4 kills exactly the H_3(T\v) generator. | opus-S56 |
| HYP-405 | Cancellation mechanism: old-face proj of new d_4 almost always in im(d_4) | n=7: 26/28 BAD vertices, old-face proj stays in im(d_4). 2/28 cases: old-face reaches H_3 direction but new-face cancels exactly. Full boundary always in im(d_4) (28/28). | opus-S56 |
| HYP-406 | b4(T)=0 for ALL b3=1 tournaments at n=7 AND n=8 (500 samples each) | Seesaw b3*b4=0 holds in the b3=1 stratum even at n=8 (where coexistence exists at 0.15%). The 500-sample test misses the rare coexistence. | opus-S56 |
| HYP-407 | chi(T) ∈ {0, 1, 7} at n=7. chi=0 iff b1>0 or b3>0. chi=7 iff Paley (b4=6) | 5000 samples. chi=1 (87.6%), chi=0 (12.3%), chi=7 (0.02%). At n=8: chi ∈ {0,1,2,3}. chi>1 from b4>0 only. | opus-S56 |
| HYP-408 | Codim-1 universality: codim(im(d_4^T)\_old\_proj, ker(d_3^T)\_old\_proj) = 1 for ALL b3=1 tournaments and ALL BAD vertices | **PARTIALLY REFUTED at n=8** (opus-S59): holds 300/300 at n=7. At n=8: 397/400 codim=1, but 3/400 codim=0 (HYP-422). Codim-0 cases = ψ escaping im(d_4^Tv). UNIVERSAL at n≤7 only. | opus-S57, updated opus-S59 |
| HYP-409 | im(d_4^{T\v}) ⊂ im(d_4^T)\_old\_proj universally | 200/200 at n=7, 122/122 at n=8. Boundaries of T\v embed into old-coordinate boundaries of T. | opus-S57 |
| HYP-410 | rank(i_*)=0 iff old-projection of embedded H_3(T\v) gen ∈ im(d_4^T)\_old\_proj | Confirmed: all failures have emb\_old\_in\_im\_old=True; all successes have False. Combined with HYP-408 (codim=1), i_*-injectivity reduces to a single linear condition. | opus-S57 |
| HYP-411 | Failure vertices have H_3(T) generator concentrated on through-v paths | Failures: 12-18 through-v path support vs success mean 9.3. Jaccard similarity of gen supports: fail 0.000-0.038 vs success mean 0.048. When H_3(T) gen relies heavily on v-paths, the embedded gen (which avoids v-paths) misaligns. | opus-S57 |
| HYP-412 | Ghost Cycle Theorem: every through-v-only cycle in ker(d_3) is in im(d_4), for b3=1 tournaments at n≤7 | n=6: trivially (no tv-only cycles). n=7: 450/450 (100%). n=8: 399/400 (99.75%). The ~0.25% failure at n=8 = rank(i_*)=0 mechanism. | opus-S58 |
| HYP-413 | Ghost Cycle ⟺ HYP-408: Ghost Cycle Theorem is EQUIVALENT to codim-1 universality, given beta_3=1 | PROVED: dim(K_tv)-dim(B_tv) = beta_3 - codim_old. Since B_tv ⊂ K_tv, equal dims ⟹ equal subspaces. Verified algebraic identity 1408/1408 across n=6,7,8. | opus-S58 |
| HYP-414 | rank(i_*: H_3(T\v)→H_3(T)) = β₃(T\v) exactly, for β₃(T)=1 tournaments | **PARTIALLY REFUTED at n=8** (opus-S59): 602/602 at n=7 ✓, but 400/400 at n=8 was SAMPLING ARTIFACT. Safe-arithmetic verification: 1199/1200 at n=8 with 2-method cross-validation (both methods agree 1200/1200). ~0.3% failure rate at n=8. n=7 confirmed universal. | kind-pasteur-S50, updated opus-S59 |
| HYP-415 | Ghost Cycle holds universally (WITH safe arithmetic): K_tv = B_tv for ALL β₃(T)=1 pairs | **PARTIALLY REFUTED at n=8**: 602/602 n=7 ✓. 400/400 at n=8 was insufficient sample (expected ~1.2 failures in 400). Safe arithmetic opus-S59: 199/200 at n=8 (full_relative_homology.py), 497/500 (relative_exactness_safe.py). Opus-S59's n=7 "failures" WERE overflow artifacts; n=8 failures ARE GENUINE. | kind-pasteur-S50, updated opus-S59 |
| HYP-416 | Block triangularity: D_{old→tv} = 0 and P_{old→tv} = 0 universally in boundary and constraint matrices | D_to=0: old p-paths have only old (p-1)-faces. P_to=0: old (p-1)-face non-allowedness only involves old p-paths. 301/301 n=7, 104/104 n=8. | kind-pasteur-S50 |
| HYP-417 | H_3 generator old component is NEVER zero for β₃(T)=1 tournaments | 0/403 at n=7 (Case 2: β₃(T\v)=0). Every H_3 class has nontrivial old-path support. D_{3,tv→tv}(z_tv)=0 always; old component NOT a cycle in T\v. | kind-pasteur-S50 |
| HYP-418 | Paley T_7 Fourier structure: each non-trivial eigenspace k=1,...,6 contributes exactly 1 to β_4=6 | Verified by Fourier eigenspace decomposition. Trivial eigenspace (k=0) has β=[1,0,0,0,0,0,0]. All 6 non-trivial eigenspaces identical β=[0,0,0,0,1,0,0]. | kind-pasteur-S50 |
| HYP-419 | T_3 Fourier: only trivial eigenspace (k=0) has nonzero Betti: β=[1,1,0] | Verified. k=1,2 contribute nothing. | kind-pasteur-S50 |
| HYP-420 | ~~β_{n-2} generically nonzero at n=8~~ | **REFUTED** (kind-pasteur-S50): opus used max_deg=6, missing im(d_7). With max_deg=7: β_6=0 for ALL 50/50 (rk(d_7)=ker(d_6) exactly). β_{n-2}=0 for ALL n<=8. | kind-pasteur-S50 |
| HYP-421 | H_4^rel = 0 universally at n=7, FAILS at n=8 | **CONFIRMED with safe arithmetic**: n=7: 200/200. n=8: 1/200 failure (full_relative_homology.py), 2/300 (relative_exactness.py). Relative d_5 surjects onto ker(d_4^rel) at n=7 (slack=0, 200/200). From LES: H_4^rel = ker(i_*: H_3(T\v)→H_3(T)). | opus-S59 |
| HYP-422 | ψ(ker) escapes im(d_4^Tv) at n=8, never at n=7 | **CONFIRMED with safe arithmetic**: n=7: 200/200 containment, codim=0 (200/200). n=8: 3/500 escapes (safe). Verified by 2-method cross-validation: both ψ-method and direct-embedding agree 1200/1200. | opus-S59 |
| HYP-423 | TOP INJECTIVITY: d_{n-1} injective on Omega_{n-1} for ALL tournaments | 1772/1772 (n=3-8). rk(d_{n-1})=dim(Omega_{n-1}) always. Consequence: β_{n-1}=β_{n-2}=0 universally. Also ker(d_{n-2})=dim(Omega_{n-1}). | kind-pasteur-S50 |
| HYP-424 | TOP VANISHING: β_{n-1}=β_{n-2}=0 for ALL tournaments on n>=3 vertices | Follows from HYP-423. Verified n=3-8 (1772 samples, 0 failures). Sharp: β_{n-3} CAN be nonzero (β_3 at n=6, β_4 at n=7 Paley). | kind-pasteur-S50 |
| HYP-425 | β_4(T) = 0 for ALL tournaments at n ≤ 7 | 1000/1000 at n=5, 1000/1000 at n=6, 1000/1000 at n=7. At n=8: β_4=1 in 5/500 (1%). ALL β_4>0 cases at n=8 have β_3=0. β_4=0 for all b3=1 tournaments even at n=8 (83/83). Simplifies LES: 0 → H_4^rel → H_3(T\v) → H_3(T). | opus-S59 |
| HYP-426 | ψ(ker(d_4^{tv→tv})) = im(d_4^{T\v}) exactly at n=7 | codim = 0 in 200/200 cases with safe arithmetic. The old-projection of relative-cycle boundaries spans EXACTLY the boundary image of T\v. Not a proper subspace — full equality. This is stronger than just ⊂. | opus-S59 |
| HYP-427 | Relative complex C_*(T, T\v) acyclic at n=7 for b3=b3_Tv=1 | All β_p^rel = 0 for all p. 100/100. Inclusion T\v ↪ T is a homology equivalence. Equivalent to i_*-injectivity since only β_3 is nonzero. | opus-S59 |
| HYP-428 | χ(T, T\v) = 5 constant for all b3=1 n=7 tournaments | 100/100 (including C_0, C_1 terms). Full χ^rel = 0 when including Ω_0, Ω_1. At n=8, χ^rel varies (7-21). | opus-S59 |
| HYP-429 | Neither d_4^{old→old} nor d_4^{tv→old} alone surjects onto im(d_4^{T\v}) | 0/200 at n=7, 0/300 at n=8 for either block alone. The surjection requires interplay between both blocks with massive cancellation (rk_oo≈9, rk_to≈9, but rk_combined≈4 at n=7). | opus-S59 |
| HYP-430 | Paley T_11: β_0-4 = (1,0,0,0,0) exact, β_5 ≤ 3400 | Direct computation via full_chain_complex_modp with increasing max_deg. max_deg=5 took 420s. Omega dims=[11,55,220,770,2255,5060]. First nonzero β beyond β_0 is at degree ≥ 5. | kind-pasteur-S50 |
| HYP-431 | Paley T_p: Omega_k = Omega_{p-k} for k=1,...,p-1 (palindrome) | **PARTIALLY-TRUE**: T_7 confirmed [7,21,42,63,63,42,21]. **REFUTED for T_11**: Omega_5^(eig)=460 but Omega_6^(eig)=700 (not 460). Per-eigenspace chi=-309 through degree 5; chi=1 per eigenspace forces sum_{m=7}^{10} (-1)^m Omega_m^(eig) = -390. For T_3: FAILS trivially. Palindrome holds ONLY for T_7 among p=3,7,11. Palindrome was motivated by self-converse symmetry, but that only implies Omega_m = Omega_{p-1-m} not Omega_m = Omega_{p-m}. | kind-pasteur-S50, REFUTED for T_11 kind-pasteur-S50 |
| HYP-437 | ALL eigenspaces of T_p have IDENTICAL Omega dims for k=0,...,p-1 | T_7: [1,3,6,9,9,6,3] same for all 7 eigenspaces. T_11: [1,5,20,70,205,460,700,690,450,180,30] same for ALL 11 eigenspaces (ALL degrees 0-10, confirmed). Follows from Paley multiplicative symmetry. | kind-pasteur-S50 |
| HYP-432 | chi(T_p) = p for Paley primes p≥7 | **PROVED for T_11**: Omega dims per eigenspace [1,5,20,70,205,460,700,690,450,180,30] sum to chi=1 per eigenspace; total chi=11=p. Palindrome NOT needed! T_7 PROVED (palindromic dims). T_11 PROVED (non-palindromic dims, but chi still = 1 per eigenspace). chi(T_p)=p is now confirmed for all computed Paley T_p with p≥7. | kind-pasteur-S50 |
| HYP-433 | T_p Betti: β_{(p+1)/2} = p-1, all other β_k=0 for k>=1 | **SUBSTANTIALLY WRONG for T_11.** T_7: β_4=6=p-1 ✓. T_11: β_6=15 NOT 10. Breakdown: β_6^(k=0)=5 (trivial eigenspace NOT contractible unlike T_7), β_6^(k≠0)=1 (k=1,2 confirmed, k=3..10 predicted by symmetry). Total β_6 = 5 + 10*1 = 15. And β_7 > 0 (predicted 5). CORRECT PATTERN: T_11 has β_0=1, β_6=15, β_7=5, chi=11=p. The "p-1=10" prediction was wrong for k=0 eigenspace. | kind-pasteur-S50 |
| HYP-438 | T_11 full Omega dims per eigenspace: [1,5,20,70,205,460,700,690,450,180,30] | CONFIRMED for k=0 and k=1 (both give identical dims and chi=1). chi from deg 8-10: 450-180+30=300 ✓ (matches the -299 deficit from deg 0-7 to reach chi=1). All dims verified with both PRIME=2^31-1 and PRIME=89. | kind-pasteur-S50 |
| HYP-439 | T_11 palindrome REFUTED: Omega_k ≠ Omega_{10-k} for all k!=5 | Omega_0=1 vs Omega_10=30; Omega_1=5 vs Omega_9=180; Omega_2=20 vs Omega_8=450; Omega_3=70 vs Omega_7=690; Omega_4=205 vs Omega_6=700. Only Omega_5=460=Omega_5 (trivially). Palindrome fails at ALL non-middle degrees. CONFIRMED. | kind-pasteur-S50 |
| HYP-440 | T_11 beta_6 eigenspace structure: β_6^(k=0)=5, β_6^(k≠0)=1 | k=0: rk(d_7)=390, ker(d_6)=395, β_6=5 (large prime). k=1..10: rk(d_7)=390, β_6 verified for ALL k via small prime 89. Total β_6=15. Unlike T_7 where k=0 was contractible, T_11 k=0 has β_6=5 (palindrome failure removes contractibility). CONFIRMED (t11_d7_all_eigenspaces.out). | kind-pasteur-S50 |
| HYP-441 | ~~T_11 predicted full Betti: β_0=1, β_6=15, β_7=5, all others 0~~ | WRONG: beta_7=5 prediction was incorrect. chi argument was garbled. Correct picture: β_5=5 (not β_7), and β_7-10=0. See HYP-443 for the correct verified result. MOVED TO REFUTED. | kind-pasteur-S50 |
| HYP-442 | Paley tournaments only exist for p≡3 (mod 4) | For p≡1 (mod 4): -1 is a QR mod p, so QR∩(-QR)=QR≠∅, making circulant digraph with S=QR have bidirectional arcs (not a tournament). Valid Paley sequence: p=3,7,11,19,23,... (primes ≡3 mod 4). p=5,13,17 are INVALID. Verified computationally. | kind-pasteur-S50 |
| HYP-443 | T_11 FULL BETTI: β=(1,0,0,0,0,5,15,0,0,0,0), chi=11=p | CONFIRMED by direct computation. All boundary ranks verified with small prime 89: rk(d_5^(k=0))=150→ker=310→β_5=5; rk(d_5^(k≠0))=151→ker=309→β_5=0; rk(d_6^(k=0))=305, rk(d_6^(k≠0))=309 (large prime, cross-checked); rk(d_7)=390 all k→ker=300; rk(d_8)=300 all k→β_7=0,ker=150; rk(d_9)=150 all k→β_8=0,ker=30; rk(d_10)=30 all k→β_9=β_10=0 (injective). Total: β_0=1, β_5=5, β_6=15, all others 0. chi=1-5+15=11=p ✓. Scripts: t11_beta5_verify.py, t11_d7_all_eigenspaces.py, t11_higher_betti.py. | kind-pasteur-S50 |
| HYP-444 | Tang-Yau symbol matrix M_m(t) for Paley T_p is CONSTANT (no t-dependence); Q+(QR_p) is EMPTY | Active t-powers = [0] for ALL degrees and ALL p tested (T_7 deg 2-5, T_11 deg 2-5). Q+(QR_p) empty: 0 exceptional t values scanned over all nonzero t in F_89. This is STRONGER than stability theorem: rank(M_m(t)) is the same for ALL t, not just generic. Explains eigenspace identity trivially. | kind-pasteur-S52 |
| HYP-445 | Constant symbol matrix is a UNIVERSAL property of ALL circulant tournaments | PROVED algebraically: face_idx=0 on D∈A_m always gives (d_2,...,d_m)∈A_{m-1}, because partial sums of (d_2,...,d_m) = partial sums of D shifted by -d_1, all distinct. Hence the t-power offset from face_0 is NEVER used. Eigenspace identity (all Omega_m^(k) equal) holds for ALL circulant tournaments, not just Paley. | kind-pasteur-S52 |
| HYP-446 | H-spectrum density → 1; only {7,21} are permanent gaps in odd integers | Exhaustive n<=6: spectrum is all odd except {7,21,35,39} at n=6 (35,39 fill at n=7). n=7: 77/95 density=0.81. Numerically d(N)→1 as N→∞ (confirmed to N=1000). Two proved permanent gaps {7,21}; all other missing values fill eventually. | kind-pasteur-S52 |
| HYP-447 | Omega(T_11) is NOT claw-free — explicit witness found | 793 odd cycles, 310167 conflict edges. Claw center=cycle #0, leaves=cycles #5,#19,#30 (pairwise non-adjacent in Omega, all adjacent to #0). Omega(T_7) IS claw-free (36 cycles, 623 edges, no claw). This confirms claw-free fails at n=11 > threshold n=9. | kind-pasteur-S52 |
| HYP-448 | Eigenspace identity for T_p verified for ALL degrees m=0,...,p-1: T_3, T_7, T_11 | T_7: [1,3,6,9,9,6] same at ALL 7 eigenspaces. T_11: [1,5,20,70,205] same at ALL 11 eigenspaces (up to deg 4). Column-scaling structure verified: C_m^(k) = C_m^(0) * diag(col_scales) at T_7 deg 3 (checked 100%). Proof: QR_p orbits + constant symbol matrix. | kind-pasteur-S52 |
| HYP-434 | chi(T) distribution at n=8: chi=1 (81.75%), chi=0 (17.55%), chi=2 (0.65%), chi=3 (0.05%) | 2000 random tournaments, full complex max_p=7. chi=2 from β_4=1, chi=3 from β_4=2. chi=3 ONLY from near-regular score (3,3,3,3,4,4,4,4) with c3=20 (maximum). No chi<0 or chi>3 seen. | kind-pasteur-S50 |
| HYP-435 | β_4>0 requires c3 >= 16 at n=8 (out of max 20) | 2000 samples: c3 range for β_4=1 is 16-20, most common c3=18,19. β_4=2 only at c3=20. Highly cyclic structure required. Contrast: β_3>0 allows c3 as low as 7. | kind-pasteur-S50 |
| HYP-436 | Unique chi=3 tournament has β_4=2, score (3,3,3,3,4,4,4,4), c3=20 | In 2000 samples, exactly 1 tournament had chi=3. Regular-adjacent score (differs by 1 from regular) with max possible c3=20. Consistent with earlier HYP-324 (β_4=5 at c3=20 from regular score). | kind-pasteur-S50 |

### REFUTED
| ID | Statement | Why it fails | First failure | Source |
|----|-----------|-------------|---------------|--------|
| HYP-441 | T_11 full Betti: β_7=5, β_6=15 | Wrong prediction: β_5=5 (not β_7=5). Higher betti β_7-10=0 by direct computation (rk(d_8)=300=ker(d_7) for all k, full exactness). | T_11 | kind-pasteur-S50 |
| HYP-101 | Per-path identity holds for all n | 3-cycle-only formula misses longer cycles | n=6 | MISTAKE |
| HYP-420 | β_{n-2} generically nonzero at n=8 | Opus used max_deg=6 (missing im(d_7)). With max_deg=7: β_6=0 for ALL 50/50. rk(d_7)=ker(d_6) exactly (top exactness). Artifact of truncated complex. | n=8 (all) | kind-pasteur-S50 |
| HYP-398 | New boundaries target only new cycles | REFUTED at n=8: fails for GOOD vertices (high rate) and ~0.65% of BAD vertices (exactly rank(i_*)=0 cases). Holds universally at n≤7. See HYP-404 for exact equivalence. | n=8 | opus-S55 (stated), kind-pasteur-S49 (refuted), opus-S56 (reconciled) |
| HYP-380 | i_* injective when b3(T)=b3(T\v)=1 | rank(i_*)=0 at n=8 even with b4=0 (trials 613,627,877). Holds n<=7 only | n=8 | kind-pasteur-S47 (stated), kind-pasteur-S48 (refuted) |
| HYP-394 | Consecutive seesaw: beta_k * beta_{k+1} = 0 for ALL k>=1 | beta_3=1 AND beta_4=1 at n=8, confirmed mod-p exact (3/2000). Holds n<=7 only | n=8 | opus-S54 (stated), kind-pasteur-S48 (refuted) |
| HYP-317 | Even Betti numbers vanish: beta_{2k}=0 for k>=1 | beta_4>0 at n=8 (~0.5%), values 1 and 5 | n=8 | kind-pasteur-S45 |
| HYP-319 | Tournament path homology simplicity: one odd hole max | beta_1+beta_5 coexist at n=8, chi=-1 | n=8 | kind-pasteur-S45 |
| HYP-321 | chi(T) in {0,1} for all tournaments | chi=-1 at n=8 (beta_1=beta_5=1), chi=2,6 (beta_4>0) | n=8 | kind-pasteur-S45 |
| HYP-102 | M(T) = M(T^op) | Fails for generic tournaments | n=4 | reversal_proof_attempt.py |
| HYP-361 | V-projection RC proof: mixed TTs span ≤2D in bad-edge subspace V | Mixed TTs span ALL 3D of V (rank=3). Each bad-edge has ≥1 mixed TT. Star-cocycle does NOT constrain to 2D. | n=5 (120/120 fail) | beta2_grok_rc_proof.py |
| HYP-371 | C₁/B₁ quotient bounds #bad: alt sum ∈ B₁ ⟹ dep mod B₁ ⟹ #bad ≤ dim(C₁/B₁) | ALL lifts have rank 0 in C₁/B₁ (every z_v ∈ B₁!). Quotient sees nothing. Lifts are independent WITHIN B₁. Even if it worked, bound would be n-1, not 3. | n=5,6 (exhaustive) | beta2_quotient_rank.py |
| HYP-374 | R = B₁/<Euler, HP 2-cycles, star annihilators> has dim(R)=3 (Grok claim) | Star generators = ALL TTs (every TT in some star), so dim(R)=0 for 97%, =1 for 3%. HP-only: varies 3-8. Coboundaries ⊥ B₁ (dim(B₁∩im(δ))=0). No "3D quotient" exists. | n=5 exhaustive, n=6 (500) | beta2_quotient_R.py |
| HYP-374b | R = B₁/<Euler + ALL HP telescopes> has dim=3 (corrected Grok) | dim(R)∈{0,1,2,4} at n=5. #bad=3→dim(R)=0, #bad=0→dim(R)=2 or 4. INVERTED correlation. | n=5 exhaustive | beta2_euler_hp_quotient.py |
| HYP-374c | i_*-injectivity: ζ_bad→3D ker(i_*), α=±1/√2 forces basis, |bad|≤3 pigeonhole | ker(i_*) per vertex is 0 or 1 dim. No 3D global space. Universal lift span=6=dim(B₁). Residual rank=3 (not 2). This IS HYP-282 restated. | n=5 exhaustive | beta2_istar_ker_test.py |
| HYP-374d | Global LES: ker(δ₃) has rank 3 exact from "Euler+HP-acyclic" | Tested: res: Z₁→⊕H₁(T\v) has rank=#bad (not fixed), ker=dim(Z₁)-#bad (varies: 3-9). coker=0 ALWAYS (surjective). No fixed 3D kernel exists. The "3" is just 6-3=dim(Z₁)-#bad at n=5; at n=6 it's 7. | n=5,6 exhaustive | beta2_global_les_test.py |
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
| HYP-207 | β₂(T) = 0 for ALL tournaments T | 0 counterexamples: exhaustive n≤6, sampled n=7-10 (5000+ tests). Cone-from-T' construction + unfiltered B always solvable | THM-100, beta2_vanishing.py, beta2_filtered_cone.py |
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
| HYP-227 | delta_\|A_3\| = (n-3)*delta_\|A_2\| under arc flip (THM-121) | PROVED algebraically; verified 0 violations n=4-9 | beta2_delta_ratio_proof.py |
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
| HYP-267 | χ(Ω_*) = 1 - β₁ for ALL tournaments | **REFUTED at n=7**: Paley T₇ has χ(Ω)=7, β₁=0, so χ≠1-β₁. Only holds when all β_p=0 for p≥2 (i.e. n≤5). At n≥7, β₃,β₄,β₅ can be nonzero and don't cancel. | beta_paley_verify.py |
| HYP-268 | dim(A_2) = C(n,3) + 2*c₃ for ALL tournaments | CONFIRMED n=4,5 exhaustive. Each transitive triple gives 1 path, each 3-cycle gives 3 | beta2_omega_poincare.py |
| HYP-269 | Ω₂ constraint rank = #{non-allowed (a,c) with mediators} | CONFIRMED n=5. Constraints from different (a,c) pairs use disjoint 2-path sets => automatically linearly independent | beta2_omega_poincare.py |
| HYP-270 | Position cone (all 4 insertion positions) fills ker(d₂) from A₃ | CONFIRMED n=4,5. sum_{v,pos} w_{v,pos} s_v^{pos} = id on ker(d₂\|Ω₂). But image NOT in Ω₃! | beta2_homotopy_proof.py |
| HYP-271 | β₂=0 extended: confirmed n=7 (2000 samples), n=8 (500), n=9 (100). Zero failures. | Extends HYP-255. β₃≠0 at n≥7 (8.2% at n=7, 19.2% at n=8). β₂ is dimension-specific! | beta2_large_n_sample.py |
| HYP-272 | Paley T₇ has β₄=6, Betti=[1,0,0,0,6,0]. Most extreme non-vanishing at n=7 | CONFIRMED. Random n=7 never has β₄≠0 (0/2000). Paley structure creates 4-dimensional holes | beta2_large_n_sample.py |
| HYP-273 | im(d₃\|Ω₃) ≠ im(d₃\|A₃) for most tournaments | CONFIRMED: 904/1024 at n=5, 32048/32768 at n=6. Ω₃ restriction loses most of A₃ image. Yet β₂=0 persists. | beta2_omega3_filling.py |
| HYP-274 | Cone-from-T' fills all swap cycles: w = Σα[(v,a,b,c)+(a,b,c,v)] | Exhaustive n=5,6 (0 failures). Filtered version fails 1/1000 at n=8. UNFILTERED works 500/500 at n=7,8,9. Multi-vertex always works | beta2_filtered_cone.py, beta2_cone_failure.py |
| HYP-275 | Cone filling is AUTOMATICALLY in Ω₃ | CONFIRMED exhaustive n=5 (120/120), n=6 (21120/21120). Breaks at n≥7: 418/425 at n=7, 539/577 at n=8. Single-cone has exposed non-allowed faces that cancel in linear combination | beta2_omega_membership.py |
| HYP-276 | β₂(T) = 0 for ALL tournaments at ALL n | CONFIRMED n≤6 exhaustive, n=7-10 sampled (total ~5000+ tests, 0 failures). Extends HYP-249/271. rank(d₃) = ker(d₂) EXACTLY for every tournament tested. Rank surplus grows with n (min 6→11→15 at n=7→8→9) | beta2_unfiltered_large.py |
| HYP-277 | Unfiltered single-vertex cone: B·alpha=z always solvable for swap cycles | CONFIRMED 500/500 at n=7,8; 200/200 at n=9. Uses ALL T' 2-paths. rank(B)-swap_dim min 6,11,15. Growing surplus means increasingly overdetermined | beta2_unfiltered_large.py |
| HYP-278 | Every tournament has vertex v with b1(T\v) <= b1(T) | PARTIALLY PROVED: 3/4 cases algebraically proved (b1=1; NOT SC; SC+kappa=1). Only kappa>=2+SC+b1=0 remains (empty at n=5, 1680/1680 at n=6, 378/378 sampled n=7). Sampled n=7-12 100%. | beta2_relative_induction.py, beta2_sc_connectivity.py |
| HYP-279 | b1(T) <= 1 for all tournaments | CONFIRMED through n=20 (0 counterexamples). b1=0 rate: 70%(n=5), 95%(n=7), 99.8%(n=9), 100%(n>=10 sampled). When b1=1, ALL vertices are monotone-good | beta2_b1_bound.py |
| HYP-280 | proj_{ker(C)} B_swap has rank = ker_dim exactly | CONFIRMED n=5-8 exhaustive/sampled. Proof: same-endpoint differences from PPQ/PQQ T' paths generate ker(C) cycle space via fundamental cycle decomposition + tournament completeness | beta2_ker_projection.py |
| HYP-281 | Swap cycles from all vertices do NOT span Z_2 | CONFIRMED: 0/1024 at n=5, 0/32768 at n=6. swap_dim << z2_dim. Cone approach fundamentally limited | beta2_z2_decomposition.py |
| HYP-282 | Sum_v b1(T\v) <= 3 for all tournaments | CONFIRMED n=5-10 (exhaustive n=5,6; sampled n=7-10). Bound INDEPENDENT of n. At most 3 bad vertices. Would give >=n-3 good vertices | beta2_sum_b1_bound.py |
| HYP-283 | codim(Z_1(T\v) + Z_1(T\w)) = 1 in Z_1(T) always | CONFIRMED n=5,6,7. The complement L_{vw} visits ALL n vertices. W_v+W_w+W_x = V for ALL triples. ALGEBRAICALLY: follows from dim formula (n-1)(n-2)/2 - n(n-3)/2 = 1 | beta2_codimension_analysis.py |
| HYP-284 | Hidden cycles at bad vertices always linearly independent | CONFIRMED exhaustive n=5 (720 tours), n=6 (27968 tours). rank = #bad always. No deficiency | beta2_hidden_cycles.py |
| HYP-285 | b1=1 implies strongly connected | CONFIRMED exhaustive n=3-6. 0 counterexamples. NOT-SC implies b1=0 | beta2_b1_characterization.py |
| HYP-286 | kappa>=2 + SC + b1=0 always near-regular scores | CONFIRMED: n=6 all (2,2,2,3,3,3); n=7 all near-regular. Regular n=5 all b1=1; regular n=7 70% b1=0 | beta2_kappa2_analysis.py |
| HYP-287 | Bad vertices always form transitive tournament (never 3-cycle) | CONFIRMED: n=5 (0/120 cycles), n=6 (0/4320), n=7 (0/1064 all transitive), n=8 (0/138). |BAD|<=3 and always transitive ordering | beta2_bad_triple_structure.py, beta2_transitive_bad.py |
| HYP-288 | f(v) >= n-3 necessary for v to be bad (f = freed cycle count) | CONFIRMED: n=5-8, 0 violations. Strict separation at n=5,6 (min_bad_f > max_good_f). Overlap at n>=7. | beta2_f_threshold.py |
| HYP-289 | Domination direction in transitive bad triple: top=(above,0), mid=(0,0), bot=(0,below) | CONFIRMED 100% at n=6 kappa>=2 (1440/1440). Top uniquely dominates from above only, bot from below only, mid has no unique domination | beta2_kappa2_good_vertex.py |
| HYP-290 | Coverage(u) < n-1 for some u in every b1=0 tournament | TRUE at n=5,6 (exhaustive). FALSE at n>=7: 17/4746 have all full coverage. Approach insufficient. | beta2_coverage.py |
| HYP-291 | β₃>0 at n=6 iff H = max in score class | CONFIRMED exhaustive 32768/32768. β₃=1 at (1,1,1,4,4,4) H=9 and (2,2,2,3,3,3) H=45 — both are score-class maxima. 0 violations. | beta3_n6_characterization.py |
| HYP-292 | β₄>0 at n=7 iff H=189 (Paley/BIBD class) among regular | CONFIRMED: 6/6 H=189 have β₄=6, 52/52 others β₄=0 (exhaustive regular sampling). β₄>0 EXACTLY picks out H-maximizers within regular tournaments | beta4_characterization.py |
| HYP-293 | β_{n-3}>0 ↔ H-maximizer within score class (general n) | REFUTED. n=6: β₃↔H-max (confirmed, but specific to n=6). n=7: β₃>0 does NOT correlate with H-max (only 5.8% at max!). β₄>0 only for Paley among regular. n=8: β₅=0 for ALL H=661. β₄>0 at high H but not exclusive. The n=6 correspondence is an anomaly, not a general pattern. | beta3_n7_hmax.py, beta5_n8_multi.py, beta4_n8_correlation.py |
| HYP-294 | Paley T_7 deletion has β₃=1 at n=6 | CONFIRMED: T_7-v has scores (2,2,2,3,3,3), H=45, β₃=1. Exactly half of H=45 tours have β₃=1 (240/480); the other 240 have β₁=1 instead. The β₃=1 ones ARE the Paley deletions | beta3_paley_deletion.py |
| HYP-295 | β₁ and β₃ mutually exclusive at n=6 AND n=7 | CONFIRMED exhaustive n=6 (0/32768), sampled n=7 (0/500). Also β₁,β₄ and β₃,β₄ mutually exclusive at n=7. Only 3 Betti vectors at n=7: [1,0,...,0], [1,1,0,...,0], [1,0,0,1,0,...,0] (plus Paley [1,0,0,0,0,0,6]) | betti_n7_full.py |
| HYP-296 | chi ∈ {0,1} for non-Paley tournaments at n ≤ 7 | CONFIRMED: n=5,6 exhaustive, n=7 sampled (500). chi=0 iff β₁=1 or β₃=1. Paley T_7: chi=7=n (the exception). n=8: chi ∈ {1,2} sampled (some H=661 have chi=2 from β₄=1) | betti_n7_full.py, beta5_n8_multi.py |
| HYP-297 | β₃>0 does NOT imply H-max within score class at n=7 | CONFIRMED: only 5.8% (22/378) of β₃>0 at n=7 are at score-class max H. β₃>0 spreads across multiple H values. The n=6 correspondence (HYP-291) is anomalous. | beta3_n7_hmax.py |
| HYP-298 | Three canonical Betti vectors at n=6,7,8: trivial [1,0,...], hole [1,1,0,...], high-β₃ [1,0,0,1,0,...] | CONFIRMED: n=6 (84.4/14.6/1.0%), n=7 (87.8/5.4/6.8%), n=8 (~82/0.7/17.7%). β₁ dominance flips to β₃ dominance at n=8. β₂=0 always. | betti_n7_full.py, beta1_beta3_n8.py |
| HYP-299 | β₁ * β₃ = 0 for ALL tournaments (mutual exclusivity) | CONFIRMED: n=6 (0/32768 exhaustive), n=7 (0/500 sampled), n=8 (0/300 sampled). β₁ and β₃ never both positive. Equivalently chi = 1 - β₁ - β₃ >= 0 for generic tournaments. | beta1_beta3_exclusion.py, beta1_beta3_n8.py |
| HYP-300 | β₁, β₃ each bounded by 1 for all tournaments | CONFIRMED: n=5-8 sampled/exhaustive. max(β₁)=1, max(β₃)=1 always. β₄ can be > 1 (up to 6 for Paley T₇, up to 5 at n=8). | betti_n7_full.py, beta1_beta3_n8.py |
| HYP-301 | β₁>0 requires strong connectivity | CONFIRMED exhaustive n=6: 100% of β₁>0 are SC. β₃>0 does NOT require SC (75% SC at n=6). β₁ detects "1-holes" only in SC tournaments. | beta1_beta3_exclusion.py |
| HYP-302 | Transitive tournament path complex = simplex: dim(Omega_p) = C(n,p+1), chi=1 | CONFIRMED n=3-8. Filling ratio f_p = 1.0 at all p for transitive. All Betti vanish except beta_0=1. | simplex_filling_analysis.py |
| HYP-303 | Filling ratio f_p = dim(Omega_p)/C(n,p+1) > 1 for p >= 3 at n >= 6 | CONFIRMED: f_3 = 1.048 (n=6), 1.157 (n=7), 1.283 (n=8). Cyclic content inflates beyond simplex. | filling_ratio_formula.py |
| HYP-304 | H(T_4) = 2*c3 + 1 for ALL 4-vertex tournaments | CONFIRMED exhaustive (64/64). Clean identity: H(transitive)=1, H(cyclic)=3, H(regular)=5. | chi_A_identity.py |
| HYP-305 | excess_4 = 2*c3*(n-3) for ALL tournaments on n vertices | CONFIRMED n=4-7 exhaustive/sampled. Each 3-cycle in (n-3) 4-subsets; H(T_4) identity gives 2 per cycle. | chi_A_identity.py |
| HYP-306 | chi_A = sum(-1)^p \|A_p\| is ALWAYS ODD | CONFIRMED n=3-7. All tournaments at all n have odd chi_A. | chi_A_identity.py |
| HYP-307 | chi_A = 1 at n=4 for ALL tournaments (tournament-independent) | CONFIRMED exhaustive (64/64). Perfect cancellation: 2*c3 - 2*c3 = 0. | chi_A_identity.py |
| HYP-308 | dim(Omega_2) is NOT determined by c3 alone at n >= 5 | CONFIRMED: (c3=4, score=(1,2,2,2,3)) at n=5 gives Omega_2 in {8,9,10}. Geometric arrangement matters. | omega2_formula.py |
| HYP-309 | dim(Omega_2) is NOT determined by (c3, score) at n >= 5 | CONFIRMED: same (c3, score) at n=5,6 gives different Omega_2 values. Path complex encodes more than cycle/score data. | omega2_formula.py |
| HYP-310 | Dimensional crossover: P(beta_1>0) peaks then decays; P(beta_3>0) grows monotonically from d=5 | CONFIRMED: P(b1>0) = 25%, 37.5%, 25%, 15.8%, 5.4%, 1.7% for n=3-8. P(b3>0) = 0, 0.98%, 8.7%, 18.7% for n=5-8. Crossover at d=6. | dimensional_crossover.py |
| HYP-311 | beta_3 fragility: beta_3>0 disappears under ANY vertex deletion at n=6 | CONFIRMED exhaustive. All 320 beta_3>0 tournaments at n=6 have trivial betti for all 6 deletions. beta_1 can survive deletion (robust). | simplex_face_restriction.py |
| HYP-312 | chi(T) in {0,1} for ALL tournaments at n <= 7 | CONFIRMED: n=3-7 exhaustive/sampled. chi=0 iff beta_1=1 or beta_3=1. Breaks at n=8 (chi up to 6 from beta_4). | euler_char_scaling.py |
| HYP-313 | Non-SC tournaments always have chi=1 at n <= 6 | CONFIRMED exhaustive. At n=7, first non-SC with chi=0 appears. | euler_char_scaling.py |
| HYP-314 | surplus = excess_paths - rank(constraints) exactly for all p | CONFIRMED n=5-7 for p=2,3. Identity is algebraic: dim(Omega_p) = |A_p| - rank(P_p). | local_redei_investigation.py |
| HYP-315 | corr(surplus, excess_H) increases with n at fixed p | CONFIRMED: p=2 correlation goes 0.53 (n=6) -> 0.87 (n=7). Cyclic content increasingly predicts Omega inflation. | local_redei_investigation.py |
| HYP-329 | dim(Z₁) = C(n,2)-(n-1) universal for all tournaments | CONFIRMED exhaustive n=5,6,7. rank(∂₁)=n-1 always (complete underlying graph). β₁ depends entirely on rank(∂₂). Note: equivalent to HYP-318 | beta2_rank_critical.py |
| HYP-330 | rank(∂₂\|Ω₂) = C(n,2)-n+1-β₁ | CONFIRMED exhaustive n=5 (1024/1024), n=6 (32768/32768). Zero violations. | beta2_rank_critical.py |
| HYP-331 | dim(ker(∂₂^T)) = n-1+β₁ | CONFIRMED exhaustive n=5,6. Clean dual formula. | beta2_bad_tt_basis.py |
| HYP-332 | Bad-vertex TT is rank-critical (removing drops rank(∂₂)) | CONFIRMED exhaustive n=5,6: n=5 (120/120), n=6 (4320/4320). At n=5: zero redundancy (all TTs RC). At n=6: redundancy∈{2,3} but bad-TT STILL RC. Bad verts always form TT, never 3-cycle. | beta2_badtt_rc_n6.py |
| HYP-333 | #RC=0 when redundancy (#TTs-rank) ≥ threshold | CONFIRMED: threshold=3 at n=5, =8 at n=6. Note: #bad=3 at n=5→redundancy=0, at n=6→redundancy∈{2,3} (nonzero). | beta2_rank_critical.py |
| HYP-334 | Flip obstruction: 3-cycle among ALL bad vertices forces β₁≥1 | CONFIRMED: 0 counterexamples at n=5 (720), n=6 (27968), n=7 (479). Flipping bad TT→3-cycle: 100% force β₁=1 | beta2_flip_deep.py |
| HYP-335 | Minimum #bad when β₁=1 grows with n | CONFIRMED: min_bad=3(n=5), 4(n=6), 5(n=7). β₁=1 bad sub-tournament often SC (55% at n=6) | beta2_flip_deep.py |

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
HYP-207, HYP-208, HYP-209, HYP-210, HYP-211, HYP-212, HYP-213, HYP-214, HYP-215, HYP-216, HYP-217, HYP-218, HYP-219, HYP-220, HYP-221, HYP-222, HYP-223, HYP-224, HYP-225, HYP-226, HYP-227, HYP-228, HYP-229, HYP-230, HYP-231, HYP-232, HYP-233, HYP-234, HYP-235, HYP-236, HYP-237, HYP-238, HYP-239, HYP-240, HYP-241, HYP-242, HYP-243, HYP-244, HYP-245, HYP-246, HYP-247, HYP-248, HYP-249, HYP-250, HYP-257, HYP-258, HYP-259, HYP-260, HYP-261, HYP-262, HYP-263, HYP-278, HYP-279, HYP-282, HYP-287, HYP-288, HYP-289, HYP-290, HYP-291, HYP-292, HYP-293, HYP-294, HYP-295, HYP-296, HYP-297, HYP-316, HYP-318, HYP-320, HYP-322, HYP-323, HYP-324, HYP-325, HYP-326, HYP-327, HYP-328, HYP-329, HYP-330, HYP-331, HYP-332, HYP-333, HYP-334, HYP-335, HYP-349, HYP-350, HYP-351, HYP-352, HYP-353, HYP-354, HYP-355

### Self-complementary / blueself
HYP-009, HYP-010, HYP-112

### Worpitzky / Forward-edge polynomial
HYP-011, HYP-012, HYP-119, HYP-120, HYP-204, HYP-205, HYP-206

### Spectral
HYP-201, HYP-202

### Dimensional meta-patterns (simplex perspective)
HYP-302, HYP-303, HYP-304, HYP-305, HYP-306, HYP-307, HYP-308, HYP-309, HYP-310, HYP-311, HYP-312, HYP-313, HYP-314, HYP-315

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
