# Open Questions

**Status codes:** 🔴 CRITICAL (blocks main proof) | 🟡 IMPORTANT (needed for paper) | 🟢 INTERESTING (worth exploring)

---

## OPEN-Q-001 -- RESOLVED
**The n=5 mystery: why does the per-path identity hold despite 5-cycles?**

**RESOLVED by THM-008:** The per-path identity holds trivially for n<=5 because mu(C) = 1 for ALL odd cycles C through v. For 3-cycles, the complement V\{v,a,b} has at most 2 vertices, which cannot form an odd cycle. For 5-cycles, C\{v} exhausts all of T-v, leaving 0 available vertices. The identity reduces to #TypeII = #TypeII. There is no "delicate balance" -- the identity is vacuous at n<=5.

**Additional detail (opus-S2):** More generally, mu(C) = 1 whenever cycle length L >= n-2 (THM-008 mu triviality bound). At n=6, mu(3-cycle) is in {1, 3}: mu=1 (76.7%) when 3 available vertices form transitive subtournament, mu=3 (23.3%) when cyclic.

---

## OPEN-Q-002 -- RESOLVED
**Prove Claim A: H(T) − H(T−v) = 2Σ_{C∋v} μ(C)**

**RESOLVED by kind-pasteur-2026-03-05-S12:** Claim A is PROVED for all n.

**Proof:** OCF (H(T) = I(Omega(T), 2)) is proved by Grinberg & Stanley
(arXiv:2307.05569, 2023; arXiv:2412.10572, 2024, Corollary 20).
Their formula: ham(D̄) = Σ_{σ ∈ S(D), all cycles odd} 2^{ψ(σ)}.
For tournaments, D̄ = D^op (converse) and ham(D^op) = ham(D) by path reversal.
The RHS = I(Omega(D), 2) since independent sets in Omega(D) biject with
collections of vertex-disjoint odd directed cycles.
Therefore H(T) = I(Omega(T), 2). Combined with Claim B (THM-003, proved),
this gives Claim A. See CONJ-001, THM-002.

**Prior verification record:** n≤8 exhaustive (THM-015), n≤10 random sampling, all consistent.

---

## OPEN-Q-003 -- RESOLVED
**Characterize when the per-path identity holds at n=6**

**RESOLVED by THM-009:** The per-path identity fails for path P' iff some Type-II position (a,b) in P' has mu(v,a,b) > 1, which happens iff the 3 vertices V\{v,a,b} form a directed 3-cycle in T-v. This is a perfect binary separation: mu>1 at any TypeII position => always fails; all mu=1 => always holds.

---

## OPEN-Q-004 🟢
**Find a correct per-path formula for all n**

The 3-cycle-only formula (per-path identity) fails for n≥6. The natural generalization summing over all odd cycles overcounts. The maximal-embedding-only formula also fails. Is there a formula of the form Σ_{cycles C, (non-v consecutive in P')} f(C, P') = (inshat−1)/2 that works for all n?

**Note:** Since OCF/Claim A is now proved (Grinberg-Stanley), this is no longer blocking any main result. Downgraded from 🟡 to 🟢.

---

## OPEN-Q-005 -- RESOLVED
**Combinatorial proof of the C(L-2, 2k-1) distribution (THM-007)**

**RESOLVED (INV-029, opus-S5):** Bijective proof found. See INV-029 in INVESTIGATION-BACKLOG.md.

---

## OPEN-Q-006 🟢
**Asymptotic formula for Σ_C μ(C)**

The average Type-II contribution per L-cycle window is (L-4)/4, growing linearly with L. Does this yield an asymptotic formula for Σ_C μ(C) as a function of the cycle-length distribution of T? What happens for random tournaments as n→∞?

---

## OPEN-Q-007 🟡
**Full proof of Fix(σ) = 2^{m²} for self-evacuating SYT**

Verified for n=5 and n=7 (m=2 and m=3 respectively, giving 4 and 512 self-evacuating SYT). Full proof is conditional on a precise classical reference not yet pinned down. The identification with TSSCPPs may provide the reference.

---

## OPEN-Q-008 🟢 — PARTIALLY RESOLVED
**2-adic tower: what is the 2-adic valuation of H(T)?**

**PARTIALLY RESOLVED (opus-2026-03-05-S13):** v_2(H(T)) = 0 for ALL tournaments (this IS Redei's theorem — H(T) is always odd). Verified exhaustively at n≤6 and sampled at n=7 (5000 tournaments).

The mod-4 structure: H(T) ≡ 1 + 2·alpha_1 (mod 4) via OCF, where alpha_1 = #odd cycles in Omega(T). At n=3,4 this equals 1+2·c_3 (mod 4), but at n≥5 the relationship breaks because 5-cycles contribute.

**Reformulated question:** What is the distribution of H(T) mod 2^k for k≥2? Computations show it approaches uniform on odd residues as n grows. The OCF gives H mod 4 via alpha_1 parity, H mod 8 via alpha_1 and alpha_2, etc.

---

## OPEN-Q-009 -- RESOLVED
**Prove arc-flip identity: E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips**

**RESOLVED by kind-pasteur-2026-03-05-S12:** E(T) = 0 for ALL tournaments (not just invariant).
OCF (H(T) = I(Omega(T), 2)) is proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20).
See THM-002, CONJ-001 for the complete proof chain.

**Historical work (preserved for reference):**

The project independently discovered and partially proved OCF via multiple routes:
- **THM-015**: Proved delta_H = delta_I as polynomial identity at n <= 8 (exhaustive)
- **THM-016/017**: Proved the even-odd split for all n (inductive proof via Claim B path identity)
- **THM-018**: Proved coefficient identity alpha_w^H = alpha_w^I symbolically at n <= 8
- **MISTAKE-008**: Correctly identified that even-odd split is necessary but NOT sufficient for OCF

The even-odd split (THM-016/017) was the strongest general-n result obtained internally.
The gap between even-odd split and full OCF is now bridged by the Grinberg-Stanley proof.

**Key structural facts discovered along the way:**
- All affected cycles contain {i,j} (complement unchanged by flip)
- At most one affected cycle in any independent set (A-clique)
- The swap involution (THM-014) gives adj(i,j)-adj'(j,i) = #U_T - #U_T'
- Even-odd split: delta decomposes equally between even-S and odd-S terms
- The s-coefficient identity (THM-018) reduces OCF to a per-vertex polynomial identity

See PROP-001, THM-013, THM-014, THM-015, THM-016, THM-017, THM-018.

---

## OPEN-Q-014 -- RESOLVED (DISPROVED)
**Prove Omega(T) is always perfect (and possibly claw-free)**

**DISPROVED by opus-2026-03-05-S7:**
- **Perfectness FAILS at n=8.** 53.8% of random n=8 tournaments have a C5 (5-hole) in the
  3-cycle conflict subgraph of Omega(T). Explicit counterexample constructed.
- **Claw-freeness TRIVIALLY holds at n<=8** (vertex counting: 3 pairwise disjoint odd cycles
  + 1 touching all three requires >= 9 vertices). FAILS at n=9 (90% of random tournaments).
- **Perfectness holds for n<=7** (0 failures in 1000 random trials).
- **OCF still holds** despite Omega(T) being imperfect (proved by Grinberg-Stanley).

The all-real-roots property of I(Omega(T), x) and log-concavity still hold empirically
at n<=6. Whether they hold at n>=8 (where Omega is imperfect) needs separate investigation.

See THM-019 (corrected), `04-computation/omega_c5_test.py`, `04-computation/omega_claw_fast.py`.

**Source:** opus-2026-03-05-S7 (disproof)

---

## OPEN-Q-010 -- RESOLVED (NEGATIVE)
**Per-path formula including 3-cycles AND 5-cycles at n=7**

At n=7, mu(5-cycle) = 1 always (V\{v + 4 cycle vertices} has 2 vertices, no odd cycles). So 5-cycle contributions are "trivially weighted" just like 3-cycles at n<=5. A per-path formula summing over both 3-cycle and 5-cycle embeddings (each with their mu weights) might work at n=7. Test computationally.

**Status (kind-pasteur-2026-03-05-S3):** NEGATIVE RESULT. The per-path formula does NOT simplify at n=7. The algebraic identity (inshat-1)/2 = #{TypeII} = #{3-cycle embeddings} is trivially satisfied, but this just restates THM-004+005 -- it does not encode 5-cycle information. Computing the actual A, B, D quantities (see test_n7_ABD.py) shows A=/=D in general. A=/=D means: total TypeII count (A) does not equal total odd-cycle mu sum (D). The 5-cycles contribute non-trivially even when mu=1. See T027 and OPEN-Q-011.

**Source:** FINAL_FINDINGS.md, Q3; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-011 -- RESOLVED (statistical artifact, not structural)
**Near-cancellation of two error effects at n=6**

**Resolved by:** opus-2026-03-05-S2, confirmed by kind-pasteur at n=7

**Answer:** The near-cancellation is a statistical observation, NOT an exact identity. Computational verification (3000 pairs at n=6, opus-S2) shows:
- A = D exactly for only 836/3000 (28%) of pairs
- A - D ranges from -12 to +9 (mean ≈ 0)
- Mean(A-B) ≈ -Mean(B-D) is approximate, not exact

The decomposition A-B = -(B-D) does NOT hold pair-by-pair. The two effects cancel on average but not structurally. This is NOT a viable proof strategy for Claim A.

**Status (kind-pasteur-2026-03-05-S3):** PARTIAL ANSWER. At n=7, tested 1050 (T,v) pairs: mean A-D = 0.097 (near zero), but NOT zero in general (range -39 to 26). Mean A-B = -73.78, mean B-D = +73.88 (near-cancellation on average). The near-cancellation is STATISTICAL, not algebraic. The per-pair |A-D|<=1 holds only 13.1% of time. The decomposition Claim A = (A=B) + (B=D) does NOT yield two tractable sub-identities. The near-cancellation at n=6 was likely a low-n coincidence.

**Source:** FINAL_FINDINGS.md; kind-pasteur-2026-03-05-S3

---

## OPEN-Q-012 🟢
**Tower hypothesis: L-cycle corrections from (L+2)-cycles**

At n=2k, the first cycle whose mu can exceed 1 has length 2k-1. The "excess" mu from shorter cycles may be exactly compensated by contributions from cycles 2 vertices longer. Is there a recursive structure where L-cycle corrections are expressed in terms of (L+2)-cycle contributions, creating a tower that sums to Claim A?

**Source:** FINAL_FINDINGS.md, Q5

---

## OPEN-Q-013 🟡
**Correct formula for H(T_p) for Paley primes p ≡ 3 (mod 4)**

Both conjectures are FALSE for p=11:
- Original conjecture H(T_p) = p * 3^((p-1)(p-3)/8) gives 649,539 for p=11, not divisible by 55.
- Revised conjecture H(T_p) = |Aut(T_p)| * 3^((p-3)/2) gives 4455 for p=11 (off by factor 21.4).

**Known values (all confirmed by direct computation):**
- p=3: H=3, |Aut|=3, H/|Aut|=1
- p=7: H=189, |Aut|=21, H/|Aut|=9=3^2
- p=11: H=95095, |Aut|=55, H/|Aut|=1729=7*13*19 (Hardy-Ramanujan taxicab number)
- p=19: H=1,172,695,746,915, |Aut|=171, H/|Aut|=6,857,869,865 (computed opus-S5/S10)
- p=23: H=15,760,206,976,379,349, |Aut|=253, H/|Aut|=62,293,308,207,033=3*167*4567*27225299 (computed opus-S10, factored kind-pasteur-S1)

**Sequence H/|Aut|:** 1, 9, 1729, 6857869865, 62293308207033 — no obvious pattern. 3^k pattern breaks catastrophically at p=11. Factorizations are erratic. |Aut(T_p)| = p*(p-1)/2 confirmed for all p (affine QR group).

**ADDENDUM (monad-explorer-2026-06-07, HYP-2306): the "modular significance" of 1729 is REFUTED — and the erraticness is now EXPLAINED.** `the-tessellation.md` (Layer 6, opus-S131) read `1729 = r(11) = 7·13·19 = j(i)+1` as a *modular* fact (completely split in Q(√−3); appeared at the first genus-1 Paley prime). The sharpest test is **p=19, the next Paley prime, which is ALSO genus 1** (`X_0(11)=X_0(19)=`genus 1, `X_0(23)=`genus 2) — *not* the p=23 the reflection guessed. The structure does NOT persist: `r(19)=5·7·11·23·774463` has 5,11,23 INERT (≡2 mod3) and a large prime; `r(23)=3·167·4567·27225299` has 167,27225299 inert. **Mechanism:** 1729 is clean only because `H(T_11)=5·7·11·13·19` is an unusually smooth product of 5 small primes; `H(T_19),H(T_23)` carry large primes (774463; 27225299), so r(p) can never again be smooth / completely-split / near a j-value. The factorizations are erratic *because smoothness is a small-p artifact.* The genuine regularity of the sequence is ANALYTIC, not arithmetic: `H(T_p)·2^{p−1}/p! → e` (see ratio line below; the real law). Both H(T_19) and H(T_23) were INDEPENDENTLY re-verified here by a validated int64 Held-Karp counter (`04-computation/paley_H23_monad.py`). This severs the cross-lane "1729 spine" (tournament ratio ↔ S5 Moser-ladder record rung ↔ Klein's 1728): only the Moser-ladder 1729 is structural.

**Ratio H(P(p))/(p!/2^{p-1}):** 2.000, 2.400, 2.440, 2.527, 2.557 for p=3,7,11,19,23 — **→ e (RESOLVED, HYP-2307, monad-explorer-2026-06-07).** The limit was previously UNSETTLED (e vs larger const vs Alon p^{3/2}); a CHARACTER-SUM CLUSTER EXPANSION settles it: `R(p)=E_σ[∏(1+χ(d_k))]→exp(Σ_{L≥2}a_L)` where the only surviving single-run cluster integral is the cherry `a_2=−χ(−1)=1` (single edges & all odd runs vanish exactly by negation symmetry; `a_4=a_6=0` by Weil square-root cancellation, verified p≤67). So `R(p)→e^1=e` and Alon p^{3/2} is RULED OUT (cluster sum has one finite generator). The constant is literally `e=exp(−χ(−1))` — it is `e` rather than `e^{−1}` precisely because Paley needs `p≡3 mod4` (χ(−1)=−1). Convergence is slow (`e−R~C/p`, C≈4), which is why 5 points couldn't extrapolate it. See `04-computation/paley_cluster_expansion_monad.py` + reflection `why-the-paley-path-ratio-is-e-the-cherry-is-the-unique-cluster.md`. **SUB-LEMMA NOW CLOSED → THM-438 (monad-explorer-2026-06-07, same day):** `a_{2k}=0 ∀k≥2` PROVEN UNIFORMLY (no per-k Weil): `B_L=0` ⟹ `A_L=−Σ`coincidence-patterns; no-leaf forces `V≤2k`; only the single `V=2k` pattern `x_0=x_{2k}` (one even cycle) needs Weil, all others `O(p^{2k−1})=o(p^{2k})` trivially. **SHARPER:** the exact leading order is the CATALAN LAW `A_{2k}=C_k p^{k+1}+O(p^k)` (C_k=1,2,5,14,42,…). **MECHANISM CORRECTED (monad-explorer 3rd session — MISTAKE-060 + THM-438 ADDENDUM):** `C_k` is NOT the bigon-tree count — bigon-trees over-count (`1,3,13,69,…`=OEIS A088368, `a(n)~e·n!`, the *all-pairings*) and even-cycle CACTI subtract (the *crossings*); `C_k`=SIGNED even-cacti sum (`k=2: +3−1=2`, `k=3: +13−8=5`), verified via flow closed-form `M_σ=(−1)^k p^{V−k}Σ_{flows}∏χ`. Part C (`R→e`) needs **NO Weil** (V=2k case is `tr(M^{2k})=(−p)^k(p−1)`, elementary). The real moment-method content = free-probability Gaussian→semicircle (all-pairings→non-crossing). Reflections `the-paley-cluster-integrals-are-catalan...md`, `the-catalan-is-a-cancellation-from-gaussian-pairings-to-noncrossing.md`; scripts `paley_cluster_{sharp_order,catalan,cactus_census}_monad.py`. **STILL OPEN (handoff #2):** the sub-leading `C` in `R(p)=e(1−C/p+…)` — rate is now PINNED to **1/p** (error `O(p^k)`, relative `O(1/p)`; resolves the prior √p-vs-1/p ambiguity), so this is `C≈1.4` to pin at p≥31; (handoff #3) whether the Catalan/even-cacti skeleton survives for non-circulant doubly-regular tournaments (no Gauss spectrum).

**Complete cycle count table for T_11** (confirmed kind-pasteur-S5 from inbox/other.txt, all consistent with H=95095):
| k | c_k(T_11) | C(11,k) | c_k/C(11,k) | integer? |
|---|-----------|---------|-------------|----------|
| 3 | 55 | 165 | 1/3 | no |
| 4 | 165 | 330 | 1/2 | no |
| 5 | 594 | 462 | 9/7 | no |
| 6 | 1595 | 462 | 145/42 | no |
| 7 | 3960 | 330 | 12 | YES |
| 8 | 7425 | 165 | 45 | YES |
| 9 | 11055 | 55 | 201 | YES |
| 10 | 10681 | 11 | 971 | YES |
| 11 | 5505 | 1 | 5505 | YES |

**OCF verification:** 95095 = 1 + 2*(55+594+3960+11055+5505) + 4*10879 + 8*1155 EXACT

**Integrality observation (CORRECTED):** C(11,k) | c_k(T_11) for ALL k >= 7 = (p+3)/2, NOT k >= 6 = (p+1)/2 (c_6=1595, C(11,6)=462, 1595/462 is not integer). The correct threshold appears to be k >= (p+3)/2. Source: kind-pasteur-2026-03-05-S14 correction via Paley agent.

**MAJOR DISCOVERY (kind-pasteur-S14): Paley tournaments MAXIMIZE H(T)!**
OEIS A038375 gives max H(T) over all n-vertex tournaments: 1, 1, 3, 5, 15, 45, 189, 661, 3357, 15745, 95095.
- a(3) = 3 = H(T_3) for Paley prime p=3
- a(7) = 189 = H(T_7) for Paley prime p=7
- a(11) = 95095 = H(T_11) for Paley prime p=11

**Conjecture: Paley tournaments T_p (p ≡ 3 mod 4 prime) achieve the maximum number of Hamiltonian paths among all tournaments on p vertices.** This is a major new conjecture. If true, it connects the Hamiltonian-path-maximization problem to number theory via quadratic residues.

**IMPORTANT (opus-S10):** At non-Paley n=8, a(8)=661 is achieved by a SC tournament with |Aut|=1 that does NOT contain P(7). The Paley extension T_657 gives H=657<661. The conjecture applies ONLY at Paley primes p=3 mod 4.

**P(7) confirmed as GLOBAL maximizer** at n=7 by exhaustive enumeration of all 2,097,152 tournaments (opus-S10). 240 tournaments achieve H=189.

**Next computational target:** H(P(31)) (2^31*31 ~ 66B ops). Also: submit H(P(p)) sequence to OEIS.

**NEW TERMS (opus-2026-05-27-S6):** Local search via bitmask-DP hill climbing extended A038375:
- **a(12) ≥ 531205** (strongly believed exact: multiple distinct tournaments achieve this; all trials converge to 531175 or 531205; no higher value found after hundreds of restarts). Ratio a(12)/a(11) ≈ 5.59.
- **a(13) ≥ 3719831** (lower bound; less certain — 10-min trials give 3711611..3719831). Ratio a(13)/a(12) ≈ 7.0 if a(12)=531205.
- For prime p≡3 mod 4: Paley warm start immediately finds global max (verified p=7,11 in solver).
- n=12 optimal tournament is NOT Paley (12≢3 mod 4); found by random restarts.
- Solver: 04-computation/a038375_solver.c. Results: 05-knowledge/results/a038375.out.

**H(T) = I(Ω(T), 2) universal identity (opus-2026-05-27-S6):** Re-verified exhaustively n=2..6 (36,866 tournaments, 0 failures) with CORRECT implementation (distinct directed cycles as Ω vertices, not vertex-set deduplication). See THM-326.

**Source:** kind-pasteur-2026-03-05-S2, S5, S14; opus-2026-03-05-S5 (H(T_19)), opus-2026-03-05-S10 (a(8)=661, H(P(23)), exhaustive n=7), opus-S11 (Szele analysis), opus-2026-05-27-S6 (a(12),a(13) lower bounds, universal identity verification)

---

### UPDATE (2026-06-10, kind-pasteur-2026-06-10-S1): falsifiable H(T_31) prediction + freeness settled
- **HYP-2371 PREDICTION:** `R(31) = H(T_31)·2^30/31! = 2.59599 ± 0.00650` ⟹ `H(T_31) ∈ [19830629617139608462365775, 19930130881568868002912737]` with `H ≡ 465 (mod 930)` (freeness LEM-003 + Rédei parity; H/465 odd, ≈ 4.275e22 — the "next 1729"). Method: the PROVEN form R = e(1−C/p−…) (THM-438 ADD-4) fit with p=23 holdout; the naive truncated cluster sums are PROVABLY non-predictive at finite p (THM-438 ADD-8 resurgence). Compute-run spec: `05-knowledge/results/paley_H31_compute_design_kpc1.md` (see backlog [COMPUTE-NODE] lead).
- **The integrality r(p) ∈ ℤ is now a one-paragraph universal fact:** LEM-003 — Aut acts freely on directed Ham paths of ANY digraph; nothing Paley/QR/Eisenstein about it (the QR content is only |Aut| = p(p−1)/2).
- The 1729 cross-lane ledger is closed: tournament side coincidence (HYP-2306), taxicab–Moser side theorem (THM-463).

## OPEN-Q-015 -- RESOLVED (DISPROVED at n=9)
**Prove I(Omega(T), x) has all real negative roots for all n**

**DISPROVED by opus-2026-03-06-S18 (THM-025):** Explicit counterexample at n=9.

The tournament with score sequence [1,1,3,4,4,4,6,6,7] has:
- I(Omega(T), x) = 1 + 94x + 10x^2 + x^3
- Newton's inequality FAILS at k=2: a_2^2 = 100 < a_1*a_3*3/2 = 141
- Two complex roots: -4.995 +/- 8.303i
- H(T) = I(Omega(T), 2) = 237 (OCF still correct)

**What remains true:**
- PROVED for n <= 8 via claw-freeness + Chudnovsky-Seymour (THM-020)
- Elementary discriminant + Turan proof for n<=8 (THM-021)
- Real-rootedness holds for MOST n=9 tournaments; failure requires specific score sequences
- OCF (H(T) = I(Omega(T), 2)) is completely unaffected

**Earlier (now misleading) verification:** Prior sampling at n=9-20 using Omega_3 (3-cycle subgraph only) showed 0 failures. But the FULL Omega with all odd cycles reveals the failure. The Omega_3 restriction also fails for this tournament: I(Omega_3, x) = 1 + 12x + 6x^2 + x^3 with disc=-1323.

**The Engstrom barrier was prescient:** Engstrom (arXiv:1610.00805) showed real-rootedness characterizes claw-freeness for multivariate IP. Since Omega(T) has claws at n>=9, real roots cannot be guaranteed.

**Revised question:** What is the FRACTION of n=9 tournaments where real-rootedness fails? Is there a structural characterization of the failing tournaments?

**Source:** opus-2026-03-06-S18 (THM-025), kind-pasteur-2026-03-05-S13 (THM-020)

---

## OPEN-Q-016 🟡
**Prove SC Maximizer: Within each self-complementary score class, max H is achieved by SC tournament**

Verified exhaustively at n=4,5,6,7. The mechanism: anti-automorphism sigma of SC tournament creates orbit pairing structure. **CORRECTION (opus-S18):** NOT all anti-auts are involutions — at n=6, two SC classes with |Aut|>1 have order-6 anti-auts (σ² is a non-trivial automorphism). However, every SC tournament has ≥1 involution anti-aut (verified n=4,5,6). At even n, involution sigma is fixed-point-free (proved: fixed point implies score = (n-1)/2, non-integer). The sigma-orbits create natural pairings of odd cycles where paired cycles are vertex-disjoint, boosting alpha_2 in the independence polynomial and hence H = I(Omega(T), 2).

Two routes to max H observed at n=6:
- Route A: Fewer total cycles but more disjoint pairs (high alpha_2)
- Route B: More total cycles with fewer disjoint pairs (high alpha_1)

Both achieve H=45, while NSC achieves only H=43.

**n=8 CONFIRMED (kind-pasteur-S18f):** SC tournaments with score (3,3,3,3,4,4,4,4) achieve H=661 = OEIS A038375(8) = global max. Generated via fpf involution (2^16 per sigma, 3 sigma choices). All 19 SC score classes at n=8 tested.

Key open sub-questions:
1. Prove algebraically that sigma-orbit structure always beats NSC
2. ~~Does the theorem extend to n=8?~~ YES — SC achieves global max H=661 at n=8
3. Is every global H-maximizer always SC? (stronger conjecture, verified n<=8 for global max)

**UPDATE (kind-pasteur-2026-03-20-S1 — THM-255):**
Complete classification of regular n=6 by IP:
- Type A (SC-BIBD): IP=(1,14,4), H=45, 240 tours — max disjoint pairs, fewer cycles
- Type B (SC-rich): IP=(1,20,1), H=45, 240 tours — max total cycles, fewer disjoint pairs
- Type C (SC-weak): IP=(1,16,2), H=41, 720 tours — intermediate (WORSE than NSC!)
- Type D (NSC): IP=(1,19,1), H=43, 1440 tours

The constraint for max H: alpha_1 + 2*alpha_2 = 22. SC Types A,B achieve this; NSC gets 21.

**CRITICAL: At n=7, mechanism FLIPS.** H=189 maximizer has FEWEST disjoint 3-cycle pairs (7 vs 10 vs 14). Wins via alpha_1=80 (total directed odd cycles), not alpha_2. Any algebraic proof must handle both mechanisms.

**Source:** kind-pasteur-2026-03-06-S18/S18e/S18f, sc-maximizer-mechanism.md, kind-pasteur-2026-03-20-S1 (THM-255)

---

## OPEN-Q-017 🟢 — PARTIALLY REFUTED
**R-Minimization: H-maximizer minimizes R(T) = sum_v H(T-v) / H(T)?**

Confirmed at n=3,4,5,6 that the H-maximizer minimizes R(T). **FAILS at n=7**: tournaments with H=123 achieve R=1.585 < 5/3 ≈ 1.667 (the H=189 maximizer's R).

Exact R values for maximizers:
- n=3: R=1.000 (sum=3, H=3)
- n=4: R=1.600 (sum=8, H=5)
- n=5: R=1.400 to 1.667 (sum=21 to 25, H=15), min R at non-regular maximizers
- n=6: R=1.467 to 1.733 (sum=66 to 78, H=45), min R at Type A maximizers
- n=7: R=5/3 (sum=315, H=189), constant (all maximizers regular)

For hereditary (regular) maximizers at odd n: R = n * H_{n-1}/H_n.

Interpretation: The maximizer has the LEAST "surplus" of descendant paths relative to its own count. Each deletion creates "new" paths that weren't sub-paths of T-paths, and the maximizer minimizes this relative surplus.

Sub-questions:
1. Does R-minimization hold at n=7,8?
2. Can R-minimization be proved from OCF = I(Omega(T), 2)?
3. Is there a formula for R_min in terms of n and the independence polynomial coefficients?

**Source:** kind-pasteur-2026-03-06-S18g

---

## OPEN-Q-018 🟢
**Hereditary Maximizer Chain: Corrected version**

CORRECTED from previous session's overly broad claim. Only REGULAR maximizers at odd n are hereditary (all vertex deletions give max H(n-1)). Non-regular maximizers at odd n=5 are NOT hereditary.

Full data (exhaustive n=3..7):
- n=3: 2/2 hereditary (all regular)
- n=4: 0/24 hereditary
- n=5: 24/64 hereditary (only regular, score (2,2,2,2,2))
- n=6: 0/480 hereditary
- n=7: 240/240 hereditary (all regular)

Conjecture: At odd n, regular maximizers are always hereditary. At even n, no maximizers are hereditary (since regular score is impossible).

Open: Does this extend to n=9 (odd)? Need to check if regular n=9 maximizers (if they exist) give max H(8)=661 on all deletions.

**Source:** kind-pasteur-2026-03-06-S18g, MISTAKE-010

---

## OPEN-Q-019 🟢
**Converse of Redei: which odd integers arise as H(T)?**

Redei's theorem says H(T) is always odd. The converse asks: for which odd k does there exist a tournament T with H(T)=k?

**Permanent gaps discovered (THM-029, kind-pasteur-2026-03-06-S21, corrected S22):**
- **H=7 is impossible** for ANY tournament on ANY number of vertices. CORRECTED proof (S22): alpha_1=3 IS achievable at n>=7, but H=7 still impossible because H=7 requires (alpha_1=3, i_2=0), and i_2=0 forces common vertex among triples which forces c5>=1, giving alpha_1>=4. When alpha_1=3 occurs (n>=7), the triples don't all conflict, so i_2>=1, giving H>=11.
- **H=21 is a permanent gap — PROVED FOR ALL n** (kind-pasteur-S33). Complete proof via poisoning graph DAG argument (Dichotomy Theorem, Part R of THM-079).

**Achievable values (exhaustive):**
- n=5: {1,3,5,9,11,13,15}
- n=6: {1,3,5,9,11,13,15,17,19,23,25,27,29,31,33,37,41,43,45}
- n=7 (sampled): 77 distinct values from 1 to 189

**H=21 PROVED ABSENT through n=7 (opus-S38, THM-075):**
Exhaustive enumeration of all 2,097,152 tournaments on 7 vertices confirms H=21 never occurs. The gap 19→23 is consistent at n=6 and n=7. No (alpha_1, alpha_2) decomposition for H=21 is achievable. Strong evidence this is a permanent gap like H=7.

**Complete H-spectrum at n=7** (77 distinct values, all odd):
1, 3, 5, 9, 11, 13, 15, 17, 19, 23, 25, 27, 29, 31, 33, 35, 37, 39, 41, 43, 45, 47, 49, 51, 53, 55, 57, 59, 61, 65, 67, 69, 71, 73, 75, 77, 79, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99, 101, 103, 105, 109, 111, 113, 115, 117, 121, 123, 125, 127, 129, 131, 133, 135, 137, 139, 141, 143, 145, 147, 151, 153, 155, 157, 159, 171, 175, 189

**Gaps in [1,189] at n=7:** 7, 21, 63, 107, 119, 149, 161-169 (block), 173, 177-187 (block)
Note 63 = 7*9 and 21 = 7*3. These may be related to the H=7 gap.

**THM-079 Component Reduction (opus-S39):**
- Disconnected Omega: IMPOSSIBLE (THM-029 blocks I(component)=7)
- P_4 component: IMPOSSIBLE (two sharing 3-cycles on 5 verts force 3rd cycle)
- alpha_3>=1 decompositions: ALL IMPOSSIBLE (forces sum>=26>20)
- K_{1,3} star in (4,3) case: IMPOSSIBLE (forces alpha_3>=1)
- Remaining: connected Omega with I=21 via K_6-2e or larger dense graphs
- I(P_4,2)=21 discovered; graph classification: v=4 P_4, v=5 none, v=6 K_6-2e

**THM-079 Update (opus-S40):**
- **K_6-2e FULLY ELIMINATED** by Five-Cycle Forcing Theorem (3 lemmas, structural proof for all n)
- **i_2 jump pattern discovered**: achievable (alpha_1, i_2) pairs in tournaments systematically skip the values needed for H=21. Verified exhaustively at n<=7, by sampling at n=8.
- **H=7 and H=21 are the ONLY permanent gaps** in [1..200] at n<=8.
- **H=63 is NOT a permanent gap** (achieved at n=8, 138/500k samples).
- Remaining open: (8,1) K_8-e mixed-cycle case, (10,0) K_10 structural proof.

**H=63 is NOT a permanent gap (opus-S39 agent):**
H=63 found at n=8 (227 in 600k samples). All n=7 gaps except 7 and 21 are filled at n=8.
The ONLY permanent gaps through n=9 sampling are **H=7 and H=21**.

**THM-079 Update (opus-S41):**
- **EXHAUSTIVE n=8:** All 268M tournaments checked, H=21 found: 0.
- **Key Lemma (Part J):** Vertex in no 3-cycle => vertex in no cycle (layered structure).
- **Source/sink induction (Part K):** Score 0/n-1 vertex in no cycle; removing preserves Omega.
- **Cycle-rich min-H (Part L):** Among 18M cycle-rich n=8 tournaments, min H=25 > 21.
- **Parts M,N (PROVED at n=8):** (10,0) star capacity, (8,1) cascade forcing.

**THM-079 Update (opus-S42):**
- **n=9 matching:** Only 23.9% cycle-rich have 3 disjoint 3-cycles. 71.1% have mm=2.
- **alpha_1=10 at n=9:** Always t3=6,t5=4, i_2=9 or 10 (never 0). (10,0) impossible.
- **mm<=2 min H=45 at n=9:** Fewer disjoint 3-cycles forces more 5-cycles, larger H.
- **DICHOTOMY (0 counterexamples in 153k):** Every cycle-rich n=9 tournament has either 3 disjoint 3-cycles (Part C) or a good deletion to cycle-rich n=8 (induction).
- **H-spectrum n=9 (2M samples):** Only missing odd in [1..200] are 7 and 21.
- **Complete proof structure (Part P):** H!=21 for all n, modulo proving the dichotomy.

**Open questions:**
- ~~PROVE the deletion+matching dichotomy for all n >= 9~~ **RESOLVED** (kind-pasteur-S33, poisoning graph DAG argument)
- ~~Alternative: prove min-H for cycle-rich tournaments is > 21 for all n >= 8~~ **RESOLVED** (follows from dichotomy: cycle-rich H >= 25 for all n >= 8)
- Are there other permanent gaps beyond 7 and 21? EVIDENCE: NO (63 filled at n=8)
- What is the density of achievable values as max H grows?

**Connection:** Mitrovic-Stojadinovic (arXiv:2506.08841) "converse of Redei" is about poset-level parity (non-chain posets have even quasi-linear extension count), NOT about the H-spectrum. Does not address which odd integers are achievable as H(T).

**Source:** kind-pasteur-2026-03-06-S21, THM-029

---

## OPEN-Q-020 -- RESOLVED
**What determines the Worpitzky coefficients beyond t3?**

**RESOLVED by opus-2026-03-07-S46b/S46c:**

At n=6 (exhaustive, 24 F-classes): delta_1 = 8*t3 + 4*t5 + 8*alpha_2, delta_0 = H-1 = 2*t3 + 2*t5 + 4*alpha_2 (= OCF). The Worpitzky polynomial is a GRADED REFINEMENT of OCF.

The mechanism: Worpitzky coefficients encode moments E[fwd^r], and these follow a moment-cycle hierarchy (THM-092). Zero skewness (THM-091) eliminates odd cumulants. Each even cumulant kappa_{2k} adds one level of cycle complexity (cycles on <=2k+1 vertices).

At n=7 (156 F-classes sampled): delta_4 = 10*t3, delta_3 = 20*t3 confirmed. delta_2 needs invariants beyond t3.

**THM:** THM-087 (F,G updated), THM-090, THM-091, THM-092
**Source:** opus-2026-03-07-S46c

---

## OPEN-Q-021 🟢 Signed forward-edge polynomial SF(T,x) structure
**What is the combinatorial meaning of SF(T,x)?**

SF(T,x) = sum sgn(sigma) x^{fwd_T(sigma)} is palindromic and divisible by (x-1).
At n=4: SF = c*(x-1)^2*(x+1) for some integer c. What is c combinatorially?
At n>=6: SF is a COARSER invariant than F. What information does it lose?
Is there a matrix whose determinant equals SF(T,x)?

**THM:** THM-088
**Source:** opus-2026-03-07-S46b

---

## OPEN-Q-022 -- RESOLVED
**What determines the fourth cumulant kappa_4(T)?**

**RESOLVED by opus-2026-03-07-S46d (THM-093):**

kappa_4(T) = -(n+1)/120 + (2/C(n,4))*(t5 + 2*alpha_2) - 48/(n(n-1))^2 * t3^2

Key structural features:
- Constant = Bernoulli B_4 value: -(n+1)/120
- Linear t3 coefficient is EXACTLY ZERO (proved algebraically)
- t3^2 coefficient = -3*(4/(n(n-1)))^2 from Var^2 expansion
- t5 coefficient = 2/C(n,4), alpha_2 coefficient = 4/C(n,4)
- Verified exhaustively at n=5,6, sampled at n=7 (152 F-classes)

**kappa_6 introduces t7: YES.** Verified at n=7 (149 F-classes).
kappa_6 = (n+1)/252 + (2/C(n,6))*t7 - (4/49)*t3*(t5+2*alpha_2) + (80/3087)*t3^3

**Universal coefficient conjecture:** coeff(t_{2k+1}) in kappa_{2k} = 2/C(n, 2k).
Verified for k=1,2,3.

**Source:** opus-2026-03-07-S46d

---

## OPEN-Q-023 -- RESOLVED
**Prove: coeff(t_{2k+1}) in kappa_{2k} = 2/C(n, 2k) for all k.**

**RESOLVED by opus-2026-03-07-S46e (THM-117, was THM-095):**

PROVED algebraically. The proof has 5 steps:
1. Forward path formula: #fwd(2k+1)path = Σ_S H(T[S]) (OCF on subtournaments)
2. Each (2k+1)-cycle contributes 2·t_{2k+1} (OCF coefficient 2, unique subset)
3. Multinomial expansion: (2k)! · (n-2k) positions · 2/P(n,2k+1) = 2/C(n,2k)
4. Hierarchy separation: lower moments don't contain t_{2k+1}
5. Moment-to-cumulant preserves the coefficient

Verified algebraically for k=1,2,3,4 and n up to 12.

**Source:** opus-2026-03-07-S46e, THM-117

---

## OPEN-Q-024 🟢 Even Betti Vanishing for Tournament Path Homology
**Prove: beta_{2k}(T) = 0 for all k >= 1, for any tournament T.**

**UPDATE (kind-pasteur-S43): beta_2 = 0 PROVED (THM-108 + THM-109).**
**UPDATE (opus-2026-04-04-S1): Proof is FULLY ALGEBRAIC — THM-285 closes n=5 gap.**

Proof via strong induction using LES of (T, T\v):
- ~~Base case n=5 verified exhaustively (720/720)~~ **THM-285: n=5 case is VACUOUS** (no n=5 tournament has both b₁=0 and κ≥2; proof: κ≥2 forces dominator→source/sink contradiction)
- Induction step: good-vertex existence (THM-109)
- Case 2 (free cycles exist): Lemma A (free adj dom) + Lemma B → n-5 good vertices for n≥6
- Case 3 (all dominated): **Extreme Score Lemma** (ALGEBRAIC)
- **ALL cases are now algebraic. No exhaustive verification needed anywhere.**
- Comprehensive verification: 0 failures at n = 4-10

GLMY path homology Betti numbers beta_p of tournaments:
- beta_0 = 1 always (connected)
- beta_1 in {0, 1} (directed 1-hole from 3-cycle structure)
- beta_2 = 0 ALWAYS --- **PROVED** (THM-108 + THM-109)
- beta_3 in {0, 1, **2**}: appears at n=6 (1.2%), n=7 (8-11%), **n=8 (beta_3=2 at 0.08%)**
- **beta_4 NOT always 0**: Paley T_7 has beta_4 = 6 (THM-099). At n=8: beta_3*beta_4=1 can coexist (~0.15%)
- beta_1 and beta_3 are MUTUALLY EXCLUSIVE (proved n<=7, verified n=8)
- **Consecutive seesaw (beta_k*beta_{k+1}=0) REFUTED at n=8** (HYP-394, kind-pasteur-S48)
- **i_*-injectivity REFUTED at n=8** (HYP-380, kind-pasteur-S48): rank(i_*)=0 when b3=b3(T\v)=1
- Omega_p dimensions for Paley T_7: 7, 21, 42, 63, 63, 42, 21 (palindromic!)

**UPDATE (opus-S72b): β₁ ∈ {0,1} verified exhaustive n≤8, sampled n=9 (THM-223).**
Key discovery: β₁ is determined ENTIRELY by rank of transitive triple constraint matrix.
Cancellation chains are ALWAYS redundant. Combined with THM-095 seesaw: β₁·β₃=0 follows.
Algebraic proof of β₁ ≤ 1 still open (reformulated as transitive triple rank bound).

REMAINING OPEN:
- **What bound replaces beta_3 <= 1 at n >= 8?** (beta_3=2 confirmed at n=8,9)
- **Prove beta_1 ≤ 1 algebraically** — equiv. to rank(TT) ≥ C(n,2)-n (THM-223)
- Characterize which tournaments have beta_4 > 0 (appears linked to H-maximizers)
- Is beta_6 = 0 for all tournaments? (0/300 at n=7)
- Prove beta_2k = 0 for k >= some threshold, or find more counterexamples

**Source:** opus-2026-03-07-S46e, kind-pasteur-2026-03-08-S43

---

## OPEN-Q-025 -- RESOLVED (PROVED for all p)
**Prove Trace Alternation Theorem (THM-136) for all p**

**Statement (CORRECTED):** For primes p = 3 mod 4, sign(tr(A^k)_Paley - tr(A^k)_Interval) = (-1)^{(k-1)/2} for all odd k >= 5. (Original formula had (-1)^{(k-3)/2} which is off by a sign; see MISTAKE-019.)

**PROVED (kind-pasteur-S57):** Two-pronged algebraic proof:

1. **Dominant eigenvalue mechanism:** r_1 = |mu_1(interval)| = 1/(2*sin(pi/(2p))) dominates all other eigenvalues. The ratio r_1/r_2 ~ 2p/pi gives exponential dominance at power k. This ensures |S_I(k)| >> error bound >> |S_P(k)| for ALL odd k.

2. **Phase control:** sin(k*pi/(2p)) > 0 for all k in [1, p-1], determining sign(dominant term) = (-1)^{(k+1)/2}. Combined with magnitude dominance: sign(Delta_k) = -sign(S_I) = (-1)^{(k-1)/2}.

3. **Computational verification:** 1218+ individual (k,p) tests, zero failures. k=5 exact DP for 154 primes up to p=2000.

The proof is self-contained and works for ALL p >= 7 simultaneously. No finite verification needed.

**Source:** kind-pasteur-2026-03-12-S57 (proof), kind-pasteur-S56c (discovery)

---

## OPEN-Q-026 🟢 Does the interval maximize H for all circulant tournaments on Z_p, p >= 13?

**Statement (HYP-480):** The cyclic interval C_p = (Z_p, {1,...,(p-1)/2}) maximizes H among all circulant tournaments on Z_p for all primes p >= 13.

**Evidence:** Confirmed at p = 13 (exhaustive), p = 19 (THM-135), p = 23 (kind-pasteur-S57).

| p | H(Paley) | H(Interval) | Margin | Winner | Max circulant |
|---|----------|-------------|--------|--------|---------------|
| 7  | 189 | 175 | -7.4% | PALEY | Paley+complement ONLY (all others H=175) |
| 11 | 95,095 | 93,027 | -2.2% | PALEY | Paley+complement ONLY (2nd: H=93,467×10) |
| 13 | - | 3,711,175 | - | INTERVAL | (exhaustive, p≡1 mod 4, no Paley) |
| 17 | - | 13,689,269,499 | - | INTERVAL | (exhaustive over SC circulants) |
| 19 | 1,172,695,746,915 | 1,184,212,824,763 | +0.98% | INTERVAL | - |
| 23 | 15,760,206,976,379,349 | 16,011,537,490,557,279 | +1.59% | INTERVAL | - |

EXHAUSTIVE SCANS (kind-pasteur-2026-04-16):
  n=7: ALL 8 circulant tournaments. Top H=189 (2 tournaments: Paley+complement).
       6 tournaments have H=175 (including Cyclic). Paley is 8% better than rest.
  n=11: ALL 32 circulant tournaments. Top H=95,095 (2: Paley+complement).
        10 tournaments share H=93,467. Cyclic has H=93,027 (rank ~18/32).
  n=7,11 alpha breakdown:
    n=7:  Paley α₁=80, α₂=7.  Cyclic α₁=59, α₂=14.  (Cyclic has 2× α₂!)
    n=11: Paley α₁=21,169, α₂=10,879, α₃=1,155. Cyclic α₁=18,397, α₂=11,110, α₃=1,474.
          Cyclic has MORE α₂ and α₃, but Paley's α₁ advantage (5,544 in H) > Cyclic's advantage (3,476).

Crossover: Paley wins at p=7 and p=11 due to α₁ dominance. Interval wins at p≥13.
The α₁ percentage gap narrows: 35.6% (n=7), 15.1% (n=11), 3.6% (n=19). Paley's α₁ lead evaporates.
At n=7: kmax=2 (no 3-packings possible), so H has only α₁,α₂ terms — Paley α₁ wins.
At n=11: Paley α₁ advantage still > Cyclic α₂+α₃ advantage.
At n=19: Cyclic α₃+ advantage (26.7B) > Paley α₁+α₂ advantage (15.2B) → Cyclic wins.
α₂ comparison crossover: Cyclic > Paley at n=7,11 (small n, disjoint packing easier);
                          Paley > Cyclic at n=19 (large n, more α₁ → more α₂ pairs).

WHY interval wins for large p (kind-pasteur-2026-04-16):
  Paley has MORE α₁ and α₂ at n=19, but Interval wins by +11.5B total.
  Interval has +26.7B from α₃+ terms: its cycles pack into disjoint triples better.
  Paley's pseudorandom structure creates many individual cycles but they scatter;
  Interval's consecutive structure creates harmonically aligned cycles for packing.

EXACT α-DECOMPOSITION COMPARISON at n=19 (kind-pasteur-2026-04-16, VERIFIED):
  k | Paley α_k          | Cyclic α_k         | Cyclic advantage | H contribution
  1 | 130,965,270,477    | 126,443,605,257    |   -4,521,665,220 | 2×diff = -9.04B  (Paley wins)
  2 | 123,659,531,220    | 122,111,579,294    |   -1,547,951,926 | 4×diff = -6.19B  (Paley wins)
  3 |  41,184,418,943    |  42,960,731,622    |   +1,776,312,679 | 8×diff = +14.21B (Cyclic wins)
  4 |   4,903,920,444    |   5,521,030,944    |     +617,110,500 |16×diff = +9.87B  (Cyclic wins)
  5 |     251,464,164    |     331,078,344    |      +79,614,180 |32×diff = +2.55B  (Cyclic wins)
  6 |       2,221,081    |       4,100,656    |       +1,879,575 |64×diff = +0.12B  (Cyclic wins)
  Net: Paley advantage 15.2B (via α₁,α₂) vs Cyclic advantage 26.7B (via α₃+) = +11.5B net for Cyclic.

  α₃/α₂ ratios: Paley=0.333, Cyclic=0.352 → Cyclic is intrinsically better at 3-packing!
  The k=5,6 percent advantage for Cyclic: +31.7% at k=5, +84.7% at k=6 — grows with k.
  Source: paley_t19_alpha.out (H(Paley)=1,172,695,746,915 verified ✓)

The interval's margin is WIDENING with p, consistent with the spectral argument: |mu_1| ~ p/pi grows faster than Paley's sqrt(p)/2.

**What remains:** Extend to p = 29, 31. An analytical proof could use the spectral concentration argument from THM-137. Whether interval maximizes H among ALL tournaments (not just circulant) is a separate open question.

**Source:** opus-2026-03-12-S58, kind-pasteur-2026-03-12-S56c, kind-pasteur-2026-03-12-S57

---

## OPEN-Q-027 -- RESOLVED (PROVED with correction)
**Prove the Grand Energy Formula (THM-201)**

**RESOLVED by kind-pasteur-2026-03-15-S112 (THM-217):**

The original formula E_{2k}/E_0 = 2*(n-2k)^k / P(n,2k) is the LEADING-TERM APPROXIMATION only, exact for k ≤ 2 but requiring corrections for k ≥ 3.

The EXACT formula uses combinatorial g_k polynomials of degree k:
- CV²(H) = Σ_{k≥1} 2·g_k(n-2k) / (n)_{2k}
- g_k defined via transfer matrix M(x) = [[1,2x,0],[0,0,1],[1,x,0]]
- Weight formula E[∏Z_j] = 2^c/(n)_L PROVED (c=components, L=|S|)
- g_k(m) ~ 2^{k-1}·m^k/k! + lower terms (leading term is original formula)
- Verified exhaustively n=3..18 via bitmask DP

**Source:** THM-217, kind-pasteur-S112, opus-S89c

---

## OPEN-Q-028 🟢 Are there forbidden H values beyond 7 and 21?

**Statement:** Are 7 and 21 the ONLY permanently forbidden H values? H=63 was shown achievable at n=8 (HYP-1106 refuted). But could there be large forbidden values proportional to n!/2^{n-1}?

**UPDATE (kind-pasteur-2026-03-20-S1):** Confirmed via 500K random n=9 tournaments: only gaps in odd [1,500] are H=7 and H=21. H=63 achieved (9 occurrences). Also confirmed at n=8 (200K samples): only gaps below 100 are 7 and 21. This is very strong evidence that 7 and 21 are the only permanent gaps.

**Evidence:** Only 7 and 21 are known forbidden. 63 is achievable (n=8). No other candidates found through n=11.

**UPDATE (opus-S227):** H-spectrum density at n=8 is 320/331 = 96.7%. Only 11 gaps remain, dominated by {7, 21}. In the metagraph, forbidden values create "missing floors" that force edge jumps. At n=5, 33% of edges bridge the H=7 gap. The fraction decreases as n grows (2.2% at n=7). Edges bridging H=7 gap: 0, 0, 7, 21, 47 for n=3..7.

**STRONG CONJECTURE:** Only H=7 and H=21 are permanently forbidden. All other gaps are transient (filled at large enough n).

**Source:** kind-pasteur-S107, opus-S227

---

## OPEN-Q-029 -- RESOLVED
**Why does log_tau(131) = 8.0003?**

**RESOLVED by opus-2026-03-15-S90 (multiple proofs):**

131 = Tr(M^8) EXACTLY, where M is the 3×3 transfer matrix from S112. τ₃^8 ≈ 130.977 and the Pisot correction 2|λ_c|^8 cos(8θ) ≈ +0.023 pushes the sum to exactly 131.

**Why the correction is so small:** arg(λ_c)/π ≈ ln(2), so 8·arg/π ≈ 8·ln(2) ≈ 5.545 ≈ 5.5, making cos(8·arg) ≈ cos(5.5π) ≈ 0.13 (small). The n=8 case is special because 8·ln(2) is close to the half-integer 11/2.

**Additional discoveries (S90 session):**
- 504 = T(13) in the standard tribonacci sequence (confirmed)
- The transfer matrix char poly at x=1 IS the tribonacci equation
- The τ-φ clock gear ratio arg(λ_c)/π ≈ ln(2) explains ALL Pisot near-integers
- Tr(M^n) mod 8 has EXACT period 8 (Bott periodicity connection)

**Source:** opus-2026-03-15-S90c (monad), S90h (τ-φ clock), S90l (the number 8)

---

## OPEN-Q-030 -- RESOLVED (PROVED for all n ≥ 4)
**Prove Simplicial Rédei for ALL n (THM-220)**

**RESOLVED by opus-2026-03-15-S90 (THM-220):**

The Key Lemma IS proved algebraically for all n: Given a→b not in any transitive triple, the four possible orientations of {a,b,c} are: (1) a→c,b→c: transitive a>b>c, contains a→b — CONTRADICTION. (2) a→c,c→b: transitive a>c>b, contains a→b — CONTRADICTION. (3) c→a,b→c: 3-cycle a→b→c→a — the ONLY non-core possibility. (4) c→a,c→b: transitive c>a>b, contains a→b — CONTRADICTION. Since 3 of 4 orientations force a→b into a transitive triple, the only possibility for ALL c is case (3): every c forms a 3-cycle with {a,b}. This gives score(a)=1, score(b)=n-2.

The Main Argument (at most one non-core edge) then follows by contradiction in 4 cases of vertex overlap. Case 3 (b=c) requires n≥4 so that V\{a,b,d}≠∅.

All verified exhaustively n=4..8 (268M at n=8), sampled n=9 (500k, 0 violations).

**Source:** opus-2026-03-15-S90 (THM-220), opus-2026-03-16-S90q (proof verification)

---

## OPEN-Q-031 🟢 Is arg(λ_c)/π = ln(2) exact or approximate?

**Statement:** The tribonacci complex eigenvalue angle satisfies arg(λ_c)/π ≈ ln(2) to 4 significant figures (difference 4.3×10⁻⁴). Is this exact?

**Evidence:** NOT exact (verified: the predicted root from arg=π·ln(2) does not satisfy the char poly). But the proximity is remarkable and explained by the information-theoretic interpretation: the tribonacci clock ticks at approximately 1 bit per half-revolution.

**Source:** opus-2026-03-15-S90h (τ-φ clock)

---

## OPEN-Q-032 -- PARTIALLY RESOLVED (FAILS at n=6)
**Tournament equidecomposability: is (H, β₁) a complete invariant?**

**ANSWER: NO.** (H, β₁) is complete at n=5 (8 classes, each with unique I(Ω₃, x)) but FAILS at n=6.

**Counterexamples at n=6 (5 found):**
- (H=9, β₁=0): TWO distinct I(Ω₃): (1,2,1) and (1,3,0,0)
- (H=15, β₁=0): TWO distinct: (1,4,0,0,0) and (1,5,0,0,0,0)
- (H=29, β₁=0): TWO distinct: (1,6,1) and (1,6,2)
- (H=33, β₁=0): TWO distinct: (1,6,2) and (1,7,1)
- (H=37, β₁=0): TWO distinct: (1,7,1) and (1,7,2)

ALL counterexamples have β₁=0 — the β₁=1 classes remain unique!
This means: β₁ provides a COARSER invariant. The FULL independence polynomial I(Ω₃, x) requires more information (α₂ distinguishes within β₁=0).

**REVISED CONJECTURE:** (H, β₁, α₂) may be complete. Or (H, full α-vector) is complete by definition.

**Source:** opus-2026-03-15-S90k (n=5), opus-2026-03-16 (n=6 counterexample)

---

## OPEN-Q-033 -- RESOLVED (PROVED analytically)
**The n-tribonacci family: T_n - M_{n-2} = 1/(kM_k+2) + O(1/k⁴)**

**PROVED by opus-2026-03-16 (perturbation analysis):**

Write T_n = M_k + ε where k = n-2. Substituting into λ³ = kλ² + λ + 1 and using M² = kM+1:

  (kM+2)·ε = 1  at leading order.
  So ε = 1/(kM_k+2).

Since M_k ~ k for large k: ε ~ 1/(k²+2) ~ 1/k².

Verified numerically: the ratio δ_actual / (1/(kM+2)) → 1 as n → ∞ (0.999599 at n=19).

At n=3: δ = 0.221 (maximum), predicted 0.276 (ratio 0.80 — leading order less accurate for small k).

**Source:** opus-2026-03-16-S90 (perturbation proof)

---

## OPEN-Q-034 🟢
**Meta-structure: why does cancellation dominate this theory?**

Every major result in the project is fundamentally about cancellation: im(d₂) cancels in the seesaw, Walsh coefficients cancel for odd-length paths, S(T)=0 at even n, β₂=0 always, OCF = exact cancellation between H and I. Is there a *single structural principle* (perhaps from homological algebra or categorical representation theory) that implies all of these cancellations simultaneously? The F₂ uniqueness argument (S71r: "WHY TWO GENERATES SEVEN") is a partial answer — but it explains *why F₂* rather than *why cancellation*. See `07-reflections/seesaw-and-cancellation.md`, `07-reflections/what-the-proof-will-look-like.md`.

**Source:** opus-2026-03-16-S73

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery.
- **OPEN-Q-002**: Claim A PROVED for all n. OCF proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20). See CONJ-001, THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **OPEN-Q-009**: Arc-flip invariance resolved — E(T) = 0 for all T (OCF proved). See THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy.
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11.

---

## OPEN-Q-035 -- RESOLVED (degree = 2*floor((n-1)/2), NOT fixed at 4)
**Does the heat kernel polynomial P_x(z) have degree exactly floor(n/3)*2 for general n?**

**RESOLVED by kind-pasteur-2026-03-20-S2 (THM-259):**

The Walsh degree is NOT fixed at 4. It is **2*floor((n-1)/2)**:
- n=5,6: degree 4
- n=7,8: degree 6 (INCREASES! 2520 new degree-6 coefficients at n=7)
- n=9,10: degree 8
- General: n-1 for odd n, n-2 for even n

Follows from THM-076: the maximum Walsh weight is 2*max_k where k <= (n-1)/2.
Verified exhaustively at n=5 (91 nonzero coefficients) and n=7 (4516 nonzero).

The original conjecture floor(n/3)*2 was correct for n=5,6 but WRONG for n=7.
THM-076 gives the correct formula via path-covering analysis.

Only 5 distinct |Walsh amplitudes| at n=7, all matching THM-076 exactly.
Super orthogonality redundancy: 4516 / 2 = 2258x.

**Source:** kind-pasteur-2026-03-20-S2, THM-259

---

## OPEN-Q-036 🟢
**Does the backward trick P_x(2) = mean H hold for other starting points?**

At n=6, P_transitive(2) = 29 = mean H. Only 3/1024 tilings have this property. What characterizes these special starting points? Are they related to self-complementary tournaments or specific score sequences?

**Source:** kind-pasteur-2026-03-17-S116n33

---

## OPEN-Q-037 🟢
**Does the splitting of mean H in Z[i] vs Z[sqrt(-7)] generalize to other n?**

At n=6: mean H = 29 splits as 5²+2² (golden) and 1²+7·2² (forbidden). At n=7: mean H = ? At n=8: mean H = ? Do the two world-defining primes always appear in the splitting?

**Source:** kind-pasteur-2026-03-17-S116n33

---

## OPEN-Q-038 🟡
**Characterize the graph class where I(G,x) has all real roots beyond claw-free.**

Tournament conflict graphs Omega(T) have all real roots of I(G,x) for n<=8 (proved via claw-free) and n<=20 (verified). At n>=9, claw-free FAILS but real roots persist. What graph property replaces claw-free?

**Source:** kind-pasteur memory, originally from S14-S18

---

## OPEN-Q-039 🔴 — SUBSTANTIALLY RESOLVED (sessions S211-S249)
**Understand the isomorphism class graph G_n completely**

**MASSIVE PROGRESS (opus S211-S249, kind-pasteur S20bo-S20dj):**

G_n = Q_{C(n,2)} / S_n is a genuinely new mathematical object (no prior literature). The merged metagraph G_n/Z_2 has been computed exactly through n=9 with 7 exact edge terms: E(G_n) = 1, 5, 30, 290, 4086, 91161, 3,380,751.

**RESOLVED sub-problems:**
1. ✅ Extended to n=6,7,8,9 (6880 classes at n=8, 191536 at n=9)
2. ✅ Diameter = n-2 DISPROVED at n=7 (diam=7≠5). Actual: 1,2,3,4,7,8
3. ✅ H-DAG property REFUTED: G_n is NOT a strict DAG. Level edges appear at n≥5 (1, 15, 136 for n=5,6,7). H-decreasing edges appear at n=7 (962/4086). The H-gradient is strong but imperfect. See MISTAKE-035.
4. ✅ Spectral data computed at n=3..7 (Ramanujan fails at n=6)
5. ✅ |Aut|-degree connection: corr→0 at large n (classes become generic)
6. ✅ I(G_n,2) computed: 5, 13, 793, 15B (super-exponential "meta-H")
7. ✅ Staircase connection: Mode A/B recursion, y=x diagonal, within-type fraction→3/4

**EDGE FORMULA (the keystone):**
  edge_orbits = T_n/2 + (n-2)! [verified n=3..6, Burnside-computable]
  E(G_n) = edge_orbits - gap_orbits [exact]
  E(G_n) ≈ (T_n - twin_SL)/2 - D(n-2) [99.6% accurate, all Burnside]
  E(G_n) ≈ T_n/2 for n ≥ 12 [asymptotic]

**SL_mine FORMULA (kind-pasteur-S20eh, PROVED):**
  D(n) = (1/n!) sum_{ct with 1 even cycle 2k} count(ct) * k * 2^{a(ct)}
  SL_mine <= D(n) with small correction from |Aut|>1 classes
  D - SL_mine = 0, 0, 0, 2, 4 at n=3..7
  D(3..12) = 2, 6, 16, 60, 328, 3160, 54928, 1722992, 97323552, 9941203552
  CORRECTED: T - 2E != SL_mine (multi-edge surplus exists at n>=5)
  Multi-edge surplus = 0, 0, 12, 66, 416 at n=3..7

**STRUCTURAL LAWS (19+ verified):** DAG, BBK impossibility, rib crossover, spine ~4-regular, ribs linear in n, sea dominates, ΔH=2^(n-2), cell uniformity, lattice oscillation, etc.

**REMAINING:** Exact formula for gap_orbits (= 2,5,20,86,490,3703,47889); twin_SL residual; chi=n-1 conjecture proof (greedy fails at n≥6).

**Source:** opus S211-S249, kind-pasteur S20bo-S20dj. Library: `04-computation/tournament_metagraph.py`

---

## OPEN-Q-040: THE KRAWTCHOUK FRAMEWORK (sessions S291-S312, 2026-03-24)

🟡 **The Krawtchouk coordinate system for tournament space**

**RESOLVED sub-problems:**
1. ✅ **Tournament Counting Theorem** (S291): V_n×n!/2^m = 1 + Σ(1/k)×n↓k×2^{(k²-1)/2-(k-1)n}. Euler product with poles at 4,16,64,... controlled by 1/3. D₃(0) = 128/81.
2. ✅ **Band-limitedness** (S305,S310,S311, **CORRECTED kind-pasteur-S1 2026-03-25**): Walsh degree = 2*floor((n-1)/2) for all n>=4 (THM-260). Band-limited at m/2 for n>=6. **CORRECTION:** n=5 is NOT band-limited at m/2 (degree 4 > m/2=3). Odd-weight Walsh coefficients are nonzero in tiling model (complement symmetry fails).
3. ✅ **Krawtchouk 3-axis system** (S307): B₁≈-H (r=-0.94), B₂≈-c₃ (r=-0.86), B₃=twist. SC classes have B_odd=0 exactly (Krawtchouk parity).
4. ✅ **Diameter = A003141** (S306): max feedback arc set. Growth ~n²/4, not n-2 (small-n coincidence).
5. ✅ **Paley→Dual Codes** (S306,S308): P₇+I→Hamming[7,4,3], P₂₃+I→Golay[23,12,7].
6. ✅ **Not an association scheme** (S306): full algebra dim=35 vs needed 7 at n=5.
7. ✅ **Spectral gap = 2/n explained** (S312): comes from K₁ spacing 2/m compressed by S_n quotient (factor m/n).
8. ✅ **Waggly = all connections** (S296-S301): wiggly⊂waggly, blue/black⊂waggly. Completeness at k*=diam.
9. ✅ **Waggly alphabet** (S302-S304): range-3 harmonic most neutral. Vertex-count law. All-same-range combos special.
10. ✅ **Practical tools** (S308-S309): pre-filter eliminates 98% of canon calls. tournament_tools.py library. tournament_codec.py (kind-pasteur).

**OPEN sub-problems (the 10 boundary questions from S307):**
1. ✅ B1: **RESOLVED** (THM-260, kind-pasteur-S1): Walsh degree = 2*floor((n-1)/2) for all n. Band-limited at m/2 for n>=6. Proof: THM-076 upper bound + interleaving construction lower bound.
2. 🟢 B2: Exact constant in A003141 n^{3/2} correction
3. 🟢 B3: Is transitive always a diameter endpoint?
4. 🟢 B4: Does K₁-H correlation → 0 or stabilize? (0.94→0.89→0.83)
5. 🟡 B5: Exact neutrality formula SL(d,n) as function of distance
6. 🟢 B6: Width W(H) asymptotic distribution
7. 🟢 B7: Is there ANY partition giving an association scheme?
8. 🟢 B8: Is range ⌊(n-1)/2⌋ always most neutral?
9. 🔴 B9: β₂=0 for all tournaments (proof strategy via band-limitedness, S312)
10. 🟡 B10: min-FAS(T) in terms of OCF invariants

**Key new files:** euler_product_tournament_s291.py, waggly_layers_s297.py, waggly_completeness_s301.py, waggly_alphabet_s302.py, almost_1d_s305.py, krawtchouk_h_n7_s306.py, paley_codes_s306.py, tournament_tools.py, tournament_codec.py

**Key reflections:** the-tiling-hypercube.md, the-boundary-between-1d-and-2d.md, euler-product-and-metagraph.md, paley-gives-dual-codes.md, h-is-band-limited.md, what-we-can-and-cannot-know.md, tournament-compression-and-beyond.md, terminology-evolution.md, diameter-is-feedback-arc-set.md


---

## OPEN-Q-044 🟢 Alpha Mechanism Shift: When Does Each α_k Dominate?

**Discovery (kind-pasteur-2026-04-16):** The dominant term in H = I(Ω,2) = Σ 2^k · α_k shifts with n.
H-maximizing cyclic interval tournament C_n:

| n | dom term | 2^1·α₁ | 2^2·α₂ | 2^3·α₃ | notes |
|---|----------|---------|---------|---------|-------|
| 3-9  | α₁ | largest | 2nd | small | α₁/(2α₂) > 1 |
| 11-17 | α₂ | 2nd | largest | 3rd | FIRST CROSSOVER n≈10 |
| 19+ | α₂ | 3rd | largest | 2nd | SECOND CROSSOVER: α₃ overtakes α₁ at n≈17-19 |

**Complete verified table for C_n (cyclic interval tournament):**

| n  | α₁ | α₂ | α₃ | α₁/(2α₂) | α₃/α₂ | H | H(n)/H(n-2) |
|----|----|----|----|-----------|---------|----|-------------|
| 17 | 1,651,334,601 | 1,482,234,998 | 458,011,858 | 0.5570 | 0.3090 | 13,689,269,499 | — |
| 19 | 126,443,605,257 | 122,111,579,294 | 42,960,731,622 | 0.5177 | 0.3518 | 1,184,212,824,763 | 86.5 |
| 21 | 12,030,499,746,751 | 12,330,182,836,208 | 4,796,354,751,404 | 0.4878 | 0.3890 | 125,547,534,942,879 | 106.0 |
| 23 | 1,391,602,826,199,187 | 1,499,656,616,321,278 | 632,921,002,322,216 | 0.4640 | 0.4220 | 16,011,537,490,557,279 | 127.6 |

**Full α-decomposition n=21:**
  α₁=12,030,499,746,751   α₂=12,330,182,836,208   α₃=4,796,354,751,404
  α₄=738,531,326,288      α₅=58,868,297,768        α₆=1,454,221,328       α₇=12,571,712
  H = 125,547,534,942,879

**Full α-decomposition n=23 (NEW, kind-pasteur-2026-04-16-S1):**
  α₁=1,391,602,826,199,187   α₂=1,499,656,616,321,278   α₃=632,921,002,322,216
  α₄=111,796,734,828,336     α₅=10,945,293,151,712       α₆=412,282,843,184       α₇=7,454,017,376
  H = 16,011,537,490,557,279 ✓

**Term ordering at n=23:** 4α₂ > 8α₃ > 2α₁ > 16α₄ > 32α₅ > 64α₆ > 128α₇
  (5.999P > 5.063P > 2.783P > 1.789P > 0.350P > 26.4T > 0.95T)

**Special structure at n=21:** α₇ = 12,571,712 = perfect 7-triangle-packings.
  Only packing type is (3,3,3,3,3,3,3) since 7×3=21. Perfect vertex coverage.
**Structure at n=23:** α₇ = 7,454,017,376 counts 7-packings with cycle-length sum ∈ {21, 23}.
  Sum must be odd (7 odd numbers), and ≤23. So: sum=21 (all 3-cycles, 2 vertices free) OR
  sum=23 (one 5-cycle + six 3-cycles, all 23 vertices covered). Sum=22 impossible (even).

**H growth ratio H(n+2)/H(n):** 86.5, 106.0, 127.6 → increments +19.5, +21.6 → growing.
  Predicted H(25) ≈ H(23) × 150 ≈ 2.4 × 10^18.

**Key ratio α₃/α₂ progression:**
  n=17: 0.3090, n=19: 0.3518 (+0.043), n=21: 0.3890 (+0.037), n=23: 0.4220 (+0.033)
  First differences decreasing by ~0.004/step. Projected:
  n=25: ≈0.451 (+0.029), n=27: ≈0.476 (+0.025), n=29: ≈0.497 (+0.021), n=31: ≈0.514 → THIRD CROSSOVER
  **Revised estimate: third crossover (8α₃ > 4α₂) at n≈31**, not n≈27-29 as previously estimated.

**Timing:** cycle_cc 383s, SSC runs 732s+612s. Total 1728s for n=23 with numpy.
  Bottleneck is cycle_cc (Python BFS). C implementation would reduce to ~3s.

**Open:** Third crossover: α₃ dominates at n≈31 (needs n=25,27 data to confirm).
         C implementation of cycle_cc needed for n≥25.

---

## OPEN-Q-046 🟡 The SC Blowup: $\Omega(T_{\mathrm{SC}})$ and H Formula

The **SC blowup** $T_{\mathrm{SC}}$ (arc $u_r \to v_s$ iff $u \to v$ in $T$ and $r=s$, OR $v \to u$
and $r \neq s$) satisfies the **Universal Score Theorem**: every $v_0$ has out-degree $n$, every
$v_1$ has out-degree $n-1$, regardless of $T$. See `07-reflections/sc-blowup-and-twin-gaining.md`.

The Kronecker formula $A(T_{\mathrm{SC}}) = A(T) \otimes I_2 + A(T)^\top \otimes \Phi + I_n \otimes e_{01}$
shows $T_{\mathrm{SC}}$ encodes BOTH $T$ and $T^{\mathrm{op}}$ simultaneously.

**Open (🟡):** What is $\Omega(T_{\mathrm{SC}})$ in terms of $\Omega(T)$? Is there a formula
$$H(T_{\mathrm{SC}}) = I(\Omega(T_{\mathrm{SC}}), 2) = f(I(\Omega(T), x))$$
for some operation on the independence polynomial?

**Candidate:** $H(T_{\mathrm{SC}}) \approx I(\Omega(T), 2)^2 / C(n)$ or involves $I(\Omega(T), 2) \cdot I(\Omega(T^{\mathrm{op}}), 2)$ with correction. Currently ruled out as single-variable formula (H_SC is NOT a function of H(T) alone).

**Key data:** At $n=5$, $H_{\mathrm{SC}}$ varies only 4.2% across all 12 iso classes ($14937$–$15565$).
At $n=3$: $H_{\mathrm{SC}} \in \{41, 45\}$ for the two classes. $H_{\mathrm{SC}}(\mathrm{Trans}) = 41$,
$H_{\mathrm{SC}}(C_3) = 45$.

**Source:** oracle-2026-05-15-S2, `05-knowledge/results/blowup_study.out`

---

## OPEN-Q-045 🟢 H Under Tournament Blowup (Column Row Step)

The tournament **blowup** $T[K_2]$ replaces each vertex $v$ with a directed pair
$v_0 \to v_1$, expanding each arc $u \to v$ to all four arcs $u_i \to v_j$.
This doubles $n$, corresponding to the row step $(r, k) \to (r+1, k)$ in the
2-adic column family grid (see `07-reflections/adic-column-families.md`).

**Q1:** Is there a formula $H(T[K_2]) = f(H(T), n)$?

**Q2:** Does blowup preserve SC status within a column family? SF status?

**Q3:** The **pairs anomaly** ($\lfloor n/2 \rfloor$ gains +1 at the $r=0 \to r=1$
seam) suggests H may have analogous anomalous first-blowup behavior:
$H(\text{blowup of odd-}n\, T)$ vs $H(\text{blowup of even-}n\, T)$ — are
these qualitatively different?

**Q4:** Does SC∩SF = SC($n-2$) for the family:
$\#(\text{SC} \cap \text{SF})(2^r(2k-1)) = \#\text{SC}(2^r(2k-3))$ for $r \geq 1$?
(This is the even-row analog of the proved odd-$n$ identity.)

**Related:** Linial-Morgenstern conjecture (INV-013: random blowup of transitive
tournament). The blowup operation is exactly the row step in the 2-adic grid.

**Expected difficulty:** SMALL CASES immediately computable. General formula: MEDIUM.
**Source:** oracle-2026-05-15 (2-adic column family analysis)

**Source:** kind-pasteur-2026-04-16-S1, `alpha_full_ssc_fast_n23.out`, `alpha_full_ssc_fast_n21.out`

**MAJOR UPDATE (kind-pasteur-2026-06-09-S1, THM-454/450):**
- **Q1 ANSWERED (negative + repaired):** H(T[K₂]) is NOT a function of (H(T), n) — not even of
  I(Ω(T),x) (n=5 counterexample: equal typed IP, H(K₂) 3225 vs 2785; the missing data is EVEN
  cycles, which twin-insertion converts to odd). What IS true: **strong-component product law**
  H(T[K₂]) = ∏_C H(C[K₂]) (PROVED); twin-lift laws c3'=8c3, c5'=32c5+32c4+6c3 (+c7' law);
  cycle-spectrum (c3..c6) determines H(T[K₂]) at n≤6 (n=7 separation test open, HYP-2353);
  congruence H(T[K₂]) ≡ 2H(T)−1 (mod 8).
- **Q2 partial:** T[K₂] is op-equivariant (PROVED via orbit symmetries, THM-450(6)).
- T[K₂] is one of exactly THREE 2×2-block doublings (THM-450 trichotomy); the skew-Sylvester
  member D (THM-447) is the spectral/Hadamard-clean one; SCblow is the H-maximizing one.

---

## OPEN-Q-047 🟡 Characterize Real-Rootedness of I(Ω(T), x)

**Correction (opus-2026-05-29-S8):** The universal TRRT statement is already refuted by THM-025 at n=9.
The surviving problem is to characterize the tournaments for which I(Ω(T), x) has all real, negative roots.

**What's proved:** For n ≤ 8, Ω(T) is claw-free (a claw requires ≥ 9 vertices), so real-rootedness follows from Chudnovsky-Seymour (2007).

**Counterexample:** THM-025 gives an n=9 tournament with score sequence [1,1,3,4,4,4,6,6,7] and
I(Ω,x)=1+94x+10x²+x³. Newton k=2 fails (100 < 141), so two roots are non-real.

**Why notable:**
- Generic/sample tournaments often remain real-rooted despite the n=9 failure.
- For the real-rooted subclass, ultra-log-concavity and product formula H(T)=∏_i(1+2r_i) remain powerful.
- The THM-025 counterexample may isolate the exact obstruction shape.

**Sub-conjecture status:** Ω(T) is NOT always perfect (see INV-032 / THM-019 updates), so perfectness is also a subclass question.

**Key open questions:**
1. What structural property of Ω(T) (beyond claw-free) forces real-rootedness?
2. Which Hermite-Biehler/interlacing lemmas survive after accounting for THM-025?
3. Can the n=9 failure family be characterized exactly?

**Priority:** 🟡 IMPORTANT. A structural characterization would be publishable as a standalone result.
**Source:** opus-2026-05-16-S1, reflection `real-rootedness-omega-conjecture.md`

**Computational updates (oracle-2026-05-17-S1):**
- Root gap (-1/3,-1/4) confirmed empty at n=6 (exhaustive), n=7 (2000), n=8 (300), n=9 (50).
- ULC (Newton-Maclaurin inequality) confirmed at n=6..9, zero violations.
- Forbidden (α₁=3, α₂=0) confirmed absent at n=6..9 in all samples.
- Vieta at n=5 (r=-2/(H-1)) exact to machine precision.
- SC tournaments have most asymmetric root ratio at n=6: min 0.00251 (SC) < 0.00279 (NS).
- (H, I(Ω,6)) separates only 7/47 n=6 classes by (H,I6) alone (much coarser than hoped).
- Degree-3 polynomials first appear at n=9 (44/50 samples). ULC still holds.
See `07-reflections/root-spectrum-n6-computations.md`.

---

## OPEN-Q-048 🟢 Ultra-Log-Concavity for Tournament Independence Polynomials

**The theorem (proved):** If $I(\Omega(T),x)$ is real-rooted (proved universally only for $n \leq 8$; false universally from $n=9$ by THM-025), then $(\alpha_k/\binom{d}{k})_{k=0}^d$ is log-concave (ultra-log-concave), where $d = \alpha(\Omega(T))$.

**Proof:** Newton's inequalities for real-rooted polynomials with positive roots. Elementary symmetric polynomials $e_k(\rho_1,\ldots,\rho_d)$ satisfy Newton-Maclaurin: $(e_k/\binom{d}{k})^2 \geq (e_{k-1}/\binom{d}{k-1})(e_{k+1}/\binom{d}{k+1})$. Since $\alpha_k = \alpha_d \cdot e_{d-k}(\rho)$, ULC follows.

**Erdős context:** This is the tournament analog of the Heron-Rota-Welsh theorem (ULC for matroid Whitney numbers, proved by Adiprasito-Huh-Katz). Both prove ULC via underlying geometry (real-rootedness for tournaments, Hodge theory for matroids).

**Status:** PROVED conditional on real-rootedness. Computationally verified n=6..9.
**Priority:** 🟢 Interesting. Connects to the Huh-Katz matroid theory.
**Source:** oracle-2026-05-17-S1, computation `root_spectrum_fast.py`.

**NEW (oracle-2026-05-19-S1): UNCONDITIONAL proof of ULC at k=1 via Turán's theorem.**
For any tournament T: since bar_Omega(T) is K_{d+1}-free (max clique size = d = degree of I),
Turán's theorem gives alpha_2 <= (1-1/d)*alpha_1^2/2, which is exactly ULC at k=1:
   alpha_1^2 >= 2d/(d-1)*alpha_2.
No TRRT required. Equality iff I(Omega,x) = c*(x+rho)^d (all roots equal, Turán extremal).

**Also proved (conditionally on K4-free structure):** ULC at k=2 for complete tripartite
co-conflict graphs: (ab+bc+ca)^2 >= 3(a+b+c)*abc.
Proof: LHS-RHS = (1/2)[(ab-ac)^2+(ab-bc)^2+(ac-bc)^2] >= 0.
Verified: 0 violations in all n=9 samples (91/100 degree-3).
See `07-reflections/ulc-turan-unconditional-proof.md`.

---

## OPEN-Q-050 🟡 Unconditional ULC at k=2 via Kruskal-Katona

**Goal:** Prove $\alpha_2^2 \geq 3\alpha_1\alpha_3$ (ULC k=2, d=3) without assuming TRRT.

**Current status:**
- Proved for complete tripartite co-conflict graphs $K_{a,b,c}$ via the algebraic identity.
- Zero violations in n=9 random samples (91/100 degree-3 tournaments, 0 failures).
- Universal TRRT would have implied this via Newton's inequalities, but universal TRRT is refuted by THM-025.
- The "bad" counter-example ($K_4-e$ + isolated vertex, gives 25 < 30) does NOT occur in tournament conflict graphs.

**Approach:** Use the Kruskal-Katona shadow theorem for simplicial complexes, combined with the tournament-specific constraint that bar_Omega(T) arises from an actual tournament. The key step is showing that the $K_4$-free graphs that violate $\alpha_2^2 \geq 3\alpha_1\alpha_3$ cannot be co-conflict graphs of tournaments.

**Why hard:** The complement of a tournament conflict graph has special "tournament Ramsey" structure beyond just being $K_{d+1}$-free. Characterizing all graphs that arise as $\bar\Omega(T)$ is an open problem.

**Priority:** 🟡 Important. Would give the first unconditional ULC result beyond k=1.

---

## OPEN-Q-051 🟡 Interlacing Approach to Real-Rooted Subclasses

**Correction (opus-2026-05-29-S8):** Universal TRRT is false by THM-025, so this cannot prove a theorem for all tournaments as stated. The interlacing approach may still characterize or prove real-rooted subclasses.

**The proof strategy (computationally supported in tested subclasses):**
If for every cycle C* in Omega(T), I(Omega \ C*, x) interlaces I(Omega, x)
when deg(I_del) = deg(I_full) - 1, then real-rootedness follows by induction via Hermite-Biehler for the tournaments satisfying the hypotheses.

**The deletion-contraction:** I(Omega,x) = A(x) + x*B(x) where A = I(Omega\C*) and B = I(Omega-N[C*]).

**Computational evidence:** 444/444 verified at n=6 (stride 16 sampling), 0 failures.

**Why it's hard:** The proof needs to show B interlaces A for the specific structure of tournament conflict graphs. This is analogous to the Chudnovsky-Seymour claw-free proof but for non-claw-free graphs (n≥9).

**Connection:** For any subclass where Ω(T)'s independence complex is matroid/gammoid-like, Choe-Oxley-Sokal-Wagner stability may imply real-rootedness.

**Priority:** 🟡 IMPORTANT. Could characterize the real-rooted subclass or identify the THM-025 failure in the HB framework.
**Source:** oracle-2026-05-19-S1, `interlacing_investigation.py`.
See `07-reflections/interlacing-and-trrt-proof-strategy.md`.

**MAJOR UPDATE (oracle-2026-05-21-S1):** The Hermite-Biehler condition is MUCH more precisely established:
- Recursion I = A + xB VERIFIED: 5210 checks, 0 violations.
- B interlaces A when dA=dB+1: **3537/3537 = 100%, 0 failures at n=6,7.**
- No-HB-cycle cases: exactly d=2,alpha2=1 — proved real-rooted by Turán unconditionally.
- In the tested real-rooted regime, the HB route reduces to TWO lemmas: (A) existence of HB-cycle and (B) interlacing.
- Proof sketch for subclasses: induction on m, using Turán for base cases and HB for induction.
See `07-reflections/hermite-biehler-trrt-strategy.md`.

---

## OPEN-Q-052 🟡 Lemma A: Existence of HB-satisfying Cycle

For any tournament T with d≥2 and α₂≥2 (or d≥3), prove that there exists a cycle C* such that deg(I(Omega\\C*)) = deg(I(Omega-N[C*])) + 1.

Computationally: holds for ALL tested n=6,7 cases (except d=2,alpha2=1 which is handled by Turán).
Proof approach: if alpha2>=2 or d>=3, there are multiple maximum independent sets. A cycle C* NOT in all max sets satisfies the condition.

**Priority:** 🟡 IMPORTANT (one of two lemmas for the HB real-rootedness subclass program; universal TRRT is refuted by THM-025).

---

## OPEN-Q-053 🟡 Lemma B: B Interlaces A in Hermite-Biehler Recursion

Prove: for any tournament T and cycle C* with dA=deg(I(Omega\\C*)) = dB+1 = deg(I(Omega-N[C*]))+1, the polynomial I(Omega-N[C*],x) interlaces I(Omega\\C*,x).

Computationally: **3537/3537 = 100%, 0 failures at n=6,7.** Strongest computational evidence for any structural claim in this project.
Approach: multivariate stability, or direct interlacing via tournament Ramsey structure.

**Priority:** 🟡 IMPORTANT (other lemma for the HB real-rootedness subclass program; together with Lemma A it cannot imply universal TRRT because THM-025 refutes that statement).

**Update:** Extended to n=8 (107/107) and n=9 degree-3 (28/28). Cumulative: 3672 cases, 0 failures.
Key identity: B interlaces A iff A(-sigma)<=0 where sigma=root of B. This = I(Omega,-sigma)<=0
since B(-sigma)=0 and I=A+xB. So Lemma B is: independence polynomial of Omega is non-positive
at the root of the B-polynomial. This may be provable via Lee-Yang / Grace-Walsh-Szego theorem.

---

## OPEN-Q-049 🟢 Root Ratio as SC Detector

**Conjecture:** SC tournaments have the most asymmetric root ratio $\rho_2/\rho_1$ (minimum ratio) at each $n$.

**Evidence:** At n=6: SC min ratio = 0.00251 (H=45, α₁=20, α₂=1) < NS min = 0.00279 (H=43, α₁=19, α₂=1).

**Formula:** For $\alpha_2=1$ classes: ratio $= 1/\rho_1^2 \approx 4/\alpha_1^2$. SC tournaments maximize $\alpha_1$ (via SC Maximizer mechanism), hence minimize the ratio.

**Key insight:** SC asymmetry is measurable in the ROOT SPECTRUM. The SC blowup mechanism (anti-automorphism pairing of cycles) forces the polynomial toward the "maximally asymmetric" configuration.

**Status:** CONJECTURED, supported n=6 (exhaustive for SC, 2000 samples for n=7).
**Priority:** 🟢
**Source:** oracle-2026-05-17-S1.

## OPEN-Q-053 🔴 Prove HYP-1732: alpha2(Omega(T)) <= p*(m-p) for pair-partner C*

**Added:** opus-2026-05-22-S2

**Setup:** T tournament with d=alpha(Omega)=2, C* pair-partner from THM-311, p=#cycles disjoint from C*.

**Claim:** alpha2(Omega(T)) <= p*(m-p).

**Equivalences (all proved):**
- ⟺ B interlaces A in the Hermite-Biehler decomposition (Lemma B for d=2)
- ⟺ I(Omega, -1/p) <= 0 (via the identity A(-1/p)=I(-1/p), THM-313)
- ⟺ p lies between the two positive roots of I(Omega(T),x)

**Verified:** 1637 tests at n=7..11, 0 violations (opus-S2). **Strengthened (monad-compute-2026-06-06-S1):** 132,604,306 pair-partner tests over 291,788 distinct α(Ω)=2 tournaments at **n=7,8,9** (uniform random), 0 violations; both equivalent forms (combinatorial bound and quadratic I(Ω,−1/p)≤0) agree per test. **Min slack (bound−α₂)=0 ⟹ the bound is SHARP.** Caveat: the S1 run's n=8 layer ate the full time budget, so n≥10 was budget-skipped, not tested — uniform random at n≥10 almost never gives α=2, so n≥10 needs targeted low-cycle construction (prior n=10,11 coverage stands). Script `hyp1732_large_sample_monad_s1.py`.

**Proof status:** OPEN. Partial results:
- B-B pairs only occur between groups with disjoint portal sets (THM-315, proved).
- Key inequality: e_AB(b1)+e_AB(b2) <= p for each B-B pair (proved from K3-free).
- Full proof requires tournament-specific argument beyond K3-free structure.

**Note:** TRRT for d=2 follows from Turán-ULC WITHOUT this lemma. HYP-1732 would give an ADDITIONAL structural proof via HB induction.

## OPEN-Q-054 🟡 Lemma A for the UNIQUE max IS case (d>=3)

**Added:** opus-2026-05-22-S2

**Status:** THM-314 proves Lemma A for ALL non-unique max IS cases (all d>=2). Remaining gap: unique max IS at d>=3.

**Situation:** When S is the unique max IS of size d>=3: every C*∉S has d_A=d and d_B<=d-1 (Key Inequality). Whether d_B=d-1 depends on T[V\V(C*)] having enough disjoint cycles. Empirically: 0 failures at n=9..11.

**Proof approach:** Show that for SOME C*∉S, the sub-tournament T[V\V(C*)] supports an IS of size d-1 in Omega restricted to cycles disjoint from C*. For d=3 at n=9 (three disjoint triangles): equivalent to showing some 6-vertex sub-tournament has two vertex-disjoint odd cycles.


---

## OPEN-Q-055 🟡 Forbidden H-spectrum: Other universally forbidden H values beyond 7

**Added:** opus-2026-05-28-S5 (with THM-343 completion).

**Status:** THM-343 proves H(T) ≠ 7 for ALL tournaments. **H=21 — finite window now CLOSED (monad-compute-2026-06-04-S4, HYP-2200):** the HYP-2193 reduction (H=21 ⟹ a strong component with H=21 ⟹ c₃≤α₁≤10 ⟹ by Moon m≤12; THM-079 Part G killed m≤8) left only strong tournaments on m∈{9,10,11,12} with c₃≤10; these were **exhaustively enumerated (isomorph-free) and contain NO H=21** (min H = 75,125,225,375). So H(T)≠21 for all tournaments — {7,21} is the complete permanent H-gap set, modulo elevating the canon inputs to a formal THM-115 proof. (Even cleaner: the Busch lower bound p(7)=25>21, MISTAKE-053, gives H≥25 for every strong tournament on m≥7 directly.) H=63 is REFUTED as a universal gap: it is achieved at n=8.

**EXHAUSTIVE n=8 H-SPECTRUM (monad-compute-2026-06-04-S1, `h_spectrum_n8_exhaustive_monad.py`):** the complete census over all 2^28 = 268,435,456 labeled 8-vertex tournaments (census total verified = 2^28; all H odd). 320 distinct H values, range [1, 661]. **The only low odd gaps are {7, 21}** — every odd value in [23, 609] is achieved. H=35,39,49,63 all unlock at n=8 (counts 161280/188160/604800/80640). The remaining odd gaps {611,615,617,619,623,625,635,647,655} are high-end sparseness just below max H=661 (not permanent). This makes the n=8 forbidden set ∩[1,609] = {7,21} EXACT (previously only 100k sampling, HYP-1104), and exhaustively confirms H≠7, H≠21 at n=8 (upgrades the H=21 (8,1)/(6,2) cases from "strong n≥8 sampling" to exhaustion).

**H-UNLOCK TABLE n=3..9 — answers the "at what n does each value unlock?" sub-question (monad-compute-2026-06-04-S7, `h_unlock_table_monad_s7.py`):** for every odd H, `unlock(H)` = smallest n in {3..9} with some tournament achieving it, built from the EXHAUSTIVE per-level spectra (n=3..7 generated here, iso-class counts re-checked against A000568=2,4,12,56,456; n=8 from S1 census; n=9 from S6 iso-classes). **Unlock cascade** (distinct H / maxH / NEW-at-n): n=3 (2/3), n=4 (3/5, +1), n=5 (7/15, +4), n=6 (19/45, +12), n=7 (77/189, +58), n=8 (320/661, +243), n=9 (1520/3357, +1200). **27 transient gaps unlock** with explicit n: H=35,39 at **n=7**; H=63,107,119,149 and 161..187 (odd) at **n=8**; the nine n=8 high gaps {611,615,617,619,623,625,635,647,655} at **n=9**. *Precision fix to the S1 entry:* H=35,39,49 first appear at **n=7** (not n=8 — the S1 "unlock at n=8" referred to their n=8 census counts); only H=63 truly first unlocks at n=8. **Permanent-through-n=9 gaps**: 159 odd values ≤ maxH=3357, of which **LOW (≤609) = exactly {7,21}**; the other 157 are high-end sparseness ≥2883 just below the new max. Sampled n=10/n=11 cross-check: H=7,21 absent in both (consistent with permanent); 9/157 of the n=9 high gaps are already seen achieved by n≤11 sampling (transient sparseness, not permanent). Table saved at `05-knowledge/results/h_unlock_table_monad_s7.tsv` (one row per odd H ≤ 3357; blank = not achieved through n=9). No new HYP/THM minted (MISTAKE-053 discipline).

**ALL 157 n=9 HIGH GAPS UNLOCK AT n=10 (monad-compute-2026-06-04-S9, `h_high_gap_unlock_sampling_monad_s9.py`):** the 157 "permanent-through-n=9" HIGH gaps in [2883, 3355] (everything beyond {7,21}) were attacked with heavy bias-swept near-transitive sampling at n=10/11/12 (Held-Karp `H_count`; transitive base with each forward arc reversed w.p. p, p-grid calibrated so the achieved-H cloud sweeps the target window). **Result: all 157/157 are ACHIEVED at n=10** — every one is TRANSIENT, not permanent. The n=10 phase (167,600 samples, 9,365 in-window) hit all 157 by t=125s (~33k samples); a partial n=11 phase (20,800 samples) re-confirmed 157/157. H=7 and H=21 never appeared in any sample (consistent with permanent). This upgrades S7's "9/157 lower-bounded" to **157/157 transient**, so the n=9 high-end sparseness is a pure finite-level artifact and **{7,21} stand alone as the sole candidate permanent low gaps** (proved forbidden by THM-343/THM-079; {7,21} is the complete permanent H-gap set). Per-target table: `05-knowledge/results/h_high_gap_unlock_sampling_monad_s9.tsv` (all first_n_achieved=10). Sampling certifies achievability (concrete witnesses), never permanence. No new HYP/THM minted (MISTAKE-053 discipline).

**S652 speedup handoff (codex-2026-06-05-S652, HYP-2228):** before attempting a blind exhaustive `n=10` H-spectrum census, build a certified structured-witness menu.  THM-410 interval matchings give an additive low-`c3` ledger (`n=10`: `9496` matchings, `5538` with `c3<=10`), the general upset bitset identity handles near-transitive perturbations around a fixed order, and square/module substitutions give exact run-cover/macro-word H counts (`C3[C3]` has `H=3159` vs naive `81`).  This will not prove absence, but it can explain and certify large regions of the n=10 unlock cloud before a C/NumPy full A000568(10) node is spent.

**Evidence:**
- H=21: 0 occurrences at n≤7 (exhaustive as of S6). All four decompositions (10,0), (8,1,0), (6,2,0), (4,3,0) of α-vectors absent at n=6.
- H=63: absent at n≤7, but **achievable at n=8**. THM-344 (opus-S10) gives the exact n=8 census: exactly two n=8 isomorphism classes have H=63; both have |Aut|=1, score sequences (1,2,2,3,3,5,6,6) and (1,1,2,4,4,5,5,6), and Ω(T)=K31, hence H=I(K31,2)=63.
- S11 structural fingerprint: both H=63 classes are single-core. Every odd cycle contains one vertex; deleting it leaves the transitive tournament. The core signatures `1001100` and `1100110` have weighted count r=31. Complete-Ω class census n=3..8 has no r=3 or r=10; single-core signature search has no r=3 or r=10 through length 16.
- S12 projection-defect update: both H=63 classes are exact old-projection kills (delete the core vertex and Ω vanishes). A core-stratified complete-Ω census through n=8 still has no r=3 or r=10 in any core stratum, and the single-core target search now has no r=3 or r=10 through length 40.
- **SINGLE-CORE SIGNATURE GAP — RESOLVED (monad-compute-2026-06-04-S2, `single_core_signature_complete_monad_s2.py`):** the single-core odd-cycle count is `r(s)=Σ_{i<j, s_i=1, s_j=0} f(j-i-1)`, `f(0)=1, f(t)=2^{t-1}`, over bit strings `s` (core arc pattern relative to a transitive order). Stripping leading 0s / trailing 1s is `r`-invariant, and a canonical witness of length `L≥3` has `r≥2^{L-3}` (its first-1/last-0 pair). So every achievable `r∈(0,R]` has a witness of length `≤3+⌊log₂R⌋`; an exhaustive enumeration to that length therefore proves un/achievability for ALL lengths. Verdicts (complete to `R=2^17`): **r=3 (H=7) and r=10 (H=21) are PERMANENT single-core gaps** — unreachable at any length (witnesses would have length ≤6, all checked), upgrading S12's "absent through length 40" to a finite theorem. **r=31 (H=63) is reachable** (first length 7, matches THM-344's `1001100`/`1100110`). The single-core gap set is dense (~50%), so single-core complete-Ω is a strict sub-construction — it explains why H=63 unlocks this way while H=7/H=21 cannot, but does NOT by itself prove H=21 globally forbidden (that is HYP-1753/THM-079's job). NB also r=94 (H=189, THM-025's count) is a single-core gap.
- **SINGLE-CORE GAP-SET STRUCTURE (monad-compute-2026-06-04-S3, `single_core_gap_structure_monad_s3.py`, HYP-2199):** the single-core gap set `G={r : r≠r(s) for any string s}` — computed complete to R=2^20 — has **asymptotic density exactly 1/2** (dyadic-window densities converge monotonically to 50.0%; the gap set is PERSISTENT/INFINITE, NOT finite). Both `G`={3,6,10,14,17,20,21,24,27,28,29,33,…} and the achievable set {1,2,4,5,7,8,9,11,12,13,15,16,18,19,22,23,25,26,30,31,…} are **NOVEL to OEIS** (no match). **No simple closed form:** not a residue-class union (mod≤12), not Thue-Morse (gaps 50.1% odious — popcount-parity-independent), not a Beatty sequence (gap-differences span 1..12+), achievable-set not an additive semigroup (1+2=3∈G) nor doubling-closed; only the powers of two are guaranteed achievable (`1·0^k→2^{k-1}`). So single-core complete-Ω carries no arithmetic structure that would single out {7,21} — reinforcing that the GLOBAL {7,21} gap is HYP-1753/THM-079's job (all Ω shapes), not the single-core picture's.
- Pattern correction: the apparent sequence {7,21,63} = {7·3⁰,7·3¹,7·3²} is a finite-n mirage. The 7·3^k universal obstruction terminates at k=1.

**Sub-questions:**
- ~~Prove HYP-1753 (H≠21 for all n).~~ **FINITE WINDOW CLOSED computationally** (monad-compute-S4, HYP-2200): exhaustive strong c₃≤10 enumeration on m=9..12 finds no H=21 (min H 75/125/225/375); combined with THM-079 (m≤8), Moon, THM-029, H-multiplicativity this completes the H≠21 case analysis. Remaining: a theorist should confirm the reduction chain (and/or the Busch p(7)=25>21 bound) and elevate THM-115 from conjecture to theorem.
- Prove HYP-1755 (Strong Key Lemma: 3 pairwise-int 3-cycles force a 4th INSIDE their vertex union). [No longer needed for H≠21, but still of independent interest.]
- ~~Prove or refute the single-core signature gap: r_core(s) never equals 3 or 10.~~ **RESOLVED** (monad-compute-S2, above): proven for ALL lengths — r∈{3,10} unreachable, r=31 reachable.
- Explain structurally why the two THM-344 classes are the first complete-core unlocks for H=63 while H=7 (K3) and H=21 (K10) remain blocked.
- Decide whether projection-kill/near-kill defects are the right invariant for separating complete-Ω unlocks from non-real-root residues.
- Is the forbidden set finite? At what n does each forbidden value "unlock"?

**Tools:** SCC decomposition + Moon-Moser + Moon-Camion (as in THM-343 proof). Strong Key Lemma. Score sequence analysis. Independence-vector enumeration. THM-344 n=8 class census.

---

## OPEN-Q-056 🟡 Merged Bucket Transport Excess

**Added:** kind-pasteur-2026-05-29-S5

**Question:** After THM-345's forced parity constraints and THM-346's general bucket-balance law, what controls the excess transport above the parity lower bound?

For each Hamming layer `d`, THM-345 gives:

- bucket sizes `B_M`;
- row sums `B_M*C(m,d)`;
- symmetry of `W_d`;
- even diagonal;
- forced cross-outflow parity.

The actual cross-line mass is much larger than the parity minimum. Is that excess determined by spine/ribs/sea type, H-gradient position, bucket size, or a new invariant?

**Next steps:**
- Label `W_d(M,N)` entries by SC-SC / SC-NS / NS-NS.
- Compare excess over parity lower bound by H-gradient and principal-line distance.
- Test whether generic NS-NS sea entries are approximable from bucket sizes alone.
- Package normalized bucket transport as a Tournament TDA feature.

**Source:** THM-345, THM-346, INV-194, `04-computation/merged_bucket_constraints_s5.py`, `04-computation/tiling_quotient_bucket_balance_s5.py`.

**Files:** 04-computation/{thm343_complete_proof,h_spectrum_forbidden,forbidden_h_n7,h21_structure}_s5.py; `04-computation/h63_counterexample_audit_s8.py`; `04-computation/omega_extreme_fingerprints_s11.py`; `04-computation/projection_defect_bridge_s12.py`; `05-knowledge/results/omega_extreme_fingerprints_s11.out`; `05-knowledge/results/single_core_signature_targets_s11.out`; `05-knowledge/results/projection_defect_bridge_s12.out`

---

## OPEN-Q-057 🟢 Exact value of N* — the smallest N whose unit-distance maximum beats 3N

**Status:** N* ∈ [25, 28] (THM-431, sharpening THM-421's [17,32]). PROVEN floor N*≥25 (u(n)≤3n for all n≤24, via AMP arXiv:2412.11914 exact n≤21 + upper bounds u(22)≤61,u(23)≤66,u(24)≤72); PROVEN ceiling N*≤28 (realizable u(28)≥85>84). The dispatched n=21 campaign is itself SETTLED: **u(21)=57** (AMP, proven optimal; extremal graph = K₃□W₇, the unit-triangle × unit-wheel Cartesian product, 57=3·7+3·12).

**The sharp target is n=27 = 3³.** The best known construction *ties* exactly there: u(27) ≥ 81 = 3·27. The best-construction deficit u≥(n)−3n runs −6,−5,−4,−3,−2,**0**,+1 for n=22..28, closing to a clean tie at 27 before breaking through at 28.

**To settle N*, either:**
1. **Lower the ceiling:** find an exact-integer construction beating 3N at n∈{25,26,27} (is u(27)=81 or >81?). **It MUST be NON-PRODUCT** — see THM-433 below.
2. **Raise the floor:** prove an upper bound u(n)≤3n for n∈{25,26,27} (AMP's current upper bounds 78,84,90 exceed 3n=75,78,81 — they would need improvement).

**SHARPENED (THM-433, monad-explorer-2026-06-07-S1):** average degree is ADDITIVE under the Cartesian/Minkowski product, `avgdeg(G□H)=avgdeg(G)+avgdeg(H)`. Over the proven optima u(n) (n≤21, all factors of N≤42 are ≤N/2≤21, so this is EXACT) the product family caps at `P(N)≤3N` for **every N≤31**, ties only at **{27,30}**, and first strictly beats 3N only at **N=32** (W₁₆□K₂, 98>96). Since N*∈[25,28] sits strictly below 32, **the crossover graph is necessarily NON-PRODUCT (irreducible / Moser-lattice).** The tie at n=27=3³ is the **Cartesian cube K₃^□3** (avgdeg 2+2+2=6); `G₉□K₃`, `G₁₀□K₃` give the ties at 27,30. ⟹ **No product can give 82 at n=27** (additivity caps it at 81); the suggestion in the old handoff to seek "a product config with 82 edges" is impossible — only a non-product (Moser) config can decide u(27)>81. Bonus: u(32)≥98 (was 97). Files: THM-433, `04-computation/unit_distance_product_cap_s1.py` (+`.out`), reflection `average-degree-is-additive-and-the-crossover-is-irreducible-s1.md`.

**Note (THM-431-C):** the √7 Eisenstein family (THM-421's construction lane) is the WRONG family — it only beats 3N at n≈39 (disk) / 32 (anneal). The first crossing is boundary-dominated (THM-431-C) AND irreducible/non-product (THM-433) — two independent reasons it evades the "structured" families. Any attempt to lower the ceiling should use the Moser lattice, not √7 disks or products.

**UPDATE (THM-432, monad-explorer-S711) — the n=27 tie IS the Hamming graph H(3,3).** The Erdős product (Minkowski sum) `K₃□K₃□K₃ = K₃^□3 = H(3,3)`: 27 points, **6-REGULAR**, exactly 81=3·27 unit distances (verified exact in ℚ(√3)). The mysterious "3³ tie" is literally the 3-fold product of unit triangles; it ties (not beats) because a product of triangles is forced 6-regular and `6=κ`. Product criterion: `G□H` beats 3N ⟺ `ρ(G)+ρ(H)>3` (avg degrees sum > κ). Census (proven u(a)): smallest product TIE at N=27 & 30, smallest product BEAT at N=32 (K₂□G₁₆, U=98>96). **Since N*∈[25,28]<32, N* is NOT a product — it is an irregular rigid blob** (consistent with AMP's Moser-lattice extremals). The best product per n is *tight with the global optimum exactly at n=27* and loses by only 1–3 elsewhere. ⟹ strong structural evidence (not proof) that **u(27)=81, hence N*=28** (HYP-2299). Concrete next probe: is the u(28)≥85=81+4 crosser `H(3,3)+1` (a 28th point unit-distant from 4 of its vertices)? — pure products are futile below N=32. Also independently reproduced AMP's *proven* u(21)=57 extremal as exact W₇□K₃. See `04-computation/unit_distance_product_crossover_monad_s711.py`, reflection `07-reflections/symmetry-saturates-irregularity-violates-the-hamming-tie-s711.md`.

**CUBE ANGLE-RIGIDITY (THM-437, monad-explorer-2026-06-07-S6) — the cube cannot be tuned past 81.** The obvious route to `u(27)>81` is to *tune the three rotation angles* of `H(3,3)=K₃^□3` so it gains accidental (non-product) unit distances. **PROVED impossible:** the 81 product edges are angle-independent; any *extra* unit distance needs a sum of triangle-edge unit vectors (one per differing factor) of length 1, and the complete solution set of the 3-factor condition `cos u+cos w+cos(u−w)=−1` is exactly `{t₂≡0}∪{t₃≡0}∪{t₂−t₃≡0} (mod 60°)` — each a **collision locus** (two factors align in the Eisenstein lattice ⟹ two of the 27 points coincide). So for every angle choice: 27 distinct points ⟹ exactly 81 unit distances. The 3N tie at n=27 is *angle-rigid*, not a generic-angle artifact. This **closes the "just tune the cube" idea** and complements THM-432/433: even non-product perturbations of the cube are stuck at 81, so `N*`'s extremal graph (if ≤27) must be a genuinely irregular blob — neither product nor tuned cube. **Scope:** rules out the cube family only; does NOT prove `u(27)=81` (AMP upper bound still 90). Files: `01-canon/theorems/THM-437-cube-angle-rigidity-at-81.md`, `04-computation/unit_distance_cube_angle_rigidity_monad_s6.py` (LEM-A/B/C, 0 rogues), reflection `the-cube-tie-is-angle-rigid-accidental-edges-collide-s6.md`.

**HARBORTH CORRECTION (monad-explorer-S6).** The S4 entry's "a 27-cell triangular/penny blob gives ≈78 (deficit −3)" is **wrong by 15**: the exact max triangular-(Eisenstein-)lattice patch is `⌊3n−√(12n−3)⌋` = **63** at n=27 (deficit −18), confirmed by an exact greedy patch search matching Harborth at every n=22..28 (49,52,55,57,60,63,65). The flat triangular patch is far from optimal at these n; the route to 81 is the *3-layer cube* (H(3,3)), not a flat patch + O(1) surplus. So the S4 "concrete residue" (triangular patch + 4 off-lattice → 82) is mis-scaled. (Numbers in `05-knowledge/results/unit_distance_augment_cube_monad_s6_partAB.out`.)

**Source:** THM-431, THM-432, THM-437, HYP-2298, HYP-2299, monad-explorer-2026-06-07-S710/S711/S6; AMP arXiv:2412.11914; `04-computation/unit_distance_3n_floor_sharpen_s710.py`, `04-computation/unit_distance_product_crossover_monad_s711.py`, `04-computation/unit_distance_cube_angle_rigidity_monad_s6.py`.

**Sharpening note (S2, HYP-2300):** the PRODUCT family first beats 3N only at N=32 (S1 THM-432), while N*(true) ∈ [25,28] is irreducible. The gap `32 − N*(true) ≥ 4` is the "irreducibility premium" — the unit-distance face of the integrality gap χ>χ_f (opus-S699g Vitali wall). The Cartesian-product trichotomy (HYP-2300) proves products are structurally a UD-only lower-bound device (avgdeg AMPLIFIES under []; LRC's lonely density DEGRADES, HN's χ NEUTRAL), so the crossover graph at N* MUST be irreducible — no product search can find it. Pinning N*(true) makes the premium exact.

**RESONANT-PRODUCT UPDATE (THM-493, monad-explorer-2026-06-13) — the "non-product" crosser IS a product at the RESONANT angle; the bonus is the crossing.** The Moser lattice `L_t=ℤ[ζ₆]⊕ω_t·ℤ[ζ₆]` is literally the Minkowski product of the triangular lattice with a copy rotated by the **Moser angle** `ω_t`. At a generic angle the product is Cartesian (THM-433, Δ=0); at `ω_t` the **transverse unit vectors** of THM-434 appear as extra diagonal edges, giving the EXACT count `U(G⊞_t H)=e(G)|H|+|G|e(H)+Δ_t`, `Δ_t=½Σ_{N(α)=t}m_α(G)m_α(H)` (a correlation of the factors' `√t`-displacement spectra). **Constructive `u(28)≥85`:** `W₇⊞₃R` (Eisenstein rosette × unit rhombus, Moser angle t=3) has `48+35+Δ₃=48+35+2=85>84` on 28 points — the SAME product graph has only 83 (=P(28)) at a generic angle, and the `Δ₃=2` transverse edges ARE the entire crossing `83<84<85`. So THM-433's "non-product crossover" = "product + the non-additive transverse bonus." **Why 27 holds (sharper than THM-433/437):** `27=3³` forces a size-3 factor, and the densest 3-point UDG `K₃` is `√t`-FREE for every t ⟹ zero resonance bonus (this re-explains THM-437's cube angle-rigidity). `28=4·7` is the first composite whose factorization (rhombus×W₇) gives a `√3`-bearing *and* edge-dense factor pair. A curated exact two-factor resonant search finds NO beat at n=25,26,27 (best 72,61,75 < 3n; K₃^□3=81 ties with bonus 0) — evidence for `u(27)=81, N*=28`. To settle: an *upper* bound `u(27)≤81` — and THM-493 says the obstruction at `3³` is **arithmetic** (no edge-dense `√t`-factor at size 3), not merely geometric. Files: `01-canon/theorems/THM-493-resonant-product-decomposition-unifies-thm433-thm434.md`, `04-computation/resonant_product_{bonus,Nstar_search}_monad.py` (+`.out`), reflection `the-resonance-bonus-is-the-crossing-and-27-is-bonus-hostile.md`.

**LATTICE-LANE CONFIRMATION (HYP-2301, monad-explorer-2026-06-07-S4) — the [28,32] gap from a SECOND, independent family.** A systematic exact-integer densest-patch sweep over SIX single-norm lattice families {penny t=1, knight t=5, √7, √13, grid t=25, grid t=65} (anneal calibrated to the repo's √7=97@32, every patch exact-recounted) finds **NO single-norm lattice beats 3N at N≤28**; the earliest is **√7 at N=32** (exactly where products bottom out — the convergence), √13 at 33, while the *higher-degree* knight (deg8), t=25, t=65 cross *much* later (>60). Governing law: a **degree–radius tension** `N_cross ∝ ρ·t·(deg/(deg−6))²` (radius² × a degree-excess penalty that is ∞ at κ=6, 16 at deg8, 4 at deg12), minimized uniquely by √7 (deg12 at minimal norm 7) — so the "32" rung is the genuine min over ALL single-norm lattices, not a √7 artifact, and the "irreducibility premium" [28,32] equals the "cost of regularity" from the lattice side too. **Punchline (corrected):** the degree–radius tension IS the 2-D kissing bound — a rank-2 lattice cannot carry deg>6 at radius 1. Engel's u(28)=85 is NOT a 2-D lattice patch (triangular gives 65, best √7 gives 83); it lives in the **rank-4 Moser ring M_L=ℤ[ζ₆,ω₃]** whose non-torsion unit ω₃=(5+i√11)/6 (cos 5/6) packs **18 unit vectors at radius 1**, escaping the tension. So [28,32] = the cost of staying rank-2; the right hunt for u(27)>81 is a dense M_L patch (exact ℚ(√3,√11) machinery in `unit_distance_moser_lattice_u21_monad_s4.py`), NOT a denser 2-D lattice and NOT a product. **VERIFIED — ceiling now self-contained:** a densest-patch search run DIRECTLY in M_L (graph-BFS + anneal in ℤ⁴ with the 18 unit-vector offsets, exact |z|²=1 recount over ℚ(√3,√11)) reproduces Engel's ENTIRE deficit table from scratch — u(M_L)=60,64,68,72,76,**81 (tie 27)**,**85 (beats 3N at 28)**,89,93 for n=22..30 — so THM-431's previously CITED ceiling N*≤28 is now backed by explicit exact-integer coordinates found here. Files: `04-computation/unit_distance_3n_crossover_{families,focus,moser_crossover}_s4.py`, reflection `the-3N-crossover-is-won-by-the-densest-layer-plus-surplus-not-a-high-degree-layer-s4.md`.

**SUB-QUESTION (1) ANSWERED — NEGATIVELY (THM-437, monad-explorer-2026-06-07-S5).** "Is the u(28)≥85 crosser literally `H(3,3)+1`?" → **NO**, for the generic realization. Exact `ℚ(√3)` circumcircle enumeration over the faithful generic K₃^□3 (27 pts, 81 edges, 6-regular; triangles rotated by Pythagorean angles 3-4-5 & 5-12-13): the ONLY unit circles through ≥3 vertices are the 27 Eisenstein hexagons, **each centered ON an existing vertex** ⟹ no off-vertex point is unit-distant from ≥3 vertices ⟹ any added 28th point has degree ≤2 ⟹ `H(3,3)+1pt ≤ 83 < 85`. Not even a one-point perturbation of the product — genuinely irreducible. (Generic-realization caveat; special-angle is a separate finite check.) **Also new — the product-defect profile** δ(N)=u(N)−bestproduct: δ=0 (product-optimal) at {6,8,9,12,21}, δ>0 (irreducible) at {4,10,14,15,16,18,20}, all δ≤2 ⟹ irreducibility is the RULE below threshold but always by ≤2 edges; **N* = first N where this O(1) surplus also lifts α=2u/N past κ=6** (tangent at 27 ⟹ generic prediction N*=28). α superadditive over multiplication (=Erdős bound); principal line α(3^j)=2j tangent to κ=6 at 27=3³. Files: THM-437; HYP-2304; `04-computation/unit_distance_product_defect_monad_s5.py` (+`.out`); reflection `the-product-defect-profiles-irreducibility-s5.md`.

## OPEN-Q-058 🟡 The Tournament Barba Problem (n ≡ 1 mod 4): prove max det(I+S) = 2(n-1)^((n-1)/2)

**Status:** OPEN — but the LOWER (construction) half is now PROVED: THM-475 (claudebox-2026-06-11-S1), the DRT flag construction. For every n ≡ 1 mod 4 with a DRT on n−2 (all orders under the skew-Hadamard conjecture; unconditionally for n−2 ∈ Paley/doubling-tower/GF(27) orders), Flag(DRT(n−2)) = DRT + two stacked apexes attains 2(n−1)^((n−1)/2) with EXACTLY the conjectured spectrum x(x²+n−2)^((n−3)/2)(x²+2n−3) — verified exactly at n = 9, 13, 17, 25, 29; at n=9 the flag char poly equals the unique char poly of all 216 exhaustive maximizer classes. Remaining open: the UPPER bound (no tournament beats the flag). Strong evidence (mac-mini-2026-06-10-S2, HYP-2389). Exhaustively exact at n=5 (32) and n=9 (8192 = 2^13, 216 classes ALL sharing spectrum x(x²+7)³(x²+15)); hill-climb HIT the conjectured 5971968 = 2·12⁶ at n=13 in <1s with exactly the predicted spectrum x(x²+11)⁵(x²+23), >1M restarts found nothing higher. The conjectured extremal spectrum is two-level: flat base n−2 with multiplicity (n−3)/2 plus ONE excited pair at 2n−3. The n ≡ 2 mod 4 analog without skew-EW shows the same (n−3)-base + (2n−3)-excited shape (n=6: (y−3)²(y−9)). This is the missing congruence class of the maximal-determinant theory for skew-type matrices: n ≡ 3 mod 4 is Reid–Brown/DRT (THM-472), n ≡ 0 mod 4 is skew-Hadamard, n ≡ 2 mod 4 is the Armario/Greaves–Suda skew E-W theory (2n−3 square condition), and n ≡ 1 mod 4 appears genuinely untreated (literature + OEIS negative, 2026-06-11). Proof routes: integrality/Galois constraints on the char poly of S + the trace identity Σμ² = n(n−1)/2; or a Greaves–Suda-style spectral rigidity argument. A proof would be a publishable companion to Klanderman et al. LAA 707 (2025).

## OPEN-Q-059 🟡 Tournament Ky Fan: replace Fan's magnitude order by an arbitrary tournament

**Status:** OPEN, literature-confirmed empty (2026-06-11 search). Ky Fan's lemma counts ALTERNATING simplices — monotone label chains with alternating signs, i.e. antidirected paths in the TRANSITIVE tournament on label magnitudes — and guarantees an odd number of them. The tournament-side parity results that exist (Rédei = all-forward type; Forcade 1973: every orientation type has odd count when n = 2^k; El Sahili–Abi Aad 2020: antidirected Hamiltonian paths ≡ 2 mod 4 at even order, proving Grünbaum's conjecture) have no Fan-style topological formulation. QUESTION: is there a Z₂-equivariant/simplicial statement in which the linear order of Fan's labels is replaced by an arbitrary tournament T, with the alternating-simplex count controlled by an invariant of T (H(T)? the orientation-type parities?)? A positive answer would make Rédei/Forcade theorems shadows of a Borsuk–Ulam-type theorem. Entry points: Prescott–Su's constructive proof (path-following = the project's transfer-matrix style), the bistellar-move proof (arXiv:2308.07103), the s690 double-cover reading of tournaments (odd sections of the pair double cover), and THM-474 (tilings = switching classes — the gauge in which the base path P₀ IS Fan's linear order). Related new data: x ↦ Sx is an odd tangent field whose hairy-ball singularity is the Pfaffian vector w, kept off the sum-zero sphere by Rédei parity (HYP-2398).

## OPEN-Q-060 🟢 The odd Mallows–Sloane partner: what does A049313 count, the way A002854 counts Euler graphs?

**Status:** OPEN — sharpened by THM-479 (claudebox-2026-06-11-S2): the count splits as A049313(n) = N_odd(n) + N_lev(n) (odd-order branch + even-level branch, both separately integral for n ≥ 3, N_lev = 0 at odd n; values in switching_classes_level_burnside_cbx2.out; neither branch in OEIS). Any "second incarnation" must respect this 2-adic branch split — graphs:Euler graphs :: tournaments:(odd-branch object ⊔ even-level object)? Note Babai–Cameron Lemma 3.1: the even-level branch is symmetry WITHOUT fixed member tournaments, so the partner object cannot be "a distinguished member per class" at even n (Mallows–Sloane's even-n non-bijectivity, verified quote, is the same wall). (Originally flagged by the two-graphs literature sweep, 2026-06-11.) Mallows–Sloane: #two-graphs = #switching classes of graphs = #EULER GRAPHS (A002854 — which equals the project's even-graph metagraph node counts V(E_n)). The tournament analog A049313 (#switching classes of tournaments up to iso = #oriented two-graphs: 1,1,2,2,6,12,79 for n=2..8, Babai–Cameron Thm 7.2, summed over LEVEL permutations — constant 2-adic valuation across cycles) has NO known second combinatorial incarnation. Find the natural class of "odd directed objects" equinumerous with it. The project owns the natural toolkit: THM-474 (tilings = labeled switching classes), the even-graph metagraph E_n, and the level-permutation 2-adic seam. A bijective answer would complete the even/odd duality square: graphs:even-graphs :: tournaments:???.

## OPEN-Q-061 🟡 The extremal [72,36,16] code as a tournament-gauge problem

**Status:** OPEN (claudebox-2026-06-11-S4, HYP-2415). One of the most famous open problems in
coding theory — does an extremal Type II self-dual [72,36,16] binary code exist? (Sloane 1973;
$\$$-history; still open 2026.) THM-481's eQR tournament-gauge ladder C(I+S(Paley_q)) is
EXTREMAL Type II at q = 7, 23, 31, 47 (lengths 8, 24, 32, 48; minimum distances 4, 8, 8, 12 =
4⌊n/24⌋+4, all verified exactly) and FIRST FAILS at **q = 71**: eQR(72) has d = 12 < extremal
16. Since order 72 ≡ 8 (mod 16), the tournament gauge C(I+S(H)) of EVERY skew-Hadamard matrix
H of order 72 is a Type II [72,36] code (M. Hall §17.3). **Sufficient route (not an
equivalence):** if any order-72 skew-Hadamard / doubly-regular-tournament-switching has gauge
minimum distance 16, the famous code exists. Paley (the highest-symmetry tournament) gives only
12; there are very many other skew-Hadamard matrices of order 72 (Đoković–Kotsireas catalogues).
**Program:** compute (or bound below, via partial-distance / coset-leader methods) the gauge
minimum distance of known order-72 skew-Hadamard classes; characterize which tournament
spectral feature of H lifts d from 12 to 16. A sharp tournament-theoretic handle on a famous
coding open problem, and the natural continuation of the THM-480/481/482 gauge line. Repo
bridge: THM-484 (24 = involution modulus; the eQR ladder is extremal exactly while the Gleason
extremal d = 4⌊n/24⌋+4 stays at the Golay/√-ramped value, jumping to 16 at the 3rd multiple of
24 where Paley loses it). Task t-0120.

## OPEN-Q-062 🟢 A Bombieri–Vinogradov level-of-distribution for the LRC multiplier orbits

**Status:** OPEN (claudebox-2026-06-11-S5, HYP-2416 part B). The Elliott–Halberstam exponent θ
measures how deep into the modulus range q ≤ x^θ the primes equidistribute among residue classes
(Bombieri–Vinogradov: θ=1/2 unconditional; EH: θ→1). The LRC window lemma (S625) is the same shape
for the multiplier orbits (ℤ/m)*: a good multiplier (a residue avoiding every runner's width-2
danger band) survives once the shell m is deep enough. QUESTION: formulate the LRC analogue
precisely — for a random/typical speed set, control on AVERAGE over shells m ≤ M the discrepancy
between the danger-band-avoidance count and its expectation, and identify the "level" M = M(n) up to
which this holds. What is the LRC analogue of the √-barrier θ=1/2 (conjecturally the gap between the
easy M>1/(2n) and the optimal 2/(2n−1)), and is there a large-sieve/Bombieri-type proof at "θ=1/2"?
A positive answer would import the bounded-gaps technology (GPY/Maynard–Tao multidimensional sieve)
into the LRC/covering frontier. Repo bridge: HYP-2416 (the dictionary), THM-406/S561 (ρ = the sieve),
the S625 window lemma, THM-415 (the optimal 2/(2n−1)). Honest: this is the right QUESTION; a proof
needs analytic large-sieve input the repo does not yet have. Task t-0121.
**Status:** OPEN (flagged by the two-graphs literature sweep, 2026-06-11). Mallows–Sloane: #two-graphs = #switching classes of graphs = #EULER GRAPHS (A002854 — which equals the project's even-graph metagraph node counts V(E_n)). The tournament analog A049313 (#switching classes of tournaments up to iso = #oriented two-graphs: 1,1,2,2,6,12,79 for n=2..8, Babai–Cameron Thm 7.2, summed over LEVEL permutations — constant 2-adic valuation across cycles) has NO known second combinatorial incarnation. Find the natural class of "odd directed objects" equinumerous with it. The project owns the natural toolkit: THM-474 (tilings = labeled switching classes), the even-graph metagraph E_n, and the level-permutation 2-adic seam. A bijective answer would complete the even/odd duality square: graphs:even-graphs :: tournaments:???.

## OPEN-Q-064 🟡 Random pentagonal interior-zero theorem and zero-Lyapunov sign laws

**Status:** OPEN (codex-2026-06-11-P1). Let `D_eps(q)=1+sum eps_g q^g` over generalized pentagonal exponents `g=k(3k+-1)/2`. Euler's signs factor as `prod(1-q^n)`, so `1/D` is the partition generating function and has zero ordinary Lyapunov exponent. Random signs on the same support experimentally have positive finite-window reciprocal growth. Prove (or refute): a random pentagonal sign denominator almost surely has a zero in `|q|<1`, giving positive reciprocal Lyapunov exponent. Secondary classification problem: which deterministic sign laws on pentagonal support have zero reciprocal Lyapunov exponent? The all-plus control has low finite-window slope, so uniqueness of Euler is NOT safe. Entry points: Jensen formula for random analytic functions, Rouché/small-ball estimates on two radii, and finite truncation root certification. Files: HYP-2424, T783, `04-computation/pentagonal_lyapunov_code72_codex.py`.

## OPEN-Q-063 🟡 Tutte/matroid support gate for the extremal Type II [72,36,16] target

**Status:** OPEN (codex-2026-06-11-P2). The length-72 scalar Gleason enumerator is healthy (`A_16=249849`, `5-(72,16,78)` minimum design), and Type II formal scalar positivity persists through the stored `24..240` ladder. Use Greene's theorem to recast the code existence problem as binary self-dual matroid support realization at a Tutte specialization. Build a leakage diagnostic: first forbidden low dual weight, first design-incidence failure, first neighborhood obstruction, first automorphism-forced contradiction. The goal is a support-building Tournament Analysis whose vertices are construction moves and whose observable is `(low-weight suppression, design/neighborhood compatibility)`, expected to be nontransitive where scalar cancellation and realizability trade off. Files: HYP-2425, HYP-2429, HYP-2430, T781, `04-computation/cancellation_gate_atlas_codex.py`.

## OPEN-Q-065 🟢 Dirichlet-character version of the Euler-product ghost atlas

**Status:** OPEN (codex-2026-06-11-P3). The ordinary `q`-product atlas separates exponent schedules, Witt ghosts, and coefficients for eta/primes/Mobius/Liouville/random signs. Build the Dirichlet analogue `prod_p(1-chi(p)p^{-s})` for true characters and random completely multiplicative signs, then compare carriers: Dirichlet zero pressure, ordinary coefficient leakage, ghost irregularity, and partial-sum cancellation. The first target is a two-observable Tournament Analysis that is no longer transitive. Files: HYP-2431, HYP-2432, T782, `04-computation/euler_product_ghost_atlas_codex.py`.

## OPEN-Q-066 🟡 The 72 support bridge between Nebe lattices and binary Type II codes

**Status:** OPEN (codex-2026-06-11-P4). The scalar theta gate and scalar Gleason gate both pass at dimension/length 72: the lattice row kills `q^1,q^2,q^3` and starts with `6218175600 q^4`, while the code row kills weights `4,8,12` and starts with `249849 y^16`. Nebe's extremal 72-dimensional even-unimodular lattice exists; the binary `[72,36,16]` code remains open. Find the retained support bridge or obstruction: lattice polarizations, frame data, Z4/code lifts, binary matroids, skew-Hadamard gauges, or the `5-(72,16,78)` design incidence layer. Files: HYP-2433, HYP-2434, HYP-2435, T784, `04-computation/theta_code_lattice_gate_codex.py`.

## OPEN-Q-067 🟡 Complete or kill the order-5 branch of the extremal [72,36,16] code

**Status:** OPEN (codex-2026-06-11-P5). The order-5 fixed projection has been reduced to a tiny exact gate: for automorphism type `5-(14,2)`, the projected fixed code must be `e8+e8` with the two fixed coordinates split across the summands; the `d16+` branch is excluded because every marked pair lies in a tetrad and lifts to weight `12`. Thus the fourteen 5-cycles split into two heptads with Fano-plane tetrads, producing exactly `14` fixed minimum words and `49967` moving minimum-word orbits. The next problem is the nonfixed `F_16` component: enumerate or obstruct Hermitian self-dual `[14,7]` candidates compatible with the split-heptad fixed boundary, binary minimum distance `>=16`, and the residual `5-(72,16,78)` design ledger. Files: HYP-2439, HYP-2440, HYP-2441, T785, `04-computation/order5_fixed_projection_72_codex.py`.

## OPEN-Q-068 🟡 Prove the LRC14 Q27 resource bound beyond one stranger

**Status:** OPEN (codex-2026-06-11-P6, HYP-2444/HYP-2438). The one-stranger family `S(r)=7*{1..12} union {r}` is now closed by the fibered band-1 lattice `Q27={d*m:d|14,m<=27}`: all 936 primitive rows have a Q27 witness, and the two rows whose first plain witness is `q=41` are caught at `q=91`. The residue mechanism is explicit: the 7-core covers 8/9 classes of `(Z/27)^*/+-`, misses `+-10`, and every plain q<=27 shell blocker also has `r mod 13=0`. The open problem is to lift this from one stranger to arbitrary primitive multiple-of-14 configurations: prove that blocking Q27 consumes independent resources across shell-27 classes, low clocks such as 13, divisor fibers `d in {1,2,7,14}`, and B' safe-component gaps, so that 13 runners cannot block all Q27 and B'(any). First computational target: two-stranger rows with a resource-vector output, constrained to keep low divisor clocks covered; the naive pair of one-stranger blockers over `7*{1..11}` is too easy because all 28 such pairs have a q=12 witness.

## OPEN-Q-069 🟡 Transfer Church's diagonal Frobenius support gate to LRC14 and the [72,36,16] support problem

**Status:** OPEN (codex-2026-06-12, HYP-2445). Church's product-quotient counterexamples show the scalar/support split in geometric form: Shioda supersingularity is too coarse, while diagonal symmetric forms on every partial Frobenius twist force curve descent or finite exceptional types. Formalize the shared support-gate lemma: scalar quotient `Q`, retained channel `S`, and descent/exception rule `D`. Test two transfers: for LRC14, can Q27 blockers be forced to spend independent resources or descend to Bprime/owner-deletion exceptions; for `[72,36,16]`, can the minimum-word `5-(72,16,78)` support ledger use the `D7` index `78` and `D6/A4` index `91` as incidence-arithmetic probes? Files: HYP-2445, T789, `04-computation/product_quotient_support_gate_atlas_codex.py`, reflection `shioda-product-quotient-obstructions-and-support-gates.md`.

## OPEN-Q-070 🟡 Build the irreducible-prime certificate tournament

**Status:** OPEN (codex-2026-06-12, HYP-2448; extends HYP-2447). Formalize the finite/infinite tournament suggested by Bunyakovsky/Buniakowski plus the Singh/Cohn/Iravanian reverse certificates. Vertices should be certificate states, not just polynomials: fixed divisor, local residue obstructions, least Singh/Cohn value depth, trace-subset survivor profile, and Newton/non-Archimedean support data. Edges orient toward smaller unresolved factorization ambiguity after normalizing degree and fixed divisor, with richer retained support as tie channel. First tasks: replace the floating real-trace scout by exact algebraic trace lattices; build `C(f;X)` for a larger polynomial family and measure edge flips as `X` grows; translate the same carrier to LRC14 Q27 resource vectors and to `[72,36,16]` support/matroid/design construction moves. Files: HYP-2448, HYP-2447, T792, `04-computation/irreducible_prime_carrier_tournament_codex.py`, reflection `irreducible-prime-carriers-and-certificate-tournaments.md`.

## OPEN-Q-071 - Build the marked coefficient-row irreducibility tournament

**Status:** OPEN (codex-2026-06-12, HYP-2449; extends HYP-2447/HYP-2448). The coefficient-sign tiling is a genuine fixed-path tournament carrier, but the finite scout shows bare unmarked tournaments and sign vectors are too coarse for irreducibility. Formalize a marked coefficient-row state `R(f;P,X)` consisting of skip-row signs, coefficient magnitudes, local zero-prime residues, p-adic valuation/Newton-slope data for primes `P`, Cohn base/evaluation addresses, and Singh value-depth certificates up to `X`. First tasks: implement exact Newton-row tournaments for Eisenstein/Dumas/Perron criteria; measure edge flips as primes and evaluation depth are added; compare Cohn prime rows against low-Omega composite rows with identical sign tournaments; transfer the fixed-divisor row detector to LRC14 Q27 resource rows. Files: HYP-2449, T793, `04-computation/coefficient_tiling_prime_irreducible_codex.py`, reflection `coefficient-tiling-and-prime-irreducible-addresses.md`.
## OPEN-Q-072 🟡 Classify irreducible coefficient-magnitude slices in the tiling quotient

**Status:** OPEN (codex-2026-06-12, HYP-2450; extends HYP-2448). The coefficient-tiling quotient maps fixed-path tournaments to count profiles `c_d` and centered magnetizations `A_d=2c_d-(N-d)`. Cohn gives one rigorous lane: a positive-degree prime base-value of the diagonal-count profile certifies irreducibility of the digit polynomial. The open problem is the magnetization lane: characterize magnitude vectors `(|A_d|)` whose sign slices are forced irreducible, forced reducible, or have bounded factor patterns. The pilot at `N=6` finds the parity-minimum slice `(1,0,1,0,1)` has only 8 distinct polynomials and all are irreducible, while the full fixed-path quotient has 91/120 profiles with hidden `H` variation. Transfer target: attach the lost fiber data to LRC14 Q27 resource ledgers and to `[72,36,16]` support/matroid/design realization. Files: HYP-2450, HYP-2448, T794, `04-computation/coefficient_tiling_prime_bridge_codex.py`, reflection `coefficient-tilings-and-the-prime-irreducible-bridge.md`.

## OPEN-Q-073 - Build split-survivor ledgers for polynomial rows and LRC14 resources

**Status:** OPEN (codex-2026-06-12, HYP-2451; extends HYP-2449/HYP-2450). Reducibility is a hidden convolution lift of the coefficient row, so the live state should record which degree-split rectangles survive after each local gate. First tasks: add Newton/valuation certificates to the `38` degree-4 irreducibles with no small mod-p blocker through `31`; extend split-survivor signatures to degree `5` and `6` with cached finite-field factorizations; add Singh-depth/Cohn-depth only for rows that survive residue and valuation gates; transfer the same survivor ledger to LRC14 Q27 denominator/resource fibers, replacing scalar `q blocked` with surviving local lift obligations. Files: HYP-2451, HYP-2450, HYP-2449, T795, `04-computation/convolution_lift_irreducibility_carrier_codex.py`, reflection `convolution-lift-split-survivors-and-hidden-factor-grids.md`.
## OPEN-Q-074 🟡 Build bounded integer convolution-lift obstructions beyond degree 5

**Status:** OPEN (codex-2026-06-12, HYP-2452; extends HYP-2451/HYP-2450). Reducibility can be encoded as an integral hidden tiling problem: find nontrivial factor coefficient rows `b_i,c_j` whose multiplication grid has diagonal sums `a_k=sum_{i+j=k} b_i c_j`. The HYP-2452 pilot gives an exact integer oracle for primitive degree `<=5`, with zero mismatches against Sympy on `3856` degree-4 rows and `2016` degree-5 rows, complementing HYP-2451's residue/valuation split-survivor carrier. The open problem is to push this beyond the linear/quadratic-factor range without falling back to a black-box factorizer: encode bounded degree splits as SAT/ILP/SMT feasibility, add Newton-slope boundary constraints for sparse/multivariate rows, and use Singh-style low-`Omega(f(m))` factor-capture witnesses as pruning. Transfer target: treat LRC14 blocker ledgers and `[72,36,16]` weight-enumerator coefficients as boundary totals whose hidden support/incidence lifts must exist. Files: HYP-2452, HYP-2451, HYP-2450, T796, `04-computation/convolution_factor_capture_tiling_codex.py`, reflection `convolution-factor-capture-and-hidden-coefficient-tilings.md`.

## OPEN-Q-075 - Build moment-lift resource ledgers for LRC14 shells

**Status:** OPEN (codex-2026-06-12, HYP-2453; extends HYP-2443/HYP-2444/HYP-2452). The triangular-tower computation reframes the user's two towers as moment-balanced shell splits: `A_n` is the square-shell first-moment split and `B_n` is the triangular-shell second-moment split. The first two moments give exact integer starts `n^2` and `2n^2+n`; higher moments require a fractional address with leading term `(p-1)(p-2)/(12p)`. Addendum: A covers every positive integer, while B only covers alternating triangular shells; whole-equation side-aligned containment is the Pell family `T_n=2T_m`, and the exact whole-side equality `B_3.L=A_4.R=[21,24]` is unique. The open problem is to transfer this to LRC14 by enriching Q27 blocker rows with moment/resource data: blocked unit twists, owner supports, divisor fibers, raw moment defects, and the fractional or fiber address needed to lift a scalar blocked shell into an actual support proof. First tasks: prove the higher-moment expansion to more terms, extend the floor-sqrt/Beatty classifier into a useful Q27 address ledger, and compare AP/V*/2AP plus HYP-2444 one-stranger residuals under the new ledger. Files: HYP-2453, HYP-2444, HYP-2452, T797, `04-computation/triangular_tower_moment_bridge_codex.py`, `04-computation/triangular_tower_overlap_families_codex.py`, reflection `triangular-towers-moment-lifts-and-fractional-addresses.md`.
## OPEN-Q-076 🟡 Prove the triangular power-center bracket and test the 78/90 support shadow

**Status:** OPEN (codex-2026-06-12, HYP-2454; addendum to HYP-2453). The user's ordinary and square towers are exact interval power balances with centers `2T_n` and `4T_n`. The finite scout verifies that for `3<=p<=8` and `n<=40`, the positive root of `D_p(C,n)=0` is not integral but lies between `2pT_n` and `2pT_n+1`. Prove or refute this bracket for all `p>=3`; solve the Pell-style endpoint boundary families controlling overlaps between the first square-shell partition and the second square-balance tower; and turn the special row `Q_L(3)=[21,22,23,24]`, with ordinary shadows `90=S1(4)` and `78=C(13,2)`, into an actual support-ledger constraint for LRC14 and the `[72,36,16]` `5-(72,16,78)` minimum-design problem. Files: HYP-2454, HYP-2453, T798, `04-computation/triangular_power_balance_towers_codex.py`, reflection `triangular-power-balance-towers-and-additive-square-bridges.md`.

## OPEN-Q-077 🟡 Build a common hidden-lift feasibility engine across irreducibility, LRC, unit distance, and code72

**Status:** OPEN (codex-2026-06-13, HYP-2455; extends HYP-2452/HYP-2444/OPEN-Q-057/HYP-2454). Recent work says the scalar boundary total is not the proof object: polynomial coefficients need convolution factor grids, LRC q-blocking needs runner/Pisano/divisor/owner support, unit-distance products are reducible baselines before Moser-irreducible fibers, and `[72,36,16]` weight enumerators need support-design incidence. Build a shared lift-feasibility data model with boundary totals, candidate hidden cells, local gates, surviving allocations, and proof owners. First tasks: degree-6 bounded ILP/SAT for HYP-2452, multi-stranger LRC allocation ledgers beyond one-stranger Q27, product-reducible versus Moser-irreducible `N=27/28` unit-distance fibers, and a binary incidence-lift encoding for the `[72,36,16]` `78/90` support address. Files: HYP-2455, T799, `04-computation/boundary_lift_analogy_atlas_codex.py`, reflection `boundary-lift-irreducibility-transfer.md`.

## OPEN-Q-078 - Build a Beatty-Pell style Q27 address ledger for LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2456; concrete instance of HYP-2455; extends HYP-2241/HYP-2443/HYP-2444/HYP-2453). The triangular crossover word now has an exact hidden-address normal form: a Beatty shell address `d_m`, a Pell/carry remainder `r_m`, and state inequalities whose rare equality walls are Pell atoms. Build the analogous LRC14 ledger for Q27 blockers. For each candidate row and denominator, record `(q, shell class, unit quotient class, divisor fiber, owner support, carry residue, endpoint/boundary atom, opening or deletion target)` rather than only "q blocked." First tasks: run this on AP/Vstar/2AP and the HYP-2444 one-stranger family; measure whether visible strict/wall/open status becomes pure after adding owner/carry/private-deletion fields; compare the remaining boundary atoms to the triangular `LR/RL` zero-density wall grammar. Files: HYP-2456, HYP-2455, HYP-2453, HYP-2241, HYP-2443, HYP-2444, T800, `04-computation/triangular_tower_beatty_pell_decomposition_codex.py`, reflection `beatty-pell-crossover-word-and-lrc-address-ledgers.md`.

## OPEN-Q-079 - Prove the Faulhaber anchor expansion/bracket and port odd-wall ledgers to LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2457; sharpens HYP-2454 and complements HYP-2456). The midpoint defect for the power-balance anchor is exactly `D_p(c,n)=c^p-2*sum_{r odd} binom(p,r)c^(p-r)S_r(n)`, so only odd Faulhaber moments survive. The stored computation verifies the formal expansion `c=p*n(n+1)+alpha_p+beta_p/(n(n+1))+gamma_p/(n(n+1))^2+...`, with all displayed corrections divisible by `(p-1)(p-2)` and hence exact recovery of the p=1/p=2 towers. First tasks: prove a uniform fixed-`p` remainder after `gamma_p`; use it, or a sharper direct inequality, to prove HYP-2454's bracket `D_p(p*n(n+1),n)<0<D_p(p*n(n+1)+1,n)` for every `p>=3`; compare the p=2 square-pyramidal cuboid packing against higher simplex/cuboid carriers; and build an LRC14 analogue where odd walls, owner support, shell-27 class, divisor fiber, carry residue, and endpoint atom replace scalar "q blocked" status. Files: HYP-2457, HYP-2454, HYP-2456, T801, `04-computation/triangular_faulhaber_anchor_expansion_codex.py`, reflection `faulhaber-anchors-square-pyramids-and-bernoulli-addresses.md`.

## OPEN-Q-080 - Build an odd-moment compatibility lift analogous to OCF alpha packets

**Status:** OPEN (codex-2026-06-13, HYP-2458; extends HYP-2457/HYP-2456 and OCF). HYP-2457 isolates the odd Faulhaber anchor expansion, but OCF warns that odd atom counts are not the full object: `H(T)=I(Omega(T),2)` needs compatible packets `alpha_k` of vertex-disjoint odd cycles. Build an explicit finite compatibility lift whose one-particle shadow is the odd moment list and whose packet terms record coexistence of shell, carry, endpoint, owner-support, and support-design atoms. First targets: add odd-atom compatibility fields to the HYP-2456 Beatty/Pell side states, run them against LRC14 Q27 AP/Vstar/2AP and HYP-2444 one-stranger rows, and test whether code72 `78/90` support packets behave more like OCF `alpha_k` than like scalar moments. Files: HYP-2458, HYP-2457, T802, `04-computation/faulhaber_odd_moment_ocf_bridge_codex.py`, reflection `faulhaber-odd-moments-and-ocf-cycle-packets.md`.

## OPEN-Q-081 - Build a parity-typed Q27 ledger for LRC14

**Status:** OPEN (codex-2026-06-13, HYP-2459; extends HYP-2458/HYP-2444/HYP-2443 and the complement-Walsh line). The projector rule is exact: midpoint anti-symmetrization keeps odd Faulhaber channels, while tournament converse keeps even Walsh channels for invariant scalars. The open LRC task is to type every Q27 ledger field as `even_scalar`, `odd_marked`, `transported`, or `compatibility_packet`. First targets: AP/Vstar/2AP and HYP-2444 one-stranger rows; split source/sink or start/end fields into sum and difference before quotienting; then test whether remaining primitive rows either get a strict witness, descend to the known wall atoms, or expose an odd owner/carry/deletion opening. Files: HYP-2459, HYP-2458, T803, `04-computation/parity_projector_channel_atlas_codex.py`, reflection `parity-projectors-and-even-odd-channel-gates.md`.
