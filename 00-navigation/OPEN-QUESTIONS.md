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
- p=23: H=15,760,206,976,379,349 (computed opus-S10)

**Sequence H/|Aut|:** 1, 9, 1729, 6857869865 — no obvious pattern. 3^k pattern breaks catastrophically at p=11.

**Ratio H(P(p))/(p!/2^{p-1}):** 2.000, 2.400, 2.440, 2.527, 2.557 for p=3,7,11,19,23 — converges toward e=2.718, consistent with Szele-Alon-Friedland asymptotic theory (opus-S10/S11).

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

**Source:** kind-pasteur-2026-03-05-S2, S5, S14; opus-2026-03-05-S5 (H(T_19)), opus-2026-03-05-S10 (a(8)=661, H(P(23)), exhaustive n=7), opus-S11 (Szele analysis)

---

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

**Source:** kind-pasteur-2026-03-06-S18/S18e/S18f, sc-maximizer-mechanism.md

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

Proof via strong induction using LES of (T, T\v):
- Base case n=5 verified exhaustively (720/720)
- Induction step: good-vertex existence (THM-109)
- Case 2 (free cycles exist): Lemma A (free adj dom) + Lemma B → n-5 good vertices for n≥6
- Case 3 (all dominated): **Extreme Score Lemma** (ALGEBRAIC) — bad vertex has score 0
  or n-1, hence in no 3-cycle, hence blocks other bad vertices. At most 1 bad.
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

REMAINING OPEN:
- **What bound replaces beta_3 <= 1 at n >= 8?** (beta_3=2 confirmed at n=8,9)
- **Prove beta_1 * beta_3 = 0 for all n** (currently n <= 7 only)
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

| p | H(Paley) | H(Interval) | Margin | Winner |
|---|----------|-------------|--------|--------|
| 13 | - | 3,711,175 | - | INTERVAL (p=1 mod 4, no Paley) |
| 19 | 1,172,695,746,915 | 1,184,212,824,763 | -1.0% | INTERVAL |
| 23 | 15,760,206,976,379,349 | 16,011,537,490,557,279 | +1.59% | INTERVAL |

The interval's margin is WIDENING with p, consistent with the spectral argument: |mu_1| ~ p/pi grows faster than Paley's sqrt(p)/2.

**What remains:** Extend to p = 29, 31. An analytical proof could use the spectral concentration argument from THM-137. Whether interval maximizes H among ALL tournaments (not just circulant) is a separate open question.

**Source:** opus-2026-03-12-S58, kind-pasteur-2026-03-12-S56c, kind-pasteur-2026-03-12-S57

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery.
- **OPEN-Q-002**: Claim A PROVED for all n. OCF proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20). See CONJ-001, THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **OPEN-Q-009**: Arc-flip invariance resolved — E(T) = 0 for all T (OCF proved). See THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy.
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11.
