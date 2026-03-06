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

## OPEN-Q-004 🟡
**Find a correct per-path formula for all n**

The 3-cycle-only formula (per-path identity) fails for n≥6. The natural generalization summing over all odd cycles overcounts. The maximal-embedding-only formula also fails. Is there a formula of the form Σ_{cycles C, (non-v consecutive in P')} f(C, P') = (inshat−1)/2 that works for all n?

---

## OPEN-Q-005 🟡
**Combinatorial proof of the C(L-2, 2k-1) distribution (THM-007)**

The number of internal signature patterns giving exactly k Type-II positions within an L-cycle embedding window is C(L-2, 2k-1). This is proved (verified computationally). A bijective proof connecting this to ballot sequences or Dyck paths is missing. See TANGENT T001.

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

## OPEN-Q-010 🟡
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

## OPEN-Q-015 🟡
**Prove I(Omega(T), x) has all real negative roots for all n**

PROVED for n<=8 via claw-freeness + Chudnovsky-Seymour (THM-020). Computationally verified:
- n=5: exhaustive (1024), 0 failures
- n=6: 500 random, 0 failures
- n=9: 50 random (3-cycle subgraph), 0 failures
- n=10: 30 random (3-cycle subgraph), 0 failures
- **FULL Omega (all odd cycles, not just 3-cycles):** n=5..9, 185 samples, 0 failures (opus-S13)
- Paley p=7: full Omega has 80 cycles {3:14, 5:42, 7:24}, I.P.=[1,80,7], roots=-11.42,-0.013 (real)
- 100% log-concave and unimodal across ALL 206 tested polynomials

At n>=9, Omega(T) has claws (86-90% of random n=9 tournaments), so Chudnovsky-Seymour doesn't apply.

**Structural hierarchy (kind-pasteur-S14b):** Each property explains real roots up to the n where it fails:
- Line graph: FAILS at n=6 (K5-e appears, 53%)
- Quasi-line: FAILS at n=8 (49%)
- Claw-free: FAILS at n=9 (86%)
- S_{2,1,1}-free (subdivided claw): FAILS at n=10 (92%)
- **Real roots: STILL HOLD at n=10,12,15,20 (0 failures)**

**Refuted approaches:** Line graph (Heilmann-Lieb), quasi-line, all standard forbidden-subgraph classes.

**BARRIER (Engstrom, arXiv:1610.00805):** Real-rootedness of the multivariate independence polynomial CHARACTERIZES claw-freeness. No broader graph class guarantees real roots via the same mechanism. The explanation MUST be tournament-specific and algebraic.

**Remaining approaches:**
- Jerrum-Patel (2026, JLMS): zero-free regions for bounded-degree H-free (H = subdivided claw). Applies at n<=9, chain grows with n.
- Irving-Omar (arXiv:2412.10572): matrix-algebraic tools for Redei-Berge symmetric function. Schur-basis resolution may connect to parity.
- Prakash-Sharma (arXiv:2510.09197): root gap quantification — could show gap stays positive for tournament conflict graphs.
- Direct structural argument from Omega(T) combinatorics (alpha <= floor(n/3)).
- Grinberg-Stanley / Irving-Omar symmetric function framework (most promising).

If true, gives: alternative proof of Redei, log-concavity of I.P. coefficients, H(T) > 0 for all x > 0.

**Source:** kind-pasteur-2026-03-05-S13 (THM-020)

---

## Resolved Questions (moved here when answered)

- **OPEN-Q-001**: Per-path identity at n=5 is trivially true (THM-008). No mystery.
- **OPEN-Q-002**: Claim A PROVED for all n. OCF proved by Grinberg-Stanley (arXiv:2412.10572, Corollary 20). See CONJ-001, THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-003**: Per-path failure at n=6 iff some TypeII position has mu>1 (THM-009).
- **OPEN-Q-009**: Arc-flip invariance resolved — E(T) = 0 for all T (OCF proved). See THM-002. (kind-pasteur-2026-03-05-S12)
- **OPEN-Q-011**: Near-cancellation is statistical, not structural. Not a viable proof strategy.
- **Paley computation (kind-pasteur)**: h_QR=h_NQR=201, c_9(T_11)=11055, H(T_11)=95095. CONJ-002 refuted for p=11.
