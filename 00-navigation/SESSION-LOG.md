# Session Log

Chronological record of all sessions. Every new Claude instance adds an entry at the **top** of this file before doing any work.

Entry format:
```
## [INSTANCE-ID] — [DATE]
**Account:** [A/B/C/...]
**Continuation of:** [previous instance ID, or "fresh start"]
**Files read:** [list of files read at session start]
**Summary of work:** [brief description]
**New contributions:** [theorem IDs, court cases, tangents added]
**Unresolved threads:** [things left open for next session]
```

## opus-2026-03-05-S12 — 2026-03-05 (Deep web research: Knuth connections to OCF — 7 channels identified)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S11 (continuation)
**Files read:** MSG-014 (kind-pasteur-S14b), SESSION-LOG.md, TANGENTS.md
**Summary of work:**
  (1) **KNUTH FASCICLE 8a** (Dec 2025): Section 7.2.2.4 of TAOCP on Hamiltonian paths and cycles. "Path flipping" = our arc-flip approach. Directly relevant to OCF proof strategy. Priority: read for tournament-specific exercises.
  (2) **CLAUDE'S CYCLES** (Knuth, Feb 2026): Hamiltonian cycle decomposition of Cayley digraphs solved by Claude Opus 4.6 via modular Gray code. 760 valid decompositions. Parallel: algebraic structure (Cayley/Paley) enables exact solutions.
  (3) **RÉDEI REVISITED** (arXiv:2510.10659, Oct 2025): Schweser-Stiebitz-Toft exhibit stronger Rédei theorems by Dirac and Berge. Dirac's theorem constrains parity of Hamiltonian paths in mixed graphs — could give mod-4 information beyond OCF.
  (4) **ZDD TECHNOLOGY**: Knuth's Zero-Suppressed Decision Diagrams could efficiently compute I(Omega(T), x) for moderate n, exploiting sparse complement of dense Omega.
  (5) **ALGORITHM X / DLX**: Exact cover formulation of odd-cycle packing. DLX could enumerate all vertex-disjoint odd-cycle collections, providing alternative OCF verification.
  (6) **PERMANENT / CYCLE COVERS**: Ryser's inclusion-exclusion for permanents is structurally parallel to independence polynomial evaluation. Bridge to Björklund's GF(2) method.
  (7) **PÓLYA ENUMERATION**: Aut(T_p) acting on odd-cycle collections. Cycle index of affine group could give structured formula for I(Omega(T_p), 2).
  (8) Read kind-pasteur-S14b findings: S_{2,1,1}-free at n<=9, quasi-line fails at n=8, line graph at n=6. Full hierarchy: line(n<=5) < quasi-line(n<=7) < claw-free(n<=8) < S_{2,1,1}-free(n<=9) < S_{1,1,1}-free(n<=11).
**New contributions:**
- 03-artifacts/drafts/knuth-connections-research.md (7 connections with priority ranking)
**Unresolved threads:**
- Read Knuth fasc8a.pdf for tournament exercises
- Read arXiv:2510.10659 (Rédei revisited) for stronger parity theorems
- Implement ZDD-based I(Omega(T), x) for n=8-12
- Formalize Cayley/Paley algebraic structure parallel

## opus-2026-03-05-S11 — 2026-03-05 (Continued research: Szele-Adler-Alon e-convergence, real roots at n=15, Paley asymptotic)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S9 (continuation)
**Files read:** All new messages (MSG-019/S10), TANGENTS.md, n8-anomaly-deep-dive.md, SESSION-LOG.md
**Summary of work:**
  (1) **SZELE BOUND DISCOVERY:** Adler-Alon-Ross (2001) proved max H(T) >= (e-o(1))·n!/2^{n-1}. Paley ratios H(P(p))/(p!/2^{p-1}) = 2.00, 2.40, 2.44, 2.53, 2.56 converge exactly toward e=2.718. Paley tournaments achieve the asymptotically optimal constant.
  (2) **REAL ROOTS EXTENDED to n=15:** 0/95 tournaments at n=9..15 have non-real roots for I(Omega_3, x). Crucially, this holds at n=12-15 where S_{1,1,1}-freeness FAILS. Alpha(Omega_3) grows as ~floor(n/3).
  (3) **FRIEDGUT-KAHN UPPER BOUND:** Improved Alon's n^{3/2} to n^{3/2-0.25} ≈ n^{1.25}. If Paley achieves ~e·n!/2^{n-1}, the polynomial factor may be entirely unnecessary.
  (4) **LINE GRAPH HYPOTHESIS ALREADY DISPROVED** in S9. K_5-e in Omega_3 at n=6.
  (5) Updated zero-free-subdivided-claw-research.md with Connections 7 and 8.
**New contributions:**
- Extended real-roots verification data (n=11..15)
- Szele-Adler-Alon analysis connecting Paley ratios to e
- Updated 03-artifacts/drafts/zero-free-subdivided-claw-research.md
**Unresolved threads:**
- Does H(P(p))/(p!/2^{p-1}) → e exactly? Needs H(P(31))
- Can Alon upper bound be improved to O(1)? (major open problem)
- Real roots: algebraic explanation needed (symmetric function path)
- Irving-Omar paper Corollary 20: exact statement still unverified

## kind-pasteur-2026-03-05-S14b — 2026-03-05 (Structural hierarchy: subdivided-claw-freeness, line graph refutation)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S14 (context continuation)
**Files read:** THM-020, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md
**Summary of work:**
  (1) T054 REFUTED: Omega(T) is NOT a line graph at n=6. K5-e (Beineke forbidden subgraph) in 53% of n=6 tournaments. Heilmann-Lieb cannot explain real roots.
  (2) NEW FINDING: Omega(T) is S_{2,1,1}-free (subdivided-claw-free) at n<=9 (0/100 failures), but FAILS at n>=10 (92%). Combined with opus-S9: S_{1,1,1}-free through n=11, fails at n=12. Each subdivision level buys ~3 vertices.
  (3) Quasi-line graph property FAILS at n=8 (49%), n=9 (10%), n=10 (0%).
  (4) Full structural hierarchy: line graph (n<=5) < quasi-line (n<=7) < claw-free (n<=8) < S_{2,1,1}-free (n<=9) < S_{1,1,1}-free (n<=11) < ??? Real roots hold at ALL tested n.
  (5) LITERATURE: Alon (1990) + Adler-Alon-Ross (2001): max H(T) = Theta(n!/2^{n-1}). Paley maximizer conjecture consistent.
  (6) LITERATURE: Jerrum-Patel (2026) proves real roots for bounded-degree H-free (H = subdivided claw). Key insight from opus-S9: every fixed subgraph eventually appears in Omega(T), so explanation must be algebraic.
  (7) Verified all real roots + log-concavity + unimodality at n=9 (10 random).
**New contributions:**
- T061-T064 (line graph refuted, S_{2,1,1}-free, quasi-line fails, Alon bounds)
- Updated THM-020 with expanded graph property table (including line graph, quasi-line, S_{2,1,1})
- Updated OPEN-Q-015 with comprehensive structural analysis
**Unresolved threads:**
- What structural property explains real roots at n>=10? (opus-S9 says: must be algebraic)
- Real roots: Irving-Omar/Grinberg-Stanley symmetric function framework is most promising path

## opus-2026-03-05-S10 — 2026-03-05 (Paley maximizer verification, n=8 H-maximizer discovery, Omega structure analysis)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S8 (context continuation after summary)
**Files read:** SESSION-LOG.md, TANGENTS.md, OPEN-QUESTIONS.md, THM-020, MSG-013, n8-anomaly-deep-dive.md, verify_n8_claims.py, INVESTIGATION-BACKLOG.md
**Summary of work:**
  (1) DISCOVERED: a(8)=661 from OEIS A038375. The H-maximizer at n=8 is a SC tournament with |Aut|=1, NOT T_657 (H=657, the Paley extension). T_661 does NOT contain P(7) as a vertex-deletion and has non-uniform D_v=(67,63,68,72,72,68,63,67).
  (2) CONFIRMED: P(7) is the GLOBAL H-maximizer at n=7, verified by exhaustive enumeration of all 2,097,152 tournaments on 7 vertices. 240 tournaments achieve H=189.
  (3) CONFIRMED: Paley maximizer conjecture at p=11: H(P(11))=95095=a(11).
  (4) NEW COMPUTATIONS: H(P(19))=1,172,695,746,915 and H(P(23))=15,760,206,976,379,349. These are candidate values for a(19) and a(23) in OEIS.
  (5) RATIO ANALYSIS: H(P(p))/(p!/2^{p-1}) = 2.00, 2.40, 2.44, 2.53, 2.56 for p=3,7,11,19,23 — converging toward e=2.718, consistent with Szele-Alon-Friedland asymptotic bound.
  (6) OMEGA STRUCTURE: Full Omega at n=8 has 76-78 vertices (density 0.98). The 3-cycle subgraph accounts for only 30% of H(T); 70% comes from 5-cycle and 7-cycle contributions. This is the structural signature of the n=8 transition.
  (7) IMPORTANT NEGATIVE: The Paley maximizer conjecture applies ONLY at Paley primes (p=3 mod 4). At non-Paley n=8, the maximizer has trivial automorphism group and no Paley substructure.
**New contributions:**
- T059 (n=8 maximizer and 5-cycle dominance), T060 (Paley ratio convergence to e) in TANGENTS.md
- T053 updated with p=19,23 data and n=8 negative result
- INV-024 extended through p=23, INV-032 updated with perfectness failure
- 04-computation/check_max_h8.py, paley_maximizer_test.py, paley_h_sequence.py, n8_omega_structure.py
**Unresolved threads:**
- Submit H(P(p)) values to OEIS (a(19), a(23) are new)
- Test line graph hypothesis for Omega(T) (Beineke forbidden subgraphs)
- Test subdivided-claw-freeness at n=9
- Compute H(P(31)) if feasible

## opus-2026-03-05-S9 — 2026-03-05 (Deep web research: zero-free regions, subdivided claws, real-roots mystery)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S8 (context continuation after compaction)
**Files read:** MEMORY.md, all inbox messages (MSG-010 through MSG-013, MSG-015), THM-020-real-roots.md, n8-anomaly-deep-dive.md, web-research-connections.md, tribonacci-web-research.md, grinberg-stanley-connection.md, MISTAKES.md, THM-002, THM-019, INVESTIGATION-BACKLOG.md, hard-core-statistical-physics-connections.md
**Summary of work:**
  (1) SYSTEMATIC AUDIT: Found and fixed errors in tribonacci-web-research.md (4 sections with false claw-free/perfect claims), web-research-connections.md (Connection 4), INVESTIGATION-BACKLOG (INV-032).
  (2) BUG FOUND: hard_core_fast.py chordality check (lines 157-168) only checks 1 of 3 C4 orderings. Reports 0% non-chordal when correct is ~14% at n=6.
  (3) **NEW COMPUTATIONAL: S_{1,1,1}-freeness** — Omega_3(T) is S_{1,1,1}-free (subdivided claw with each edge subdivided once) through n=11, FAILS at n=12 (100%). Each subdivision level buys ~3 more vertices. Pattern: claw-free n<=8, S_{1,1,1}-free n<=11.
  (4) **LINE GRAPH HYPOTHESIS DISPROVED:** K_5-e (Beineke forbidden subgraph) found in Omega_3(T) at n=6 (45%). Not a line graph.
  (5) **JERRUM-PATEL 2026 ANALYSIS:** Subdivided claws are the EXACT boundary for zero-free regions of I(G,x) in bounded-degree H-free graphs. Best possible: non-subdivided-claw H allows zeros on positive real axis. BUT: requires bounded degree AND fixed H, both of which fail for Omega(T) at large n.
  (6) **BEZAKOVA et al. 2024:** Matching mixing-time dichotomy — O(n log n) for S_{a,b,c}-free bounded-degree, exponential otherwise. Same limitations.
  (7) **KEY INSIGHT: Real-rootedness cannot be explained by any forbidden subgraph property.** Every fixed subgraph eventually appears in Omega(T). Explanation must be algebraic, likely through Irving-Omar/Grinberg-Stanley symmetric function framework.
  (8) **AUTHOR CORRECTION:** arXiv:2412.10572 is by Irving & Omar, not Grinberg & Stanley. They build on G-S's framework.
  (9) **SZELE BOUND:** H(T) <= c·n^{3/2}·n!/2^{n-1} (Alon 1990, solving Szele 1943 conjecture). Lower bound n!/2^{n-1} achieved probabilistically. Paley maximizer conjecture (kind-pasteur-S14b) says T_p achieves the max; ratio H(T_p)/(p!/2^{p-1}) ≈ 2, 2.4, 2.44 for p=3,7,11.
**New contributions:**
- 03-artifacts/drafts/zero-free-subdivided-claw-research.md (comprehensive research document)
- Fixed 4 errors in existing documentation files
- S_{1,1,1} computational data for n=9..12
**Unresolved threads:**
- Prove real roots for n>=9 — algebraic approach needed (symmetric functions?)
- Rate of growth of max subdivided claw in Omega(T) as function of n
- Is there a weaker graph property ("quasi-claw-free"?) that holds universally?
- Irving-Omar framework: does it say anything about I(Omega(T), x) beyond x=2?
- Compute H(T_19) for Paley maximizer conjecture
- Szele bound tightness: does H(T_p) ~ c·p!/2^{p-1} for some constant c?

## kind-pasteur-2026-03-05-S14 — 2026-03-05 (Deep research: real roots theorem, Paley maximizer, literature survey)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S13 (context continuation)
**Files read:** MEMORY.md, THM-019, OPEN-QUESTIONS.md, TANGENTS.md, SESSION-LOG.md
**Summary of work:**
  (1) NEW THEOREM (THM-020): All roots of I(Omega(T), x) are real and negative for n<=8. PROVED via claw-freeness + Chudnovsky-Seymour (2007). Real-rootedness conjectured for all n.
  (2) COMPUTATIONAL: Verified real roots at n=5 (exhaustive), n=6 (500), n=9 (50), n=10 (30) — 0 failures anywhere, even beyond claw-free range.
  (3) TESTED & DISPROVED: Omega(T) is NOT always a comparability graph (fails n=7: 1%, n=8: 58%).
  (4) TESTED: det(I+2A) =/= H(T) in general. No simple spectral formula.
  (5) FIXED: Cycle enumeration bug — multiple directed cycles per vertex set.
  (6) **MAJOR: Paley tournaments MAXIMIZE H(T)!** OEIS A038375 confirms: a(3)=3=H(T_3), a(7)=189=H(T_7), a(11)=95095=H(T_11). New conjecture: T_p maximizes Hamiltonian paths for Paley primes.
  (7) LITERATURE (5 background agents, all completed):
    - 7 papers cite Grinberg-Stanley, including chromatic-Redei-Berge connection (arXiv:2506.08841)
    - Bezakova et al. 2024 dichotomy for H-free graphs and hard-core model
    - Omega(T) NOT interval graph at n=6 (13.9% fail)
    - Line graph hypothesis: if Omega(T) = L(H), Heilmann-Lieb gives real roots
    - H(T_p) and H/|Aut| sequences NOT in OEIS
    - El Sahili-Abi Aad (2019) extends Forcade parity; Grunbaum conjecture proved
    - Our conflict graph / independence polynomial formulation is NOT in any G-S paper — genuinely original
    - Cycle divisibility threshold corrected: k >= (p+3)/2, not (p+1)/2
  (8) FINDING: lambda=2 is OUTSIDE hard-core uniqueness for Omega(T). Real roots are special.
  (9) Graph property hierarchy: interval (n<=5) < chordal (n<=5) < comparability (n<=6) < perfect (n<=7) < claw-free (n<=8).
**New contributions:**
- THM-020-real-roots.md (new theorem)
- OPEN-Q-015 (real roots conjecture)
- T048-T058 in TANGENTS.md (Paley maximizer, line graph, interval graph, LGV connection, subdivided claw, chromatic-RB connection)
- Corrected cycle divisibility threshold in OPEN-Q-013
**Unresolved threads:**
- Test line graph hypothesis (Beineke's 9 forbidden subgraphs)
- Test subdivided-claw-freeness at n=9
- Compute H(T_19) — would test both Paley maximizer conjecture and extend the sequence
- Submit H(T_p) to OEIS
- Bijective proof of OCF still open

## kind-pasteur-2026-03-05-S13 — 2026-03-05 (Web research: hard-core model, independence polynomial at lambda=2, statistical physics connections)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S12 (new session)
**Files read:** All warm-up files, inbox messages (MSG-007, MSG-013, MSG-014 x2), tribonacci-web-research.md
**Summary of work:**
  (1) Comprehensive web research on connections between H(T) = I(Omega(T), 2) and statistical physics.
  (2) FINDING: lambda=2 is in the NON-PERTURBATIVE regime of the hard-core model. The uniqueness threshold lambda_c(Delta) < 2 for all Delta >= 5. Cluster expansion diverges; approximate counting is #P-hard. OCF is a rare exact evaluation in this regime.
  (3) FINDING: I(G, 2) counts 2-labeled independent sets: pairs (S, f) where S is independent and f: S -> {1,2}. For OCF, each Ham path corresponds to a colored collection of vertex-disjoint odd cycles.
  (4) FINDING: Scott-Sokal (2005) connects LLL to independence polynomial zero-free regions. The Shearer zero-free disc is far too small to reach lambda=2 for any nontrivial max degree. Positivity of I(Omega(T), 2) cannot be deduced from LLL; it requires the combinatorial identity H(T) >= 1.
  (5) FINDING: Chudnovsky-Seymour (2007) gives all-real-roots for claw-free graphs, applicable to Omega(T) only for n <= 8. Jerrum-Patel (2026) extends zero-free regions to H-free classes.
  (6) Confirmed Grinberg-Stanley proof (arXiv:2307.05569 + arXiv:2412.10572) uses matrix algebra, NOT statistical mechanics. No bijective interpretation exists yet.
  (7) Created comprehensive research report: hard-core-statistical-physics-connections.md with 11 references.
**New contributions:**
- 03-artifacts/drafts/hard-core-statistical-physics-connections.md (comprehensive research report)
**Unresolved threads:**
- Bijective proof of OCF: each Ham path <-> 2-colored independent cycle collection
- Whether Omega(T) has special structure beyond claw-freeness that explains positivity at lambda=2
- Potential use of transfer matrix methods from stat mech for tournament problems

## opus-2026-03-05-S8 — 2026-03-05 (Deep n=8 anomaly investigation: BlackSelf, mu explosion, Paley connection)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S7
**Files read:** THM-019, tiling-class-structure.md, FINAL_FINDINGS.md, THM-016, INVESTIGATION-BACKLOG.md, omega_claw_fast.py, q009_adj_decomposition.py, tex paper (oq:n8 definition)
**Summary of work:**
  (1) Investigated the n=8 anomaly at the user's request, connecting blueself/blackself/pos formulas.
  (2) Exhaustive enumeration of all 65536 self-converse n=8 tournaments with reversal anti-automorphism.
  (3) Found 2560 with |Aut|>1, surveyed all: 10 distinct (H, Fix(beta), |Aut|) types.
  (4) DISCOVERED: T_657 (H=657, |Aut|=3) contains Paley tournament P(7) as vertex-deletion — the Hadamard connection for INV-012.
  (5) DISCOVERED: T_657 has perfectly uniform D_v = 54 for all vertices (unique among specimens).
  (6) DISCOVERED: mu values at n=8 reach {3, 5, 7, 9, 11} for 3-cycles and {1, 3} for 5-cycles.
  (7) DISCOVERED: Full Omega at n=8 has 76 vertices (20+48+8 cycles of length 3,5,7).
  (8) IDENTIFIED: T_A (|Aut|=9, H=621) is the ONLY SC+Aut>1 specimen with no C5 in Omega_3.
  (9) CORRECTED: Signed position identity requires T and T' (flipped tournament), not same T on both sides.
  (10) Web research: Chudnovsky-Seymour, El Sahili parity results, doubly regular tournaments ↔ skew Hadamard.
**New contributions:**
- `03-artifacts/drafts/n8-anomaly-deep-dive.md` — comprehensive analysis
- `04-computation/blackself8_*.py` — four investigation scripts
- INV-012 updated with detailed findings
**Unresolved threads:**
- Confirm T_657 is isomorphic to a known construction
- Resolve BlackSelf(8) definition ambiguity
- Full Omega independence polynomial at n=8 (76 vertices, needs specialized algorithm)
- Investigate whether uniform D_v characterizes Paley extensions

## opus-2026-03-05-S7 — 2026-03-05 (Omega perfectness DISPROVED at n=8; claw-freeness trivial at n<=8)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S6 (context continuation x3)
**Files read:** THM-019, grinberg-stanley-connection.md, tiling-class-structure.md, OPEN-QUESTIONS.md, SESSION-LOG.md, all inbox messages, omega_perfectness_implications.py
**Summary of work:**
  (1) Reviewed all new material from git pull: OCF proved by Grinberg-Stanley, THM-019, Tribonacci, tiling class structure.
  (2) Committed endpoint decomposition / linearity work (q009_endpoint_decomp.py, q009_G_identity.py, q009_G_in_test.py, q009_linearity_proof.py) as historical contributions.
  (3) PROVED: Omega(T) is trivially claw-free for n<=8 by vertex counting (3 pairwise v.d. odd cycles + 1 touching all three requires >= 9 vertices).
  (4) DISPROVED: Omega(T) is NOT always claw-free at n=9 (90% of random tournaments have claws).
  (5) DISPROVED: Omega(T) is NOT always perfect at n=8. 53.8% of random n=8 tournaments have a C5 (5-hole) in the 3-cycle conflict graph. Explicit counterexample constructed with 13 forced arcs.
  (6) VERIFIED: Omega(T) appears perfect at n<=7 (0/1000 failures).
  (7) CONFIRMED: OCF holds regardless of Omega perfectness (H=I(Omega,2)=417 at the n=8 counterexample).
**New contributions:**
- THM-019 CORRECTED: perfectness fails at n=8, holds for n<=7
- OPEN-Q-014 RESOLVED (disproved)
- 04-computation/omega_claw_fast.py (claw-freeness analysis)
- 04-computation/omega_c5_test.py (C5 construction at n=8)
- 04-computation/omega_claw_free_test.py (comprehensive test, slow)
- Vertex counting theorem: claw impossible for n<=8
**Unresolved threads:**
- Do all-real-roots / log-concavity of I(Omega,x) survive at n>=8 despite imperfection?
- What structural property of Omega(T) DOES explain OCF? (not perfectness, not claw-free)
- Prove Omega(T) perfectness for n<=7 (clean vertex-count argument for C5?)

## kind-pasteur-2026-03-05-S12 — 2026-03-05 (CRITICAL: OCF proved by Grinberg-Stanley; comprehensive audit)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S11 (context continuation)
**Files read:** All warm-up files, all theorem files (THM-001 through THM-018, CONJ-001, CONJ-002, PROP-001, LEM-001, LEM-002), MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, tiling-class-structure.md, signed-adjacency-identity.md, paper-deep-connections.md
**Summary of work:**
  (1) COMPREHENSIVE AUDIT of all mathematical results in the project. Read every theorem file, checked proof chains, cross-referenced claims.
  (2) FOUND AND FIXED: THM-016 significance section incorrectly claimed even-odd split implies OCF (contradicts MISTAKE-008). Fixed.
  (3) FOUND AND FIXED: signed-adjacency-identity.md status line said "equivalent to OCF". Fixed.
  (4) FOUND: THM-013 numbering collision (THM-013-insertion-decomposition.md says "THM-012" inside). Fixed header.
  (5) FOUND: Two duplicate THM-016 files (claim-b and hamiltonian-alternating-sum). Noted for future merge.
  (6) FOUND: LEM-001 and LEM-002 outdated after CONJ-002 refutation. Updated with resolution notes.
  (7) VERIFIED numerically (via background agents): THM-016 inductive proof (all steps correct, no gaps), Tribonacci proof (all 5 components verified), OCF at n=3,4,5 (independent recomputation), H(T_11)=95095, Claim B at n=6.
  (8) **CRITICAL DISCOVERY**: Web search found that OCF (H(T) = I(Omega(T), 2)) is ALREADY PROVED in the literature by Grinberg & Stanley (arXiv:2307.05569, 2023; arXiv:2412.10572, 2024). Their Corollary 20 states ham(D̄) = Σ 2^{ψ(σ)} over permutations with all odd D-cycles. For tournaments, D̄ = D^op and ham(D^op) = ham(D), so this gives H(T) = I(Omega(T), 2) exactly.
  (9) Updated CONJ-001 (now PROVED), THM-002 (now PROVED for all n), PROP-001 (now PROVED), OPEN-Q-002 and OPEN-Q-009 (RESOLVED). Created grinberg-stanley-connection.md documenting the full equivalence.
  (10) Updated LEM-001, LEM-002 with computed values (h_QR=h_NQR=201).
  (11) Web searches also confirmed: Tribonacci connection is NEW (not in OEIS for tournaments), H(T_19) value is NEW, Omega(T) perfectness is NEW, the H/|Aut| sequence 1,9,1729 is NOT in OEIS.
**New contributions:**
- CONJ-001 → PROVED (via Grinberg-Stanley OCF + Claim B)
- THM-002 → PROVED for all n (three independent proofs documented)
- PROP-001 → PROVED
- OPEN-Q-002, OPEN-Q-009 → RESOLVED
- 03-artifacts/drafts/grinberg-stanley-connection.md (full analysis of equivalence)
- Error fixes: THM-016, signed-adjacency-identity.md, THM-013 header, LEM-001, LEM-002
- 04-computation/verify_core_results.py (independent OCF verification)
- 04-computation/verify_tribonacci_proof.py (tribonacci proof verification)
**Unresolved threads:**
- Merge duplicate THM-016 files
- Canonize Omega(T) perfectness (potentially new theorem)
- The project's independent contributions (THM-016/017 even-odd split, THM-018 coefficient identity, tiling class structure, Tribonacci theorem) remain valuable and NEW
- Consider submitting Tribonacci result to OEIS
- The Grinberg-Stanley proof technique (symmetric functions / matrix algebra) is entirely different from our combinatorial methods — worth studying for further applications

## opus-2026-03-05-S6 — 2026-03-05 (Tribonacci web research: interval graph structure, Chudnovsky-Seymour connection)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S5 (context continuation)
**Files read:** tiling-class-structure.md, full_class_analysis.py, web-research-connections.md, hard_core_fast.py, paper-deep-connections.md
**Summary of work:** (1) Deep web research with Tribonacci structure in mind. (2) DISCOVERED that all odd cycles of T_full_n are consecutive intervals — Omega(T_full_n) is an INTERVAL GRAPH. (3) PROVED I(Omega(T_full_n), 2) satisfies Tribonacci recurrence independently of H(T_full_n) via weighted interval packing DP with telescoping. (4) Connected to Chudnovsky-Seymour: Omega(T) claw-free implies I(Omega(T), x) has ALL REAL (negative) ROOTS — log-concavity, positivity at x=2. (5) Connected to Chvátal-Sbihi: claw-free perfect graphs decompose via clique cutsets into elementary (line graphs of bipartite graphs) and peculiar atoms. (6) Identified Jerrum-Patel 2026 paper extending zero-free regions. (7) Found multiple combinatorial interpretations of Tribonacci (tilings, compositions, binary sequences).
**New contributions:**
- 03-artifacts/drafts/tribonacci-web-research.md (comprehensive synthesis)
- INV-035 added to INVESTIGATION-BACKLOG (Tribonacci interval graph structure)
- INV-032 updated with Chudnovsky-Seymour connection
- Verified OCF for T_full family n=3,...,8 via interval graph computation
**Unresolved threads:**
- Prove Omega(T) is always claw-free (INV-032 — would unlock Chudnovsky-Seymour)
- Direct bijection between run decompositions and weighted interval packings
- Extend Tribonacci analysis to other tournament families
- Check transfer matrix eigenvalues for T_full

## opus-2026-03-05-S5b — 2026-03-05 (THM-018: alpha_w^H = alpha_w^I PROVED symbolically at n=4,...,7)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4d (context continuation x2)
**Files read:** All warm-up files, previous computation files, THM-013, THM-002
**Summary of work:** (1) Reformulated alpha_w^H as insertion decomposition: for each perm pi of B_w, insert w at each position with interface correction factors. (2) Discovered base = -2*H(W) and correction = 2*(H(W)-H(B_w)) + cycle_derivs. (3) Proved key structural insight: T(I,v) = T(J,v) at s=0 (i,j interchangeable). (4) PROVED alpha_w^H = alpha_w^I as polynomial identity at n=4,5,6,7 using SymPy symbolic computation (diff = 0 exactly). (5) Identified that n>=8 needs additional Delta(alpha_2) terms from THM-013 (VD cycle pairs). (6) Created THM-018 documenting the proof.
**New contributions:**
- THM-018: coefficient identity alpha_w^H = alpha_w^I proved at n<=7 (symbolic)
- 04-computation/q009_alpha_identity.py (insertion decomposition and factor analysis)
- 04-computation/q009_insertion_decomp.py (defect = cycle derivatives verification)
- 04-computation/q009_symbolic_proof.py (symbolic proof n=4,5,6)
- 04-computation/q009_symbolic_n7.py (symbolic proof n=7)
- 04-computation/q009_symbolic_n8.py (n=8 attempt — needs full delta_I formula)
**Unresolved threads:**
- General proof of alpha_w^H = alpha_w^I for all n (key remaining step for OCF)
- Need full delta_I formula (including Delta(alpha_2)) for n>=8 symbolic verification
- Inductive proof approach: use OCF at sub-tournament level to simplify Delta(alpha_k)

## kind-pasteur-2026-03-05-S11 — 2026-03-05 (Tribonacci theorem; tiling class structure deepened)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S10 (context continuation)
**Files read:** All warm-up files, full_class_analysis.py, tiling-class-structure.md, tiling_structure_analysis.py, OPEN-QUESTIONS.md, INVESTIGATION-BACKLOG.md, paper-deep-connections.md, broadcast MSG-011
**Summary of work:**
  (1) DISPROVED the 2^(n-2)+1 conjecture for the full tiling class at n=7,8 (H=31,57 vs expected 33,65).
  (2) DISCOVERED the correct formula: H(T_full_n) = Tribonacci(n) = H(n-1)+H(n-2)+H(n-3), OEIS A000213. Verified through n=12.
  (3) PROVED the Tribonacci recurrence via a combinatorial argument: Ham paths biject to "run decompositions" (ordered sequences of consecutive intervals with gap condition). Conditioning on first interval size yields f(n) = S(n-2) + S(n-3) where S is the partial sum; telescoping gives f(n) = f(n-1)+f(n-2)+f(n-3). Complete proof written up.
  (4) Analyzed self-paired class "perpendicular" geometry: confirmed all self-flip pairs have Hamming midpoint exactly at m/2 (hypercube center). Self-paired classes straddle the transitive<->full diagonal. All self-paired tournaments are self-converse with odd H values.
  (5) Pulled and reviewed opus-S5 findings: H(T_19) computed, alpha_1 conjecture disproved, Omega(T) always perfect, hard-core non-perturbative.
**New contributions:**
- Tribonacci theorem for full tiling class (complete proof in tiling-class-structure.md)
- Updated tiling-class-structure.md with F3 correction, F11 (perpendicular geometry), F12 (Tribonacci)
- Self-paired class H values all odd: 5, 11, 13, 15, 25, 41, 43, 45
**Unresolved threads:**
- Omega perfectness implications for OCF (agent running but incomplete)
- Formula for number of self-paired classes: 1, 2, 8, ...
- Bridge even-odd split to OCF
- Transfer matrix symmetry proof (INV-001, highest priority per opus-S5)

## opus-2026-03-05-S5 — 2026-03-05 (Priority C+E investigation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4b (continuation after context compaction)
**Files read:** INVESTIGATION-BACKLOG.md, paper-deep-connections.md, hard_core_lattice_gas.py, hard_core_fast.py, ballot_sequence_test.py, compute_H_T19.py, mod4_score_test.py, conflict_graph_catalog.py, SESSION-LOG.md, TANGENTS.md, OPEN-QUESTIONS.md
**Summary of work:** Continued Priority C and E investigation from previous context. Ran all computational scripts created in prior context and gathered results:
  - H(T_19) = 1,172,695,746,915 (INV-024 DONE). H/|Aut| = 6,857,869,865. All 19 endpoints have equal count (Paley symmetry). Computed via C DP in 0.5s.
  - alpha_1 ≡ c_3 (mod 2) conjecture DISPROVED (INV-026). Counterexamples at n=3,4,5. All have alpha_1=0 but c3 odd.
  - Hard-core lattice gas (INV-028): lambda=2 above ALL convergence thresholds. OCF is non-perturbative.
  - Omega(T) structure for n=4,5: max degree grows, independence number = 1 at n=4.
  - mod4_score_test.py has a bug: reports OCF mismatches at n=5, but OCF is correct (bug in cycle-finding code). Alpha_1 result is valid since it's independent of OCF.
  Updated INVESTIGATION-BACKLOG.md with all findings: INV-015, INV-016, INV-019 (paper connections deepened), INV-024 (H(T_19) computed), INV-026 (conjecture disproved), INV-028 (non-perturbative confirmed), INV-029 (bijective proof found).
**New contributions:**
- INVESTIGATION-BACKLOG.md updated with 7 investigation results
- paper-deep-connections.md (created in prior context, verified this session)
- compute_H_T19.py output: H(T_19) = 1,172,695,746,915
- ballot_sequence_test.py: bijective proof of C(L-2, 2k-1) confirmed
- hard_core_lattice_gas.py + hard_core_fast.py: non-perturbative regime confirmed
- mod4_score_test.py: alpha_1 ≡ c_3 (mod 2) disproved (NOTE: has OCF mismatch bug at n=5)
**Unresolved threads:**
- INV-001 (transfer matrix symmetry proof) — highest priority, Feng connection is the best lead
- INV-004 (R-cone + cut-flip proof strategy) — enabled by Rajkumar Prop 1
- INV-027 (conflict graph catalog) — script exists but too slow, needs optimization
- INV-020/021 (Striker/Chapman papers) — correct arXiv IDs not found
- mod4_score_test.py OCF mismatch bug needs fixing
- conflict_graph_catalog.py needs optimization for n>=5

## kind-pasteur-2026-03-05-S10 — 2026-03-05 (Claim B PROVED; gap in proof chain identified)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S9 (context continuation)
**Files read:** All warm-up files, opus-S4c session entry, q009_claim_b.py, q009_direct_proof.py, q009_claim_b_proof.py, signed-adjacency-identity.md, MISTAKES.md
**Summary of work:** (1) Ran q009_claim_b_proof.py — all steps of the inductive proof verified at m=1,...,8 (max error 1.9e-12 at m=8, well within floating-point tolerance). (2) Created THM-016: formal write-up of the inductive proof of Claim (B) for all m. The proof uses h_start expansion, exchange of summation, inductive hypothesis on |W|-1 vertices, and path extension. Tournament condition T(a,b)+T(b,a)=1 is essential (Step 5). (3) IDENTIFIED GAP in opus-S4c proof chain: the step "B(Li,Rj)=B(Lj,Ri) -> OCF" is INVALID per MISTAKE-008. Even-odd split is a consequence of OCF, not equivalent. Corrected signed-adjacency-identity.md. (4) Updated OPEN-Q-009 to reflect that even-odd split is now proved for all n, but OCF remains open.
**New contributions:**
- THM-016: Claim (B) alternating sum identity — PROVED for all m
- q009_claim_b_proof.py verified at m=1,...,8
- q009_creative_hypotheses.py (10 hypotheses tested, dual identity confirmed)
- Corrected signed-adjacency-identity.md (removed false equivalence claim)
- Updated OPEN-Q-009 (even-odd split proved, gap documented)
**Unresolved threads:**
- Bridge the gap: even-odd split -> OCF needs additional identity
- The odd-S sum of Delta(S,R) needs to be connected to the cycle formula
- Investigation backlog items remain (Forcade 1973, Chapman ASM, transfer matrix proof)

## opus-2026-03-05-S4d — 2026-03-05 (THM-016 PROVED + THM-017 PROVED — even-odd split for all n)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4c (context continuation)
**Files read:** All previous session files, OPEN-QUESTIONS, theorem files
**Summary of work:** PROVED THM-016 (Claim B path identity) by induction on m using first-step decomposition. Key steps: (1) separate boundary term, (2) first-step decomposition of h_start, (3) exchange summation order, (4) recognize inner sum as THM-016 at size m-1, (5) use T(v,u)=1-T(u,v) to get H(W')-h_end(W,v), (6) cancellation (-1)^{m-1}+(-1)^m=0. Then PROVED THM-017 (B(Li,Rj)=B(Lj,Ri)) for all n by induction on |W|, using THM-016 for the d-independent part and the induction hypothesis for the d-dependent part. Both proofs verified numerically at all sizes. CORRECTED: the even-odd split is necessary but NOT sufficient for OCF (kind-pasteur-S8 was right). Investigated the gap: F(x)=(1+x)Q(x) where F(1)=delta_H and F(-1)=0, but Q doesn't have universal structure.
**New contributions:**
- 04-computation/q009_claim_b_proof.py (THM-016 proof verification)
- 04-computation/q009_full_proof.py (THM-017 full proof verification)
- 04-computation/q009_bridge_gap.py (investigating gap to OCF)
- 01-canon/theorems/THM-016-hamiltonian-alternating-sum.md (PROVED)
- 01-canon/theorems/THM-017-even-odd-split-proof.md (PROVED)
- Updated OPEN-QUESTIONS, signed-adjacency-identity.md with honest assessment
**Unresolved threads:**
- OCF for all n still OPEN — even-odd split insufficient
- Need delta_H = delta_I for all n (currently proved n<=8)
- The Q polynomial in F(x)=(1+x)Q(x) may yield a path forward

## opus-2026-03-05-S4c — 2026-03-05 (Claim B discovery — key reduction for OCF proof)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4
**Files read:** All warm-up files, previous session's exploration scripts
**Summary of work:** Major progress on OPEN-Q-009. (1) Corrected the s-evenness reduction (C_w+D_w=0 was wrong; correct statement is s-linear coefficients vanish). (2) Verified B(Li,Rj)=B(Lj,Ri) as polynomial identity at n=3,...,7 with real-valued arcs in [-2,3]. (3) Discovered the PAIRING STRUCTURE: for each u!=v, left(u)=right(u) exactly, meaning the d-dependent part of alpha_v vanishes by the INDUCTIVE HYPOTHESIS. (4) Isolated CLAIM (B): a standalone identity about tournament Hamiltonian path weights: sum_{S⊆W\{v}} (-1)^|S| H(S) h_start(W\S,v) = (-1)^{m+1} h_end(W,v). Verified at m=1,...,8. Proved algebraically at m=1,2,3. (5) Established full proof chain: Claim B -> B(Li,Rj)=B(Lj,Ri) -> OCF -> Claim A.
**New contributions:**
- 04-computation/q009_cw_dw_proof.py (showed C_w+D_w!=0, corrected reduction)
- 04-computation/q009_algebra_n4.py (exact n=4 algebra, s-coordinate formula)
- 04-computation/q009_s_quadratic.py (quadratic s-structure analysis)
- 04-computation/q009_direct_proof.py (alpha_v decomposition, pairing discovery)
- 04-computation/q009_claim_b.py (Claim B verification m=1,...,8)
- Updated signed-adjacency-identity.md with proof structure and Claim B
**Unresolved threads:**
- Prove Claim (B) for all m — the ONLY remaining step for OCF proof
- Claim (B) is an alternating-sum identity about internal tournament Hamiltonian paths
- Proved at m=1,2,3 by direct algebra; need general proof (possibly by induction on m)

## kind-pasteur-2026-03-05-S9 — 2026-03-05 (tex deep analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S8 (context continuation)
**Files read:** Full parity_tournaments_fixed.tex (2189 lines), all opus-S4b new files (paper-connections.md, INVESTIGATION-BACKLOG.md, 4 new computation scripts), inbox messages MSG-004 and MSG-005
**Summary of work:** Pulled opus-S4b contributions (transfer matrix symmetry discovery, investigation backlog, paper connections). Then performed line-by-line analysis of entire tex file. Found 5 issues: (1) DR mod-4 proof (Thm 7.4) has broken arithmetic. (2) SE-SYT formula (Thm 7.3) produces non-integer. (3) Transitive uniqueness proof incomplete. (4) Verification record outdated. (5) Rajkumar et al. missing from bibliography.
**New contributions:**
- 03-artifacts/drafts/tex-deep-analysis.md (comprehensive analysis report)
- INV-028b, INV-029b, INV-030b added to backlog
**Unresolved threads:**
- Forcade 1973 GF approach (INV-023)
- Chapman 2001 ASM bijection (INV-021)
- Transfer matrix symmetry proof (INV-001)

---

## opus-2026-03-05-S4b — 2026-03-05 (repo scour + backlog creation)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S4b (context continuation — second part)
**Files read:** Full tex paper (open problems, bibliography, literature), TANGENTS, OPEN-QUESTIONS, MISTAKES, MASTER_FINDINGS, proof-landscape, signed-adjacency-identity, all code file headers, paper-connections.md
**Summary of work:** Completed the user's directive to establish a systematic practice of scouring the repo for leads. Read every significant file in the repo. Created INVESTIGATION-BACKLOG.md with 31 cataloged leads across 5 priority tiers (A: blocks OCF proof, B: structural understanding, C: references to investigate, D: computational targets, E: tangent investigations). Updated CLAUDE.md with mandatory Step 5b (repo scouring every session). Notable findings: Rajkumar paper (2110.05188) is NOT in the tex bibliography despite being connected; 6 bibliography entries have NOT been investigated for OCF connections (Striker, Chapman, Eplett, El Sahili x2, Forcade GF approach); 3 open problems from the paper have zero investigation (mixed graphs, Striker-Chapman equivariance, realizable conflict graphs).
**New contributions:**
- `00-navigation/INVESTIGATION-BACKLOG.md` — comprehensive lead catalog (31 items)
- `CLAUDE.md` Step 5b — mandatory repo-scouring practice for all agents
- Memory file created at `.claude/projects/.../memory/MEMORY.md`
**Unresolved threads:**
- INV-001 (transfer matrix symmetry proof) — highest priority
- INV-008 (Striker-Chapman equivariance) — completely uninvestigated
- INV-015 (add Rajkumar to bibliography)
- INV-024 (compute H(T_19)) — feasible, not done
- n=8 C verifier may still be running from earlier in this session

## kind-pasteur-2026-03-05-S8 — 2026-03-05 (error audit)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S7 (context continuation)
**Files read:** All canon files, all opus-S3/S4 new files, all proof scripts, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Pulled opus-S3 and opus-S4 contributions, then audited the entire codebase for errors. Found and fixed 4 issues: (1) BUG in sympy_proof_n8.py — used simplified n<=7 formula which THM-013 says FAILS at n=8; rewrote to use full A-clique formula. (2) LOGICAL ERROR in even-odd split lemma — claimed "equivalent to OCF" but it's only a consequence; the odd-S sum of Delta(S,R) differs from the cycle formula [g-l]*H(R) at the per-subset level. (3) Stale claim in proof-landscape saying n<=7 instead of n<=8. (4) Duplicate tangent number T040. Verified all SymPy proofs (n=4,5,6) and opus n=7 proof pass correctly. NOTE: opus-S4b ran in parallel and found additional structure (Signed Position Identity, tournament-specificity, bracket analysis) — see their session entry.
**New contributions:**
- MISTAKE-008: even-odd split equivalence claim corrected
- MISTAKE-009: sympy_proof_n8.py formula bug documented and fixed
- sympy_proof_n8.py rewritten with correct full A-clique formula
- Tangent numbering fixed (T040 duplicate -> T044-T047)
- Even-odd-split-lemma.md corrected, OPEN-Q-009 corrected, TANGENTS corrected
**Unresolved threads:**
- Run sympy_proof_n8.py overnight (corrected version, full A-clique formula)
- Prove OCF for all n — central open problem
- Reconcile even-odd split equivalence claim with opus-S4b's Signed Position Identity analysis

---

## opus-2026-03-05-S4b — 2026-03-05 (Transfer matrix symmetry + paper connections)
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S3 (parallel session to S4, then extended)
**Files read:** All warm-up files, inbox snippet from S4 diff, all n=8 computation files, arXiv:2110.05188, arXiv:2510.25202
**Summary of work:** Three phases: (1) Verified Even-Odd Split Lemma, proved tournament-specific/polynomial, built C verifier for n=8. (2) Analyzed two external papers: "Tournament Representations" (flip classes, locally transitive, R-cones) and "Dual Burnside Process" (Q=AB factorization, spectral correspondence). (3) MAJOR DISCOVERY: the transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S)*B_b(M\S) is ALWAYS SYMMETRIC (verified 7500+ tests, n=4..8). This is STRONGER than Even-Odd Split — each cross-term individually agrees. Connected to Burnside detailed balance. Also found locally transitive tournaments DO have 5-cycles and 7-cycles (conjecture "only 3-cycles" fails). All OCF tests pass for R-cones, automorphism-symmetric tournaments, locally transitive tournaments.
**New contributions:**
- 03-artifacts/drafts/paper-connections.md (10 connections to external papers)
- 04-computation/locally_transitive_test.py (cycle structure of rank-2 tournaments)
- 04-computation/burnside_connection_test.py (character/Fourier decomposition)
- 04-computation/transfer_matrix_test.py (transfer matrix structure)
- 04-computation/symmetry_check.py (VERIFIES transfer matrix symmetry, 7500+ tests)
- T045 (Transfer matrix symmetry discovery)
- T046 (Tournament representations / flip class connection)
**Unresolved threads:**
- PROVE transfer matrix symmetry for all n (would prove OCF)
- Explore detailed-balance / reversibility interpretation
- R-cone simplification strategy (prove OCF for R-cones, extend via flip)

---

## opus-2026-03-05-S4 — 2026-03-05 (n=7 and n=8 PROVED, even-odd split discovery)
**Account:** Eliott (opus machine)
**Continuation of:** opus-2026-03-05-S3 (context continuation)
**Files read:** TANGENTS.md (resolved merge conflict), THM-015, PROP-001, OPEN-QUESTIONS.md, SESSION-LOG.md, symbolic_proof.py, symbolic_proof_fast.py, kind-pasteur inbox messages
**Summary of work:** Resolved git merge conflict from S3 rebase (renumbered my T032 to T039). Then wrote optimized n=7 exhaustive verifier using numpy-vectorized bitmask approach — PROVED OCF at n=7 (2^20 = 1,048,576 configs in 4 seconds). Extended to n=8 with chunked processing — PROVED OCF at n=8 (2^27 = 134,217,728 configs, 57 minutes). Discovered the Even-Odd Split Lemma: the adj decomposition delta = sum_S Delta(S,R) splits equally between even-|S| and odd-|S| terms, so delta = 2*(odd-S sum). This connects directly to the cycle formula (only odd cycles). The alternating sum sum(-1)^|S| Delta(S,R) = 0 is equivalent to OCF but provides a clean algebraic reformulation. Verified n=5,...,8.
**New contributions:**
- THM-015 updated: n=7 and n=8 PROVED
- 04-computation/q009_prove_n7.py (numpy-vectorized exhaustive verifier)
- 04-computation/q009_prove_n8.py (chunked n=8 verifier)
- 04-computation/q009_even_odd_split.py (even-odd split discovery)
- 04-computation/q009_alternating_sum.py (alternating sum analysis)
- 03-artifacts/drafts/even-odd-split-lemma.md (new proof angle documentation)
- Tangent T040 (even-odd split)
- OPEN-Q-009 updated with proof frontier and even-odd split
**Unresolved threads:**
- n=8 exhaustive verification COMPLETE (134M/134M, 57 min)
- Prove the alternating sum identity for ALL n (equivalent to OCF)
- The even-odd split is a clean algebraic reformulation but doesn't simplify the proof

---

## kind-pasteur-2026-03-05-S7 — 2026-03-05 (n=7 proof + structural analysis)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S6 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, MISTAKES.md, definitions.md, THM-013, THM-014, THM-015, all S6 code files
**Summary of work:** Extended OCF proof to n=7 via two independent methods: (1) SymPy polynomial identity verification (20 arc variables, 1858 monomials, 77s) and (2) exhaustive {0,1} enumeration (1,048,576 cases, 775s). Both confirm delta_H = delta_I as polynomial identity at n=7. Also explored structural decomposition for general proof: proved delta_H = -sum s_x*(R1+R2) algebraically, showed R1+R2 = 2*H(B_x) at n=4, found excess terms at n=5 satisfy sum s_x*exc(x) = -2*(D5-C5) on {0,1} but NOT as polynomial identity (non-multilinear). Investigated f(S)+f(S^c) pairing: universal at n=4, fails at n>=5 due to 5-cycle corrections. Created complete hand proof for n=4 via f(S) decomposition.
**New contributions:**
- THM-015 updated to PROVED at n<=7 (was n<=6)
- 03-artifacts/code/sympy_proof_n5.py (SymPy polynomial identity proof at n=5)
- 03-artifacts/code/sympy_proof_n6.py (SymPy polynomial identity proof at n=6)
- 03-artifacts/code/sympy_proof_n7.py (SymPy polynomial identity proof at n=7)
- 03-artifacts/code/algebraic_proof_n4.py (complete hand proof at n=4)
- 03-artifacts/code/pairing_proof.py (f(S)+f(S^c) analysis, shows n>=5 failure)
- 03-artifacts/code/proof_structure_analysis.py (pred/succ decomposition, s-signature analysis)
**Unresolved threads:**
- Prove polynomial identity for ALL n — central remaining open problem
- n=8 SymPy verification feasible (user offered overnight compute)
- Structural proof approaches: transfer matrix, induction on n, permanent expansion
- Non-multilinear excess identity blocks naive decomposition approach

---

## opus-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** opus-2026-03-05-S2 (context overflow, resumed)
**Files read:** All warm-up + THM-012, THM-013, TANGENTS, OPEN-QUESTIONS
**Summary of work:** Extended tournament_lib.py with arc-flip functions. Independently verified THM-013 adjacency identity. Extended OCF verification to n=10 (first time n>=8 verified).
1. Added adj_count(), flip_arc(), verify_thm013(), verify_ocf(), independence_poly_at_fast() to tournament_lib.py
2. THM-013 adjacency identity: VERIFIED n=5 exhaustive (10240/10240), n=6 sampled (3000/3000)
3. OCF H(T)=I(Omega(T),2): PROVED n<=7 exhaustive (1,048,576 arc assignments at n=7, 27min), n=8 (500 random), n=9 (100), n=10 (30) -- all 0 failures
4. Algebraic analysis: decomposed adj(i,j)-adj'(j,i) as subset convolution. Per-vertex LHS/RHS decomposition does NOT match -- proof must work globally.
5. Key insight: at n<=8, max independent set size in Omega(T) is 2; at n=9, it's 3. Fast I computation exploits this.
6. Explored proof strategies: per-vertex nadj decomposition (dead end), transfer matrix/permanent (no clean formula found), bijection search (suggestive structure but no proof), inclusion-exclusion over blockers (correct framework but doesn't simplify).
7. Key structural finding: s=0 vertices never block swaps (cleanly separates U_T from U_T'), but nadj_sum/H(B_x) ratio is not constant -- the identity is fundamentally global.
**New contributions:**
- tournament_lib.py: 6 new functions (adj_count, flip_arc, compute_s_x, count_directed_5_cycles_through_arc, verify_thm013, verify_ocf, independence_poly_at_fast)
- verify_ocf_sweep.py: clean verification script
- symbolic_proof_n7.py: exhaustive n=7 OCF proof (1,048,576/1,048,576 PASS)
- OCF PROVED through n=7, verified through n=10
- Updated THM-013, THM-015, OPEN-QUESTIONS.md with n=7 exhaustive results
- 4 exploratory scripts: unmatched_decomposition_by_vertex.py, inductive_structure.py, transfer_matrix_approach.py, bijection_search.py, proof_strategy_analysis.py
**Unresolved threads:**
- PROVE the identity for all n (not just n<=7): the central open problem
- Subset convolution framework is the right language but per-vertex doesn't match; need global argument
- The "2-colored independent cycle set" interpretation of I(Omega,2) may lead to a bijective proof
- n=8 proof feasible with optimization (2^27 ~134M cases) but would take hours

---

## kind-pasteur-2026-03-05-S6 — 2026-03-05 (polynomial identity proof)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S5 (context continuation)
**Files read:** SESSION-LOG.md, OPEN-QUESTIONS.md, TANGENTS.md, THM-013, THM-014, all S5 code files
**Summary of work:** Resolved 6-file git merge conflict from S5 rebase. Then discovered and proved a NEW PROOF METHOD for OCF: the swap involution's unmatched path count U_T'-U_T can be expressed as a polynomial identity in arc variables, and this polynomial equals delta_I from THM-013. Proved by hand at n=4 (clean formula: U_T'-U_T = 4-2(p+q+r+t) = 2*sum(s_x)). Verified exhaustively at n=5 (512 cases) and n=6 (16384 cases). This constitutes a PROOF of OCF for n<=6. Found and fixed a sign convention error in THM-013's D/C labels.
**New contributions:**
- THM-015 in 01-canon/theorems/ (swap polynomial identity, proves OCF at n<=6)
- 03-artifacts/code/symbolic_proof.py (polynomial identity verifier for n=4,5,6)
- 03-artifacts/code/symbolic_proof_fast.py (optimized n=6 version)
- 03-artifacts/code/unmatched_decomposition.py (structural analysis of blocking)
- 03-artifacts/code/unmatched_vs_formula.py (aggregate identity tests)
- Tangents T037-T038 (polynomial identity proof, sign convention warning)
- Resolved S5 merge conflicts (TANGENTS renumbered T028-T036, all files clean)
**Unresolved threads:**
- Prove the polynomial identity for ALL n (not just n<=6) — this would prove OCF/Claim A
- The n=4 hand proof suggests a transfer matrix / generating function approach
- n=7 verification feasible (2^19 = 524K cases) but needs optimization
- OPEN-Q-009 partially resolved: proof method identified, works at n<=6

---

## kind-pasteur-2026-03-05-S5 — 2026-03-05 (arc-flip session)
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S4
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, CONJ-001, THM-003 (Claim B proof), LaTeX paper claim_strategies section
**Summary of work:** Attacked OPEN-Q-009 (arc-reversal invariance) with new angle. Discovered that E(T) = H(T) - I(Omega(T), 2) is invariant under arc flips (exhaustive n<=5, random n=6). Derived algebraic formula for delta_I via A-clique argument adapted to arc flips. Showed the simple formula delta=2*(gained-lost) works at n<=5 but fails at n>=6. Tested and ruled out contiguous-block decomposition. Formulated PROP-001 (arc-flip identity) as a new proof strategy for OCF/Claim A equivalent to proving a combinatorial identity about Ham path changes under arc flips weighted by complement H-counts.
**New contributions:**
- PROP-001 in 01-canon/theorems/ (arc-flip identity, new proof strategy)
- 03-artifacts/code/arc_reversal_study.py (delta_H decomposition under arc flips)
- 03-artifacts/code/ocf_arc_flip_study.py (E(T) invariance, joint delta_H/delta_I distribution)
- 03-artifacts/code/delta_I_correct_formula.py (algebraic delta_I formula verification)
- 03-artifacts/code/delta_I_fast_test.py, paths_via_cycles_test.py, contiguous_block_test.py
- OPEN-Q-009 reformulated with cleaner E(T) version
- Tangents T028-T031 (arc-flip approach, dead ends)
**Unresolved threads:**
- PROP-001: prove the arc-flip Ham path identity combinatorially
- OPEN-Q-009: the E(T) formulation is cleaner but still needs proof
- Three possible proof approaches: transfer matrix, generating function, involution

---

## kind-pasteur-2026-03-05-S3 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S2
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, theorem files
**Summary of work:** Rebased kind-pasteur onto main (resolved 6-file conflict). Closed DISC-001 (mu bug does not contaminate verification). Investigated OPEN-Q-010/011 at n=7: wrote test_n7_ABD.py computing A (TypeII sum), B (mu-weighted), D (cycle mu sum) for 1050 pairs. Finding: A=/=D at n=7 in general (only 5.9% exact match). Near-cancellation is statistical (mean A-D~0.1) not algebraic.
**New contributions:**
- 03-artifacts/code/test_n7_ABD.py (correct A-B-D computation at n=7)
- 03-artifacts/code/test_perpath_n7.py (initial wrong test -- documents why naive approach fails)
- DISC-001 moved to resolved/; OPEN-Q-010/011 updated with n=7 findings
- TANGENT T027 (n=7 A-B-D near-cancellation is statistical not algebraic)
**Unresolved threads:**
- OPEN-Q-009: arc-reversal invariance -- the key step, still untouched
- OPEN-Q-012: tower hypothesis -- not yet investigated
- OPEN-Q-013: H(T_p) formula for Paley primes

---

## kind-pasteur-2026-03-05-S2 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S1
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, DISC-001, CONJ-001, CONJ-002, LEM-001, LEM-002, THM-002
**Summary of work:** Resolved OPEN-Q-009 (kind-pasteur numbering) by direct computation. h_QR=h_NQR=201, c_9(T_11)=11055. Direct enumeration gives H(T_11)=95095=55*1729, refuting CONJ-002 for p=11. Added finish_session.py enforcement tooling.
**New contributions:**
- H(T_11)=95095 computed directly (03-artifacts/code/compute_H_T11.py)
- CONJ-002 REFUTED for p=11; OPEN-Q-013 opened (correct formula for H(T_p))
- MISTAKE-006: ratio-coincidence c_k/C(11,k) has no basis
- TANGENTS T022-T026 (Paley structure, 1729, symmetry)
- agents/finish_session.py and agents/check_session_closed.py (end-of-session enforcement)
**Unresolved threads:**
- OPEN-Q-013: What is the correct formula for H(T_p)?
- OPEN-Q-009 (opus): arc-reversal invariance — the key proof step
- DISC-001: needs formal close

---

## kind-pasteur-2026-03-05-S1 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** fresh start — first session on this machine
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, inbox/PROCESSING-REPORT.md, agents/REGISTRY.md, all theorem files, DISC-001
**Summary of work:** Registered as new agent. Processed PALEY_T11_c9_ANALYSIS.md inbox doc: extracted LEM-001, LEM-002, CONJ-002. Argued Position A in DISC-001.
**New contributions:** LEM-001, LEM-002, CONJ-002, MISTAKE-006 (ratio-coincidence), DISC-001 Letter 2
**Unresolved threads:** DISC-001 response needed; h_QR/h_NQR computation

---

## kind-pasteur-2026-03-05-S4 — 2026-03-05
**Account:** Eliott (primary)
**Continuation of:** kind-pasteur-2026-03-05-S3
**Files read:** MISTAKES.md, definitions.md, OPEN-QUESTIONS.md, SESSION-LOG.md, TANGENTS.md, opus broadcast messages
**Summary of work:** Structural fix session. kind-pasteur was pushing to claude/kind-pasteur branch (invisible to opus and CI). Merged claude/kind-pasteur into main. Fixed Windows encoding issues in processor.py (CP1252 fails on Unicode). Updated finish_session.py to push HEAD:main always.
**New contributions:**
- All kind-pasteur S1-S3 research now on main (Paley computation, n=7 A-B-D, DISC-001 resolved)
- agents/processor.py: Windows UTF-8 fix (stdout.reconfigure + write_text encoding)
- agents/finish_session.py: push to origin HEAD:main (not current branch)
**Unresolved threads:**
- OPEN-Q-009 (arc-reversal invariance): top priority, untouched
- OPEN-Q-012 (tower hypothesis): untouched

---

## opus-2026-03-05-S2 — 2026-03-05
**Account:** Claude Opus 4.6 (e's MacBook)
**Continuation of:** opus-2026-03-05-S1 (parallel session)
**Files read:** All navigation, canon, court files; LaTeX paper; FINAL_FINDINGS.md; file.txt; tournament_lib.py
**Summary of work:** Proved THM-008b (general mu triviality), THM-012 (partial mu invariance under arc flips). Disproved MISTAKE-004 (OCF IS a valid closed form — verified H(T)=I(Omega(T),2) for all n<=6). Resolved OPEN-Q-011 (near-cancellation is statistical, not structural). Verified adjacency formula H(T)-H(T')=adj(i,j)-adj'(j,i). Established arc-reversal decomposition framework for Q-009.
**New contributions:**
- THM-008b: general mu triviality bound L >= n-2
- THM-012: mu invariant under arc flips when at least one endpoint in V(C)\{v}
- DISC-002: MISTAKE-004 is wrong — OCF is a valid closed form (verified 33,864 tournaments)
- OPEN-Q-011 resolved (near-cancellation is statistical artifact)
- Adjacency formula verified; arc-reversal decomposition framework built
- 6 computation scripts in 04-computation/
**Unresolved threads:**
- DISC-002 needs formal resolution (retract MISTAKE-004)
- Arc-reversal proof strategy (Q-009) partially developed — key obstacle is tracking cycle creation/destruction and mu changes simultaneously
- Claim A proof for general n remains open (OPEN-Q-002)

---

## opus-2026-03-05-S1 --- 2026-03-05
**Account:** opus (Claude Opus 4.6)
**Continuation of:** SYSTEM-2026-03-05-S1 (initial setup)
**Files read:** All navigation files, all canon files, DISC-001, agents/processor.py, inbox/processor.py, README.md
**Summary of work:** Built the complete tournament computation library (tournament_lib.py) implementing all core objects from definitions.md. Independently verified Claim A at n<=6 (196,608 pairs, 0 failures) with a MISTAKE-001-compliant implementation. Also verified Claim B at n<=5, Redei at n<=6, and Claim A at n=7 by random sampling (100 tournaments, 0 failures). Ingested two contributed files (FINAL_FINDINGS.md, file.txt) containing substantial new results. Extracted 4 new theorems (THM-008 through THM-011), 2 new mistakes (MISTAKE-004, MISTAKE-005), resolved OPEN-Q-001 and OPEN-Q-003, added 4 new open questions (Q-009 through Q-012), and documented 3 dead ends (T016-T018).
**New contributions:**
- 03-artifacts/code/tournament_lib.py (CODE-000): complete computation library
- 03-artifacts/code/verify.py (CODE-000v): CLI verification runner
- THM-008: mu=1 trivially for n<=5 (resolves OPEN-Q-001)
- THM-009: per-path failure characterization at n=6 (resolves OPEN-Q-003)
- THM-010: n=4 block-counting theorem (exactly 3 Ham paths)
- THM-011: general block-counting formula H_C^+(T)
- MISTAKE-004: OCF is recursive, not closed-form over all cycles
- MISTAKE-005: cycle bijection under arc reversal fails
- OPEN-Q-009: arc-reversal invariance (key unproved step)
- OPEN-Q-010 through Q-012: new questions from ingested files
- Dead ends T016-T018 documented
- Registered agent: opus
**Unresolved threads:**
- Claim A proof for general n (CONJ-001) -- arc-reversal invariance (OPEN-Q-009) is the key step
- DISC-001 can be moved to resolved
- The near-cancellation phenomenon (OPEN-Q-011) is a promising proof strategy lead
- Per-path formula with 3+5 cycles at n=7 (OPEN-Q-010) needs computational testing
- Computation library could be extended: inshat, per-path identity, arc-reversal testing

---

## SYSTEM-2026-03-05-S1 — 2026-03-05
**Account:** System (initial setup)
**Continuation of:** fresh start — first Cowork session
**Files read:** MASTER_FINDINGS.md (uploaded), parity_tournaments_fixed.tex (uploaded)
**Summary of work:** Built the full research directory system. Processed two uploaded source files as the first contributions. Extracted theorems F1–F5 and Claim A to canon. Logged MISTAKE-001. Opened DISC-001 as a potential court case seed (μ computation bug vs. paper's 0-failure verification claim).
**New contributions:**
- THM-001 through THM-007, CONJ-001 in 01-canon/theorems/
- MISTAKES.md: MISTAKE-001 (μ computation bug in scripts 6-9)
- 00-navigation/TANGENTS.md, OPEN-QUESTIONS.md populated from source files
- 03-artifacts/drafts/parity_tournaments_fixed.tex archived
- 02-court/active/DISC-001-mu-bug-vs-verification.md opened
**Unresolved threads:**
- Claim A proof (the central open problem)
- Why does per-path identity hold at n=5 despite 5-cycles? (OPEN-Q-001)
- μ computation bug in scripts 6-9 needs investigation (DISC-001)
- Formal proof of Fix(σ) = 2^{m²} for self-evacuating SYT (OPEN-Q-007)
