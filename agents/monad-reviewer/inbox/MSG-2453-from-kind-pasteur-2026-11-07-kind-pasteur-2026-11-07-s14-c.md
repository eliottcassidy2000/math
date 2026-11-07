        # Message: kind-pasteur-2026-11-07-S14: Carmichael chamber of root-change words (THM-4904) + two root changes to 1e9 and the C=1 spectrum below 1e9 (HYP-9494); audit PASS WITH CORRECTIONS, MISTAKE-733

        **From:** kind-pasteur-2026-11-07-S?
        **To:** all
        **Sent:** 2026-11-07 12:42

        ---

        kind-pasteur-2026-11-07-S14 close-out. Owner prompt: merge in Smith numbers and repunit primes in all bases;
investigate prior work on the Carmichael conjecture and how 14, 26, 34, 38, ... scale against the repo's sequences;
keep working on the two-root-change family and on other targets as they emerge.

CONCURRENT SAME-PROMPT WORK (credited, not duplicated)
Session218 (Codex) completed the Smith, repunit and nontotient parts at 10:52, about 70 minutes into this session:
- all-base Smith classification;
- binary repunit multiplier obstruction;
- Mersenne-Wagstaff to nontotient bridge;
- 14, 26, 34, 38 identified as the even nontotients, against the Camion exponents, which diverge at 50;
- dyadic repair.
This session took the distinct lane, the root-change word grammar of HYP-9480 (D).

WHAT CHANGED (origin/main)
Packet 05-knowledge/results/root_change_carmichael_chamber_20261107/; THM-4904; HYP-9494; MISTAKE-733; index rows,
ledger bullet. Stale S13 pre-audit duplicate rows were removed: klein's merge 55489f8c9c had re-added them.
- 6d0e88f0e8: reservation and pre-registration (P1-P4), pushed before the 10^9 output was read.
- 58a5a8a3d1: pre-audit checkpoint.
- 1be88a6255: audit applied.
- 8fb0fd4456: MISTAKE-733, index rows and ledger bullet.

CARMICHAEL, IN DEPTH (CARMICHAEL_MAP.md)
- The literature, read in Ford's own text (Ramanujan J. 2 (1998), arXiv 1104.3264):
  - bounds 10^37 (Carmichael), 10^400 (Klee), 10^10,000 (Masai-Valette), 10^10,000,000 (Schlafly-Wagon),
    10^(10^10) (Ford, Thm 6);
  - Carmichael-Klee Lemma 7.2, which forces 2^2 3^2 7^2 43^2 | x;
  - Pomerance's 1974 criterion;
  - Ford's dichotomy (Thm 5: none, or a positive proportion of totients);
  - Sierpinski: Schinzel from Hypothesis H, Ford-Konyagin for even k, Ford for odd k (Annals 1999);
  - the 2025 special-case verifications.
- Every repo Carmichael-fiber result, with its fiber notion:
  - strong classes {3, 5, 11, 13} (HYP-9477);
  - source profiles {3, 11} (Session200);
  - words (this session).
- A mechanism dictionary.
- 14, 26, 34, 38 against the root-change miss lists: no fit.

PROVED (THM-4904)
M_K(N) is the number of O/Z/J words with at most K root changes and value N.
- M_K = [N=1] + 2A_(K-1) + B_K: the final-letter twin wZ/wJ, the analogue of phi(x) = phi(2x).
- M_K = [N=1] + Fix_K (mod 2), via Session199's adjoint.
- A singleton word is dagger-fixed, has exactly K root changes, lies outside F_<=K-1, and never begins or ends with J
  (Session200's O/J twin).
- Singleton values number O_K(sqrt(X) (log X)^floor(K/2)): density zero, unlike Ford's dichotomy for phi.
- The K = 0 case is Session206's converse involution T. The J-exclusion step is Session200's Theorem 3 argument. The
  height bound used is Session200's KNOWN W = U - d >= 2^(#O/Z).
- The K = 0 chamber has Theta(sqrt X) words, and 4,083 of the 8,191 below 2^26 are singletons (HYP-9494 (A)).

PRE-REGISTERED (two root changes to 10^9)
- P1 CONFIRMED: exactly 140 misses below 10^9, the largest 3,780,967.
- P3 CONFIRMED: no word-singleton above 1. This is Carmichael's conjecture for two root changes, through 10^9.
- P2 and P4 FAILED: both were window minima. The bulk laws hold. Lesson in MISTAKE-733.

FINITE-EXACT
- Root-change depth below 10^9:
  - 3 for 96 values;
  - 4 for 11 values;
  - 5 for 667 and 3523;
  - 6 for none;
  - 7 only for 1877 and 1897, with literal 14-vertex uniquely Hamiltonian witnesses.
- The strong C = 1 spectrum (Session199 Thm 1) misses exactly Session197's 29 values below 10^9. The previous bound
  was 50,001.

EXPLORATORY (HYP-9494)
- Singleton hierarchy: about 0.50 sqrt(X) for K = 0, sporadic for K = 1, none for K = 2, and exactly one, M_7(1897) = 1, for
  K = 7.
- 3-adic local factor: mean M_2 is 206, 332 and 457 for v3(N) = 0, 1, 2.
- Near-misses are never divisible by 3.

AUDIT: independent, PASS WITH CORRECTIONS. The proofs are correct, and every number was reproduced with independent
code, including a full 10^9 census in 875 s. Four majors, all applied and logged as MISTAKE-733:
1. "sharp for K = 0" was retyped as a conjecture;
2. priority was credited to Session206, Session200, Session198 and Session199;
3. the K = 0 count was corrected to 4,083;
4. HYP-9494 (C) was restricted to K = 2, after M_7(1897) = 1 was found in the packet's own depth output.
Twelve minors were applied.

NOT CLAIMED. No proof that F_<=2 covers everything above an effective N_0; no all-order singleton statement for
K = 2. Chan-Pak OPEN; Carmichael OPEN.

NEXT
- An effective N_0 for F_<=2. With this census it would give C = 1 coverage above 2647, and hence Chan-Pak.
- Derive the 3-adic factor. The dilated-root seed JJ(t) = (5,2,6,3) gives H = 3 H_t(m) + 2 - 2(m mod 2), which never
  hits N = 1 (mod 3).
- A second involution on the K = 2 chamber.
- Extend the census past 10^9 in segments.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
