# Independent audit: THM-4853 and the groups-of-five PRIMEGAME note (commit 8bc3d1c905)

**Verdict: SOUND WITH CORRECTIONS.** THM-4853 is correct (checked beyond its claims); the metagraph replication is
exact; pre-registered outcomes are reported truthfully. Problems: one false statement, one STRUCTURAL over-typing, two
credit gaps, wording.

**Abbreviations.**
- NOTE = `05-knowledge/results/groups_of_five_primegame_merged_metagraph_20261031.md`.
- THM = `01-canon/theorems/THM-4853-...md`.

**Method.**
- **Own code, run from the repo root.** Each script writes the `.out` file beside it.
  - `audit_primegame.py`: checks A1-A5.
  - `audit_metagraph.py`: checks B1-B3.
  - `audit_correspondences.py`: checks C(a)-C(f).
- **OEIS data, fetched 2026-10-31.**
  - JSON records for A007542, A007546, A007547 and A032742.
  - The b-files of A007546 (5000 terms) and A007547 (1000 terms), stored here.
  - The A007542 b-file (terms 0..8102), checked by hash.
- **`git fetch`.** The only newer commit, 426e639fe7 (session127), is unrelated; nothing further to credit.

## Findings

1. **ERROR** — the Fermat "exactly".
   - *Where:* NOTE §3, "Added here": "the strong points of symmetry are exactly the Fermat-order circulants plus K1".
   - *Evidence:* this is false. 17 is a Fermat order, yet all 256 circulant tournaments on Z_17 have a non-palindromic
     Hamiltonian path (C(a)). THM-4730 (B) excludes every odd order from 7 upward.
   - *Fix:* write "exactly the circulants of order at most 5". That 3 and 5 are Fermat primes is coincidental.

2. **OVERCLAIM (typing)** — the Fermat reading typed STRUCTURAL.
   - *Where:* NOTE §3 row 2: "THM-871 already explains their uniqueness: at Fermat orders every rotational tournament
     has Aut = Z_n. STRUCTURAL (inherited)". The same reading appears in §2.6 and §6.4.
   - *Evidence:*
     - Order 17 is also rigid, yet it has 16 rotational classes (C(a)) and none of them is pure-blue.
     - Uniqueness at orders 3 and 5 is just the count 2^((n-1)/2) = n - 1, which holds only for n = 3 and 5.
     - No predicate-preserving map exists. By the note's own definitions this is NUMEROLOGY.
   - *Fix:* retype it as NUMEROLOGY and delete "explains".

3. **OVERCLAIM (wording)** — "points of symmetry" versus "strong points of symmetry".
   - *Where:* NOTE §2.6: "The points of symmetry are the circulants of order 1, 3, 5".
   - *Evidence:* the points of symmetry are TT5, 1->C3->1 and T5, and only T5 is a circulant.
   - *Fix:* use §6.4's wording: the *strong* points of symmetry, K1, C3 and T5.

4. **CREDIT** — the black-line pairing is already canon.
   - *Where:* NOTE §5: "Why black lines come in pairs (a two-line corollary of THM-643's (sigma, tau) frame)".
   - *Evidence:* this is THM-4731 (3), the sigma-twin theorem (proved 2026-10-20). The argument is correct; I
     re-verified n = 3..6 (B1).
   - *Fix:* cite THM-4731 (3).

5. **CREDIT / OVERCLAIM** — novelty claims that OEIS already covers.
   - *Where:*
     - NOTE §2.1: "the OEIS attribution and the growth constants are new";
     - THM, "New here": the constants 2(1-c) and c/2.
   - *Evidence:*
     - **(a)** OEIS already makes both attributions:
       - P. von Brömssen's comment on A007546 (2026-04-30) names the 15/14 program;
       - A007542's program is the owner's list, and its 2016 comment identifies the 2^k indices as A007547;
       - M. F. Hasler (2026) gives A007542(A007547(n)) = 2^prime(n).
     - **(b)** A032742 carries the formula sum_{k<=n} a(k) ~ 0.165049... n^2 (A. Eldar, 2022). Its series is c/2, so
       s_O - s_G ~ (c/2)p^2 is an OEIS corollary. Only 2(1-c) is new.
     - **(c)** The closed forms sum session126's block lengths, so they are summed here rather than discovered.
   - *Fix:* credit A032742 in (e) and the OEIS comments in (d).
   - *Residual risk:* neither the session nor I consulted Guy 1983.

6. **MINOR** — HYP-9448's wording defect, repeated in THM.
   - *Where:* three places.
     - THM (a): "to the next state-free configuration 2^N 7^(D(N)-1)". The next state-free configuration after 2 is 15,
       and after 4 it is 30 (A1).
     - THM's proof: "this is the phase law (HYP-9448's statement)". HYP-9448's literal statement is false at step 1, and
       THM proves only its phase-entry form.
     - NOTE §9 suggests promoting HYP-9448. That needs this rewording first; NOTE §7 itself cites the counterexample 15.
   - *Fix:* write "the next configuration produced by 1/17", and "the phase-entry form of HYP-9448".

7. **MINOR** — test outcomes and cross-references in §1.
   - *Where:* NOTE §1: "P4 and P5 FAILED as stated (sections 6.2 and 6.4)".
   - *Evidence:* P5 is in §6.3, and only P4(a) failed.
   - *Fix:* write "P4(a) and P5 FAILED (§6.4, §6.3)".

8. **MINOR** — the graph non-isomorphisms are decided by edge counts.
   - *Where:* NOTE §2.3 and §6.1.
   - *Evidence:* the edge counts already differ: 19 and 28 against 18. Under a one-way reading of "preserves adjacency",
     the metagraph support does embed as a spanning subgraph of the co-occurrence graph (C(e)).
   - *Fix:* write "not isomorphic (edge counts differ)".

9. **MINOR** — the calibration is exhaustive, not random.
   - *Where:* NOTE §2.5 and §6.3: "9.6% of random targets".
   - *Evidence:* the 136 targets are enumerated exhaustively, not sampled. The post-hoc chance of 0.33 also applies the
     hit rate for (0,a,b) targets to tau = (2,4,12), which lies outside that family.

10. **MINOR (wording)**:
    - NOTE §5, "bit 0 orients the arc downward": this means a clear bit gives x -> y.
    - THM's proof, "share one 4-counter, 7-state machine": the state-free rows differ between the two lists.
    - THM's proof, "The divisors d = N-1, ... are tried": these are candidate divisors.
    - NOTE §6.5, "The one exact count-level contact found here is THM-4853": THM-4853 is not a contact with the
      metagraph.

## Independently confirmed

**THM-4853 (a).**
- Literal integer FRACTRAN on both lists, through the first 25 primes:
  - each configuration produced by 1/17 is 2^N 7^(D(N)-1);
  - the only pure powers of 2 that occur are 2^p.
- On the counter machine, in 3540 cases (list, N, d) with N <= 60:
  - the trial costs are 6N + 2q + 2 and 4N + 2N/d + 2, with the right end configurations;
  - the state-free segment costs N or N + m;
  - the reload costs 2N - 1.
- Proof checked line by line; no off-by-one.

**THM-4853 (b)-(d).**
- The closed forms agree with my accelerated simulator for all 168 primes below 1000.
- They match the 34 OEIS data terms, **all 5000 terms of the A007546 b-file** and **all 1000 terms of the A007547
  b-file**, so the theorem can claim more than 34 terms.
- The difference law holds on those 1000 terms.
- The owner's orbit equals the A007542 b-file (orbits first differ at index 131).

**THM-4853 (e).**
- c = 0.33009804, from the primes up to 10^7, with a tail below 3.5e-9.
- The quoted ratios 1.0015 and 0.9932 are correct (1.00006 and 0.99993 at p = 48611).
- The proof is correct and "~" is justified; its "+O(log log x)" is superfluous but harmless.

**Metagraph.**
- The THM-643 census holds for n = 3..6.
- I recomputed the full 12-class table.
- The merged metagraph has 8 blue lines (one cycle, 5-6-9-8) and 12 black pairs, two of which are loops.
- The S70 / THM-4729 dictionary is exact and forced: there is exactly one colour-preserving isomorphism.
- session126's least-mask numbering and mask facts (7, 17, 55) are reproduced; its gauge gives converse tournaments.

**Correspondences.**
- Circulants: all those of orders 1, 3 and 5 are pure-blue; none of order 7 is.
- Terras cycles: {0} (D = 1), {1,2,4} (D = 5), {5,10,19,20,40} and {7,11,14,22,28} (D = 23), five cycles (D = 101).
- Relabelling counts: 864, from three disjoint sources of 288; then 96 and 48.
- The battery values and the 13/136 calibration.
- Graph sizes: 24 arcs giving 19 pairs; 28 edges; 18 edges plus 2 loops.
- The residue densities, including q/14, 3q/7 and q/2.
- F = U + S holds exactly at N = p, p^3 and p^2qr, checked up to N = 20000.

**Pre-registration.**
- The test script is byte-identical in the pre-registration and final commits.
- P4(a) and P5 failed exactly as reported.
- "No HYP" is justified.
