# Independent audit of prime290

Scope: independent audit of the adjacent prime290 proof and checks. Independent test `independent.py` uses an undirected-component traversal, checks all 290 center residues, computes prime hits and enumerates every state. Both producer test modes (`python3` and `python3 -O`) produced identical receipts. The independent audit was adapted to explicit failure gates, so its optimized run enforces the same checks.

## Findings

No substantive mathematical errors found in the reviewed claims. Verified independently:

- φ(290)=112, and the directed +2 wheel has 81 edges, 31 path components: 25×4, 2×3, 2×2, 2×1, with isolates 27 and 263. The list of exceptional actual twin pairs not represented in the wheel is exactly (3,5), (5,7), (29,31). The edge (7,9) shows that local admissibility is not primality.
- The first +1 and −1 primes are 1451 and 2029; the first boundary twin pair is (6089,6091). The first union hit has eligible-prime rank 227. These are finite searches and not infinitude claims.
- For the alternating unit-filtered MRB series the periodic mean is −112/290=−56/145, so its leading coefficient is −28/145. Opposite ±1 orientations have zero mean; each positive/negative subseries diverges, so the ordinary series is conditionally but not absolutely convergent. Since f is eventually decreasing, f(290k+289)−f(290k+291)≥0 eventually and their sum telescopes by comparison to the total variation; the paired series is absolutely convergent. f(1)=0 deals with the initial +1 term.
- F(a)=(a²−29)/4 on odd numerators descends to 145 odd residue states. Direct computation independently reproduces all eight 3-cycles, one 5-state basin with distance counts (3,2,0), seven 20-state basins with (3,7,10), max preperiod 2. The cycle −7→5→−1→−7 corresponds to 283→5→289→283; 5 is not a 290-unit.
- All eight odd integers of absolute value at most 7 enter the 3-cycle. If |a|≥9 then F(a)≥13 and subsequently F(a)>a, proving divergence, so exactly those eight odd numerators give bounded x=a/4 orbits. In particular, modular periodicity does not imply bounded real orbit: 9→13→35→299 while residues return to 9.

Source verification: the stated prior work is present at origin/main commit 4c6f0ca13528f249e22184129966ef727d4db0e3 via `git show`; this checkout's HEAD is 74e81ed, so using an unqualified checkout path can misidentify it. The linked English translation of Dirichlet's 1837 paper is indeed https://arxiv.org/abs/0808.1408. For a supplemental source on cyclicity of finite-field multiplicative groups, see theorem 2 of https://www.maths.tcd.ie/pub/Maths/Courseware/FiniteFields/GF.pdf. No claim of novelty or original source for the rational 3-cycle is verified by this audit; `known exact orbit` is arithmetically true but not historically attributed.
