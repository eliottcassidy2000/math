# Independent audit: valuation compatibility and cover geometry

Verdict: PASS for the scoped elementary statements inspected. No source changes requested. Audit date 2026-10-23. This is not an audit of all prior cited repo theorems or a transfer to native tournament wiring.

## Least divisor-height realization

The proof in PROOF.md section 6 is correct, including k_b=0 and N_min=1. Every candidate valuation vector dominates the lower envelope L. Exact-height conditions are disjunctions of upper bounds. If a disjunction already fails at L, it fails at every larger vector. If none fails, L itself is a realization and divides every other realization. This is a least element, not merely a numerical minimizer. The positive-N>1 caveat is essential and properly present.

The cycle fixes e_2=5,e_5=3,e_3=7,e_7=4. Its composite condition h(4)=2 is redundant but compatible, and no other prime exponent is constrained. Adding 13->2 and 11->1 correctly recovers their exact exponents. The quotient identifying 13 with 4 is a forward congruence on this functional graph, not a congruence of the integers or an incoming-neighborhood isomorphism.

Independent extension test: all 720 two-base, four-target systems on ten bases supported at primes 2,3,5. Direct exponent-box search agrees with the lower-envelope criterion. The box contains every potential lower envelope since each base exponent is <=2 and targets <=3. Source's 2240 three-base test also replays exactly in normal and optimized Python.

## Monodromy and decomposition

The double-star polynomial's derivative and normalization are correct. Independent symbolic expansion checks derivative, critical value and power decomposition for every 1<=a,b<=7, including edge cases omitted by the source's >=2 test.

The generated group proof is sound. Detail clarifying the semidirect-product notation: in cyclic coordinates modulo n=dm, distance-a transpositions generate all permutations within each residue class modulo d. The n-cycle rotates these classes. To see an actual splitting, choose the block rotation (r+dk)->((r+1 mod d)+dk). It has order d and differs from the n-cycle only by a permutation within a block, which already lies in the base group. Thus the extension splits as S_m wreath C_d. The conclusion S_n iff gcd(a,b)=1 holds, including a=b=1. Source BFS checks replay.

The finite (5)(2) monodromy cannot be a connected unbranched cover of C*, since one generator then has two orbits. Adding the joining transposition yields a connected cover with extra branching. The explicit polynomial supplies existence independently of an abstract realization theorem.

## Sphere and orbifold distinction

The source correctly separates genus-zero source surface from the induced marked orbifold. Branch orders for (5,2) are (10,2,7), not (5,2,7). The characteristic is -9/35. Pullback cone orders 2 at z=0, 5 at z=1, and 2 at the five simple points over 1 give -9/5, equal to degree 7 times the base value. This is compatible with the underlying sphere's Euler characteristic 2. The complex-plane/logarithm discussion also correctly distinguishes a branch convention from an intrinsic singular negative real axis.

## Replay and bounds

- Both supplied scripts pass normal/-O with identical outputs.
- independent.py passes normal/-O identically: 720 compatibility systems and 49 symbolic polynomial cases.
- Exact general claims rely on proofs, not those finite bounds.
- Historical group-order origin and prior THM-2996/4076 are provenance citations, not freshly verified group-order or prior-theorem audits here.
