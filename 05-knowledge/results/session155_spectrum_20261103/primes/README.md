# Strong-tournament prime-lane gadget notes

Status: rigorous elementary bound plus finite, exactly reproducible witness computations. No infinite prime-lane completeness claim. No literature novelty claim.

## New escape certificates for two exceptional seven-prime values

Use vertex set 0,...,6. The i-th integer in a row-mask list has bit j set exactly when i beats j. Then

    H=35: [62, 60, 120, 16, 96, 8, 43]
    H=49: [62, 60, 56, 80, 96, 8, 39]

Both are strong, have exactly the displayed number of Hamilton paths, and have minimum feedback vertex number exactly 2. Deleting vertices 5 and 6 leaves a transitive tournament, whereas every single-vertex deletion contains a directed triangle. `verify.py` checks strongness by Boolean transitive closure, the deletion condition by explicit triangle checks, and H independently by (i) all 7! permutations and (ii) subset endpoint dynamic programming. The same file checks additional witnesses H=13 (order 5) and 29 (order 6). It uses standard Python only.

Neither 35 nor 49 occurs for any FVS-at-most-one tournament. Strip source/sink padding from an apex/transitive-deletion representation. The remaining nontrivial strong core has a binary orientation word b_1...b_m beginning 1 and ending 0. (A non-strong word has a leading 0 or trailing 1 and these are padded source/sink vertices; the residue becomes strong or trivial.) With state (H,A,p) initialized (1,0,0), append 1: (H,2A+1,1), append 0: (H+A+p,2A,0). As A before the final zero is at least 2^(m-2), H>=2^(m-2)+1. Therefore all strong cores of H<=49 have m<=7, and the verifier exhausts these words, yielding

    m=2: {3}; m=3: {5}; m=4: {9,11};
    m=5: {15,17,19,23};
    m=6: {25,27,31,33,37,39,45,47,51};
    m=7: {45,51,53,61,63,65,69,73,79,81,87,95,97,107}.

The recurrence and existing H=13 obstruction appear in `math/05-knowledge/results/session26_random_bridges_20261022/session26-tournament/CUT_ADDENDUM.md`. The exact 35 and 49 certificates were not located in a targeted local search; this does not establish novelty.

## Clone-growth obstruction to linear ladders

Let T be any strong tournament with at least three vertices, and v any vertex. Replace v by k>=1 transitive clones having all the same external orientations as v. Every v lies in a directed triangle: choose an out-neighbor and in-neighbor of v; if no out-neighbor beats an in-neighbor, those in-neighbors cannot be reached from v without a first crossing contradiction. The blow-up therefore contains the tournament with a transitive k-block and two singleton blocks directed cyclically. It has exactly 2^k+1 Hamilton paths: its transitive-deletion word is 1^k0 and the above recurrence gives H=2^k+1. The full blow-up has at least as many Hamilton paths, by vertex-insertion monotonicity. Thus no fixed strong host with one growing transitive clone class can have H(k) polynomially bounded, let alone form an unbounded arithmetic progression. This includes repeating one fixed port type in a two-apex transitive-word construction once the host is strong. It is not an obstruction to arbitrary FVS-2 constructions or an infinite-prime impossibility result. Adjacent prior art: THM-4145 treats homogeneous pair expansion, including strongness and exact H for two clones; the cyclic triangle calculation follows from the Session26 one-apex word recurrence. No novelty is claimed for this corollary.

For the triangle host itself, H(k)=2^k+1 never lies in the 7p lane: powers of 2 modulo 7 cycle among 1,2,4, none of which is -1. If H(k) is prime, k must be a power of 2: if k=2^a d with odd d>1, then 2^(2^a)+1 is a proper factor of 2^k+1. This leaves infinitude of prime values open.

## General clone-slot formula

This generalizes THM-4145 from two homogeneous clones to arbitrary k by the same rooted deletion/defect mechanism. Take a fixed tournament T, root r, and n-1 other vertices U. Replace r by k transitive, externally identical clones, k>=1. For each permutation pi=(u_1,...,u_m) of U, call a slot *valid* if a nonempty clone run fits its two boundary arcs: the left endpoint requires r->u_1, the right endpoint u_m->r, and internal slot i requires u_i->r->u_{i+1}. A defect is an internal adjacent pair oriented u_{i+1}->u_i. Every defect slot must be occupied. If any defect slot is invalid, pi contributes zero. Otherwise let t be its number of valid slots and q its number of defects. Its exact contribution is

    F_pi(k) = sum_{j=0}^q (-1)^j binom(q,j) (t-j)^k.

Each of the k distinguishable clones chooses one valid slot, and within each slot their transitive orientation fixes their unique order. Inclusion-exclusion forces the q distinct defect slots to be nonempty. Summing over all pi counts every Hamilton path once. If U is empty, H=1 separately. The formula extends to k=0 with 0^0=1, so then H is the count of the root-deletion tournament.

Consequently H(k)=sum_{j=0}^{n} c_j j^k for integers c_j depending on (T,r), and its ordinary k-generating function is rational with denominators among 1-jz. At k=2, the summand equals t^2, 2t-1, 2, or 0 as q=0,1,2,>2, exactly recovering THM-4145. It supplies no general prime-generation theorem. Together with the triangle lower bound, it shows every such strong-root family is exponentially bounded below and is not a polynomial in k.

`clones.py` checks the formula against an independent vertex subset-DP for all 1,140 combinations of every labeled tournament of orders 1-4, every root, and k=1..4; normal and optimized outputs agree. This is a finite control for the proof, not a substitute for it. A targeted search found THM-4145's k=2 case; no novelty or literature-priority claim for the all-k formula is made.

## Exploratory bounded census

`fvs2.cpp` enumerates all orientations on a transitive base of size m=0,...,8 and two marked apex vertices (2^(2m+1) codes each). It tests strong connectivity and counts paths by subset DP. The last line's m=8 strong-code count is 72,960. It compares results <=100,000 with the exact all-order one-apex spectrum in that numerical range. The latter exactness follows from the above length bound (m<=18 here); the implementation enumerates up to m=21 and prunes a prefix only when it cannot end in a zero without exceeding the bound. This census identifies 236 primes and 42 seven-times-primes realizable in this finite FVS-2 scan but absent from the entire one-apex spectrum. These are candidate counts from a single implementation, unlike the independent elementary small-witness checks; treat as exploratory, not audited or evidence of eventual completeness.

Reproduce: `python3 verify.py`; `g++ -O3 fvs2.cpp -o fvs2 && ./fvs2`. The file `fvs2.out` records code output through H=3000; `fvs2.progress` gives the full strong-code counts. The linear search was exploratory and is superseded by the clone-growth proof.
