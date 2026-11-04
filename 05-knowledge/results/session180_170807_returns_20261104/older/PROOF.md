# Prefix realizability closes the arbitrary-prefix loophole

## Result

For the complete suffix menu in `../returns/search.cpp` (zero, one, or two edge reversals of the pinned fourteen-vertex strong suffix, keeping only strong suffixes, and every actual marked vertex), **arbitrary nonempty tournament prefixes add no return beyond the transitive four-vertex prefix already found**. The unique surviving arithmetic row is

    suffix edge IDs 31,40 (zero-based lexicographic pairs {4,31},{6,22});
    marked index 2, actual physical name 4;
    h_B=5655, F_w=31673;
    H(A)=1, P_A(2)=16;
    H(T)=5655+16*31673=512423=3*170807+2.

There is no +4 target 512425 in this entire suffix menu, even with prefixes of unbounded order. This is a menu-specific exclusion, not an obstruction to larger suffix rewrites, nonuniform attachments, or all tournament constructions. Legality of the old-edge reversals, named joint ports, and finite-field phases are separate obligations handled by the return/port packets; this prefix proof does not infer them from scalars.

## 1. The faithful prefix observer

For a tournament A on n vertices let H(A) count directed Hamiltonian paths, taking H(empty)=1. Let P_A(2) count ordered covers by two directed paths, allowing either slot to be empty. Then

    P_A(2)=sum_{S subset V(A)} H(A[S]) H(A[V(A) minus S]).   (1)

This retains the exact two-slot interface of the zero-ready attachment. For a strong marked suffix (B,w), adding A -> B and a root r beaten only by w gives

    H(T)=H(A)H(B)+P_A(2)E_w(B,2).                        (2)

The suffix and actual named mark remain part of the state. Formula (2) is not a general graph replacement invariant.

## 2. Uniform finite reduction for an unbounded prefix

Every induced tournament has at least one Hamiltonian path. Hence every summand of (1) is positive, proving

    P_A(2) >= 2^n.                                      (3)

Let a_1(A) count permutations with exactly one backward consecutive adjacency. A permutation with no backward adjacency admits n-1 internal cuts into two nonempty paths and two placements with an empty slot. A permutation with one backward adjacency admits exactly the cut at that adjacency; a permutation with two or more admits none. Therefore

    P_A(2)=(n+1)H(A)+a_1(A) >= (n+1)H(A).              (4)

This is the degree-two specialization of the older endpoint-jet/path-cover dictionary, THM-3134, equations (3), (5), (6). It is proved directly here and does not rely on the older numerical transcript.

Rédei oddness also gives P_A(2) divisible by four when n>=2: pair each subset with its complement in (1); there are 2^(n-1) pairs, each with odd product. The n=1 exception is P_A(2)=2. H(A)=1 forces transitivity: a directed triangle has three paths, and repeated insertion of a new vertex into a path injects old paths into new ones. Thus any nontransitive tournament has H>=3. For a transitive A, each induced H is one and P_A(2)=2^n.

Consequently, if an exact arithmetic candidate specifies (h,p), every realization has

    n <= min(floor(log2(p)), floor(p/h)-1).             (5)

No arbitrary order cutoff is necessary.

## 3. Complete arithmetic candidate handoff

Writing h_B=h and F_w=f, a prefix for target N must have an odd positive x=H(A) and even y=P_A(2), with

    hx+fy=N, y>=2x.

The return search enumerates every odd x satisfying x(h+2f)<=N, tests divisibility, and retains every such (x,y). Its 4,187 suffix candidates (3,815 strong) yield exactly 53 arithmetic rows for the two targets. The largest retained y is 198. This is a maximum of the complete necessary-candidate list, not an imposed cap.

By (3), every realizable candidate therefore has n<=7, because 2^8>198. The complete prefix atlas need only cover orders one through seven.

## 4. Independent prefix census

`prefix_atlas.cpp` independently enumerates all

    1+2+8+64+1024+32768+2097152 = 2,131,019

labelled tournaments of orders one through seven. It uses endpoint subset DP to compute H on every subset, then evaluates (1). There are exactly 22 distinct pairs with P_A(2)<=198; they are listed with a realizing labelled tournament in `prefix_atlas.json`.

`verify.py` checks every one of the 22 witnesses by a separate permutation algorithm implementing the zero/one-backward-adjacency count above. It also independently enumerates every labelled tournament through order five using that permutation algorithm and checks the complete small atlas. It then tests all 53 arithmetic candidates against the atlas; the only survivor is (1,16) at the row displayed above. The return worker separately built its own full prefix census, providing another implementation.

The C++ exhaustive component proves completeness within its finite universe; the direct bound (3) supplies completeness over every prefix order. This combination is the concrete transfer of the old owner-window pattern: first prove an unbounded parameter is bounded on every residual candidate, then exhaust that finite residual. The LRC metric and ownership formulas themselves are not transferred.

## 5. Older techniques recovered, with exact scope

- `01-canon/theorems/THM-3134-tournament-endpoint-jet-and-c3-newton-profile-transform.md`: permutation cuts give the precise prefix response and bound (4). The full endpoint jet is needed when an operation sees three or more module runs; H alone is insufficient.
- `01-canon/theorems/THM-4162-rooted-pair-mixed-two-ear-tensor-and-enumeration-free-johnson-cosets.md`: full birth fields have directed-cut capacity tensors, and fixed-cardinality gcd/cosets can be computed in O(N^3). This is a concrete fast future filter for suffix birth searches; it does not assert every lattice point is realizable. The 52 rejected arithmetic prefix rows here illustrate that exact distinction.
- `01-canon/theorems/THM-4145-rooted-homogeneous-pair-expansion-two-defect-formula.md`: a nonuniform extension must retain actual rooted defect slots and the mixed two-ear tensor. Squaring the zero-defect term loses paths.
- `01-canon/theorems/THM-1169-three-four-carrier-owner-window-closure.md`: the usable methodological pattern is a proved finite residual bound plus endpoint-exact enumeration. There is no claimed map from its LRC carrier windows to tournament return counts.
- `session33-ports/` contains the context-relative continuation-row quotient: merging states is safe exactly when their responses agree for every allowed continuation. Simultaneous symmetry or matching scalar counts does not establish that equivalence.

F49/F343/81/27 observations are not inferred or identified in this packet. A valid grammar must preserve their separate carriers and named joint responses in addition to satisfying (2).

## Reproduction

    g++ -O3 -std=c++17 prefix_atlas.cpp -o prefix_atlas
    ./prefix_atlas > prefix_atlas.json
    python3 verify.py
    python3 -O verify.py

`verify.py` reads the sibling return packet's complete `hits.jsonl`. The checks use explicit exceptions, not Python assertions.
