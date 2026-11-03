# Three-port replacement and uniform-triangle obstruction

Status: exact all-order identity plus a finite, representative-specific three-port construction and finite census. No general affine transition or spectrum theorem is claimed. The values 41 and 43 were already known, e.g. the adaptive strong-tournament witnesses in Session 150 (`adaptive/check199.out`); the interest here is in contrasting uniform and nonuniform replacements of the same seed root. Replacement deletes the old root, and is **not** an extension preserving the parent as a subtournament.

## Uniform module identity

Let T be a tournament, v a vertex, U=T−v, and H=H(T). Replace v by a directed triangle x→y→z→x. Give every triangle vertex the same in/out-neighborhood in U that v had. Denote the result T*. Define P_r, for r=2,3, as the number of ordered decompositions (A_0,...,A_r) of U into directed paths, with A_0 and A_r allowed empty and each internal A_i (1≤i<r) nonempty, subject to:

- the last vertex of A_0 points to v if A_0 is nonempty;
- v points to the first vertex of A_r if A_r is nonempty;
- for each internal A_i, v points to its first vertex and its last vertex points to v.

These are ordered path decompositions on disjoint vertex sets, not necessarily partitions into contiguous pieces of an existing T-Hamiltonian path. Empty paths count once. Then

    H(T*) = 3 H + 6 P_2 + 6 P_3.                         (1)

Proof: split each T*-Hamiltonian path into maximal consecutive triangle blocks. One block uses all 3 triangle vertices; their three possible directed orders, after contraction, correspond bijectively to T-Hamiltonian paths. Two triangle blocks have sizes (1,2) or (2,1), with three choices of the directed pair for each case, and correspond to P_2. Three blocks are singletons, in any of 3!=6 orders, and correspond to P_3. These cases are disjoint and exhaustive. In particular H(T*) ≡ 3H (mod 6) and H(T*) is always divisible by 3. Thus this one uniform-module architecture cannot realize 3H+2 or 3H+4 for any parent. This does not exclude nonuniform ports, root changes, arbitrary interfaces, or adaptive constructions. The divisibility also follows from the free action of rotating the triangle on Hamiltonian paths.

## Strong nonuniform outputs for one parent

Use vertices 0,...,4, with the 10 pairs i<j in lexicographic order, least-significant bit first; bit 1 means i→j. Let the seed be code 41. Its H is 13 and it is strong. Delete root 0 and relabel retained vertices 1,2,3,4 as local indices 0,1,2,3. Add x,y,z with x→y→z→x. A port mask m encodes x→local-index i if bit i of m is 1 (similarly for y,z). Leave the induced tournament on the retained vertices unchanged. The exact masks and independently checked results are:

| masks (x,y,z) | H | strong? |
|---|---:|---|
| (15,7,1) | 41 = 3·13+2 | yes |
| (15,15,9) | 43 = 3·13+4 | yes |
| (6,1,1) | 41 | no |
| (3,4,1) | 43 | no |

Literal permutation enumeration checks each path count, and transitive closure checks strongness. This is not an operation on all H=13 representatives or on an invariant family, and does not supply the missing 3H branch. It illustrates how unequal ports break the uniform module's threefold symmetry without asserting that the port choices extend to arbitrary parents.

## Finite search and limits

For each root, `census.cpp` checks every one of 16^3=4096 ordered port-mask triples, with the same directed triangle and fixed retained tournament. It finds the following numbers of labeled mask triples giving each target; entries are all / strong:

    root 0: 39: 0/0, 41: 81/45, 43: 54/18
    root 1: 39: 0/0, 41: 24/24, 43:  0/0
    root 2: 39: 0/0, 41: 24/24, 43:  0/0
    root 3: 39: 0/0, 41: 81/45, 43: 54/18
    root 4: 39: 0/0, 41: 24/24, 43:  0/0

This finite negative result concerns precisely this seed, deletion of one original vertex, and a *cyclic* replacement triangle; it does not rule out 39 via a transitive replacement, a different seed, other attachments, or other methods. Session 150 already constructed strong 39. It does not show failure of an interface bank at larger sizes.

## Reproduction

Run in this directory:

    python3 verify.py
    python3 -O verify.py
    g++ -O3 -std=c++17 census.cpp -o /tmp/census154 && /tmp/census154
    g++ -O0 -std=c++17 census.cpp -o /tmp/census154_O0 && /tmp/census154_O0

The verifier enumerates permutations for selected examples and the five uniform replacements, and separately constructs boundary path covers for (1). The census uses a forward subset/last-vertex DP from `helper.hpp`, different from literal enumeration. Frozen outputs are `verify.out`, `verify.opt.out`, `census.out`, and `census_O0.out`. Their identical optimization-mode pairs were compared. The finite checks support the implementation but do not replace the proof of (1).
