# Complete small paired-gate spectrum and its least 7p output

Session197. **FINITE-EXACT with an all-order cutoff.** In the paired construction

    P(T)=C3[T,T^op,Tr2], T nonempty,

the least output of the form7p with p prime is

    37050461=7*5292923.

This is a minimum within this named construction family, not within all tournaments. The output is attained on16 vertices. The complete paired spectrum below73529109 has119 values, including8 primes and exactly3 values of the form7p.

## 1. Named all-order input and cutoff

The frozen `THREE_COVER_CYCLIC_GATE.md`, SHA256 `0630c4ff7f9a52e69cbce4f18cdac8678609d4565a3e42f0afec8d588afe37cf`, proves that for the first three unordered path-cover counts (h,u,v) of T,

    H(P(T))=5h^2+12hu+28u^2+48uv+72v^2.

It also proves the sharp order-specific lower bound H(P(T))>=Phi_n, where n=|T| and

    Phi_n=2*9^n-8*6^n+13*4^n+4*3^n-22*2^n+15.

Phi_n is strictly increasing and Phi_8=73529109. Therefore every paired output strictly below73529109 has an input with at most7 vertices. The enumeration below is complete at every order for that output interval, despite examining only orders1–7.

This packet is an arithmetic/computational companion. Its all-order reduction uses exactly that named gate theorem; it does not independently replace the gate theorem's proof or audit.

## 2. A complete label filter, not an isomorphism claim

For every n=1,...,7 the checker enumerates all2^(n(n-1)/2) orientations of the labelled complete graph. It retains a graph only if its outdegrees, in vertex-label order, are nondecreasing.

Every tournament is isomorphic to a retained labelled graph: simply relabel its vertices in nondecreasing order of outdegree, breaking ties arbitrarily. Relabelling preserves all path-cover counts and the paired output. Thus this filter loses no possible value. It is not a canonical-labelling algorithm; tied degrees can leave many copies of an isomorphism class. No count of isomorphism classes is asserted.

The exact retained counts for n=1,...,7 are

    1,1,3,9,61,529,11559,

for a total12163 inputs. For each retained graph, endpoint subset DP computes H(S) for every induced vertex subset S. Then

    c_2(S)=1/2 sum_(empty!=A proper subset of S) H(A)H(S-A),
    c_3(V)=sum_(A containing the least vertex) H(A)c_2(V-A).

The second formula selects the unique block containing the least vertex, so it counts each unordered three-path cover once. Empty residual two-covers contribute zero. The resulting integer triple gives the exact paired count. Values below the cutoff are deduplicated by their integer count, retaining an actual input graph and triple as a witness for each value.

All intermediate counts fit comfortably in signed64-bit integers: the input order is at most7, and every paired output before cutoff comparison is a Hamilton count on at most16 vertices, at most16!, which is below2^63. The subset counts and their component products are smaller still. No floating-point arithmetic is used.

## 3. Exact outputs

The complete119-value list, one physical input row list and cover triple for each value, and primality/compositeness data are in `paired_gate_catalog_normal.json`. The8 prime outputs are

    5,509,120293,189701,38997197,49745741,55619813,58696949.

The3 outputs of the form7p are:

- 37050461=7*5292923, input covers(23,273,611), input rows[2,48,3,7,13,29,63].
- 46937093=7*6705299, input covers(37,339,673), input rows[40,17,3,6,13,30,63].
- 56228501=7*8032643, input covers(47,387,729), input rows[2,96,3,7,15,29,61].

Every listed input has seven vertices, so the paired graph has16 vertices. Primality is determined by complete trial division: the recorded field is1 if no divisor2<=d<=floor(sqrt(N)) exists, and otherwise is the least such divisor. Each output H has this field, as does H/7 whenever7 divides H. Thus the classification of all119 values is reproducible without a probabilistic primality test.

The first7p value is globally least within the paired family. Every smaller candidate is below Phi_8 and therefore appears in this complete enumeration if attainable. Every order at least8 is bounded above this entire interval. The proof does not need to enumerate larger orders or assume that arbitrary triples (h,u,v) are physically realizable.

## 4. Independent physical controls supplied with the source

`verify_paired_gate_catalog_witnesses.py` reconstructs every retained witness input and recomputes its cover triple. For each of the8 prime and3 seven-prime output values it builds the actual paired graph, checks every anti-involution incidence and strongness, and directly counts its Hamilton paths by vertex subset DP. It also rechecks every primality/compositeness field by integer trial division. This is a complementary physical control, not a substitute for the completeness of the finite orientation loop.

The output witness for the first7p value can differ by labelling, converse or a count-equivalent input from the separately displayed graph in the gate theorem. Equality of a cover triple is used only for the exact count formula, not to infer an isomorphism.

## 5. Reproduction and exact scope

    g++ -O3 -std=c++17 verify_paired_gate_catalog.cpp -o paired-catalog
    ./paired-catalog > paired_gate_catalog_normal.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify_paired_gate_catalog.cpp -o paired-catalog-debug
    ./paired-catalog-debug > paired_gate_catalog_debug.json
    cmp paired_gate_catalog_normal.json paired_gate_catalog_debug.json
    python3 verify_paired_gate_catalog_witnesses.py
    python3 -O verify_paired_gate_catalog_witnesses.py

The C++ enumeration is standalone. The physical witness checker uses the frozen `verify.py` helper, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`. All checks use explicit exceptions and remain enabled under optimized execution.

No equality classification of input isomorphism classes, global tournament minimum, infinite prime family, or unrestricted-spectrum coverage is inferred from this restricted catalog.
