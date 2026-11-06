# Independent review of the refined rooted-profile order cap

**PASS, as a corollary of the cited sharp strong-tournament minimum.** No profile occupancy or eventual spectrum threshold follows from this order bound alone.

Reviewed proof SHA256: `41d8a8bde4d31b480a6324b4f924857f8fec60832e944317ca790d669c3a4b02`. Every producer-manifest input is frozen under `reviewed-inputs/`.

## Mathematical review

For s>0 the root is in the initial strong component C. Deleting all later singleton components preserves the physical pair (B,s). If q is the product of the retained tail counts, the exact identities h=q h(C), L=q L(C), and s=q s(C) imply B=q B(C) and q divides gcd(B,s). The tail has no singleton components and hence order at most G(q), using the inherited exact reduced-component minimum.

The boundary B=s is precisely the singleton initial-component case. In a nontrivial strong component each vertex starts a Hamilton path, so h(C)>s(C), while L(C)>=d(C)>=s(C). Thus B(C)>s(C). For B=s, the same literal reduced representative has order at most 1+G(B).

For B>s, the actual edge extension E(C), with root->a->C-root, is strong, has one more vertex and H(E(C))=B(C). This gives |C|<=S(B/q)-1. The actual q is in the stated divisor set Q, so the maximum K over that set bounds the same singleton-deleted representative.

The inequality K(B,s)<=G(B)-1 is valid even when a maximizing divisor is not physically realizable. The sharp minimum realizations at orders S(B/q) and G(q) have an ordered join of count at most B and reduced order S(B/q)+G(q). Hence that order is at most G(B). This argument compares extrema; it does not claim an exact B-witness exists. The alternative direct proof is also correct: E(U') has no singleton components and count B.

The combined elementary bound applies to the same induced representative U', so taking the minimum of the two bounds does not mix incompatible witnesses. The definition G(1)=0 and the empty-tail product 1 handle q=1 and the singleton source correctly. An empty Q proves nonexistence, whereas a nonempty Q is only a necessary numerical condition.

At base X=137, the eighteen caps are exactly

    r=1: 9,9,9,9,9,9,9,8,8
    r=2: 9,10,9,9,9,9,10,8,8.

The maximum is 10. Both (B,s)=(135,5) and (125,20) have a maximizing q=5. No sharpness of this order cap for each individual pair is asserted.

The external input remains CITED: Busch's sharp f(n), together with the previously proved reduced-component minimum g(N). Their frozen source and immutable provenance are retained in the input packet. The audit does not supply a new proof of Busch's theorem.

## Independent controls

The independent C++ checker uses literal directed-path enumeration, including every subset for h, s, and L. Its SCC routine uses reachability searches, separately from the producer's transitive-closure routine. It recomputes the reduced minima through order 40 by minimum-product partition dynamic programming from f, rather than using the proposed closed formula to define G, and then compares that result with the formula.

Both optimized and unoptimized builds agree on:

- All 33,867 labeled root-zero tournaments through order 6.
- 27,476 positive-start inputs and 1,100 inputs with B=s.
- Exact SCC scaling of h, L and s, singleton deletion, the physical edge-extension identity, the B=s boundary, both order bounds, and the elementary bound on the same representative.
- 259,649 arithmetic pairs with 1<=B<=10,000 and 1<=s<=26, B>s.
- Every entry of the eighteen-request base-137 cap table.

The producer's optimized and debug controls were also rerun and matched its frozen receipts. Validation uses explicit exceptions, not removable assertions.

Independent source: `check.cpp` plus `independent_paths.hpp`. The path-counting helper was authored for the independent small-profile audit; it imports no producer count or SCC code. The header also contains unused isomorphism helper functions from that audit.

## Reproduction

    g++ -O3 -std=c++17 check.cpp -o check
    ./check > NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 check.cpp -o check_debug
    ./check_debug > DEBUG.json
    cmp NORMAL.json DEBUG.json
    g++ -O3 -std=c++17 reviewed-inputs/verify.cpp -o producer
    ./producer > PRODUCER_NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 reviewed-inputs/verify.cpp -o producer_debug
    ./producer_debug > PRODUCER_DEBUG.json
    cmp PRODUCER_NORMAL.json PRODUCER_DEBUG.json
    cmp PRODUCER_NORMAL.json reviewed-inputs/NORMAL.json

Independent check.cpp SHA256: `649f4d41d66683240d667b475894ea089cf5afe46411536219736646740abe0e`.

Identical independent normal/debug receipt SHA256: `3a9578fc3e39dc076116834074bab2de92fdc78794dceda5ab6b5bf6a6c3d0c7`.

The manifest freezes all sources, reviewed inputs, and receipts, excluding local binaries. F(8,137) and the eventual threshold remain open.
