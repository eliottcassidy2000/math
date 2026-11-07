# Session214 forced-shuffle independent audit

**Verdict: PASS.** The all-order combinatorial proofs and the exact small-increment classifications in the frozen packet are valid. No mathematical correction was required. Chan–Pak completeness remains OPEN.

## Audited bytes

Packet: `05-knowledge/results/session214_chan_pak_combined_proof_20261107/forced-shuffle/` (author staging directory `session214-combine`).

- PROOF.md: `851b8187b36b5fe0f4e42adf4858bf3caad32417aac9fcd9f8c68bd5ab0c28b6`
- README.md: `b8e12105a7a5b69aa97bd0d50c9d94db5190cd158971454d391a6bf6047ab894`
- verify.py: `96596f895f086d37a44c2333955caeebb41658b6751f13c8206d016d6fbd9c04`
- RESULT.json and RESULT.optimized.json: `3502cb3482f3250a754b3a062444781eafbaf09d5980107fa30e7f62ce1e1468`
- SOURCES.json: `3f33f71c7d1ae39a2f835344caad12dd135cf934d350272d2305a1b0fedbe3e8`

The author manifest was rehashed successfully. Source inspection used freshly fetched `origin/main=4c63fd58b566f3e006eb87eaedcf48824741cc84`, without modifying the shared dirty checkout. HYP-9469 still called the general ordered-path inequality numerical there. THM-4745 statement 1 already credits Busch's directed-word merge and restriction recovery, and the packet correctly acknowledges that inheritance. No external novelty claim is certified.

## Mathematical audit

The maximal A index with any exit into B exists for every pair of monochromatic Hamilton words under the stated forward-arc hypothesis. Its least B outneighbor gives the predecessor join; maximality gives the final B-to-A suffix join. Greedy head comparison preserves both path orders. Empty prefixes, empty suffixes, singleton colour classes, and first/last index cases are valid. Restriction to colours is a left inverse, and the selected A-to-B step makes every output nonordered. Therefore all H(A)H(B) inputs yield distinct images disjoint from the q ordered words.

This proves h>=p+q, while ordered-word restriction proves q<=p. Rédei oddness gives 2q+1<=h. The exact slot identity is 1 plus the number of A-to-B transitions. Its injection is independent of the shuffle injection and gives H(Y)>=2h-q. Consequently H(Y)>=max(h+p,2h-p)>=(3h+1)/2; the last inequality uses integral p and odd h. Even increments give the stated ceiling. The zero-growth case includes empty cut parts and is exactly the dominance join.

The arbitrary induced-extension conclusion follows by the first positive step; it does not assume a fixed growth ratio at zero steps. The asymptotic R_n sole-exit count correctly partitions paths by v's position and the forced a,b,0 segment. The n=3 exception is explicitly separated. The ratio tends to 3/2.

For increment 2k, h<=4k-1 and H(Y)<=6k-1. Factoring out components outside the marked strong component gives k=gd with odd g. For k<=4 only g=1,3 are possible. The cited strong-order bound excludes order>=7 because the child count is at most23<25. Thus an exhaustive labeled strong-core enumeration through6, followed by the permitted factors, genuinely proves all-order exhaustion. Singleton padding does not contradict finite sets of parent counts.

## Independent computation

`audit.cpp` was written independently of the author's verifier. It generates every labeled tournament through order6 and directly enumerates all directed vertex words, not subset-DP path totals or inherited class representatives. It explicitly constructs and compares every image pair at every active cut.

- 33,867 labeled tournaments through6.
- 2,044,814 active cuts.
- 7,326,350 injected path pairs, all valid and pairwise distinct within each cut, with exact inverse subsequences and no ordered images.
- 30,606 fresh child enumerations at every active cut through parent order5.
- Every labeled strong rooted deletion through order6 independently checked.
- Parent sets at increments2,4,6,8: {1,3}; {1,5}; {3,5,9}; {1,5,9,11}.

`sharpness.py` separately uses forward subset DP to check the parent, marked endpoint, and sole-exit child formulas for R_n at n=3,...,14.

Reproduce:

    g++ -O2 -std=c++17 -Wall -Wextra -Werror -pedantic audit.cpp -o /tmp/session214-shuffle-audit
    /tmp/session214-shuffle-audit > /tmp/session214-shuffle-audit.json
    cmp receipt.json /tmp/session214-shuffle-audit.json
    python3 sharpness.py > /tmp/session214-sharpness-audit.json
    cmp sharpness.json /tmp/session214-sharpness-audit.json

All checks use explicit runtime failures. A second compiled replay reproduced receipt.json byte-for-byte. This audit certifies the stated mathematics and finite universes, not unrestricted spectrum completeness, generic-cut equality classification, or literature priority.
