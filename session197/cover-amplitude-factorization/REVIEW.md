# Independent review: minimum-cover leading-amplitude factorization

**PASS. No correction requested.** The theorem a_(ν+1)=ν! L_C R_C>0 holds for every minimum crossing cover, with the stated labelled-slot and empty-path conventions. This review applies to the following frozen source files under `session197-prime-atoms/cover-amplitude-factorization/`:

- `PROOF.md`: `2a0451a8473cd17f38a9140c18788fca50325d6f07d8a2716281252f871fbfd6`.
- `verify.cpp`: `3abdad3bcfcf2953af4f0a811fb3dac38610e9e2777d8c95db65715aaeb69581`.
- `verify.py`: `10e9c13dae637d8d1fbfa59951b09f69aaa0ba67492436af54c7a8c2443a3a08`.
- `receipt-normal.json`: `90b8d34e9d5dbad1769a1882ffa1c9b3f306e2b65948f7aba0158866be038e62`.

## Mathematical review

The maximal-run word has exactly ν disjoint nonempty internal fragments. Each must use an O→I edge and therefore a cover vertex. Since the minimum cover has exactly ν vertices, every internal fragment uses exactly one, and no cover vertex can occur in an outer fragment. This exhaustion is the crucial step; no assumption that a fragment contains only one crossing edge is needed.

After removing the cover, every I'→O' arc is uniform. Before the unique cover vertex an internal fragment cannot enter I', because that would require an uncovered O'→I' edge. After the cover vertex it cannot revisit O' and still finish in I' for the same reason. The two outer fragments consequently lie wholly in I' and O'. This gives precisely the two independent half-cover objects, including the mandatory/optional slots and their endpoint incidences. Their labels are the actual cover vertices; the ν! factor orders these distinct internal fragments without identifying any graph automorphisms.

For positivity, every edge in a maximum matching meets the size-ν cover exactly once, and every cover vertex is met. Otherwise the ν disjoint edges would require more than ν cover vertices. These matched outside vertices provide the mandatory singleton half paths. The two unused halves have Hamilton paths and serve as the outer slots. This proves positivity for each chosen minimum cover, including covers whose factors differ.

The ν=0 case, empty sides, and singleton source are treated correctly: the graph is I→r→O and the count is constant. The root-fixed anti-involution corollary is valid only when it preserves the chosen cover setwise; then the two half-cover objects are exchanged and ν is even. The source does not overextend this special square factorization to every symmetric cover.

The distinction between positive clone multiplicities and the deletion value at multiplicity zero is explicit and correct. No realization claim for arbitrary factor pairs or conclusion about lower coefficients follows from the theorem.

## Code review and independent checks

The source exterior-gap DP and the two half-cover subset convolutions count different objects and agree over every enumerated minimum cover. On the right side, reversing the induced I' tournament turns a required starting vertex into the endpoint used by its path DP; the boundary incidence is still tested in the original direction. Empty mandatory slots and the unrestricted outer slot are handled exactly.

The complete unlabelled extension recursion and degree-class canonicalizer are exact. The complete labelled universe is a separate enumeration. Larger seeded controls are correctly labelled nonexhaustive. Counts and intermediate sums are safely within unsigned 64-bit range at the stated orders.

- Source normal and Python-optimized replays reproduce receipt 90b8d34e9d5dbad1769a1882ffa1c9b3f306e2b65948f7aba0158866be038e62 exactly: 9,774 rooted tests and 17,002 cover identities.
- An assertion-enabled unoptimized C++ replay gives the same receipt.
- My separate `literal.py` uses no matching algorithm and no subset path DP. It finds every minimum cover by brute subsets, enumerates root-run words by literal exterior permutations and gap choices, and enumerates both half-cover objects by permutations and weak compositions into labelled slots. On every labelled tournament through order 5 it independently checks 5,405 roots and 8,747 factorization identities. Normal and optimized Python agree.

The independent files are in `audits/cover-factor/`:

- literal.py: b1bbdb276dfd5cb885ef7b5631ca596d5af3161dd5b8476e3c186b929fb2cc16.
- literal.json: 31ab27ebd540a2431b3871c1ad2459f6e714bbb82219ca05b8ac5928692a14d8.

No original source file was changed. The separately cited polynomial lift was not needed by this scalar proof or this review.
