# Independent review of the opposed-apex cycle packet

Verdict: **PASS after one corrected scope omission**, bound to:

- `session197-opposed-cycle/PROOF.md`: `7d4f2cc29697dcb59a2ad518a1765b36922ec909b880024d58f5901e75b27b7a`
- `session197-opposed-cycle/verify.py`: `b5a244dc80682e16ea21dfe52a362cab8af0b3e9f4fde5574127c25afc084eac`
- `session197-opposed-cycle/receipt.json`: `a96797631e747b042d5a09765a0ed36487c56f6f4ccb627fcb99c333bb1d5442`
- `session197-opposed-cycle/bounded_image.json`: `beb25ee9aa7504c98070a8680c8d545cf2f37c229faf249e592625860d4b46c0`

The explicit odd-subset cycle injection is correct, including the empty 1-slot and nonempty 2-slot. Odd subset size gives odd cycle length; different subsets give different vertex sets. The odd-subset difference is exact. The proof correctly uses only separate single-cycle terms of the inherited odd-cycle-collection identity, rather than incorrectly multiplying overlapping cycle weights.

The resulting strong-word length bound and Euler-product closure are sound. The disjoint-block forbidden-pattern estimate is an upper bound in the right direction. Length-one forbidden patterns are handled explicitly, including the all-transitive case. The Fibonacci exponent, exact comparison with 7/10, and explicit partition-function constant are correct.

The initial version of Section 3 omitted the condition lambda>=1. I requested that condition because a finite language can satisfy a growth estimate with lambda<1, whereas a negative scalar exponent cannot bound a product support containing one. The frozen final proof includes the correction.

I independently confirmed the asymptotic sharpness formula: deleting the unique 2-vertex gives the one-apex word `1 0^a 1^b 0`. The three endpoint categories contribute exactly the displayed odd-cycle counts, so H=(2^a+1)(2^b+1). At a=b=t, the ratio to 2^(m-1) tends to one. No claim of pointwise equality of the general lower bound is needed.

I reran the final verifier under normal Python and `python -O`; outputs agree. It checks 2,187 literal odd-cycle words on 120 binary input words, 682 forbidden-pattern inventories, 16 physical sharpness graphs, and the inherited 1,530 general slot/SCC controls. The bounded image gives 157 strong generator values and 414 product values through 5,001, using the independently justified length-13 cutoff.

The general theorem remains conditional only on the explicitly cited, established odd-cycle-collection identity. The rest of its all-order argument was reviewed directly; finite tests are implementation checks, not the source of the all-order conclusion. No unresolved correction remains.
