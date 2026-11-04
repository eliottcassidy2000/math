# Integrate the incoming disjoint reversal cube without duplicating it

Read directly at commit [`b6135263ba67ee2908df7625c67a88a3df8fbe03`](https://github.com/eliottcassidy2000/math/commit/b6135263ba67ee2908df7625c67a88a3df8fbe03): [PROOF, S1-S3](https://github.com/eliottcassidy2000/math/blob/b6135263ba67ee2908df7625c67a88a3df8fbe03/05-knowledge/results/simultaneous_module_reversal_20261104/PROOF.md) and its [review](https://github.com/eliottcassidy2000/math/blob/b6135263ba67ee2908df7625c67a88a3df8fbe03/05-knowledge/results/simultaneous_module_reversal_20261104/REVIEW.md).

## What is inherited from this checkpoint

S1 proves exact isolation of one disjoint block by giving both ports its quotient incidences outside. Every other block remains a module of the actual port extension, so its complete converse-invariant ordered-path-cover profile makes its orientation invisible to that query.

S2 chooses one such query per selected block with difference 2 modulo 4. These fixed coordinates decode the entire reversal vector and distinguish every pair of the 2^r states, not merely each state from the original one. This is the same disjoint mechanism proved in our concurrent packet, with its useful fixed-coordinate decoder stated explicitly. It is credited to the incoming S1-S2, rather than counted twice as a new mechanism.

S3 gives a genuine broader formulation: nuisance blocks may vary through any prescribed same-named family with identical complete q profiles. A differing distinguished block is separated whenever its core responses differ modulo four. Equality of H alone is not enough for exact isolation.

The incoming packet also retains the same covered-laminar C3[TT2,TT2,TT2] hostile as our checker. At (a,c)=(1,5), its quotient orientation pair gives 251 and 249; every block-uniform query agrees modulo four, although 24 of 64 differ as integers. Its complete sixteen-state orbit is already checked there and is not presented as another new census here.

## Interactions do not become additive

The incoming four-state interaction has the exact quartic

    [a_0 a_2 c_1 c_4](F_0-F_1-F_2+F_3)=-2.

`check.py` independently verifies it by 64 literal Hamilton counts in a sixteen-term finite difference. It also checks all 24 decoder values in the eight-state fixed-quotient cube: each is 393 or 483 according only to the selected bit.

Even smaller actual-query controls show nonadditivity modulo four without any quartic extraction. Start with transitive rows (0,1,3,7):

- Disjoint modules with masks 3 and 12, query (5,0): the four counts are 19,11,15,9, giving mixed difference 2.
- Nested modules with masks 3 and 7, query (9,0): the four counts are 15,9,5,9, giving mixed difference 10, again 2 modulo 4.

These are literal fixed-query interactions. They do not contradict exact isolation on specially lifted decoders, nor exact invariance of the averaged anchor battery. No proof here adds individual module-response differences as a linear map over F2.

## Distinct progress retained after integration

The incoming S1-S3 explicitly stop at pairwise-disjoint blocks. The current session's additional results remain:

- the exact singleton-port averaging lemma;
- the all-order odd/even anchor battery;
- exclusion of every proper nested pair;
- exclusion of exposed-anchor laminar families;
- companion crossing-pair and interval positivity arguments;
- companion arbitrary-original-module XOR reduction to scalar-preserving legal inversion classes.

Together the pair-specific arguments cover every pair of distinct nontrivial original modules. They do not yet prove that every fully covered laminar family or every arbitrary legal inversion sequence is visible. The full named cubic and quartic terms remain essential in that residual.

Normal and optimized integration receipts agree byte for byte. These checks verify selected incoming exact controls, not a fresh replay of its entire 40-state full-response packet.
