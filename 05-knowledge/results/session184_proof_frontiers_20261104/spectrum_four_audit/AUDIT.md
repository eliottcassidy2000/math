# Independent audit: SCC-order-at-most-four return obstruction

Verdict: **PASS**. No proof or certificate defect found. This is a scoped all-length, arbitrary-fixed-clock-shift obstruction, not a global Hamilton-path spectrum obstruction.

## Proof review

1. SCCs of tournaments form an ordered sequence. No strong order-two component exists; C3 exhausts strong order three, and the producer explicitly classifies all 64 labelled order-four tournaments, finding the 24 strong ones all isomorphic to Q. The bank represents every actual port subset: cardinality suffices for C3 by rotations, while all 16 Q subsets are retained.
2. The Session182 ordered-sum and attachment formulas were inspected. They apply to arbitrary suffix tournaments, so importing them for Q components does not add an unproved strongness restriction. Terminal-component marking is precisely the necessary and sufficient full-attachment strongness condition.
3. Removing a maximal initial unmarked prefix yields c=h_P*c_D and F=p_P*F_D. Since each digit target C is nonzero modulo 3, no unmarked C3 can occur in that prefix. Each unmarked Q contributes h=5 and p=32; each unmarked singleton contributes h=1 and p=2. Thus all Q-prefix counts 0 through v5(|C|) must be scanned. Their order and arbitrarily many unmarked singletons cannot change the normalized c or odd part of F.
4. All marked bank entries have v2(F)<v2(p), with v2(F)<=4. The odd h of every bank block and strictly more divisible p*F term preserve this valuation after any append. This establishes the finite five-scale normalization without a word-length assumption.
5. Search pruning is complete: F is nondecreasing; any successful longer word must eventually append a marked final block with F_Z>=1, hence final F>=current p. Products p grow by factors at least two, giving a proved finite recursion. Testing a final block before refusing further recursion correctly retains the possible one-extra-block endpoint.
6. Every normalized hit lifts to a target under an integer clock shift; sufficiently large clocks make negative shifts legal. Accepted transitions have no unmarked Q prefix, so extra initial unmarked singletons only change dyadic scale, leaving the normalized successor search invariant at the next step. The seven-state successor closure therefore proves the full address-language statement.
7. In particular C=-35 must include r=1, giving normalized constant -7, in addition to r=0. Both cases, each at all five scales, are present and rejected. Ignoring that case would have left a real coverage hole, but this implementation does not ignore it.
8. The conclusion that a missing sibling +4 return requires some suffix SCC of order at least five follows within the same architecture and clock-uniform definition only. The text properly excludes one-clock coincidences, other architectures, and claims of native operation legality.

## Replays and independent implementation

Producer verify.py was rerun normally and under -O. Both outputs byte-match its stored normal.out. It reports 24 labelled strong order-four tournaments, 1794 literal mixed-word signature controls, 1097 normalized valuation controls, 75 complete target scans, and 2445362 visited prefixes.

The independent verify_independent.py does not import producer code. It reconstructs all 22 component signatures by explicit path permutations and ordered two-slot covers, rather than subset DP or the displayed table. It performs one joint traversal for all 75 distinct target pairs, up to the mathematically derived largest F=27216, visiting 588218 prefixes and examining 285 terminal candidates at requested F values. All fourteen transition hit lists, including literal words and full signatures, agree exactly with receipt.json. All six existing edges and all eight missing edges are confirmed. The independent normal and optimized outputs byte-match.

The independent traversal still uses the proved seven-coordinate assembly law; its independent component derivation and differently organized exhaustive census corroborate implementation, while the proof review above establishes unbounded coverage. No timing or speedup claim is inferred from these counts.

No source modifications, installations, forbidden-file reads, or external writes were needed.
