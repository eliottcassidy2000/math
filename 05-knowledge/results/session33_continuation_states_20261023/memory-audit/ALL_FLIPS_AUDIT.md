# Independent audit of the all-labelled-flips addendum

Verdict: PASS. No mathematical correction requested. Only new audit files were created; earlier main-memory and composition audits were not modified. No Git mutations.

## Exact independent checks

2,134 explicit checks pass under normal and optimized Python with byte-identical receipts. The audit rebuilds actual tournaments and canonicalizes over every vertex permutation plus converse, not using producer node maps or helper functions. Directed triangles are found by an independent search over every ordered triple rather than adopting the producer's shortest-backward-edge algorithm.

Radius-d observation signatures use every XOR shift of Hamming weight at most d. They directly reconstruct continuation refinement without the producer's iterative label-refinement routine:

- n3: full flips2→2; complement only2→2.
- n4: full flips3→7→8→8; complement only3→5→5.
- n5: full flips10→64→64; complement only10→38→38.

Every unordered pair of distinct tilings is separated by the explicit anchor shift: 1+28+2016 pairs. The translated singleton indicators have exactly the distinct one-point supports, an exact permutation-basis certificate of full rank2,8,64. The converse-merged transitive fibre is singleton in each independent census. All other tilings contain a directed triangle.

## Proof review

The all-forward tiling is the unique transitive tournament compatible with the fixed directed path. A transitive tournament's topological order is forced by the path through all vertices. Converse merging cannot add another tiling because converse preserves transitivity and every transitive tournament is isomorphic to its converse.

The alternative shortest-backward-edge proof is correct: for backward j→i choose k=i+1. The path forces i→k, and k→j is either forced consecutive or shorter than the chosen backward edge and therefore forward. These three arcs give a directed cycle.

All labelled coordinate flips generate the translation group on the Boolean cube. For distinct x,y, translation by x XOR anchor sends x to the singleton fibre and y outside it. Thus equality of all continuation observations implies equality of source tilings. The proof needs pinned coordinate labels and all words; merely allowing an unlabelled choice of an arc is not this action.

The observable node-indicator space includes the singleton's delta function. Its group translates are every coordinate delta function. Therefore the closed linear observation space is all R^X. This proves the all-order linear dimension bound for all signed sources, independently of finite tests.

The positive lower bound is valid even for nonlinear probability encodings and without hidden-state intertwining. For each x≠y, an allowed common word produces distinct pure observed outputs. As the corresponding OT_word is column-stochastic, extremality forces the hidden supports of p_x and p_y to be disjoint. All2^m initial tilings need nonempty disjoint supports, so at least2^m hidden coordinates are necessary. Literal tilings attain equality. The requirement to represent every pure tiling and every labelled operation word is explicit and indispensable.

These are state distinguishability and linear/positive realization bounds, NOT static enumeration time or memory complexity bounds. An implicit algorithm or static node census is not required to expose an independently addressable continuation state for every tiling. The addendum correctly preserves that distinction, and the smaller complement-only quotient is not contradicted.

## Inheritance boundary

Reviewed the current THM-840 operation-kernel and linear-sidecar statements and THM-796's stated operation-relative scope, along with the previously reviewed boundary-support metagraph note. These supply the general criterion and source carrier, not the explicit full-labelled-flip singleton-anchor theorem. A bounded filename/text search in relevant theorem/result material did not locate that exact stronger statement. This is not a repository-wide or literature novelty certificate. The group-action singleton-anchor lemma is elementary standard reasoning and is correctly presented without historical-priority claims.

## Files and pins

New files: audit_all_flips.py, all_flips_audit_receipt.json, all_flips_audit_receipt.normal.json, all_flips_normal.out, all_flips_optimized.out, ALL_FLIPS_AUDIT.md, ALL_FLIPS_SHA256SUMS.

The receipt pins all three audited producer addendum files by SHA256. Reproduce with python3 memory-audit/audit_all_flips.py and python3 -O memory-audit/audit_all_flips.py. No producer functions are imported, and no previously frozen audit files are overwritten.
