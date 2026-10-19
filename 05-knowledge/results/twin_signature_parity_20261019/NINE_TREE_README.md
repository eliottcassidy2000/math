# Complete nine-tree parity classification and pendant-pair cancellation

Read [the nine-vertex theorem](NINE_VERTEX_GOOD_TREES.md), [two-pendant-pair theorem](TWO_CHERRY_PARITY.md), and [independent review](review/CHERRY_EXTENSION_REVIEW.md). These extend the [frozen twin-signature core](TWIN_SIGNATURE_PARITY.md).

Exactly 36 of the 47 nine-vertex tree types are good; 12 are not hereditarily symmetric. The proof gives an explicit structural criterion and critical-leaf orientation selector. Every nine-tree has affine embedding parity as a function of orientation bits, including the nongood trees whose coefficients may depend on the tournament. All eleven excluded classes have independently recomputed pairs of hosts with different parity. The all-order good-tree classification remains open.

Two disjoint pendant leaf-pairs force even embedding counts into every tournament, even when the orientation has no nontrivial automorphism. A retained isolated vertex is essential to the deletion argument.

Run from this directory, beside the unchanged core utilities:

    python3 -B verify_two_cherry.py
    python3 -B -O verify_two_cherry.py
    python3 -B verify_nine_vertex.py
    python3 -B -O verify_nine_vertex.py
    python3 -B review/independent_cherry_extensions.py
    python3 -B -O review/independent_cherry_extensions.py
    python3 -B review/independent_nine_criterion.py
    python3 -B -O review/independent_nine_criterion.py
    python3 -B review/independent_nine_witnesses.py
    python3 -B -O review/independent_nine_witnesses.py

Author and independent normal/optimized receipts match exactly. The 1,169-check independent cherry audit, independently regenerated 47-tree census, and all eleven 9!-permutation witness checks retain their separate finite scopes. The universal assertions rest on the proofs and their explicitly inherited inputs.
