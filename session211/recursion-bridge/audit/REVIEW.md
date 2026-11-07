# Independent recursion-bridge review

**Verdict: PASS.** Final proof SHA-256: `987706d4116dc6df55a7332fe65d3c4621ad483d42c40edf41da03cadf69116b`.

## Mathematical review

- RB1 correctly distinguishes order `n->n+2` from any affine path-count law. For a nonempty core, the changed apex arc is absent from every padded Hamiltonian word, so none is destroyed. The new words correspond bijectively to ordered two-path covers. The empty-core exception is explicit.
- RB2's sharp minimum and nontransitive second minimum are proved with complete equality conditions. The triangle-resolved lower bound follows by summing inherited nonnegative odd-cycle contributions over the two shores. The exact gated-class values at most 31 are `{3,5,9,15,17,27}`; this is not an unrestricted-spectrum classification.
- RB3's all-pairs feedback energy identity separates two within-sheet terms from the constant cross-sheet contribution. The reflected image is an affine diagonal section of the ordered-sheet product face, with Euclidean scale `sqrt(2)`, not generally a face. Its edges are diagonals of genuine commuting square faces.
- RB4 uses endpoint counts for the first-sheet/second-sheet sector and starting counts for the reverse sector. Both formulas and the lower-bound inequalities are valid. The uncounted remainder consists of paths with at least two sheet changes and is not assumed to vanish.
- RB5 distinguishes equal object orders from isomorphic outputs. A sink detects failure of padding/doubling commutation; an outdegree-one gate tail detects failure of gate/doubling commutation. The source remains a nonempty tournament throughout.

The low-count OCF facts, old padding/doubling laws, and global omissions 7,21 are inherited with exact source routes. No novelty or unrestricted count-surjectivity claim is supported or made. All entries of `SOURCES.json` were independently compared with the immutable reference Git snapshot and match their SHA-256 values.

## Reproduction

The final verifier was copied to an independent directory and run normally and under `python3 -O -B`. Both stdout files match frozen `RESULT.json` byte for byte:

- `verify.py`: `3be195d176f5ab491ffe883159831825a881f30306bc236070d6ae56bc5ac9b6`
- `RESULT.json`: `3ea466d9c12d018d8a522747630fde42fb94c9c083feac21e81ba1c8eb21d334`

The producer's exact universe contains all 1,099 labeled tournament cores of orders 1..5, all 124,469 reflected permutation words, and direct two-block enumeration for all 75 labeled cores of orders 1..4. Both noncommutation degree obstructions and the sharp gate bounds are checked.

## Independent geometric route

`independent_folded_squares.py` imports no producer functions. It constructs adjacency matrices directly and verifies all four corners and all four energy differences for each commuting square. Its universe is all 75 labeled tournament cores of orders 1..4, all 1,589 permutation words, and all 4,708 adjacent-interchange instances. It also verifies the paired-position equations for the affine diagonal section. Normal and optimized JSON outputs agree exactly.

The independent script, both outputs, and final copied-verifier replay outputs are approved for inclusion with this review. Finite checks are controls; the all-order claims rest on the reviewed proofs.
