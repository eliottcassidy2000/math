# Arbitrary shared-module correction

Proof candidate pending independent review. No eventual-spectrum threshold is claimed.

PROOF.md extends the audited A/B/C physical core rewrites from transitive modules to every tournament module U, with H(empty)=1:

    H(child(U)) = 2 H(parent(U)) + H(U).

The same proof gives exact cycle counts, including B/C replacement from H(U) Hamilton cycles to 3H(U) at unchanged vertex order. The correction and the bin-cover coordinates are actual coupled graph data, not free integers.

The optimized and debug C++ verifiers independently check all 76 labelled input modules through order four, 684 literal output graphs, and 228 rewrite pairs. The general proof uses module path-run decompositions and the independently audited finite core skeletons. Frozen dependencies are included.
