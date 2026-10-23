# Sharp five-arc unindexed window threshold in the affine doubled-Paley family

Status: PROVED for the explicit orbit of INFINITE_ARC_FLAG_GAP.md.

For every prime q=3 mod4 with q>=7, the minimum length of a contiguous
unindexed path window that identifies the displayed P within its affine orbit
is exactly five arcs. This differs from two arbitrary unindexed arcs and one
indexed flag.

Proof. Coordinate-changing arcs occur every four steps. Any window of at most
four arcs therefore contains at most one such arc. If it contains neither such
arc nor a center-to-copy/copy-to-center pin, it appears in all orbit paths.
If it contains exactly one coordinate-changing arc and no pin, its admissible
maps are the q-1 translations described in the incidence proof; accompanying
same-coordinate steps add no further restriction. If it contains a pin but no
coordinate-changing arc, it leaves (q-1)/2 maps. The first pin and first
coordinate-changing arc are four positions apart, as are the last changing
arc and final copy-to-center pin, so a window of four arcs cannot contain both.
Thus every such window has multiple survivors.

The initial five-arc window contains both the first pin and first changing
arc. The pin imposes b=0 and the change imposes a=1; hence it identifies P.
This proves both bounds. At q=3 a pin alone suffices, so that boundary is
intentionally excluded from the five-arc theorem.

The standalone checker exhausts every window of lengths1 through5 at
q=3,7,11,19,23,31, constructs actual source tournaments and verifies all affine
automorphisms. Normal and optimized receipts match. Nothing asserts that the
window identifies P among Hamiltonian paths outside this declared orbit.
