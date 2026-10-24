# Validation contract

The universe is all almost-regular ten-vertex tournaments and their unique
regular-eleven extensions, not an unrestricted order-ten or order-eleven atlas.
The complete source nine-parent atlas and all its full actions are inherited
from a byte-pinned independently audited packet. The small five/seven input
atlases are reconciled against the archived certified JSON.

Three complementary completeness mechanisms are retained:

1. The exact score filter and full source-group subset orbits generate every
   pointed almost-regular parent, with no duplicate pointed types.
2. The stable-cell cyclic normal form is a proved complete isomorphism invariant
   of this universe, retaining every crossing bit and every minimizing witness.
3. Labeled masses from source extensions, deduplicated parent groups, regular
   component groups, and an independent degree-sequence DP all agree.

All source rows are checked to be tournaments; all action lists are checked for
permutation validity, every arc, identity, distinctness, and closure. Full-group
completeness follows from the normal-form witness proof, not closure alone.
The lower order7 and order9 controls reproduce the known regular corrections;
order9 also matches the inherited full orbit-multiset histogram.

All three orders are reproduced under O2 and O3/NDEBUG, with byte-for-byte
matching JSON, parent atlases, state TSVs, and component TSVs. Guards are runtime
exceptions rather than disabled assertions. Normal and Python -O degree DP
results match exactly. Six malformed source inputs are rejected in both modes:
truncation, a trailing token, a self loop, a nonpermutation action, an orientation
reversing transposition, and a duplicate representative.

`verify_packet.py` checks every manifest and source hash; all D rows/components;
parent arc actions; stabilizer-derived full child group orders; all converse
pairings and D intertwining; exact degree-DP masses; lower controls; and the
combined correction histograms. It runs identically under Python -O.

`reproduce.py` accepts an external clean output directory, recompiles both modes,
and checks all output artifacts against the frozen O2 baseline. No run depends
on the original scratch path. Binaries, Python caches, and verification receipts
are outside the payload manifest. Source inputs and complete finite outputs are
inside it.

An independent audit packet may additionally canonicalize ten-parent objects
using unrestricted stable-cell permutations and search eleven-vertex full
automorphism groups. Such audit checks do not alter the production boundary:
there is no target-eleven canonicalizer in either the production or intended
audit, and no unrestricted ten-atlas is used. The audit remains separately
versioned; its receipt, rather than this statement, determines audit completion.
