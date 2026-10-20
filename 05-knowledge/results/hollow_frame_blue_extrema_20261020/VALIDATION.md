# Validation record

## Complete finite universe

For each `1<=k<=5`, all `2^(k(k-1))` labelled external hollow frames were
enumerated without filtering. This is not a sampled frame collection or a
full tournament census.

## Independent algorithms

The C++ producer enumerates potential signed Hamiltonian prefixes and adds
one to each supporting external-frame subcube. The independent C++ replay
reconstructs each frame and builds its path count by forward dynamic
programming over used pair subsets. Every array entry agrees, including all
1,048,576 order-eleven frames.

Python independently checks all 2,850 inherited classes at orders 3–9 and
the order-eleven witness. A separate direct enumeration of all 1,920
potential order-eleven prefix words produces 195 paths. Every full word is
checked against the supplied 11-vertex adjacency and anti-involution.

## Symmetry and quotient checks

All 120 pair-index permutations of the maximizing order-eleven witness
produce precisely the forty maximizing external masks from the complete
census. The frame group has order three. A one-bit matching completion has
trivial stabilizer and therefore blue degree 195.

All 32 matching vectors were quotiented explicitly. The maximizing hollow
type contains eight classes of blue degree 195 and eight of blue degree 65.
The eight maximum classes split 4–4 by triangle parity.

## Independent external audit

The separate audit worker reports agreement on the complete binary
histograms and moments, the direct word count, the frame group and
stabilizer, the exact maximum hollow orbit, and the matching quotient.
The labelled mean and Gluck regular-subset reduction were independently
checked, including the primary 1983 source.

The elementary pair-block upper bounds and the explicit cyclic three-block
family also passed independent proof review, including the rank-class proof
that its frame group is exactly C3.

## Isolated reproduction

An isolated temporary copy with no sibling census ran `verify_extrema.py`
and `build_witness.py` successfully. It explicitly reported zero inherited
class checks while retaining all core results, including maximum 195. In
the full staged workspace, the optional inherited checks still cover all
2,850 classes. Publication links and reproduction commands are relative to
the packet's eventual published location.

## Boundaries

No exhaustive search at k=6 or beyond was attempted. No all-k variance
formula is claimed in this packet. The auxiliary Q3 is not asserted to
preserve actual blue adjacency. The general extremal growth problem remains
open; no literature-priority claim is made.
