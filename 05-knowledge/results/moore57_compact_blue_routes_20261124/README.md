# Condensing Moore observer tournaments and finding their native lines

**2026-11-24.** This follows the [original 10/20/100 observer construction](../moore57_tournament_observers_20261124/README.md). A degree-57 Moore graph remains open; all its sizes below are conditional.

| Decoding contract | Pentagon | Petersen | Hoffman–Singleton | Degree 57 |
|---|---:|---:|---:|---:|
| Original separate full observer | 10 | 20 | 100 | 6500 |
| Exact same binary Gram, sharing sources/probes | 5 | 13 | 93 | 6493 |
| Compact ordinary tournament, unmarked | **5** | **10** | **50** | **3264** |
| Compact ordinary blue tournament, unmarked | **5** | **10** | **50** | **≤6097** |

- [Short arrow rules and decoders](COMPACT.md): the finite-field rule uses five-site cyclic blocks and a single affine exception. For a general unknown degree-57 graph, binary beacons and two twin markers recover all labels for only 14 extra sites. The separate [strong blue bound](moore_observer_compression.md) deletes observers and uses the Moore matrix equation to recover the missing block. The [strict Gram rank/overlap bound](moore-gram-overlap.md) explains why simply trying to fuse identical types of observers cannot cut the original 6500 in half.
- [Native metagraph paths](NATIVE-PATHS.md): the first line from transitive is unique and blue. From there we certified shortest all-blue paths of lengths 3, 6 and 7 respectively to the regular pentagon, the original doubled pentagon, and the compact Petersen tournament. Their classes are identified by minimum masks over all Hamilton roots. At 50 vertices we verified a partial 91-move route out of the compact HS (not a complete route to transitive).
- [Blue incidence and hidden counts](moore-observer-blue.md): a blue path splits over half as many labels by reflection. The eight original doubled pentagon classes have 980 physical blue incidences to 857 external classes and no internal line of either color. For the original Petersen observer one specified blue line changes the exact Hamilton count, and its open vs closed blue roots are counted independently.

A change of decoding contract is not a native line: 100 and 50 vertices lie in different fixed-size metagraphs. Likewise a chart preserving minimum edit distance and a physical line obtained by complementing every free tile are not the same adjacency. The first-neighbor score slice gives explicit equal-height, equal-score classes that take different native paths.

## Reproduce

Run compact.py, gram-small.py, native-replay.py, hs-prefix-check.py, and moore-observer-blue.py with Python 3. The partial observer decoder uses NumPy/SciPy in moore_observer_compression.py. Compile moore-gram-overlap8.cpp and moore-observer-height.cpp with C++17 and optimization for the two finite exact certificates. No script purports to construct a degree-57 Moore graph.
