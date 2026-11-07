# Square-sum count and port bridge

Read [PROOF.md](PROOF.md) for source/status distinctions, the two-port insertion theorem, the obstruction to mapping square-sum chain counts to the tournament odd index, and inherited Pythagorean scaling.

- [RESULT.json](RESULT.json): all counts through order 25; the 14 paths at onsets 15,23,25; endpoint pairs; two-path predecessor forests.
- [MACRO_RESULT.json](MACRO_RESULT.json): 22 restricted three-square path constructions from all 11 primitive triples with smaller leg s and larger leg t<=60. This is not a full path census of the large square-sum graphs.
- [SOURCES.json](SOURCES.json): repository source paths pinned to main commit 214586af6c8803d8b21959951e8f151e64bff36f; saved source copies match the recorded Git blobs exactly.

Run from this directory:

    python3 -B verify_square_bridges.py
    python3 -O -B verify_square_bridges.py
    python3 -B verify_pythagorean_scaling.py
    python3 -O -B verify_pythagorean_scaling.py

All checks use standard-library Python and explicit exceptions, so optimized execution retains the tests. Normal and optimized outputs agree. The count and macro control claims are finite; global tournament-spectrum completeness remains open.
