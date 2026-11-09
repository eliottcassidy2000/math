# Reproduce Session231

Run from a checkout containing the packet and inherited `04-computation`
inputs. Python 3, a C/C++17 compiler, NumPy and SymPy are required. Standard
library Python suffices for the positive interval certificates. No network,
credentials or new package installation is part of these commands.

## Compact theorem and audit receipts

From this packet directory:

    python3 -I -B replay.py > replay.json

The driver copies the packet to a temporary directory, compiles the small
source inverse there, runs 33 producer/audit checks in normal and optimized
Python, compares frozen receipts, and leaves tracked files untouched. Its
frozen root replay is `REPLAY.json`. This includes the 35,567-block ribbon
spectral audit. It does not rerun the full order-nine discovery census.

For just the main constructive theorem:

    python3 -I -B fibers/verify_short12_unit.py
    python3 -I -B fibers/verify_short12_independent.py
    python3 -I -B fiber_audit/short12_audit.py
    python3 -I -B rerooted/verify.py

Repeat with `-O`; explicit gates remain active. Some standalone audit scripts
write their receipt beside themselves. Use the temporary-copy driver if
preserving a clean checkout is important.

## Full independent physical controls

    python3 -I -B replay.py --full-physical > replay-full.json

This additionally compiles `rerooted_audit/physical.cpp` and independently
counts both 25-vertex original witnesses and the 9/21-vertex bare witnesses,
including every named root. It uses about 3.7 GB RAM. The C++ implementation
checks endpoint counts and uses 128-bit intermediates to reject overflow of
its retained 64-bit totals. Its independent frozen receipt is
`rerooted_audit/CHECK.json`, with a matching optimized-Python receipt.

## Complete fixed-core context language through 10^8

    python3 -I -B replay_coverage.py > replay-coverage.json

This regenerates producer and independent unions in a temporary directory,
runs both implementations and their controls normally and under `-O`, verifies
all retained witnesses, tests deliberate-failure gates, and compares the two
union bitsets byte-for-byte. `COVERAGE_REPLAY.json` is the frozen root replay.
The 6,250,000-byte bitsets are intentionally not tracked. The universe is
arbitrary normalized singleton contexts around the specified fixed core,
not arbitrary tournaments. All moving-core contexts are already absorbed
in that same universe.

## Large discovery and completeness censuses

See `fibers/README.md` for classgen, erased-interface and saturation-bank
regeneration. `classgen.c` is copied unchanged from Session230, SHA-256
cb78b6418a7754d656f628eede3cb3628a5100ad93af14d947b25e48223f562e.
Generate the class lists with prefix `classes9`, then independently certify
and replay them using:

    g++ -O3 -std=c++17 fiber_audit/orbit_mass.cpp -o /tmp/s231-mass
    g++ -O3 -std=c++17 fiber_audit/replay_census.cpp -o /tmp/s231-census
    /tmp/s231-mass 9 fibers/classes9_9.txt
    /tmp/s231-census 9 fibers/classes9_9.txt 64

The retained independent outputs are MASS9.txt and REPLAY9.txt. The audit
also checks pairwise nonisomorphism, not only orbit mass. Lower orders use
the corresponding class list and order argument. Once class lists exist,
`fiber_audit/audit.py`, `physical.py` and `weighted.py` reproduce the earlier
reachability, graph and weighted-SCC audits. Their shared BASE path resolves
to sibling `fibers`; no original scratch directory is required.

For the compact filtered minimality certificate:

    python3 -I -B fibers/check_shorter_loops.py
    python3 -I -B fiber_audit/minimal_walks.py

The independent seed rebuilder is run after generating all class lists:

    g++ -O3 -std=c++17 fiber_audit/minimal_seeds.cpp -o /tmp/s231-seeds
    /tmp/s231-seeds fibers > seeds.txt 2> seeds.stats
    cmp seeds.txt fiber_audit/MINIMAL_SEEDS.txt

The proof and full regeneration lineage are in
`fiber_audit/MINIMALITY_PLACEMENT_AUDIT.md`. Seed-only replay does not itself
prove the census complete; the independent class audit and seed regeneration
do. The all-order positive interval theorem needs none of this discovery
completeness machinery.

## Integrity and corrections

`SHA256SUMS` covers final packet files, excluding itself and generated binaries
or bytecode. Run `sha256sum -c SHA256SUMS` from this directory. Frozen manifests
named INPUT or AUDITED_INPUTS describe the submitted prepublication inputs;
editorial status changes and repository-relative path repairs can change those
files after the original audit. The final manifest covers the published copies.

Read MISTAKE-741 and the individual audit reports for empty-context typing,
phantom diagnostic entries, saturation, endpoint changes, k=0 qualifications,
optimized-gate repairs and finite/all-order boundaries. The later smaller
scalar family is not silently substituted for the common-full-interface family.
