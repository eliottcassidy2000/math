# Independent audit: the order-eight skew-conference characterization

**PASS.** The unique root normalization is an all-order matrix argument, and its specialized order-eight enumeration and complete set comparison prove exactly the claimed finite equivalence. No general-order equivalence or full-response collision is inferred.

For S'=[[0,1^t],[-1,C]], direct multiplication gives S'S'^t=7I precisely when C1=0 and CC^t=7I-J. Fixing the root diagonal switch to +1 makes the source normalization unique. Conversely the 128 switches are distinct, and different normalized blocks cannot overlap. This validates the completeness step rather than treating switch sampling as enumeration.

The independent checker imports no source code. It enumerates regular seven-name tournaments by degree-constrained row choices rather than scanning the producer's full edge-code interval. Every regular tournament appears once; exactly 2,640 are found. The independent common-outneighbor criterion leaves 240 blocks, whose full sorted code list agrees with the source. Every bank entry is checked directly by Seidel-row Hamming distances, independently normalized, and assigned to one of these 240 blocks with multiplicity 128.

Forward switching in a separate row representation reconstructs all 30,720 codes and compares the entire integer set with the reviewed h-zero bank. The audit also directly recomputes all shared-cubic columns and Hamilton endpoint parities on all 30,720 graphs using subset endpoint recurrence, without the source's seven-name lookup tables: every shared-cubic matrix is zero and every endpoint-XOR mask is 255. Completeness of the h-zero bank itself remains the explicitly cited independently audited dependency.

All 40320 relabelings of each of the four representatives are generated independently. The four disjoint orbit sizes are 1920, 13440, 13440, 1920 and cover the complete class. The middle two representatives are prime, the outer two have a universal vertex, and whole converse pairs the stated orbits. The TT4 boundary has zero shared cubics but row inner product 2, correctly preventing an unrestricted-order inference.

Independent normal and optimized audit outputs agree. Fresh normal and optimized source drivers also reproduce their frozen receipts exactly, including the driver's separate O2/O0 C++ enumeration check. Every audit condition uses explicit exceptions or a nonzero executable exit.

Reproduce with `python3 verify.py > normal.json`, `python3 -O verify.py > optimized.json`, and compare. An already installed C++17 compiler is required for the independent endpoint recurrence; no installation is performed. Frozen audit allowlist: `REVIEW.md`, `verify.py`, `endpoint.cpp`, `normal.json`, `optimized.json`, `SHA256SUMS`. Redundant producer replays are excluded. Full normalized codes, orbit data and source/dependency pins are retained in the receipts.
