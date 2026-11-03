# Audit, inherited mechanisms, and search scope

**Status:** same-agent proof review; exact computations with two languages and
independent enumerations. All checks use explicit exceptions rather than
Python `assert`, so optimized execution retains them.

Normal and optimized replays both passed with certificate SHA256
`5f2b6b42a324057868459f7c0f14c157b2e9abc923f514983e4554c0f6542d7a`.
The independent C++ complete census and repository documentation checks passed.

## Literal branch census

Vertices and bit rows are named; bit j in row i means i->j. A birth mask marks
old vertices pointing to the newborn. The complement is the actual
Hamilton-path operation retaining consecutive path arcs and reversing all
others. The fixed universe is each supplied labelled order-12 source, zero
or one such complement, followed by one of all 4096 newborn masks.
Intermediate sources may be weak. No isomorphism or phase quotient is taken.

| Source H | legal paths | all +2 hits | strong +2 hits | all +4 hits | strong +4 hits |
|---|---|---|---|---|---|
|1457|1457|10|4|10|3|
|1459|1459|12|6|14|8|
|1463|1463|10|5|24|12|
|1465|1465|4|0|10|10|

The four weak +2 hits in the last row are retained in `census.json`.
Consequently two complements are necessary and sufficient for a **strong**
4397 target in the language "complements first, then exactly one birth" from
this particular 1465 source. This is not a lower bound if births and other
operations may be interleaved, nor a requirement for attaining scalar4397.

Python `search.py` chooses an enumeration priority from the entire named
root-defect vector, beginning with zero guards. Exact suffix-completion DP
only removes path prefixes that cannot finish; it does not change the legal
universe. It stops when both strong hits exist. For the last source it
exhausts all1466 programs. `reset.py` then orders first-complement sources by
H and finds the extra reset at outer rank1, after54 second paths. These are
search costs under specified orders, not optimal-complexity bounds.

C++ `lift_search.cpp` independently enumerates every source Hamilton path,
computes all newborn values using the inherited positive cut quadratic, and
reproduces the full census. Skipping intermediates with H above both targets
is safe because induced extension cannot decrease H. The bound is exact
int64 arithmetic: every evaluated tournament has at most14 vertices, and
all intermediate counts/products are bounded by factorial-scale enumeration
well below2^63. The Python lift audit uses arbitrary-precision integers.

## Six all-size program certificates

The C++ search keeps the common ordered predecessor/successor at an extreme
port, checks that birth preserves its sole incoming/outgoing neighbour, and
tests both m=1 and m=2. All six targets are strong. These finite checks invoke
L1/L2's proofs; they are not a numerical extrapolation. `audit.py` expands both
paths literally for m=1,...,4, verifies every arc and the full target rows,
then independently checks joint responses and run bounds.

Only positive witnesses in `lift_certificates.json` are promoted. Bounded
searches for the missing uniform +2 lifts are recorded in `search_bounds.json`;
their failure, if any, excludes only the stated ordered candidate prefix.
It does not negate the already constructed finite +2 branches.

All original vertices survive the programs. Clone addresses `(old_name,index)`
make the original root2 explicit. Full target relations are retained beside
the small response coefficients. The ancestral F343 state is a source-history
sidecar; the current F49 observation is recalculated from actual joint counts.
No unsupported Frobenius action on arbitrary complement states is asserted.

## Prior and incoming work used as operations

* `ready_return_pentagon_20261103/return_program.py` already supplied ordered
  Stirling run coefficients and all-five-Q4 checks. L3 uses that mechanism;
  the new ingredient is the k+1 bound and the six concrete return programs.
* `two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md` already
  exhibited a spectral mismatch between a 1,2,3-mode observer and a native
  tribonacci output. L4 explains its arithmetic boundary for every fixed
  transitive-module context; it does not rediscover that example.
* `session137_overnight_proof_frontiers_20261102/cycle_carriers/THEORY.md`
  explains the failed abstract conflict-graph route. K2 is natural but adding
  a universal cycle-vertex gives K3, which is not. We replace that non-realized
  operation with legal paths and retain the intermediate incidence relation.
* The October `native_return_necklaces_20261026/README.md` showed why an
  apparently reset readout still needs the original source and clock. The
  frozen F343 ancestry here serves that precise role; it is not the current
  clone clock.
* Incoming Session168, commits `59c2fc6fae` through `e92fc100d5`, provides
  equal-profile marked computations and controlled cellular clocks. Its
  `turing/NOTE.md` explicitly retains program, register positions and schedule.
  Our counterpart is the named path pair plus block address. Both demonstrate
  why equality of path-cover readouts does not make the scheduled operation
  autonomous. We do not import its single-arc-flip realizations as actual
  complement programs or claim its clocks equal ours.

The inherited twelve-vertex 81/27 hostile remains pinned to
`overnight_rank_return_20261029/DEFECT_TRANSFER.md`: equal two-run mass9011898
has ranks24/22 modulo3; equal three-run mass80614764 has ranks430/435. These
controls are not rerun or newly claimed. The H53 same-F49/different-return
control and the prior 161 dyadic/tribonacci near miss likewise remain live.
They prohibit replacing our full literal source by H, its current phase, or
the three golden entries alone.

## Reproduction

From this isolated worktree, compile the search into a temporary executable:

```text
g++ -O3 -std=c++20 05-knowledge/results/adaptive_four_returns_20261103/lift_search.cpp -o <temporary-executable>
<temporary-executable> 05-knowledge/results/adaptive_four_returns_20261103/source_1465.txt 0
python -B 05-knowledge/results/adaptive_four_returns_20261103/audit.py
python -B -O 05-knowledge/results/adaptive_four_returns_20261103/audit.py
```

Repeat the census with source1457,1459,1463. A positive-search invocation takes
`source_file exclusive_outer_limit [offset=0] [start_outer_rank=0]`. Original
successful prefixes end at ranks3,28 (1457);2,24 (1459);4 (1463/+4);8
(1465/+4). The source path universe and all masks are explicit in the code.
`python -B .../search.py` and `python -B .../reset.py` reproduce the finite
branch search with their deterministic path priorities.

The replay independently enumerates unpruned Hamilton paths for the eight
finite strong targets, verifies all stored rows and masks, checks weak census
targets, all five Q4 phase readouts and normal-basis coordinates, and checks
the exact observer periods through modulus343. General all-depth statements
use the proofs in LIFTING, not those finite modular checks.

No claim is made of a full spectrum proof, a new A000568 value, an all-size
strong-atom carrier classification, or minimality under arbitrary re-marking.
