# Audit and reproducibility

Status: same-agent proof review with complete independent computational
algorithms for the principal atlas and eight-vertex boundary results.
No separate-agent audit, unrestricted spectrum theorem or priority claim.

## Exact universes

| Check | Universe / result | Independent evidence |
|---|---|---|
|Strong source charts|All 32768 spine charts at 7;30773 strong|Reachability versus all-forward interval cuts and SCC word recurrence|
|Seven-vertex matrices|Every entry of 30773 matrices|Subset DP versus full-cube literal permutation/ordered-cut expansion at 7,8,9|
|Eight-vertex boundary|All 2097152 spine charts;2032504 strong;zero 14-fringe hits|Subset DP versus full-cube permutation/cut expansion|
|65-bank saturation|8583 states at H<=1000|Independent inherited Python forward closure plus C++ parent/closure checker|
|960-bank saturation|32313 states at H<=1000|Every state has an actual predecessor; all eligible transitions close|
|67-bank saturation|30016 states at H<=2671|Same exact parent/closure certificate|
|Full through-seven bank|488230 states at H<=2671;497715 eligible distinct-matrix transitions|Independent checker reads the proposed set and checks every predecessor and transition|
|95/97 exact obstruction|74/84 full-bank states;none with both lifts|Complete fibres, including all word lengths|
|Compact actual repairs|All four full states and native programs|Literal H/end/cover subset counts and complete operation replay|
|Repeated offsets|All 30 edges through depth 4, endpoints at most 12 vertices|Every selected source independently counted;1112 native moves replayed;peak 13|
|Five Q4 phases|All five continuations of the first repair|Literal source graphs at both ends, full correction and carry retained|
|Actual clock twins|Two order-eight H293 sources|Independent literal states and q2 continuations 3557/3865|
|Additional transporter|All 31733 individual symbols|Exactly b,q2,s78 preserve tau up to sign; literal s78 seed control|
|Small local root probe|1232/1992 cases on95/97|Primary subset DP; bounded negative only. Independent permutation control on180 H13 cases|

## Why saturation is an all-order certificate at a value bound

All compiled matrices are nonnegative and dominate P_a entrywise. A word
of l letters therefore has at least h(a^l) paths. This bounds l at a fixed
H; it does not assume a finite cone inherited from a smaller alphabet.
H never decreases under extension, so no path leaving a value bound can
return. A listed predecessor proves soundness, and transition closure plus
the initial state proves completeness. Matrix duplicates can be collapsed
for reachability because their source representatives remain explicit;
this is not a quotient of automorphism groups or source multiplicities.

The eight-vertex negative result uses exactly the 14 reverse-fringe states
under 65 inverses. No C16 bound is used. It excludes only a single initial
processing atom of order 8 followed by 65 processing letters. It is not a
complete through-eight alphabet calculation and is not an obstruction for
arbitrary higher atoms.

The compiler and independent cube audit use different counting algorithms
and different strongness tests. All integers are bounded far below signed
64-bit limits: direct counts occur on at most ten vertices in the compiler,
and the saturation's matrices and small value bound give safe products.
The Python literal checks use unbounded integers.

## Native programs and guards

`native.py` implements the inherited one-auxiliary-vertex theorem. It searches
for a Hamilton path containing the chosen arc by contracting that arc to a
two-ended unit, then uses exact subset search. If the arc is untraced, it
constructs the reverse-direction gadget from the target. Both native path
guards, both full relations and the old restriction are checked. Padding
births, temporary births, native moves and named deletions have separate
clocks. The receipt retains the complete command list for every finite
branch; it is not an existence claim based on metagraph connectivity alone.

No monotonicity of intermediate H, fixed-order connectivity, strongness of
intermediates or minimum program length is claimed. The actual endpoints
are counted and strong. Names are transported in a fixed supplied gauge;
ordinary class-valued resetting is never an operation.

## Proof audit boundaries

- Three q2 probes are invertible over Q and modulo 7 because det O=96.
  This makes exact equivariance stricter than congruence equivariance.
- The four nonzero carries generate the full transverse integer plane;
  closing under q2 forces L7. Linear combinations need not be individual
  physical source differences. Smith minima are additive-observer minima.
- The growing-carry theorem uses rationality and conjugate quadratic
  eigenspaces, not a numerical eigenvalue fit. It forbids a finite fixed
  error set, not compact formula representations with an integer clock.
- The 42-state digit clock is an abstract residue controller. The current
  actual realization certificate has only four levels. The extra zero
  field state in the larger source language is verified on an actual graph.
- The source alphabet expanded beyond 65; its old C16 pruning guard is not
  inherited. Chart t0=(1,32,32) is an explicit hostile to that extension.
- The 81/27 joint-port hostile is inherited and remains a warning against
  reading full isomorphism or symmetry from three counting coordinates.

## Reproduction and sealing

`research.py` builds the five C++ programs in a dedicated temporary directory,
runs all complete atlases/certificates, independently checks selected literal
graphs, and compares the resulting compact receipt. `--write` regenerates it.
Normal and Python-optimized runs use explicit exceptions, and C++ checks
use exceptions rather than optimization-removable assertions. The primary
and audit C++ programs are algorithmically independent, not merely different
compiler optimization settings.

The receipt pins the inherited proofs, selected source states and generated
atlas/closed-set LF hashes. `seal.py` checks every local link and exact packet
file/hash set. The repository documentation checker is run before publication.
Only the packet and its navigation entry are intentionally changed.

Remaining work: an all-size closed realization family, its 3H return branch,
the complete seed interval and full prime/7p coverage. No A000568 extension,
LRC theorem, Collatz theorem or unrestricted metagraph minimum is claimed.
