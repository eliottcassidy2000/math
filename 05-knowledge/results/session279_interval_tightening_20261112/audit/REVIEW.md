# Independent new-claim review

**PASS WITH SMALL REPAIRS APPLIED.** This review covers the Session279 packet,
not a re-audit of the inherited repository. Chan–Pak remains OPEN. No global
integer-enumeration speedup is certified.

## Scope and repaired wording

1. The backbone API initially accepted repeated lengths, although its mathematics
   requires a set. Duplicate/nonodd/too-small lengths now raise an error; empty
   order, all-backbone and out-of-range-backbone controls are retained.
2. The physical procedure initially called itself an exact selector while
   allowing an unresolved short-depth ambiguity. It now explicitly says
   **certifying partial selector**, complete when its widths are below one;
   smaller depths need the stated finite-box or direct-bank fallback.
3. The marked-graph note no longer conflates the centralizer denominator of a
   recombined sector with every intermediate marked summand's denominator.
   Its fixed-depth polynomial bit-size conclusion survives.
4. Absolute coordinate widths oscillate: first crossing below one cannot justify
   all later depths. The physical proof uses four consecutive rational bounds
   and a rigorous four-step multiplier, not a monotonicity assertion.

These are domain/wording/continuation repairs. No numerical result was reversed
by this review. New inputs added during review (integer and divisor jets,
order64 marked reconstruction, general active-word spectral criterion) were
also reviewed on their final supplied source.

## A000568: backbone and moment compilers

The conditioned transfer is exact because its multiplier depends on the old
unmarked mass and the *fixed* backbone response, not on a forgotten unmarked
partition. First multiplicity moments therefore suffice for C1; second raw
moments and correction-multiplicity cross moments suffice for C2. The general
budget `ell+|alpha|<=K` is closed: adding h correction edges consumes at most h
multiplicity degree while lowering correction order by h.

The factorial/dyadic integer compilation preserves these equations exactly.
Its insertion factor is integral by labelled cycle placement; the final factor
is integral because the backbone centralizer divides its support factorial.
The primitive-binary-word identity gives the stated divisor response. The
cutoff d<=n/3 is valid only because queried lengths are distinct odd lengths;
the code treats equal-length interactions exactly before this quotient.

The marked-graph compiler independently exposes at most 2k incident cycle
instances for k edges, divides by m! after summing labelled markings, and sums
all unmarked cycles in their external field. Repeated lengths and unmarked
copies of marked lengths are allowed. The claimed fixed-B,fixed-K arithmetic
polynomial bounds are sound; they do not imply fixed depth suffices for
arbitrary integer accuracy.

The rational prototype is slower in the reported comparisons. The integer
first jet improves the matched C0/C1 workload at the larger tested sizes.
The divisor quotient has mixed timings. None is a benchmark victory over the
repository's optimized full exact counter.

## Ground moments and exactness band

The affine elementary-moment transform is exact, including its zero-coordinate
boundary. The listed four-root polynomials are exactly the feasible vertices
of the fixed-constant degree-four dual: upper roots come in adjacent pairs;
lower roots include 1 and r with an adjacent internal pair. The zero lower
bound covers nonpositive objectives. Their optimization is sharp for the
scalar nonnegative moment relaxation, not for every arithmetic constraint on
Burnside sectors.

The B={3,5,7,9,11} proof of rational exactness through69 is complete: its four-
and five-cycle case lists exhaust the possible residual cores. The core
(13,15,17,25) at70 has five nonzero shifted probabilities and a full-support
six-atom law. Its parity perturbation makes the degree-four relaxation genuinely
ambiguous. Singleton padding proves unbounded absolute ambiguity at fixed
backbone/depth within that relaxation, without claiming the perturbed laws
are alternative tournament arithmetic sectors.

The no-reference consumer at64 returns the same 1721-bit integer as the
independent marked-graph C0..C6 reconstruction. The stratified moment producer
still enumerates partitions; the graph producer is a separate fixed-parameter
construction. Their agreement is exact, but neither fact by itself is a
runtime advantage for full enumeration. The saved packet contains additional
minimum-defect and exceptional-rank information; the width gain is not from
identical information alone.

## Physical suffix recovery and spectral classification

The cone decomposition coefficients sum to the initial W, so the scale bound
N<=B a_max is correct. Both coordinate-width bounds follow from the exterior
minor/denominator identity. The three image rays have positive A after one
step. In particular, the formal ray with A=0 is a cone boundary, not a claimed
physical seed.

Independently confirmed:

- wedge matrix `[[3,3,4],[-3,1,0],[-4,0,0]]`;
- four-step infinity norm992;
- C3 growth bracket `[97/15,797/123]`;
- `992*(G/g^2)^4<1`;
- all four U_k(96), k25..28, are strictly below one.

Thus the all-k>=25 packet injectivity/unique recovery claim follows without
floating spectral estimates. The inverse retains divisibility by4, cone
membership, and actual seed ownership. It separates integer packets, not
rooted isomorphism classes. The fixed bank's O(log X) scalar supply cannot
prove all-odd coverage.

For an active word, the positive matrix gives Perron root lambda>2^length.
The other roots have positive product below one. The two endpoint signs at
+1 and -1 therefore characterize both transverse roots lying strictly inside
the unit circle. Equality is correctly excluded. The length8 hostile
10101010 has an expanding real transverse root. This classifies the stated
constant-word spectral mechanism, not arbitrary switching or every possible
decoder.

## Independent controls

Run from any working directory:

    python3 /absolute/path/to/session279_interval_tightening_20261112/audit/check_independent.py

The normal-Python run uses live assertions; an optimized `-O` run is **not**
claimed as a second assertion audit. `independent.json` records:

- 175 first-jet comparisons against a separately written occupancy sum;
- 63 second-jet comparisons, including non-prefix backbones;
- 99 marked-graph coefficient comparisons, including degree-three stars;
- 75 integer-jet and75 divisor-jet comparisons to that independent sum;
- all3002 four-contact polynomial vertices for supports5..14;
- 70 arbitrary positive moment mixtures unrelated to tournament partitions;
- independent exact rational uniform-threshold calculations;
- 44 literal physical packet checks by permutation enumeration, not subset DP;
- active-word characteristic coefficients independently calculated by minors,
  giving the same pass/fail counts through length12;
- independent residual-core search confirming first exceptional rank5 at70;
- agreement of the two order64 exact reconstruction representations.

The new all-order claims rest on the preceding proofs, not on extrapolation
from these checks. Final-source hashes are in `SHA256SUMS`.
