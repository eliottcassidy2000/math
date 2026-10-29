# Delayed activation and half-deficits in exact A000568 counting

**PROVED implementation equivalence / FINITE-EXACT measured costs.** The
Burnside formula and response quotients are inherited. This is a small
implementation improvement, not a new unrestricted combinatorial formula,
not class generation, and not a polynomial-time complexity claim.

**Concurrent precedence:** before publication, `f8846701cd` supplied the
equivalent [forced-zero agenda contraction](../metagraph_scheduled_recursion_20261029/SCHEDULE.md),
with exactly the same88915 saved transfers at200. Our parked-bucket layout
corroborates that mechanism; it is not a distinct new quotient. The useful
extension of that incoming work is now
[completion-sensitive precision pruning](PRECISION_EXTENSION.md).

## Inheritance and exact state

The comparison uses the Session87 degree-bucket engine as pinned by
Session105 and Session105's
[integer projective response](../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md).
The exact snapshots are `baseline87.py`, `baseline105.py`; original paths,
commit and hashes are in `source_pins.json`. This retains the stronger
full-prefix engine rather than treating the newest single-target engine
as an automatic speed improvement.

At a descending odd-cycle stage p, the prefix engine stores processed
degree d, integer-scaled Burnside weight, and future gcd responses
h_j=sum_a gcd(a,j). The projective single-target engine at total order N
stores remaining capacity r=N−d and deficits

    D_j=sum_a (j−gcd(a,j))=j t−h_j,

with t the number of processed cycles, and inherited weight
U=odd(d!) sum 2^(e+r t)/z. All odd cycle lengths, including composite
ones, are retained. It is not a prime-only automorphism approximation.

## Two exact changes

**Delayed activation.** At stage p, a bucket with remaining capacity r<p
cannot accept a positive number of p-cycles. Its zero-multiplicity operation
is the identity, including its response coordinates. To see the latter,
the previous cap was min(p,r) and the next cap is min(p−2,r). If r<p−1,
they coincide; if r=p−1, they retain the same odd indices because p−1 is
even. Thus keep these dictionaries in place until a smaller p becomes
feasible. Only active buckets are traversed or copied.

Positive insertions from an active bucket may enter a parked bucket. This
is safe: the old active source dictionaries are separate from their new
targets, and parked dictionaries are not read by the active loop in that
stage. New contributions therefore cannot reuse the same cycle-length
stage. This is the mutation boundary checked by the implementation.

The prefix engine uses the equivalent condition d>N−p. The projective
engine uses r<p. There is no state identification or discarded arithmetic;
the saved work is repeated identity processing of inaccessible buckets.

**Half-deficit packing.** All cycle lengths and response indices are odd,
so every j−gcd(a,j) is even. Pack D_j/2 instead of D_j, and multiply the
extracted digit by two in the transition exponent. Descending order gives
j−gcd(a,j)<=a−1, hence D_j<=d<=N. Therefore D_j/2<=floor(N/2), and
`bit_length(N//2)` bits per digit suffice (one bit for the zero case).
No carry between packed digits is possible. The response quotient and
integer weight remain identical.

The inherited exact-shift guard checks every negative projective binary
exponent. Final divisions by odd factorials remain checked. The prefix
engine also supports any positive integer number of vertex colours and
returns the whole coefficient list from 0 through N.

## Matched-output measurements and controls

`check_counts.py --bench` compares seven methods in alternating order over
three repetitions. Median seconds on this machine:

| N | inherited full prefix | parked full prefix | inherited projective | parked half-deficit projective |
|---|---:|---:|---:|---:|
| 100 | 0.02831 | 0.02633 | 0.02364 | 0.02302 |
| 200 | 1.37623 | 1.30589 | 1.97627 | 1.97098 |

At N=200 the parked prefix traverses 546833 states and skips 88915 identity
visits. The projective variant traverses 495459, skipping 83325; its packed
digit width falls from eight to seven. Despite that reduction, its measured
runtime is essentially unchanged, and it remains slower than the full-prefix
engine at N=200. Smaller keys are not by themselves a demonstrated speedup.
The modest measured prefix gain is about 5.1%; it is not an asymptotic claim.

Controls cover every N=0..30, both parking modes, both projective packing
modes, and vertex-colour counts 1,2,3. Direct odd-partition Burnside counts
give an independent small uncoloured control. A separate grouped-core kernel
checks all three coloured prefixes through30. Larger matched outputs at100
and200 are authenticated by exact equality and a decimal-value SHA-256;
no timing tolerance substitutes for equality. Receipts retain all individual
timings, not only the medians. Optimized Python runs keep explicit checks.

## What this does and does not change

The all-order formula is already exact; its large response-state cost is
the unresolved performance target. This change avoids processing states
that the current operation cannot touch. It does not bound the number of
distinct active profiles uniformly, generate the metagraph, or identify
classes from Hamiltonian counts.

The comparison with the seam observer is precise: the counting engine may
discard an expired coordinate only after **every possible continuation**
stops reading it. A seam-equivalent pair cannot be discarded before a native
complement, because the five-vertex hostile in SEAM gives a continuation
that distinguishes it. Likewise the Collatz decoder must retain a marked
source even when it can recover a longer suffix from the endpoint.
