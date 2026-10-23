# Independent audit: actual shared square-clock ports

23 October 2026. Reviewer: lattice-attachments agent, independently of the
author of [SHARED_CLOCKS.md](SHARED_CLOCKS.md) and its compiler.

**PASS, scoped all-size proof and finite implementation controls.** The
result classifies supplied ordered clock paths with actual indexed port
equalities. It does not discover a port identification from an unindexed
graph, identify arbitrary tournament wiring, or silently force every pair
of vertices in different paths to be distinct.

## Local legality and signed propagation

The recurrence `c_0=0, c_i=s_i-c_(i-1)` gives exactly
`x_i=(-1)^i*t+c_i`. Pulling back the host inequalities gives the stated
integer interval. Equal-parity repeated constants force a permanent
collision. Opposite parities give coefficient two or minus two and therefore
at most one forbidden integral start. Deduplicating collision holes is
necessary for the count and is correctly implemented.

The shared-port equation has the correct sign and offset:

    t_b=(-1)^(i+j)*t_a + (-1)^j*(c_(a,i)-c_(b,j)).

Reversing it gives offset `-eta*q`. Tree propagation preserves one actual
root coordinate with coefficient plus or minus one. Every original equation
is then tested, including parallel edges and loops. A zero coefficient
tests an exact translation; a nonzero coefficient pins twice the root and
must pass integer divisibility. Conflicting pins, nonintegral pins, and
pins outside the legal root domain are rejected. Absolute indexed marks
have coefficient plus or minus one and therefore pin an integer root
directly.

These are exact integer clocks. Passing a congruence check alone is not
sufficient. The example of sums four and `(M+2)^2` has corresponding owner
values equal modulo `M` but unequal as integers for every positive `M`.

## Complete component domain and counting

The interval and hole pullbacks are correct for both signs. The intersection
of the intervals and UNION of the pulled-back forbidden values retains
every local legality condition. A root value determines every path start
and vertex; conversely the root path's start recovers that root. Thus the
compiler establishes an actual bijection, not just a necessary condition
or an orbit count. It has no hidden multiplicity within a root component.

The initial script correctly returned zero solutions when holes covered
an entire interval, but its `empty` flag could remain false. The author
repaired this unpublished interface issue before final replay, both for
local and component domains. The independent checker now requires that
the flag agree with literal domain emptiness. In particular a sum-four
edge in host `[1,2]` has the singleton candidate start two, removed by the
half-diagonal collision, and is correctly flagged empty.

The proof's count uses a constant number of integers per interval and at
most quadratically many holes in each supplied path length. It does not
enumerate all local starts in order to count an uncoupled component. The
test script's explicit root enumeration is only a finite verification
method, not the all-size counting argument.

## Pins, collision guards, and the product boundary

For same-component inequalities the forbidden equation has coefficient
zero or plus/minus two. The four cases (identically equal, identically
unequal, integral root, nonintegral root) are handled correctly. Cross-
component inequalities remain binary forbidden affine diagonals in the
two different roots. They are not unary holes and cannot be removed by
multiplying component counts.

This restriction is substantive: singleton paths with domains `[1,N]`
and arbitrary pairwise distinctness constraints already include ordinary
proper graph-coloring counts. No general efficient counting theorem for
these extra joint factors follows from the one-root-per-component result.
The source appropriately retains them for further exact contraction or
enumeration.

The theorem requires actual index pairs and the declared local distinctness
conditions. Additional overlaps between different paths are allowed unless
the caller explicitly forbids them; this is stated in the input contract.
Square roots must be supplied correctly when the result is interpreted as
a square-sum graph realization. The signed-affine theorem itself also holds
for arbitrary prescribed integer sums.

## The useful boundary against cyclic CRT fill-in

I independently found the following literal square-clock domains:

- no edge, host `[1,2]`: starts `{1,2}`;
- one edge of sum four, host `[1,3]`: starts `{1,3}`;
- one edge of sum nine, host `[1,7]`: starts `{2,3,4,5,6,7}`.

Every pair can identify its start port; all three cannot identify one
common start. The author added this hostile to the proof and receipt.
It demonstrates precisely why the cyclic-coset Helly-two theorem in the
[lattice lane](../lattice/README.md) must not be imported after arbitrary
interval/collision-hole filters. The present compiler keeps the entire
root-domain intersection and therefore detects the obstruction.

## Executed finite controls

I read the full source and independently executed both normal and optimized
replays of `check_shared_clocks.py`. Both returned **PASS, 30,743 checks**,
including the frozen receipt comparison. Its 15,276 exhaustive cases are
exactly the three finite universes declared in SHARED_CLOCKS.md; they do
not claim an exhaustive census of arbitrary clock networks. The explicit
positive triangle `(6,19,30)`, half-integer cycle pin, integral pin hitting
a collision hole, out-of-host triangle, and the incoming order-36 twin
exercise distinct possible failures.

The separate [independent checker](audit_shared_clocks_independent.py)
calls only the source compiler. It builds its expected paths by literal
tuple enumeration and adjacent-sum comparisons, without using the source's
affine constants, reflection recurrence, or brute-force reference solver.
It tests 1,000 seeded networks on one to four actual square-clock paths,
lengths zero through two, with independently selected hosts four through
eight. They include 548 cases with loop equalities, 605 with absolute marks,
691 with explicit inequalities, and 227 retaining cross-component factors;
297 are feasible. Both normal and optimized independent replays return
**PASS, 4,620 checks**.

Reproduction:

    python -B 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/audit_shared_clocks_independent.py
    python -O -B 05-knowledge/results/attachment_lattices_square_sum_20261023/square_sum/audit_shared_clocks_independent.py

The corresponding [receipt](shared_clocks_independent.json) records seed
42319 and the exact finite universe. No network access, external package,
or removable assertion is needed.
