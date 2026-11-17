# Exact separator center and the three-tape entropy ceiling

Use charged tiles A (directed triangle, H=3, half-charge −1) and
B (H=161, half-charge +54), with the masks in export.cpp and
transfer.py. Place one before the other; the shared internal
separator V is ordered transitive. Bit 0 means losing to both
ports; bit 1 means beating both. Either tile order gives unit
half-charge, \(-161+3\cdot54=1\), independently of the separator.

## Structural obstruction to arbitrary high entropy

In a fixed A→V→B quotient, choose distinct departure singletons
A1,A2 into A from x,a, and distinct arrival singletons B0,B1
from B to x,a. Let A0 and B2 be any Hamilton paths on their
remaining nonempty vertices. These distinct choices exist for
both tile orders (two distinct nonempty proper masks suffice).
The fixed path template
\[
 A_0V_0B_0\,x\,A_1V_1B_1\,a\,A_2V_2B_2
\]
now admits every allocation of separator vertices among the
three increasing runs (each may be empty). These are \(3^m\)
distinct paths for m separator vertices, regardless of contact
word or branch. Consequently the whole fixed two-tile family
with arbitrary internal binary separator has only
\(O(X^{\log_3 2})\) address presentations with \(H\le X\).

This also holds with external shared tapes L,R, up to
\(O(\log^2(2X))\) additional multiplicity. Use the [external
gate convention](../signed_ordinal_unit_centers_20261117/TWO_SIDED.md), trim the
initial global sources and final global sinks, and add L_j
before A_j and R_j after B_j in the same template. If a tape
has at least two gates, force the two ports' adjacent endpoints
to use the prescribed gate; then there is no marked contact
constraint involving A_j or B_j on that side. Thus each high
tape contributes independently \(3^{i-1}2^{r-1}\) allocations
in the preceding proof's depth/gap variables, in addition to
\(3^m\). On a low side put its whole tape outside the template;
its independent one-gate bound gives only O(log X) possible
lengths (or the tape trims to empty). When both sides are
high, the summation is over \(m+i+j-2\) at fixed gaps, at
most O(log² X) summands per level, and the two gap sums are
geometric. The one-high and zero-high cases are bounded
identically after charging the low lengths. This proves the
all-three-tape sparse upper bound. In particular the union
over both tile orders misses infinitely many primes.

## An exact finite linear representation

Remove the two ports from a full path. Each remaining slot
contains an ordered sequence of portions of A,V,B, some empty.
Enumerate ordered assignments of A to the three slots and of B
to the three slots. The internal path counts in A and B,
with imposed starts and ends as needed, are exact subset
quantities. For each hypothetical mask of nonempty V slots,
the count conditional on V is a product of these A/B path
counts, with only the following V tests: if it follows a port
without a preceding A piece, its first bit must be 0; if it
precedes a port without a following B piece, its last bit
must be 1. An entirely empty middle slot is allowed only
in the direction of the switched port arc. Sum over both
orders of the ports. This gives a finite nonnegative integer
readout table on three triples of masks (nonempty slots,
required first-zero slots, required last-one slots), exported
by export.cpp.

While reading V left-to-right, each slot's state is empty
or has a pair (first,last) in {0,1}². For each new bit one
chooses one slot and updates its last bit (also its first
if empty). These two exact integer transition matrices
have 125 states. Slot symmetry compresses them to the
35 multisets of three states drawn from five options; the
read vector is averaged over each orbit. The center is the
minus readout plus one.

Rational row-reduction of the exact matrices (rank.py)
finds all 35 dimensions reachable on this compressed space.
The observable space has dimension 10 for either tile order,
and their joint observable space has dimension 18.
Equivalently the two separate center series have Hankel
rank 10 and the joint series rank 18. The centered
commutator section for exchanging 01 with 10 has rank 6:
the matrix with entries center(u01v)−center(u10v), over arbitrary
prefix and suffix, has that rank in either tile order. These
ranks refer to the exact signed tile data above, not universal
ranks of all such gadgets. In particular a terminal commutator
is not just a single-bit substitution: for AB at empty prefix
the swap response is 310480, whereas the one-bit difference
is −111862.

Two state encoding boundaries matter. With codes 1,2,3,4 for
00,01,10,11, a first appended 1 creates state 4, not state 2;
and a required final 1 forbids states 1 and 3 (odd), not 2
and 4. The independent literal subset DP checks the exact
counts, not just their difference.

## Replay

    g++ -std=c++17 -O3 export.cpp -o /tmp/sepex
    /tmp/sepex > weights.txt  # reference table; embedded in transfer.py
    python transfer.py
    python rank.py

The direct check runs both port orientations, both tile
orders, and every word of length through 3. The larger
transfer-only check verifies odd counts, gap two and the
symbolic lower bound through length 10. This constructs
arbitrary exact center values as a finite-matrix readout,
not an exact inverse at a prescribed target.
