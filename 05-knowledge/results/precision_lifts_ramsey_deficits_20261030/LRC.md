# Two arbitrary attached speeds, in their common integer clock

**PROVED:** after a fixed core's complete safe intervals are available, exact
simultaneous attachment of any two specified positive integer speeds can be
decided using Euclidean floor sums, without enumerating either new comb.
**FINITE-EXACT plus an analytic tail:** every pair of positive integer speeds
can be attached at radius1/14 to each of four explicit eleven-speed cores.
**OPEN:** coverage for arbitrary eleven-speed cores, and universal LRC(14).

The closest proved mechanism is the preceding packet's
[full phase cover for one attachment](../terminal_fees_phase_covers_20261030/LRC.md).
The canonical hostile is core(1,2) at radius1/4: attaching3 alone succeeds,
and attaching12 alone succeeds, but attaching both fails. The correction is
to retain their common native time interval. The least-used sidecar is the
exact integer floor sum; this is a change in the consumer of interval data,
not another endpoint-only approximation. No literature priority is claimed
for the standard Euclidean floor-sum routine.

## L1. Count actual intersecting tooth pairs

Fix a closed core interval I=[L,R], rational r=s/d in(0,1/2], and positive
integer speeds u,v. A u-tooth has interval

    J_k=[(k+r)/u,(k+1-r)/u].

The intersecting indices are ceil(uL-1+r)<=k<=floor(uR-r). The fully
contained indices are

    a=ceil(uL-r) <= k <= b=floor(uR-1+r).

There are at most two clipped teeth outside this latter range: the first
and last intersecting ones. Process these directly, retaining closed
endpoints, even when a component or tooth is a single point.

For a fully contained tooth, the v-tooth indices that meet it are exactly

    v(k+r)/u-1+r <= l <= v(k+1-r)/u-r.

Put M=du, A=dv, B=v(d-s)-su, C=vs-u(d-s)-1. Their number is

    f(k)=floor((Ak+B)/M)-floor((Ak+C)/M).                (3)

The minus one in C is essential: it preserves integer equality at a closed
boundary. Since the original interval for l is nonempty as a real interval,
f(k) is nonnegative. The exact number of tooth-pair intersections is

    sum_(k=a)^b f(k) + the two clipped-tooth counts.     (4)

Each term is a Euclidean floor sum. Summing over the complete disjoint core
intervals is an exact decision criterion. Distinct tooth pairs have disjoint
nonempty intersections because r>0, so they also count the resulting
components when the core cover was merged first. This retains the paired
tooth labels; marginal success bits cannot substitute for (4).

If the full range has positive count, binary search its nondecreasing prefix
counts to find a k with f(k)>0. Set l to the least integer in its allowed
range and return t=max((k+r)/u,(l+r)/v). On a clipped tooth also intersect
with I. This produces a rational witness in the original clock, not a
residue representative that would need a separate lifting argument.

For K core components, decision needs O(K log M) arithmetic operations;
witness extraction adds a logarithmic search over the tooth-index range.
Bit lengths of speeds and rational endpoints remain part of the cost. Core
construction, component count and denominator size are not claimed bounded
uniformly. No O(1) big-integer claim is made.

## L2. A finite stopping certificate for two infinite tails

Suppose 0<r<=1/6 and the core has an interval of width W>0. Put B=ceil(1/W).
If min(u,v)>=B, both speeds attach successfully. To prove it, order u<=v
and take a subinterval of length1/u in that core interval. The u-safe part
has total length(1-2r)/u and at most two pieces. One piece has length at
least(1-2r)/(2u)>=2r/u>=2r/v. On an interval of length2r/v, every starting
v-phase can reach its closed safe set. This proves the claim with equality.

For 1/6<r<=1/4 the same conclusion holds with
B=ceil((2-2r)/W): the interval then contains an entire u-safe tooth, of
width(1-2r)/u>=2r/v. The smaller B is not asserted by this proof in that
larger-radius range.

Only finitely many smaller speeds1<=u<B remain. Compile the **joint** core
C union{u}. If its widest interval has positive width W_u, every
v>=ceil(2r/W_u) is accepted by the preceding one-tail theorem. Check the
finite range u<=v<ceil(2r/W_u) using (4). A point-only or empty joint core
requires its own exact residue analysis; it is not waved away by this rule.
These tests and the two tails form a complete certificate for all unordered
positive pairs for the fixed core. Repeated speeds are harmless and included;
the claim in particular covers all distinct admissible speeds.

## Computation and independent controls

Run `python3 joint_phase.py` from this packet. It performs780 signed floor-
sum checks against literal sums, 8,620 literal tooth/complete-comb comparisons
including points, and42,767 finite bank queries checked independently with
the preceding one-tail compiler on the already intersected core. Every
returned witness is checked against every actual speed.

At r=1/14 all positive pairs are accepted for these eleven-speed cores:

* (1,2,...,11);
* 26*(1,2,...,11);
* (1,2,3,5,8,13,21,34,55,89,144);
* (2,3,5,7,11,13,17,19,23,29,31).

The bank stores every small-speed branch, its tail cutoff and its complete
finite rejection list, all empty. Positive pairs with101-digit and1001-digit
coordinates also produce exactly checked witnesses. The hostile(1,2)+3+12
is correctly rejected at1/4. See[joint_phase.json](joint_phase.json).

## Relation to the incoming odd-descent work

[odd_descent_local_relaxations_20261030/PROOFS.md, P3](../odd_descent_local_relaxations_20261030/PROOFS.md)
proves that odd-denominator compatible marginals descend when every constraint
is preserved by an appropriate symmetric operation. On a fixed tooth, the
phase equation y=vt-k is increasing affine, so coordinatewise medians preserve
it. Across wraps, the relation y={2t} fails this property: the three tuples
(1/8,1/4),(3/8,3/4),(5/8,1/4) have median(3/8,1/4), which is not a tuple
of the relation. A unary subset of a line, including a union of teeth, IS
median closed; the failure concerns this coupled time/phase relation.
The actual transfer used here is to retain the integer tooth labels and
test their common interval. Separate successful attachments do not supply
a common time marginal, so the descent theorem cannot combine them.

The reusable move is to compile the **joint** labelled relation before taking
an existence bit. In Collatz the analogous obligation is retaining every
compatible multiplier lift, then paying its charge. In Ramsey, matching
neighborhood inventories will likewise need a shared-edge consistency test.
