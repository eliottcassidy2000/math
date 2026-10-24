# Sharp joint eight-slot / twenty-two-clock connectivity in the complete core

Status: **FINITE-EXACT** for the complete authenticated 28-source binary menu at
**tournament order eight**. The reduction, complement transport, lower-bound
argument and primitive-move replacement below are proved. The uniform upper
bound is an exhaustive finite theorem, independently replayed with two distinct
algorithms. This is not a theorem at arbitrary tournament order.

## Concurrent priority and independent contribution

The incoming `sixslot570_20261024` packet (commit ae8d0f15654) already states
the all-margin sharp(8,22) theorem, proves six-slot connectivity at margins
6–22, and provides the optimal220 extra eight-slot types. Our work began
against the preceding1d94c1d1d9e snapshot before that packet was materialized.
The uniform bound and lower witness below are therefore an INDEPENDENT
verification of the concurrent result, not a first claim. The new additions
are the exact coordinate bridge, independent sparse-bank audit, explicit
weighted-circuit/complement distinctions, and the denser-graph depth bound.

## 1. Exact statement and supply

Keep the original 28 source types, each with its two actual pinned paths to
U=1965 and V=2234. No new helper, relaxed pin, repeated source copy or extra
normalized resource is introduced. At source s let b_s=1 mean target V, and let

    sigma_s = (X_s mod128, Y1_s mod2, ..., Y4_s mod2),
    Q = Z/128 x (Z/2)^4,
    c_s = number of differing indexed owner clocks between its two paths.

For every nonempty fibre

    F(k,q) = {b in {0,1}^28 : sum b_s=k, sum b_s sigma_s=q},

and every pair of its states, there is a route such that every event has BOTH

    changed source slots <= 8,
    actual ordinary owner-clock occurrences <= 22.

Both constants are sharp uniformly over this supplied family. The same
margin-four two-state fibre attains both lower bounds.

The exact clock distribution is twenty cost-two sources, five cost-three
sources and three cost-four sources. Costs are recomputed from literal path
orders, not inferred from charge coordinates or source labels. All events
retain the full original paired dyadic nodewise-current bank, fixed target
margins, original common-clock records and chosen ordinary per-slot X phase.
They need not retain exact integer-X current or odd-primary guards. Every
ordinary individual-port harmonic phase is a different, stronger predicate;
the authentic character 4e_U already forbids these target changes under it.

Inherited native authentication is replayed in normal and optimized Python:
11,828 actual Hamiltonian paths; trivial automorphism groups for all28 sources;
56 available normalized F/converse resources with224 distinct closed endpoint
words. Each arrival uses28 resources /112 words. An immutable remainder must
avoid the complete alternative reservoir. These are partial-packet routes, not
a rewiring of the fully occupied canonical atlas.

There is one raw common-clock path at source1923 with target1559. The specified
zero target1559 margin excludes it. Thus completeness is under pins PLUS
specified target margins, not under pins alone. Independent binary mixtures
are compatible because the alternative resource orbits are disjoint. No
additional cross-source owner identifications are silently imposed.

## 2. Complete finite upper certificate, including the middle margins

Every subset of each size0 through14 is enumerated. Both graph computations
agree that every nonempty fibre is connected under the JOINT bounds (8,22).

| k | All states | Nonempty fibres | Disconnected at (8,22) |
|---:|---:|---:|---:|
|0|1|1|0|
|1|28|28|0|
|2|378|351|0|
|3|3,276|1,646|0|
|4|20,475|2,048|0|
|5|98,280|2,048|0|
|6|376,740|2,048|0|
|7|1,184,040|2,048|0|
|8|3,108,105|2,048|0|
|9|6,906,900|2,048|0|
|10|13,123,110|2,048|0|
|11|21,474,180|2,048|0|
|12|30,421,755|2,048|0|
|13|37,442,160|2,048|0|
|14|40,116,600|2,048|0|

This is154,276,028 directly enumerated states. Complementation supplies all
remaining margins, covering all268,435,456 binary states and47,060 nonempty
fixed-margin/full-charge fibres. In particular margins9–19 are no longer an
untested gap for the joint(8,22) bound. Our independent run does NOT test support-six
connectivity at those margins; the incoming packet separately settles it.

The primary checker recursively generates subsets, groups them by exact Q
charge and computes union-find components. Its early exit is valid only when
all already processed vertices form one component: omitted edges among those
vertices cannot change that partition. If multiple components remain, it
checks every preceding vertex. Every subsequent vertex is still processed.

The independent checker uses numeric Gosper subset enumeration, independent
14+14 charge/weight lookup tables, and deletion-list breadth-first search.
It has no union-find or producer enumeration dependency. Both enumerate
entire fibres, not selected charges or sampled endpoints. A third relaxed
clock-only22 run agrees. No floating-point arithmetic is used.

## 3. Complement isometry, not a within-fibre palindrome

The full source sum is3733, hence the total charge is (21,15), with the second
coordinate packing the four parities. For C(b)=1-b,

    C: F(k,(x,y)) -> F(28-k,(21-x mod128, 15 xor y)).

Furthermore C(b) xor C(b')=b xor b', so both support and weighted work are
preserved. This proves the upper-margin transport without any extrapolation
in tournament order or source supply.

At margin14 complement NEVER fixes a charge fibre: 2x=21 modulo128 is
impossible. Thus this is a cross-fibre symmetry. The total clock weight67 is
odd, but that does not make each individual fibre's cost enumerator
palindromic or force an equal even/odd split inside that fibre. The restricted
18/37 circuit has a different, within-fibre endpoint-complement symmetry.

## 4. The sharp two-state lower witness

The ENTIRE margin-four, charge(123,11) fibre is exactly

    A={1532,1691,1923,1974},
    B={1841,2176,2837,3360}.

An independent20,475-subset enumeration checks this completeness. The states
are disjoint, so their unique nontrivial event changes eight slots. Their
literal clock costs are

    A: 4+2+4+2 = 12,
    B: 3+2+2+3 = 10,
    total: 22.

No budget21 clock route or budget7 support route exists because there is no
third feasible state. Their direct event attains (8,22). Their actual relative
owner permutations comprise four transpositions, two3-cycles and two4-cycles.
The cost is therefore4*2+2*3+2*4, not a numerical identification of prime11.

The integer source sums are507 and635. Their difference128 disappears in Q,
but the actual target-owned integer-X current changes by+256 at U and-256 at
V for the transition A to B. This is NOT an exact-integer-current witness.
The earlier charge(15,3) eight-slot hostile had equal integer sums527 and
established a different exact-X-compatible support obstruction.

## 5. A weighted binary move-bank interpretation

Write z=b'-b. A move has sum z_s=0 and sum z_s sigma_s=0, with every z_s in
{-1,0,1}; its cost is sum c_s |z_s|. Endpoint margins force an even changed-slot
count. Distinct source charges exclude nontrivial two-slot moves.

The complete unordered disjoint balanced move counts are:

| Changed slots | All equal-Q moves | Moves of clock cost<=22 | Primitive among those |
|---:|---:|---:|---:|
|4|27|27|27|
|6|1,917|1,917|1,917|
|8|53,220|51,925|51,593|

Here primitive means no proper conformal balanced subrelation. With distinct
source labels, six-slot moves are automatically primitive; an eight-slot
move is nonprimitive exactly when it contains a proper two-per-side relation.
Such a move splits into two feasible disjoint four-slot moves. Each has no
larger clock cost than the original, and the intermediate state stays binary.
Therefore the census also proves that the53,537 primitive moves in this
bounded bank connect every binary fibre. This is a deliberately redundant
binary move bank, not a claim of a minimum basis or a Markov basis for all
unbounded nonnegative-integer fibres.

The old27+1,917 four/six trades generate the full additive charge kernel but
fail to connect some binary boundary fibres. The two-state lower witness
forces an eight-slot,22-clock move inside the complete28 supply. Large
conformal circuits, such as18/37, can nevertheless be bypassed using other
slots. The distinction is between conformal primitiveness in a restricted
endpoint union, and connectivity in the full capacity-one cube.

## 6. Incoming coordinate bridge and sparse-bank audit

The incoming packet writes the sharp charge as(binary6,cyclic123), whereas
our inherited domain file writes(cyclic123,binary11). These are the SAME
fibre. The input bases are(Y0,Y1,Y2,Y3) and(Y1,Y2,Y3,Y4), respectively. On
every supplied source, Y0 = (X mod2) xor Y4. Writing our charge as(x,y),
the incoming packed parity is

    p = ((x mod2) xor y_bit3) | ((y &7)<<1).

Its inverse is y=(p>>1) | (((p&1) xor (x&1))<<3). This is an automorphism
of Q, checked on all2048 charges and their additive generators. It
preserves every fibre and sends our(123,11) to incoming(6,123). The packed
integer11 is therefore not an intrinsic prime label.

`audit_incoming.py` independently reconstructs every margin4/5 component
using our earlier native records. It verifies all220 incoming bridge types,
every occupied one-source lift at margin5, and all21,673 audit checks in
normal/optimized Python. The bridge cost histogram is

    16:35, 17:63, 18:60, 19:37, 20:18, 21:6, 22:1.

Every chosen type reduces the margin4 quotient-component excess by exactly
one; the total excess is220, proving that the incoming lower bound is
attained. With all1944 small moves this is a2164-type alphabet, far smaller
than our deliberately redundant full bank. We do not relabel this inherited
minimum-alphabet theorem as new.

## 7. A denser grammar gives a shorter uniform root route

The fully allowed joint(8,22) graph has53,869 unordered move types. A complete
numeric-subset/deletion-BFS replay additionally records the depth of every
state from the largest-mask root in its fibre. The maximum depths at
margins0 through14 are

    0,0,1,1,2,2,2,3,3,3,3,3,3,4,4.

Thus EVERY fibre has a rooted spanning tree of depth at most4. Complement
transports these trees to the remaining margins. Every pair can be joined
via the root in at most8 events and at most176 actual clock occurrences.
These are upper bounds, not a minimum-diameter or minimum-total-work claim.

This improves the incoming10-event/220-clock tree-route bound only after
allowing a larger grammar. It does not assert that the optimal2164-type
alphabet admits these shallower trees. Removing nonprimitive moves retains
the176 total-work bound because conformal splitting partitions changed
slots, but can increase the event count.

For any positive state law pi and reversible kernel giving positive
conductance Q(e) to this new tree, the same exact cut-demand argument gives
length-weighted congestion at most

    max_e 8*pi(A_e)*(1-pi(A_e))/Q(e).

The new tree edges and their actual conductances must be used. One cannot
silently replace10 by8 for a kernel supported only on the incoming sparse
alphabet. No scale-uniform spectral-gap assertion follows.

## 8. Constructive comparison with Session46's18/37 circuit

Let the inherited18-slot endpoints have nine V sources each, charge(9,7), and
integer source sum1161. Add the ten omitted core sources, initially and
finally at U. The complete28-slot fibre has3337 dyadic states; its exact-X
subfibre has1355. Restricting all ten omitted sources permanently to U leaves
exactly the original two states, at distance18 slots /37 clocks.

The saved `route18_37.json` gives independent full28 replays:

- Dyadic route: five events, clock costs12,12,12,9,12; total57; peak12.
  Only outside source2071 moves. Intermediate integer sums1033 produce a
  target-owned debt of-256 at U, invisible to the original dyadic bank.
- Exact-X route: four six-slot events, clock costs13,12,12,12; total49;
  peak13. Outside sources176 and2071 move, and all arrivals have sum1161.

The complete full28 lower-budget component separates these endpoints at
clock11 in the dyadic fibre and clock12 in the exact-X fibre. Thus the sharp
endpoint clock bottlenecks in these two full28 domains are12 and13. At these
budgets support<=6 is automatic, so these lower checks do not hide an extra
support restriction. The routes restore all ten outside slots and retain
all full native records and normalized resource constraints by the inherited
binary realization.

These numerical route peaks/totals are NOT new improvements over Session46:
its cost37/variants packet already supplied one-buffer2071 routes with the
same peaks12/13 and totals57/49 in a19-slot menu. The new calculations replay
these phenomena in the complete28 domain; the saved exact-X path happens to
choose a different two-buffer representative. No minimum total-work or event
count theorem in the enlarged28 domain is asserted here.

Under strict temporal nonreturn, restored auxiliary slots cannot move. In
that stronger model the original primitive18/37 obstruction remains intact.
The routes here use returnable capacity-one resources, never simultaneous
reuse. Every repeated event is charged.

## 9. Conditional induction interface

`SHEET_GLUING.md` proves an all-size occupied-sheet gluing lemma with the
correct deletion-relative hypotheses, and an exact component-overlap quotient
when every allowed event leaves a common occupied source. A complete native
probe finds eight margin4 fibres in seven single-deletion banks that need23
clocks even though the complete28 bank needs only22. This is a concrete
obstruction to a tempting but invalid induction, not an all-size resolution.

## 10. Reproduction and boundary

Run from any working directory:

    python /path/to/this/reproduce.py --repo /path/to/math

This compiles the two independent exhaustive methods plus the relaxed
clock-only control, reconstructs routes and the move-bank census, reruns the
inherited full native authentication, and performs the new checker normally
and with Python optimization enabled. The final checker receipts are
byte-identical. `MANIFEST.json` hashes the deliverables; `normal.json` pins
all native source inputs. No discovery cache or Git mutation is needed.

This independently verifies the incoming closure of the complete finite
core's uniform (8,22) question. It does not
supply native resources at another tournament order, compatibility with an
arbitrary occupied remainder, exact-integer-X uniform22 connectivity, odd
harmonic guard preservation, or universal per-port phase preservation. The
larger helper controller concerns a different family, allowing core charge
to change while helpers pay the residual. Its claims are neither needed nor
superseded by connectivity of each fixed-core-charge fibre proved here.
