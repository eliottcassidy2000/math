# All-order image tests, and where the construction stops

**PROVED** termination and completeness of bounded-H search; **FINITE-EXACT**
receipts for the indicated integer ranges. A missing value in this packet
means absent from the named family at every order. It is not a new gap in
the unrestricted spectrum.

## R1. Why a finite search can exclude all vertex orders

For every nonempty golden word its state satisfies

    A,B>=0,  H=A+B>=1,  H<=C<=8H.                            (R1)

The lower bound C>=H also has a direct injection: cut any Hamilton path
just before the port a. Its suffix starts at a; its prefix is either
empty or the second path. The original path is reconstructed by
concatenating prefix and suffix, so this is injective into C's objects.

For each nonsingleton matrix M, every row satisfies

    1<=M_iA+M_iB<=M_iC<=8(M_iA+M_iB).

For s, H'=A+C>=H and C'=H+C<=2C<=2H'. Thus (R1) is inductive, including
each first letter from the formal empty state. Every letter has H'>=H.
The auxiliary integer W=A+B+C strictly increases at every step: for s,
W'-W=A+C>0; for the other matrices each row sum is greater than1.
The initial unit also increases W.

Fix X. There are finitely many integer states satisfying H<=X and (R1).
No path can cross above X and later return. Breadth-first closure from
(1,0,0), discarding only H>X, therefore enumerates exactly all reachable
states with H<=X over all word lengths. W bounds word length in this
sublevel set. Parent pointers give actual positive witnesses, and
finite forward-closed state sets containing the unit certify absence.
This is stronger than enumerating tournaments only through orderN.

The exact three-coordinate observer also makes this state identification
lossless for all future H queries within the stated language. Source words
remain attached to witnesses whenever a literal construction is required.

## R2. Golden prime and 7p lanes

[reachability.json](reachability.json), reproduced by reach.py --limit100000,
contains119590 states,30194 distinct H values, and the state-set SHA256.
Dijkstra's algorithm uses the actual vertex costs to find minimum order
within this language. The largest such minimum among these states is20.
For H1 the formal empty witness may be replaced by the physical singleton.

| Required lane | Values tested | Attained in golden image |
|---|---:|---:|
| odd prime p<=100000, p!=7 |9590|4907|
|7p<=100000, p odd prime, p!=3|1674|851|

The first missing prime is37; the next are43,47,67,71,79,83,103.
The first missing7p is133=7*19. Its next missing prime factors are61,67,83,89.
The full small closed certificates have61 states at cap37,206 at cap133,
and536 at cap343. verify.py checks every in-range successor explicitly.
These certify37,63,133,343 absent from the golden image at every order.
All are known to occur in other tournament families; none is a proposed
unrestricted gap.

Examples, in processing order (reverse source-component order):

| H | A supplied word |
|---:|---|
|49|s s q5 s|
|55|s s q1 s|
|57|q1 q2|
|77|q0 q5|
|91|s q1 c|
|111|c s s q2|
|119|s s s s s q2|

The target-lane reduction itself is inherited from
[THM-4094, Hamiltonian matching deficit and two-prime-lane completeness](../../../01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md).
It uses strong realizations of all such p and7p, together with49,63,343
and ordered-join multiplication. A finite count range does not discharge
those two infinite obligations. This golden family cannot discharge them
even in principle, since37 and133 are permanently missing.

## R3. The complete strong-five block extension

Of the64 fixed-spine five-vertex charts, exactly50 are strong. Include all
of them as source atoms before applying actual K. Together with the original
seven letters, this is precisely the source language with all strong
components of size at most5. The same observer has57 explicit nonnegative
integer matrices, independently reconstructed from literal graph paths.
The cone becomes

    H<=C<=16H,

and the monotonicity and W proofs are unchanged. Thus its bounded-H
closure also classifies membership over all vertex orders.

At H<=10001 there are147655 states and4991 distinct odd values. Exactly
the following odd values are missing:

    7,21,47,83,103,491,709,833,3227,4619.

The prime-lane misses are47,83,103,491,709. The only7p-lane miss in this
range is3227=7*461. The cap47 closed certificate has171 states, proving
47 absent at every vertex order. Hence repairing47 while remaining an
actual path-complement-of-SCC-word construction requires a source SCC
of at least6 vertices. This is a block-size necessity, not a lower bound
on all possible tournament constructions or a sufficiency assertion.

A particularly small repair of the original obstructions uses one extra
chart, p5_2, with source out-neighbour masks [10,4,9,18,7]. Its matrix is

    [[5,0,14],[5,0,13],[27,5,46]].

Then s p5_2 has H37, and c p5_2 has H133. Another actual chart p5_6,
with source rows [26,4,9,18,6], has matrix

    [[10,1,13],[9,1,12],[22,5,38]],

and p5_6 c has H63. The full-bank witness for343 is p5_5 s s q0, where
p5_5 has source rows [22,4,8,19,6]. Selected literal targets, their states,
and strong-connectivity tests are retained in verification.json.

No claim of a new unrestricted finite coverage record is intended.
Incoming Session150 has a much larger general finite prefix; this result
classifies the expressive limits of specific small-alphabet recurrences.
The decisive next problem is a source-preserving rule for choosing larger
blocks, or a different operation, with an infinite prime/7p covering proof.

**Incoming completion of the finite repair stage:**
[INCOMING](INCOMING.md) transports the independently published strong6 atlas,
replays its eight repairs, and cites its fixed65-symbol bank's holdout through
30000. The next work should use that bank rather than repeat the strong6
enumeration. The arbitrary-block recurrence is now proved; all-value image
coverage still needs a separate theorem.
