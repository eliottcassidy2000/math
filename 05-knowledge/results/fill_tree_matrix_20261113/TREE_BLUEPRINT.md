# An exact arithmetic tree for the allowed odd targets

**PROVED arithmetic construction; physical tournament lift OPEN.**
The construction describes the desired label set without assuming that
every label is a Hamiltonian-path count. Its missing lift is a concrete
form of the remaining proof obligation, not a proof of Chan–Pak.

**Pure-birth lift REFUTED:** the complete calculation in
[INDUCED_BOUNDARY.md](INDUCED_BOUNDARY.md) proves that no H27 tournament
contains an induced H13 source. Five arithmetic parent edges already
fail through target73. Repairs changing old arcs or other transitions
remain outside this obstruction.

## The binary tree and its two missing vertices

Write an odd target as H=2m+1. The maps

    L(H)=2H−1,  R(H)=2H+1

are exactly m→2m and m→2m+1. Keep the root1, discard its self-loop,
and otherwise use the usual binary parent. In H coordinates,

    p(H)=(H+1)/2 if H≡1 mod4,
    p(H)=(H−1)/2 if H≡3 mod4,       H>1.

The parent is odd and strictly smaller. Thus this is a spanning rooted
tree on all positive odd integers. The forbidden tournament targets7
and21 correspond to the binary addresses11 and1010.

Remove7 and21, and contract the resulting two gaps:

    3→13 and3→15 replace3→7→{13,15};
    11→41 and11→43 replace11→21→{41,43}.

All other parent edges remain. The result is a rooted tree on exactly
the positive odd integers other than7 and21. There are no new roots:
every remaining nonroot vertex has its original parent unless that
parent was removed, in which case the displayed grandparent is used.
Every edge still strictly increases its label, so backwards iteration
terminates at1. Equivalently, remove those two vertices from any ordinary
binary ancestry chain. The two holes do not lie on one another's ancestry.

The splice arrows are4H+1 and4H+3 at H3, and4H−3 and4H−1 at H11.
These formulas are declarations of four particular arithmetic edges,
not universal tournament operations.

## A precise sufficient construction theorem

Suppose a class of states consists of actual finite tournaments with
finite auxiliary data and a designated state above count1. Suppose
every state reached by the following recursion has actual child states
for every edge of the contracted arithmetic tree, with the indicated
Hamiltonian count. Then every allowed odd target is attained.

Proof: use its finite decreasing parent chain, reverse it, and compose
the certified child operations from the root. No probability, density
limit, or exchange of quantifiers is needed. A state may carry different
ports, defects, cores and sizes at different steps. Counts alone do not
certify the next operation.

This gives a sufficient format for a proof by finite local rules **if**
one can supply an operation-closed state class and prove all its edge
obligations. It does not show such a class exists. Unrestricted binary
closure is impossible:3→7 and11→21 would contradict the all-order
omissions. A unit-ready-only state class is also impossible for this
whole tree, because its allowed vertex53084231 has no indegree-one
unit-ready representative at any order. Ordinary representatives at
that count exist.

One can instead use SCC multiplication and the inherited prime-lane
reduction, THM-4094 at
`01-canon/theorems/THM-4094-hamiltonian-matching-deficit-and-two-prime-lane-completeness.md`:
strong realization of every odd prime p≠7 and every7p with odd prime
p≠3, together with the proved carry atoms49,63,343, suffices. This is
a multiplicative forest of witnesses rather than a binary construction
tree. Neither criterion is discharged by a large finite prefix.

## The numbers5,17,7,21

There are several exact connections, with different meanings:

* 5→9→17 follows L twice; both5 and17 lie on the complete fan line
  H=1+2^s. Their least strong witnesses have orders4 and6.
* 3·5+2=17 is the arithmetic C2 relation, while3·7=21 is a factorization.
  The first requires a suitable physical interface to be a construction;
  the second explains why ordinal multiplication cannot rescue21 using7.
* Their binary count addresses are10,1000,11,1010 respectively. The two
  omissions are holes in this arithmetic tree, repaired by the four
  splices above when only labels are considered.
* Under odd-accelerated ordinary Collatz,7→11→17→13→5→1, whereas21→1.
  Thus the four numbers belong to the same known basin. Basin membership
  does not distinguish attained tournament counts from the two omissions.

The strong-floor proof of the omissions is inherited from
THM-1370 at
`01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`.
There is no implication from being a Fermat prime or lying in the Collatz
tree to being a tournament count. The actual5/17 multiplicities and their
neutral-growth geometry are recorded in the multiplicity lane.

## Decisive checks and scope

`tree_blueprint.py` checks every odd label through131071, its complete
parent chain, all splice destinations, the ordinary odd Collatz paths
above, and the arithmetic identity between binary shifts and L/R.
These finite controls support the elementary all-label proof; they do
not test the open physical lift. The binary addresses are encodings of
target integers, not tournament vertex labels or adjacency rows.
