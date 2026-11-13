# A neutral core-growth bridge reaches the strong rank-one reset

**Status:** FINITE-EXACT (literal bridge, isomorphism, response coefficients,
and two-port tests); PROVED (universal single-port inflation obstruction from
those coefficients). The bridge was found by the
[absorption lane](../absorption/corridor.json). The rank-one reset was found
by the [orbit lane](../orbits/PROOF.md). This note connects and independently
checks the two objects; it does not claim a new universal reset construction.

## The same actual seven-vertex tournament in two constructions

Use the outgoing adjacency rows

    A=(44,17,2,6,13,30,63),
    B=(110,16,2,6,13,30,62).

They are the absorption census codes 1162 and 33931. Reversing just the arcs
on {0,1} and {0,6} transforms A into B. Both have H=25, verified independently
by all 7! permutations and subset dynamic programming. A is a singleton
source 6 preceding a strong six-vertex H=25 core; B is strongly connected.
Thus this move absorbs neutral padding into a larger strong core. A split
into two H=5 SCCs cannot explain it: each such component needs at least four
vertices, exceeding the seven vertices available.

At the retained vertex 1, the orbit lane's rooted state (h,E,O) changes as

    (25,12,15) -> (25,10,15).

All simple directed paths starting at vertex 1 were independently enumerated
to check the parity-weighted complement sums. The unmarked P₂ changes from
608 to 548, and the full cover polynomials are

    pi_A(z)=z^7+12z^5+12z^3,
    pi_B(z)=z^7+10z^5+12z^3+2z.

The orbit construction M6² has rows

    (124,97,114,118,2,16,48),     root=4.

The explicit isomorphism from these labels to B's labels is

    (0,4,5,6,1,2,3).

It sends root 4 to B's vertex 1 and preserves every directed edge. Consequently
the endpoint is exactly the rooted strong H=25 reset of the orbit lane, whose
matrix is

    M6²=(5,2,3)^t (1,4,4).

This is the substantive cross-lane connection: a count-preserving core repair
lands at the same rooted object as a rank-one state reset. A rank-one state
reset need not be an SCC decomposition, since B is strong.

The two arc reversals must be treated as a simultaneous move for this exact
target. Reversing only {0,1} gives H=15; reversing only {0,6} gives H=49.
The Boolean two-edge count polynomial is therefore

    H(x,y)=25-10x+24y-14xy,

with (x,y)=(0,0) and (1,1) the two H=25 endpoints. Cancellation depends on the
interaction term; neither individual flip is neutral or a count step of ±2.

## The smallest universal module lift fails at every port

Replace a vertex p by an arbitrary nonempty tournament M, giving all its
vertices the old vertex's orientations to the six other vertices. Write
c_k(M) for ordered k-tuples of nonempty directed paths covering M. Then

    H(A[p<-M]) = sum_k a_k c_k(M),
    H(B[p<-M]) = sum_k b_k c_k(M).

These exact coefficients have a small independent compiler. Enumerate the
6! orders of the six external vertices. Each of their seven gaps is forced
empty, forced occupied by an M-path, or optional, according to the two module
boundary arcs and the direct external arc. A choice of k occupied gaps pairs
with exactly one ordered k-path cover of M. Thus the resulting coefficient
list is exact for every M, without imposing a size bound on M.

The computed differences, for nonempty M, are:

| Replaced vertex | H(B[p←M])−H(A[p←M]) |
|---:|---|
| 0 | −9c₂−6c₃ |
| 1 | −23c₂−14c₃ |
| 2 | −2c₂ |
| 3 | −2c₂ |
| 4 | 22c₂+12c₃ |
| 5 | −2c₂ |
| 6 | 10c₂ |

**PROVED consequence.** No vertex can be replaced by any tournament of at
least two vertices while preserving the H equality of this bridge. Indeed,
c₂(M)>0: split M into two nonempty induced tournaments, each of which has a
Hamiltonian path. Every displayed difference has a fixed nonzero sign. For
a singleton, c₂=c₃=0 and the original equality returns. This completely
classifies the single-port equality boundary.

The full coefficient arrays, including the k=0 deletion terms, are in the
[receipt](receipt.json). Independent literal subset DP checks the formulas
for each port with transitive M of every order 0 through 7, and also C3 and
the strong four-vertex atom Q: **140 wrapper checks**. The transitive tests
alone uniquely determine every coefficient: c_k(TT_m)=k! S(m,k), a triangular
matrix for 0≤m,k≤7 with nonzero diagonal. This provides an independent exact
coefficient check, not just sampled large-module agreement.

## Two-port hostile probe and the preserved older positive result

For each of the 21 vertex pairs, replace both vertices by the same transitive
two-vertex tournament, and separately by the same transitive three-vertex
tournament. Every pair fails equality in both cases. The smallest absolute
two-vertex difference is 6, at the pair {1,6}. All exact differences are in
the receipt. This refutes a universal identical-module lift at any two ports.
It does not classify isolated larger modules or unequal simultaneous modules.
The bounded inflation probe stops here.

This obstruction does not retract the older positive
[REPAIR23.md](../../surplus_selector_20261113/REPAIR23.md). That construction
uses two specifically matched rooted six-cores and an arbitrary *ordinal
prefix*, preserving the marked response needed by a later operation. It has
an explicit prefix guard and reusable chart entry. The failed proposal here
substitutes a module for one vertex of a different H=25 bridge. The precise
operation and the preserved response differ.

The corrected next question is to search for a **matched core family with
equal external-cover coefficient vectors at the chosen interface**, then ask
whether the required rooted capability differs. Equality at one scalar input
is insufficient; this seven-vertex example shows exactly which c₂ and c₃
coordinates a proposed lift destroys. No further size search is justified
without changing that interface condition.

Reproduce all checks with the same
[check.py](check.py) command listed in [PROOF.md](PROOF.md).
