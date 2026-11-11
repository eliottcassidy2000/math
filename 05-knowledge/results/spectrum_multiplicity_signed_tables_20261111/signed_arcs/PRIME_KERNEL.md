# The 25-prime square-sum graph has a three-bit cycle kernel

Date: 2026-11-11. **PROVED:** the kernel and its arbitrary-subdivision parity
compiler; the directed three-cycle example and its completion bound.
**FINITE-EXACT:** the 25-prime edge inventory and the literal H13 comparison.
This is a new operation on an inherited graph, not a claim about why Ogg's
prime list exists or an identification of square-sum paths with tournament HPs.

## Recovered object, without changing its vertex convention

The precise old graph is in
`05-knowledge/results/machines_automata_fractran_tournaments_dictionary_20261101/README.md`,
§3.3, and THM-4860,
`01-canon/theorems/THM-4860-the-prime-square-sum-graph-is-bipartite-by-residue-mod-4-two-primes-1-mod-3-are-never-adjacent-chebyshev-bias-obstruction-to-prime-square-sum-paths-none-below-615731-and-the-25-primes-below-100-by-ogg-layer.md`.
Its vertices are the **25 primes themselves**, not their ranks and not the
first25 integers; distinct p,q are adjacent iff p+q is an integer square.
Earlier provenance runs through `session132_fractran_tournament_recursion_20261031/prime_graphs/`
and `session32_square_monster_20261023/PRIME_BLOCK_CARDINALITY_HOSTILE.md`.
The old edge inventory and Ogg partitions are inherited, not re-audited.

The useful new view is to remove leaves and suppress vertices of degree2.
The core consists of five internally disjoint paths on the three ports

    a=3,       b=23,       c=47:

    Q:   3–13–23                           length2,
    A0:  3–97–47                           length2,
    A1:  3–61–83–17–47                     length4,
    B0:  23–2–47                           length2,
    B1:  23–41–59–5–11–53–47              length6.       (1)

The attached trees are `2–7–29–71–73`, `2–79`, `17–19`, `5–31`, and
`11–89`; primes37,43,67 are isolated. This is a complete reconstruction
certificate for all25 vertices and24 edges. The core has14 vertices and16
edges. Its suppression has three vertices and five edges: a triangle with
the a–c and b–c sides doubled. No degree2 contraction is being asserted to
preserve Hamiltonian paths of the original graph.

## Six cycles, one common port, three independent phases

**Kernel theorem.** The core has exactly six simple undirected cycles:

    A0+A1,    B0+B1,    Q+Ai+Bj  (i,j∈{0,1}).         (2)

Every cycle contains c=47, and the cycle-space dimension is3. The original
cycle lengths are6,8,6,10,8,12 respectively, so every cycle is even.

*Proof.* A cycle avoiding Q must use the two parallel a–c paths or the two
parallel b–c paths. A cycle using Q must choose one path from each doubled
side. These exhaust all possibilities and prove the common-port claim.
The connected core has16−14+1=3 independent cycles. The attached trees
contain no cycle. QED.

Now replace **any** of the five paths by an arbitrary positive-length path
with the same endpoints and fresh interior vertices. Equivalently, subdivide
edges repeatedly. Retain their five length parities q,a0,a1,b0,b1 and form

    u=a0+a1,    v=b0+b1,    w=q+a0+b0     in F2.       (3)

The six cycle parities are exactly

    u, v, w, w+u, w+v, w+u+v.                        (4)

**Finite compiler theorem.** The three-bit state (u,v,w) completely describes
the odd/even status of every cycle after arbitrarily many edge subdivisions.
Concatenating a path segment updates its path parity by XOR, and hence acts
on (3) by a fixed XOR vector. All eight states occur. There is no smaller
state space that can recover the labelled six-cycle parity vector for all
subdivisions, because u,v,w themselves are three entries of (4).

The number o of odd cycles can only be

    o=0 in one state,       o=3 in four states,
    o=4 in three states.                                 (5)

In particular no sequence of subdivisions can make all six cycles odd.
If u=v=0 then o=4w. Otherwise the last four coordinates in (4) have exactly
two ones, so o=2+u+v. This proves (5) and its equality cases, with no bound
on subdivision lengths. The exact retained observer is a linear functional
on the cycle space; five raw path phases have two redundant vertex gauges.

Every cycle passes through47, so the vertex-disjoint cycle-packing polynomial
for the odd cycles is simply `1+o t`. At t=2 the possible values are1,7,9.
This is an undirected-cycle packing statistic at this point, not a Hamiltonian
path formula for a tournament.

## Why a K6 overlap picture does not by itself produce H=13

The inherited H13 core has adjacency bit rows

    Q13=(18,8,3,5,14).

Its directed odd cycles, with a canonical first vertex, are

    (0,1,3), (0,4,2), (0,4,3), (1,3,2),
    (0,4,1,3,2), (0,4,2,1,3).

They meet pairwise, so their conflict graph is K6 and their odd-cycle packing
polynomial is1+6t, giving H=13 by the inherited tournament odd-cycle formula.
Their **total intersection is empty**. The prime graph's six simple cycles
also have conflict graph K6, but they share47, and (4) forbids making all six
odd by any subdivision. Thus the proposed map "list cycles, retain pairwise
overlap" is concrete and accurately records the conflict graph, while losing
two relevant coordinates: higher intersections and the cycle-space XOR law.
The missing coordinate is not recovered by a larger scalar path total.

The comparison suggests a constructive use: search for gluing operations
that change the cycle-space dependencies while preserving pairwise conflicts
and the rooted response. Subdivision alone is now completely classified
and cannot perform that task. No arbitrary tournament completion is assumed
to preserve the six-cycle list.

## A directed three-cycle carrier exposes what tournament completion adds

There is a particularly explicit positive specialization of the kernel.
Orient A0 from a to c, A1 from c to a, Q from a to b, and both B paths from
b to c. Insert one fresh vertex, denoted100, in A1 between3 and61. Keep the
other paths as in(1). The resulting directed graph on15 vertices has exactly
three directed cycles:

    A0+A1                  length7,
    Q+A1+B0                length9,
    Q+A1+B1                length13.

All are odd, all meet, and there are no other directed cycles. The proof is
the same port enumeration: A1 is the only route out of c; at a one returns
along A0 or uses Q, then one of the two B paths. Its directed odd-cycle
packing polynomial is **exactly1+3t**, taking the forbidden-looking value7
at t=2. This is a real directed graph carrying the statistic, not a formal
coefficient choice.

**Completion lemma.** Every tournament completion that retains these arcs
has at least four directed odd cycles, and consequently at least9 Hamiltonian
paths by the tournament odd-cycle formula.

*Proof.* It retains the displayed three odd cycles, all of length at least7.
Every cyclic tournament contains a directed triangle: take a shortest directed
cycle; a chord of a cycle of length greater than3 gives a shorter directed
cycle in one of its two orientations. The new triangle is distinct from all
three displayed cycles. Thus the sum of single-cycle contributions alone is
at least8, in addition to the empty collection's1; other packing contributions
are nonnegative. QED.

This supplies a small, exact realization boundary for the seven obstruction:
a cycle conflict pattern with polynomial1+3t is easy outside tournaments,
and completing the pair relation necessarily introduces an additional odd
cycle. It does not classify every possible tournament with hypothetical H=7;
the canon's full exclusion remains the inherited theorem. The new value is
the explicit connection from the user's25-prime graph to an arbitrary-length
three-bit gluing grammar, with the tournament completion cost visible.

## Reproduction and connection contract

`python3 -B check.py` independently generates square-sum edges, peels the
2-core, enumerates all simple cycles, checks all32 path-parity assignments,
and enumerates the directed15-vertex carrier's cycles. It also checks Q13
by literal HP subset DP and direct odd-cycle enumeration. These small checks
support the unbounded proofs; no large census is involved.

Source = the25-prime square-sum graph. Target = a three-port, five-path
multigraph. Map = prune trees, then suppress degree2 vertices. Preserved =
simple cycles and their intersections at branch ports. Lost = path lengths,
labels, Hamiltonian endpoints, and tree obligations. Sidecar = the explicit
five paths and attached-tree list in(1).

Source = arbitrary subdivisions of that kernel. Target = F2³. Map =(3).
Preserved = each cycle's parity, exactly and recursively. Lost = length and
integer labels. Decisive hostile = trying to make all six cycles odd.

Source = directed odd-cycle graph. Target = tournament completion.
Map = orient every formerly missing pair. Preserved = existing arcs/cycles.
Added information = forced chord cycles. Decisive test = the15-vertex carrier
with odd-cycle packing value7, whose every completion has HP total at least9.
