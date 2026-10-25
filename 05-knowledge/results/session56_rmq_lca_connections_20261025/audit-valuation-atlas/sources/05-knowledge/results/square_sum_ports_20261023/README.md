# Square-sum paths: births, actual endpoints, and forbidden-edge guards

**PROVED:** the endpoint-marked arrival recurrence and cut-budget inequality
below, for arbitrary finite simple graphs; a finite structural certificate
that no Hamiltonian path of \(G_{25}\) uses edge \(1--15\).
**FINITE-EXACT:** all 896 Hamiltonian paths, up to reversal, in \(G_n\) for
\(1\le n\le32\), their individual edge usages and endpoint pairs; exhaustive
cut searches in the explicitly listed universes.
**PROVED:** every edge of a finite simple Hamiltonian graph lies in some
Hamiltonian path and is omitted by another.
**COROLLARY with a CITED dependency:** every edge of \(G_n\) is optional
among its Hamiltonian paths for all \(n\ge32\), using the inherited
Gerbicz all-order Hamiltonian-cycle theorem.
**OPEN here:** exact all-size usage counts, fixed endpoint constraints,
and compatibility of several prescribed edges.
No new claim about Collatz or tournament enumeration follows from this census.

Here \(V(G_n)=\{1,\ldots,n\}\), and distinct \(a,b\) are adjacent exactly
when \(a+b\) is a square. A path and its reversal count once. The singleton
path at \(n=1\) counts once. A graph with no Hamiltonian path receives the
status **no_paths**: we do not exploit vacuous quantification to call its
edges forced.

## Inheritance and the genuine addition

The starting point is
[Session32 square-sum README](../session32_square_monster_20261023/square-sum/README.md):
the arrival degree, the marked degree-two ear correspondence, endpoint-aware
chain compression, and the corrected cut/degree arguments. Its companion
[PORT_TRANSFER.md](../session32_square_monster_20261023/ports/PORT_TRANSFER.md)
already proves the general named-port spanning-forest transfer, gives all
ten paths at size 25, and exposes the degree-correct but disconnected
path-plus-cycle hostile. The earlier
[square-sum Hamiltonicity note](../collatz_mod6_20260922_w6_square_sum_hamiltonicity.md)
already records the counts through 36. These counts are independently recovered,
not claimed as discoveries.

The closest hostile is the retracted inference that both edges at every
degree-two vertex are forced. They are forced only after that vertex is
certified not to be an endpoint. In \(G_{23}\), a path ending at 22 omits
\(14--22\); cuts.json retains a literal witness.
The least-used sidecar in this extension is the pair of *other* endpoints
left after a newcomer joins two components.

The additions are the complete edge timelines through 32, a marked
subset-convolution implementation, minimal tight-cut certificates and their
explicit failure boundary, and the hybrid proof for \(1--15\) at size 25.
The general cut inequality is an elementary refinement of the inherited
standard cut obstruction, not a priority claim about graph theory.

## The proposed groups describe births, not permanent edge roles

The new neighbors of \(n\) number
\[
d_{\rm new}(n)=\lfloor\sqrt{2n-1}\rfloor-\lfloor\sqrt n\rfloor.
\]
The first eight proposed edges are exactly the birth edges at traceable
sizes \(15,16,17,23,25\), through 25. Sizes15,23,25 are **onsets**:
the previous graph is not traceable, and a two-path join creates a path.
Sizes16,17 extend an existing path at an endpoint.
The second fourteen edges are the degree-two births at
\(13,14,19,20,21,22,24\), where the new graph has no Hamiltonian path.
This classification was already recognized in the Session32 inheritance;
the exact usage inventory shows why it cannot be promoted to a permanent
forced/forbidden partition.

| First group edge | Birth | Uses among the ten \(G_{25}\) paths |
|---|---:|---:|
| 13--23 | 23 | 9 |
| 2--23 | 23 | 10 |
| 11--25 | 25 | 10 |
| 24--25 | 25 | 10 |
| 1--15 | 15 | 0 |
| 10--15 | 15 | 10 |
| 8--17 | 17 | 10 |
| 9--16 | 16 | 10 |

| Second group edge | Birth | Uses among the ten \(G_{25}\) paths |
|---|---:|---:|
| 3--13 | 13 | 2 |
| 12--13 | 13 | 8 |
| 2--14 | 14 | 8 |
| 11--14 | 14 | 2 |
| 6--19 | 19 | 10 |
| 17--19 | 19 | 10 |
| 5--20 | 20 | 10 |
| 16--20 | 20 | 10 |
| 4--21 | 21 | 10 |
| 15--21 | 21 | 10 |
| 3--22 | 22 | 9 |
| 14--22 | 22 | 10 |
| 1--24 | 24 | 4 |
| 12--24 | 24 | 6 |

The corrected third group is \(7--18,1--3,4--12,2--7\), with respective
usages 10,7,6,1. **Input correction:** the user confirmed that the original
\(3--12\) was a typo and should be \(4--12\). The corrected pair is used
throughout the inventory; its sum is 16. This preserves the correction
provenance rather than treating the old typo as an unresolved objection.

The listed valid edges leave six edges of \(G_{25}\) unlisted:
\[
1--8,\quad3--6,\quad4--5,\quad5--11,\quad6--10,\quad7--9.
\]
Thus the supplied lists identify 26 of the 32 edges; the third group's
examples do not exhaust the remainder.

![All G25 edges, placed along one literal Hamiltonian path](G25_edge_roles.svg)

## Complete finite inventory and loss of monotonicity

There is one path at each of 1,15,16,17, and none at the other sizes below 23.
All remaining counts in this experiment are:

| \(n\) | Paths | Edges | Forced | Sometimes | Never |
|---:|---:|---:|---:|---:|---:|
| 23 | 3 | 28 | 20 | 4 | 4 |
| 24 | 0 | 30 | — | — | — |
| 25 | 10 | 32 | 14 | 17 | 1 |
| 26 | 12 | 34 | 12 | 21 | 1 |
| 27 | 35 | 36 | 7 | 29 | 0 |
| 28 | 52 | 38 | 7 | 30 | 1 |
| 29 | 19 | 40 | 15 | 22 | 3 |
| 30 | 20 | 42 | 15 | 24 | 3 |
| 31 | 349 | 44 | 4 | 39 | 1 |
| 32 | 392 | 46 | 0 | 46 | 0 |

The never-used sets are
\[
\begin{array}{c|l}
15,16,17&\{1--3\}\\
23&\{1--15,3--6,3--13,4--5\}\\
25&\{1--15\}\\
26&\{1--3\}\\
27&\varnothing\\
28&\{7--9\}\\
29,30&\{2--7,5--20,7--9\}\\
31&\{7--9\}\\
32&\varnothing .
\end{array}
\]
In particular, the edge born in the first onset, \(1--15\), is forced
there, forbidden at 25, and used again by 27.
Every edge already present at 25 has an occurrence witness by 27.
Consequently none of those edges is “never used at any later size.”
The finite census alone does not prove an eventual claim, but the following
cycle argument does, when combined with the inherited all-order cycle
theorem. The verification.json file records exact usages, not just roles.

**Cycle surgery lemma.** Every edge of a finite simple graph with a
Hamiltonian cycle is optional among its Hamiltonian paths.
For omission, cut a cycle edge \(e\) itself; a noncycle edge is omitted by
any cycle cut. For inclusion of a cycle edge, cut a different cycle edge.
For inclusion of a chord \(uv\), orient the cycle and delete
\(u\,\mathrm{succ}(u)\) and \(v\,\mathrm{succ}(v)\).
The two remaining path components have \(u,v\) on different components,
each as an endpoint. Add the chord \(uv\), obtaining one spanning path.
These constructions prove both inclusion and omission, without enumerating
Hamiltonian paths.

In particular, all 46 edges at 32 have both kinds of explicit witnesses in
cuts.json. The consumer checks every cycle-edge/chord position at sizes
3 through 64, and all 549 graphs containing the fixed labelled cycle at
sizes 3 through 6, covering every subset of available chords.
The spanning-cycle hypothesis matters: a triangle 1,2,3 with pendant
edges 1--4 and 2--5 has the unique Hamiltonian path 4--1--3--2--5; its edge
1--2 is forbidden. A Hamiltonian path alone is insufficient.
Nor does the lemma preserve prescribed endpoints: in the four-cycle,
the only path from 1 to 2 is 1--4--3--2 and cannot use edge 1--2.

**Correction lineage.** An earlier draft of this packet falsely said that
a Hamiltonian cycle did not itself imply inclusion of every noncycle edge.
It considered only deleting one cycle edge and missed the two-cut,
one-chord surgery above. That false boundary is superseded by this lemma,
its explicit witnesses, and the small chord controls. The actual finite
inventories never depended on the false statement.

The inherited Gerbicz result cited in the Session32 and earlier square-sum
notes supplies a Hamiltonian cycle for every \(n\ge32\); the primary
[OEIS A090461 comment](https://oeis.org/A090461) is Gerbicz's January 21, 2018
statement. Consequently **every edge is optional at every size \(n\ge32\)**.
Equivalently, for every distinct pair of positive integers \(a,b\) with
\(a+b\) square, the fixed edge \(ab\) is optional for **every**
\(n\ge\max(32,a,b)\). This quantifies the behavior of each edge as the graph grows.
Together with the complete finite inventory below 32, this gives the
qualitative never/forced/sometimes classification at every size, relative
to that cited theorem. The exact usage numerators remain a separate problem.

## An arrival recurrence that retains actual ports

Let \(H\) have at least two vertices, and adjoin \(v\) with neighbor set \(N\).
Write \(h(H)\) for its number of unoriented Hamiltonian paths and \(E_H(a)\)
for those ending at \(a\).
For distinct \(a,b\), let \(C_H(a,b)\) count spanning covers by **exactly
two paths**, with \(a,b\) in different components and exposed as the chosen
joining endpoints. Singletons are permitted.
Then
\[
h(H+v)=\sum_{a\in N}E_H(a)+
       \sum_{\{a,b\}\subseteq N}C_H(a,b).
\tag{1}
\]
The first term has \(v\) as an endpoint. Deleting an interior \(v\) gives
the second term. Conversely, append \(v\), or join the two named ports
through \(v\). These inverse operations are literal bijections. There
is no factor of two: the unordered neighbor pair is counted once.

For the individual new edge \(va\),
\[
U_{H+v}(va)=E_H(a)+\sum_{b\in N\setminus\{a\}}C_H(a,b).
\tag{2}
\]
For an old edge \(e\), the same partition applies with each \(E,C\)
restricted to paths/covers containing \(e\). Thus an edge-marked state
propagates the actual usage inventory.

To retain future endpoint feasibility, let \(E_{H[A]}(a,r)\) count paths
on the induced set \(A\) with named endpoints \(a,r\). A singleton
has \(a=r\) and count 1. With the pair \(a<b\) fixed,
\[
C_H(a,b;\{r,s\})=
\sum_{\substack{A\subset V(H)\\a\in A,\ b\notin A}}
\ \sum_{\{u,w\}=\{r,s\}}
E_{H[A]}(a,u)\,E_{H[V(H)\setminus A]}(b,w).
\tag{3}
\]
The inner sum assigns the other two endpoints to their actual components.
The condition \(a\in A,b\notin A\) prevents double counting the split.
The recurrence.py script enumerates fixed-start simple paths on subsets and uses
(3), independently of the new full-path catalog, to predict its entire
endpoint-pair distribution and each newcomer-edge usage.

| Arrival \(n\) | Endpoint extensions | Two-path joins |
|---:|---:|---:|
| 15 | 0 | 1 |
| 16 | 1 | 0 |
| 17 | 1 | 0 |
| 23 | 0 | 3 |
| 24 | 0 | 0 |
| 25 | 0 | 10 |
| 26 | 2 | 10 |
| 27 | 3 | 32 |
| 28 | 5 | 47 |
| 29 | 0 | 19 |
| 30 | 2 | 18 |
| 31 | 20 | 329 |
| 32 | 59 | 333 |

At29 the previous graph has 52 paths but neither new neighbor,7 or20,
can be an old endpoint. The entire count 19 comes from the two-cover
state. A recurrence using only the old path count, or even only its
endpoint counts, is missing that state.

## A cut budget, its exact failure boundary, and its repair

For a Hamiltonian path \(P\) and \(S\subsetneq V(G)\), let \(e_P(S)\)
be its selected edges with both endpoints in \(S\), and \(b_P(S)\)
the number of its two endpoints in \(S\). Direct vertex/edge counting
in the forest \(P-S\) gives
\[
c(P-S)=|S|+1-b_P(S)-e_P(S).
\]
Since \(G-S\) can only join those components,
\[
\boxed{b_P(S)+e_P(S)\le |S|+1-c(G-S).}
\tag{4}
\]
Equality holds exactly when extra graph edges join no two components of
\(P-S\). A negative right side forbids all Hamiltonian paths.
A zero right side forbids both endpoints in \(S\) and every selected
edge inside \(S\).

The following are **minimum-cardinality** tight cuts containing the
specified edge, found by exhaustive increasing-size search:

| \(n\) | Forbidden edge | Smallest tight cut \(S\) |
|---:|---|---|
| 15 | 1--3 | {1,3} |
| 23 | 1--15 | {1,4,6,7,13,14,15} |
| 23 | 3--6 | {3,4,6,7,14,15} |
| 23 | 3--13 | {3,4,7,13,14} |
| 23 | 4--5 | {3,4,5,7,14} |
| 26 | 1--3 | {1,3,5,7,10,12,14,23} |
| 29,30 | 7--9 | {7,9,20} |

There is **no** such tight cut for \(1--15\) at 25, \(7--9\) at 28,
\(2--7\) or \(5--20\) at 29 or30, or \(7--9\) at 31. These are exhaustive
negative statements in those finite graphs, not a failure to find a
small cut. A tight cut has
\(|S|+1=c(G-S)\le n-|S|\), so only
\(|S|\le\lfloor(n-1)/2\rfloor\) need be checked. The search checks
460,095,104 candidate cuts in total; all positive certificates are
independently replayed in Python.

For \(G_{25}\), even the strengthened endpoint version of (4) is
insufficient by itself. Assume \(1--15\) is selected. Vertex18 is a leaf,
so one endpoint is18. For each possible other endpoint \(t\), every cut
must satisfy
\[
\mathbf1_{18\in S}+\mathbf1_{t\in S}
 +\mathbf1_{\{1,15\}\subset S}
\le |S|+1-c(G_{25}-S).
\tag{5}
\]
Testing all 21,977,516 sets with \(|S|\le13\) eliminates21 of the 24
endpoint choices, leaving **8,10,21**. No larger set can violate (5):
its left side is at most3, while for \(|S|\ge14\) the right side is
at least \(2|S|-24\ge4\).

Now keep the actual endpoint labels and propagate the necessary degrees:
1 at 18 and \(t\), 2 everywhere else. Initially select \(1--15\).
At each vertex, exclude undecided edges once the target degree is met;
select all undecided edges when they are all needed.
This forces a degree deficit at vertex 21 when \(t=8\).
For \(t=10\) or21 it forces the proper cycle
\[
2--14--22--3--13--23--2,
\]
which cannot lie in a Hamiltonian path. Hence all 24 cases are eliminated.
The guarded_forbidden.json file records the 21 explicit cuts and the complete
three propagation traces. Its consumer needs no Hamiltonian-path search.

This is the useful common mechanism with the tournament port work:
degree totals and scalar inequalities do not determine global wiring.
Named endpoints plus the connectivity/cycle guard carry the missing
information. The connection preserves feasibility of a prescribed
spanning structure; it is not an asserted bijection between square-sum
graphs and tournaments.

## Reproduction and scope

From the repository root, with Python 3 and g++:

    python 05-knowledge/results/square_sum_ports_20261023/reproduce.py

The driver compiles into a temporary directory and reproduces the literal
path and cut universes byte-for-byte after replacing CRLF by LF. The same
explicit protocol applies to the manifest and input-path SHA-256 hashes;
no other byte normalization is performed. This makes Windows-generated
text and Git's LF checkout representation authenticate identically. A
CRLF round-trip control is checked for every hashed file.
It then replays all Python consumers
normally and with -O; they use explicit checks, not removable asserts.
The path census is independent in two ways:

1. paths.cpp: fixed-start vertex-order DFS, only sound leaf/connectivity
   continuation pruning, canonical reversal representatives.
2. verify.py: choose endpoint pairs, propagate degree 1/2 constraints,
   branch on edge subsets, and reject cycles or disconnected allowed graphs.
   It compares full path sets, not just the counts.

The consumers also check edge-incidence and endpoint-incidence sums.
The recurrence.py script is a third independently constructed route to the arrival
endpoint and new-edge counts. The manifest.json file freezes all core sources,
literal fixtures and outputs. The draw.py script produces the native SVG from the
exact size 25 path and usage data.

The earlier repository already cites the all-\(n\ge25\) path existence
result and all-\(n\ge32\) cycle result. This packet does not re-prove those
literature claims. Its all-size qualitative corollary uses the general
cycle surgery lemma and the cited cycle theorem, rather than extrapolating
the finite inventory. The next useful target is a named-port rule for
**prescribed endpoints or several prescribed edges together**, where
one-edge optionality no longer settles compatibility. Equation (4) remains
a cheap necessary filter; endpoint-guarded connectivity retains the
missing obligation.
