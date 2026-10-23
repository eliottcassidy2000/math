# Independent correction: a Hamiltonian cycle makes every edge optional

23 October 2026. Reviewer: hall_bott_recovery. **PROVED**, for finite simple
undirected graphs with at least three vertices and unrestricted Hamiltonian
path endpoints. This is an elementary consequence of a supplied Hamiltonian
cycle, not a new square-sum existence theorem.

Let C be a Hamiltonian cycle in G, and let e=uv be any graph edge.

**Inclusion.** If e lies on C, cut a different cycle edge. The resulting
Hamiltonian path contains e. If e is a chord, orient C solely to name its
successors, and delete u--succ(u) and v--succ(v). The remaining two spanning
path components have endpoint pairs {u,succ(v)} and {v,succ(u)}. Add uv to
join them. The union is connected; u,v again have degree two, their two
successors have degree one, and every other vertex has degree two. It is a
Hamiltonian path containing e.

**Omission.** If e lies on C, delete it. If e is a chord, delete any cycle
edge and use the remaining cycle path. In either case a Hamiltonian path
omits e.

Thus every graph edge is optional: its usage lies strictly between zero and
the total number of unrestricted Hamiltonian paths. Combined with the
inherited all-n>=32 square-sum Hamiltonicity theorem, this proves the complete
eventual edge-role answer: every edge of every such square-sum graph is
optional. The statement that Hamiltonicity proves only non-forcing but not
positive edge usage was false; the missing operation was inserting the
chord after cutting two cycle edges.

The endpoint qualification is essential. In the four-cycle 0--1--2--3--0
with chord 0--2, the only Hamiltonian path with prescribed endpoints 0,1 is
0--3--2--1 (up to reversing its orientation). It avoids the chord. Hence the
unrestricted theorem does not prove every edge can be used at every endpoint
port. The construction also reverses one cycle segment when expressed as an
ordered path; it is not automatically valid for directed graphs.

The producer's explicit inclusion/omission witnesses at n=32 provide finite
controls. The argument above supplies the all-size implication whenever a
Hamiltonian cycle is supplied or cited.
