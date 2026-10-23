# Independent audit of the first square-sum cycle threshold

Status: **PASS for the bounded exact graph claims.** The independent checker
`audit_cycle_threshold.py` works on the full original edge sets, propagates
vertex-degree obligations, then branches on an undecided edge. It neither
imports the producer nor repeats a large fixed-cardinality subset enumeration.

## Exact cover census

A spanning 2-factor is an edge subset in which every original vertex has
degree two. In a finite simple graph it is a disjoint union of cycles, and
it is Hamiltonian precisely when connected. The exact results are:

- G31: 44 graph edges; two 2-factors, with cycle lengths 11+20 and 14+17.
  Neither is a Hamiltonian cycle.
- G32: 46 graph edges; five 2-factors, with cycle lengths 3+6+23, 3+29,
  5+27, 9+23 and 32. Exactly one is Hamiltonian.

The complete edge sets are in `cycle-certificate.json`. Only three search
nodes were needed at 31 and nine at 32. Normal and optimized runs give
identical results, with explicit checks rather than removable assertions.

Every order 1<=n<=30 has a vertex of degree less than two. The certificate
lists one such vertex at each order. Thus these orders cannot have a cycle,
G31 is excluded by its complete two-cover classification, and the explicit
32-cycle proves that 32 is the first Hamiltonian-cycle order. This does not
by itself prove existence for every n>=32; that stronger statement remains
the separately cited Gerbicz theorem.

## Why only two covers occur at 31

Propagating degree-two obligations in G31 leaves 23 forced edges and 15
undecided edges, with six other edges excluded. The optional graph is a
subdivision of K4, whose branch vertices are 1,3,4,15. Its six branch-to-branch
chains are

    1–3,
    1–15,
    1–24–12–4,
    3–6–10–15,
    3–22–14–11–5–4,
    4–21–15.

Every internal vertex of an optional chain needs exactly one optional edge.
Thus its selected edges alternate. The first five chains have odd length,
so their two branch incidences are equal. The last chain has even length,
so its two branch incidences are complementary.

Let p,q,r,s,t be the selected incidence values on the first five chains
respectively, and let u be the incidence of the last chain at vertex 4.
All six variables are Boolean. The residual degree demands at branch
vertices 1,3,4,15 are respectively 2,1,2,2, giving

    p+q+r=2,
    p+s+t=1,
    r+t+u=2,
    q+s+(1−u)=2.

Adding all four equations gives p+q+r+s+t=3. Adding the first two and
comparing then forces p=0, hence q=r=1. The remaining conditions reduce to
s=u and t=1−u. There are exactly two choices. Alternating extension along
the chains is unique in each case. Expanding the forced edges gives precisely
the 11+20 and 14+17 covers. This is an elementary two-case proof, not just a
negative Hamiltonian search result.

## The endpoint information carried by arrival 32

Vertex 32 has exactly the neighbours 4 and 17. Every Hamiltonian cycle
therefore uses both 4–32 and 32–17. Deleting 32 is a bijection from
Hamiltonian cycles of G32 to Hamiltonian paths of G31 with prescribed
endpoints 4 and 17. Conversely, adjoining 32 closes each such path.

The independent solver finds five G31 edge sets with degree one at 4,17
and degree two elsewhere. Exactly one is connected. Its edge set is the
unique 32-cycle with the two edges incident to 32 deleted. This sharpens
the explanation of the first cycle: the degree-two arrival is useful
because of its particular endpoint pair, not merely its degree.

## Final producer and surgery review

The finalized source packet in `session32-cycle-threshold` was read in full.
Its graph edges, both forcing traces, forced/allowed edge sets, complete
cycle-cover edge sets and prescribed-endpoint path match this independent
audit. Each individual forcing inference is replayed against its actual
premises, rather than accepted solely from the final counts.

The port-loss example is valid: the bare eight-block quotient cycle uses
port 5 twice while omitting its forced-path partner 22, and uses port 11
twice while omitting partner 24. Thus its apparent spanning cycle cannot
expand to a 2-factor of the original graph. This is a false-positive
Hamiltonicity example caused specifically by discarding endpoint ports.

The stated local surgery is exact. From the 11+20 G31 cover, delete
1–24, 4–12 and 8–17, then add 1–8, 12–24, 4–32 and 17–32. The resulting
edge set is the unique 32-cycle. The prose correctly separates the
hand-checkable first-cycle threshold proof from the finite-exact uniqueness
claim and from the cited all-order persistence theorem.

The producer README, checker and both certificates are pinned in the
independent cycle certificate. No producer files were modified.
