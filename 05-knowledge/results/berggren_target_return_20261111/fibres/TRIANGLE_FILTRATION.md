# Triangle chains: complete first content and an exact all-depth filtration

**PROVED in this note:** the first extension content for every ordered chain
of directed triangles and singleton gaps; the full extension filtration
when the internal gaps vanish; general upper and lower lifetime bounds.
**FINITE-EXACT:** literal small controls and the explicitly seeded probes.
The full arbitrary-gap filtration is now **PROVED** by the independently
reviewed [half-contact theorem](../arithmetic/TRIANGLE_ALL_K_PROOF.md).
Prescribed target realization remains **OPEN**. This is a family
classification, not Chan–Pak.

## 1. Source family and exact first step

Let r>=1 and g_0,...,g_r>=0. Put

```
T(g)=TT_(g0) -> C3 -> TT_(g1) -> ... -> C3 -> TT_(gr).
```

It has H=3^r and feedback vertex number exactly r. The latter follows
from the r disjoint triangles and deletion of one vertex from each.
Use the extension-content C_k of the incoming P13 theorem, credited in
[EXTENSION_CONTENT.md](EXTENSION_CONTENT.md). Then

```
C_1(T(g)) = 3^(r-2), if r>=2 and some internal g_i=0;
             3^(r-1), otherwise.                         (1)
```

Proof uses the exact quadratic extension polynomial in that companion
note. An internal pair in one triangle has coefficient `-2*3^(r-1)`.
For endpoints i,j in distinct triangle components A_p,A_q, p<q, let G
be the total singleton count strictly between these components. Then

```
K_ij = 3^(r-2) if q=p+1 and g_p=0, and 0 otherwise;
K_ji = 3^r * 4^(q-p-1) * 2^G.                            (2)
```

Indeed all intermediate triangles have ordered two-cover total12,
endpoint triangle split sums are3 each, exterior triangles contribute3,
and intermediate singletons contribute2. Forward paths can bridge no
intermediate component. Every pair with at most one endpoint in a
triangle has a coefficient divisible by `3^(r-1)`; with neither in a
triangle it is divisible by3^r. This follows by the same forced-whole /
free-two-cover component decomposition. Start/end linear terms are also
divisible by `3^(r-1)`. Hence (2) lowers the gcd to `3^(r-2)` exactly
when an adjacent triangle pair exists; its coefficient is
`-10*3^(r-2)`. Otherwise the internal coefficient and H give
`gcd(3^r,2*3^(r-1))=3^(r-1)`. The r1 case gives1 directly.

For fixed r these source classes have generating function
`z^(3r)/(1-z)^(r+1)`. The no-adjacent-triangle subfamily has
`z^(4r-1)/(1-z)^(r+1)`. Thus (1) partitions the complete family at every
order. Across all r, including r0, the ordinary class enumerator with u
marking triangle count is `1/(1-z-u*z^3)`: this is the ordered SCC word
enumerator for alphabet {singleton,C3}. These are actual ordinary classes
contributing to A000568; no converse quotient or labelled multiplicity is
inserted.

## 2. Every k-birth extension retains a quantified factor

For every gap vector and k>=0,

```
3^max(r-2k,0) divides H(U) for EVERY k-vertex extension U,
and hence divides C_k(T(g)).                             (3)
```

This direct path action proves (3), independently of the coefficient
compiler. In a Hamilton path of U, call an old triangle touched if some
vertex of it is adjacent in the path to a new vertex. At most2k triangles
are touched. Cyclically relabel the three vertices of any untouched
triangle in the Hamilton path, leaving all positions fixed. Every used
edge incident to that triangle joins old vertices. The cyclic permutation
is an automorphism of the entire old source, because all other old SCCs
see the triangle uniformly. It therefore preserves every used edge;
unused old-new edges need not be preserved. The positions adjacent to
new vertices are unchanged, so the touched set is preserved. Rotations
of distinct untouched triangles commute and act freely. Partition paths
by the touched set;
each orbit has cardinality `3^(number of untouched triangles)`, divisible
by the factor in (3). Sum the orbits. Ends of the full path and empty
exterior gaps introduce no exception.

**Correction lineage:** the first checkpoint incorrectly said untouched
triangles must occupy consecutive positions. They may be split across
old-path segments by a new-vertex detour whose neighbours lie in other
SCCs. The automorphism-on-used-edges argument above is the repaired proof.
For example TT1->C3->TT1, with triangle1->2->3->1, extends by vertex5
beating only0; the path1,4,5,0,2,3 splits the untouched triangle.
The theorem and quantitative bound are unchanged. The independent
arithmetic-lane review found this error before final synthesis.

This gives the lifetime bound `L(T)>=ceil(r/2)`. It is stronger than
testing a few random extensions: every old-new and new-new orientation
is allowed in the quantifier.

## 3. Exact full filtration without internal singleton gaps

Suppose g_1=...=g_(r-1)=0; exterior gaps remain arbitrary. Then for ALL k,

```
C_k(T)=3^max(r-2k,0),       L(T)=ceil(r/2).               (4)
```

Here are literal constructors attaining the required gcd. One triangle
extends to S4=(2,4,9,3), of H5. Two adjacent triangles extend by one
vertex beating old vertices0 and3 (mask9), giving outgoing rows recorded
in `triangle_filtration.json` and H41. Both retain every old edge.

For k<=floor(r/2), extend k disjoint adjacent pairs in this way, uniformly
joining all new block incidences. The resulting count is
`41^k*3^(r-2k)`. A competing k-extension by sinks has count3^r.
Their gcd is exactly `3^(r-2k)`, which with (3) proves equality. If r is
odd, the remaining triangle at k=ceil(r/2) extends to S4, giving
`41^floor(r/2)*5`, coprime to3. Subsequent sinks preserve that count.
Even r needs no singleton step. All exterior TT vertices survive and
do not change any displayed count. This proves every stage of (4).

Consequently increasing structural complexity has a completely known
arithmetic escape time in this family. It supplies a useful comparator
for the older endpoint-one caps with H=3^j and fvs=ceil(j/2), recovered in
[CAP_WIDTH.md](../../compensation_returns_20261103/CAP_WIDTH.md). The caps
are different graphs; their full extension filtration is not identified
with (4) by matching the number3^j.

## 4. General gap bounds and the completed matching law

Let nu be the maximum matching size in the path on the r triangles whose
edge i is available exactly when g_i=0. Extend nu disjoint available
pairs to H41, and each remaining triangle to H5. This gives the actual
upper bound

```
ceil(r/2) <= L(T(g)) <= r-nu.                            (5)
```

More generally the same disjoint constructors imply

```
3^max(r-2k,0) | C_k(T(g)) | 3^max(r-k-min(k,nu),0).       (6)
```

For the upper divisibility bound, use min(k,nu) disjoint pair births and
as many additional singleton births as available, then sink-pad unused
births. Take the gcd with the original count supplied by all-sink births.

The right side of (6) is now known to be equality for every gap vector
and every k, by the new [half-contact proof](../arithmetic/TRIANGLE_ALL_K_PROOF.md):

```
C_k(T(g))=3^max(r-k-min(k,nu),0),    L(T(g))=r-nu.        (7)
```

The proof groups coefficients by actual endpoint constraints and charges
triangle losses to old-new contacts, available singleton assignments,
terminal endpoints and zero-gap matchings. Independent analytic reviews
and local-cover/budget controls passed. The former conjecture and its
120 seeded probes were discovery evidence; they are not substituted for
that proof. The [joint enumerator](TRIANGLE_METAGRAPH.md) now classifies
the family by graph size, triangle count, matching size and exact lifetime.

## Controls

`extension_content.py` independently checks (1) on twelve varied gap
vectors through four triangles, including the15-vertex all-separated
case. `triangle_filtration.py` directly counts H41 and the explicit
filtration witnesses of manageable size, verifies their exact gcds, and
records seeded probes with all actual matrices. The shared C++ counter
is an independent subset Hamilton-path implementation; it uses neither
(1) nor (4). Run both scripts in this directory, normally and with `-O`.

Connection contract: ordered SCC graph -> extension-count gcd filtration;
map actual births and path-rotation action; preserved actual old graph
and divisibility; scalar H loses gap adjacency and escape time; sidecar
is the ordered gap vector. The smallest decisive hostile is still
H9 at two adjacent versus separated triangles.
