# Dihedral boundary states beyond the ternary gasket

Status: PROVED IN NOTE; finite checks are separately scoped in `check.json`.
This extends the audited [ternary boundary theorem](../../hopf_sierpinski_port_recovery_20261023/REPORT.md), not its already more general arbitrary-permutation theorem. The improvement is the small, exact state count for the signed-affine subgroup, its even-modulus boundary, and its shared-sign CRT obstruction.

## D1. Complete connected classification

Let a finite connected graph have at least one named boundary vertex. Give every directed edge a bijection of Z/m of the form x -> epsilon*x+delta, epsilon in {1,-1}, with the inverse on the reversed edge. A section assigns vertex labels obeying these maps. Boundary vertices are ordered and distinct; there are b>=1 of them. Graphs and edge maps vary, while m and b are fixed.

For m>=3 the distinct boundary relations have exactly these forms:

* empty;
* the full graph {(x,g2(x),...,gb(x)):x in Z/m}, where every gi is signed-affine;
* for odd m, any singleton;
* for even m, any antipodal pair {a,a+(m/2,...,m/2)}.

Consequently their exact number is

```
odd m>=3:  1 + m^b   + (2m)^(b-1),
even m>=4: 1 + m^b/2 + (2m)^(b-1).
```

For m=2 the two sign choices are the same permutation, and the answer is instead `1+2^(b-1)`: empty and the graphs of independent translations. Each allowed boundary tuple has exactly one interior completion.

Proof. Root at the first boundary vertex. Tree transport expresses every label as a signed-affine function of the root label x. The possible x are exactly the common fixed points of the based cycle holonomies. A nonidentity translation has no fixed point. Two distinct reflections multiply to such a translation. Thus a nonempty root fixed set is either all of Z/m or the fixed set of a single reflection, `2x=delta`. This is a singleton for odd m; for even m it is empty or one antipodal pair. When m=2 every nonidentity transport is a fixed-point-free translation. Conversely tree transports realize all full graphs. A cycle with reflection holonomy realizes every stated proper nonempty set, and a cycle with nonzero translation realizes the empty set. Subdivision makes all these constructions simple graphs. Connectedness and root determination prove uniqueness.

The count has no quotient by simultaneous relabelling: boundary labels are observed. Normalize the first coordinate to x. For m>=3 all 2m signed-affine permutations are distinct. There are m^b singleton tuples, or m^b/2 antipodal pairs, respectively.

These relations remain realizable using only reflection edges `R_s(x)=s-x` with s a square modulo m. Indeed R_0 and R_1 generate all signed-affine maps: R_1 R_0 is translation by 1. Replacing an edge by a word in these reflections and padding a cycle by R_0 R_0 preserves transport and gives a simple subdivided graph. This is an abstract square-labelled constraint network. It is not a claim that all its residues lift to distinct positive integer vertices of a literal square-sum graph.

## D2. Exact triangular gluing and minimality

For three relations on `(a,x,z)`, `(x,b,y)`, `(z,y,c)`, define the parent relation by existentially eliminating x,y,z. The glued support is connected, so D1 classifies the result and every admitted outer tuple has exactly one inner completion. Thus these finite state spaces are closed under every finite-depth Sierpinski triangular recursion supplied with connected cells.

The number of three-boundary states is 64 at m=3 (recovered), 97 at m=4, 226 at m=5, 540 at m=7, 17101 at m=30, and 44101 at m=42. They are much smaller than the old arbitrary-permutation state counts because the allowed transports are restricted.

These are also minimal exact states for section-existence queries under arbitrary compatible gluing contexts. If two relations differ, choose a tuple in their symmetric difference. For odd m a context pins it by reflection holonomy and boundary transports. For even m pin its antipodal pair, which is either wholly in or wholly outside each relation. For m=2 pin its translation orbit. Fresh subdivided edges avoid duplicate-edge conflicts. Therefore every distinct pair of states is distinguished by an allowed context. Empty hidden components and disconnected cells fall outside the theorem.

## D3. What survives on the circle and under CRT

For m>=3, embed the labels as m-th roots of unity. The same signed-affine maps act on the circle as `z -> exp(2*pi*i*delta/m)*z^epsilon`. Nontrivial rotations have no fixed point; a reflection has two. For odd m every reflection has exactly one discrete fixed point. Hence discrete and circle section *existence* agree, but counts do not. For even m a reflection with odd delta has two circle fixed points and no discrete fixed point; a three-edge cycle with that holonomy is a counterexample. At m=2 the sign is lost by the permutation action itself, so a circle lift needs the original sign as additional data.

Chinese remaindering does not make the local signs independent. A global map has the same epsilon in every odd-prime coordinate. Thus the actions for m=30 and m=42 have index two in the products `C2 x D3 x D5` and `C2 x D3 x D7`, respectively (Dq has order 2q). As a hostile example, the product of a free Z/2 coordinate, a free Z/3 coordinate and one pinned Z/5 coordinate has six root values. It cannot be a global mod-30 section set, whose nonempty sizes are only 2 and 30. Similarly a free 2 and 3 coordinate with a pinned 7 coordinate is impossible globally mod42. The missing data are the shared sign/holonomy, not more local cardinalities.

## Scope and next question

This completely classifies the stated connected dihedral boundary category, not arbitrary tournaments, their Hamiltonian paths, or the whole metagraph. The old native-port audit already exhibited authentic tournament matching systems with zero sections: an empty compatibility relation is not a rejection criterion for a tournament. The useful next construction must carry the Hamiltonian endpoint and cycle-defect state alongside these compatibility states.
