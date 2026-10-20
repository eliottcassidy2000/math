# Boundary hubs and the cycle content of Hamiltonian relocation

Status: **PROVED** exact identities and sufficient detour criteria;
**FINITE-EXACT** bounded mechanism tests. No universal black compiler is
claimed. The invariant lost by a scalar count is the boundary edge or the
actual cycle support, depending on the move.

## 1. Open paths and unordered boundary-cut hubs

Let P=(p_0,...,p_(n-1)) be a directed Hamiltonian path in T, n>=3, with
p_0->p_(n-1). Its complement S_P=K_P(T) contains the directed Hamiltonian
cycle formed by P and p_(n-1)->p_0. Cut this cycle at the P edge p_i p_(i+1)
and let U_i=(p_(i+1),...,p_(n-1),p_0,...,p_i). Then

```
K_(U_i)(S_P) = T XOR {p_0 p_(n-1), p_i p_(i+1)}.          (1)
```

Here XOR reverses the two unordered pairs. Proof: twice taking the whole
tournament complement cancels; the two path edge sets differ precisely in
the closing edge and the cut edge. This is a literal labelled identity.

Consequently, two open paths P,Q whose unordered pairs
{boundary edge,cut edge} coincide have a common hub. If both cut rotations
are black in their respective partners, S_P and S_Q are joined by two black
moves. Crucially, the boundary edge of P may be the cut edge of Q: retaining
only ordered roles discards legal connections.

Construct a bipartite incidence graph with open paths on one side and
two-edge hub supports on the other, retaining a branch only when its cut
rotation is black. If a marked blue path P connects to a black open path Q
by r hub transitions, the associated blue metagraph edge has a black detour
of length at most 2r+1. The final step is the black T--S_Q move; traversing
the hub chain in reverse reaches S_P. Coincident quotient nodes cause no
problem for connectivity, but the argument does not classify individual
multigraph edge bridges.

The older [fixed-boundary hub](../defect_transitions_20261020/blue.md) is the
special case in which the two roles never exchange. The new sidecar is the
unordered pair of physical vertex pairs. This is more informative than just
the size of the Hamiltonian-path set.

`boundary_transport_probe.py` exhausts all standard-blue presentations through
order seven, choosing the open side, requiring strong mixed endpoints, and
removing merged loops. It finds 1 of 3 nonloop presentations at order five,
11 of 25 at order six, and 221 of 239 at order seven. Of the last 221,
219 use a length-three detour and two use length five. Combined with the
independent tests in [detours.md](detours.md), the unresolved mechanism cases
number 2,7,1. The order-seven residual is T_3, now repaired by the explicit
family certificate. These counts concern presentations, not isomorphism
classes. Failure of this incidence test does not imply a genuine bridge.

## 2. Interior relocation is a pair of directed triangle reversals

Suppose P=(A,x,B,C) and Q=(A,B,x,C) are directed Hamiltonian paths in the
same tournament T, with A,B,C nonempty and |B|>=2. Set

```
a=last(A), b=first(B), z=last(B), c=first(C).
```

All five displayed vertices are distinct. The complement S_P=K_P(T)
contains the two directed triangles

```
a -> x -> b -> a,       x -> z -> c -> x.                 (2)
```

Their edges are disjoint, and the triangles meet only at x. Moreover S_Q
is obtained by reversing exactly these two triangles.

Indeed, the removed P joins are a->x, x->b, z->c, while the added Q joins
are a->b, z->x, x->c. The former are retained by K_P and the latter are
reversed, giving (2). All other consecutive edges of the two paths coincide.
This proves both directedness and the exact six-edge support.

Each directed cycle reversal preserves every outdegree. Therefore S_P,S_Q
have identical labelled score vectors and equal directed-triangle counts,
using

```
c_3(S) = binom(n,3) - sum_v binom(d_S(v),2).
```

The move also changes an even number of arcs. These observations reconcile
interior word relocation with score-preserving cycle moves, including the
odd-order parity obstruction to compiling a *single* triangle reversal.
They do not supply a Hamiltonian-path-reversal factorization whose moves
are all black.

This recovers a precise image of the operation from
[Hamiltonian relocation connectivity](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt).
The source is a pair of legal HP words with shared endpoints; the target is
a pair of edge-disjoint directed triangles in K_P(T). The map preserves
the literal final tournament and score vector. Keeping only the triangle
count loses the vertex x, the two supports, and the decomposition into path
blocks. All of these are needed for a prospective black compiler.

## 3. Sharp boundary and next obligation

The two-block theorem in [detours.md](detours.md) is a working compiler for
a different, explicitly controlled rerouting. Its endpoint-preserving
internal paths are exactly the data a purely combinatorial relocation
argument leaves unproved. The [double-return family](double_return_counts.md)
shows why finding any arc-disjoint replacement is insufficient: the unique
replacement is always blue.

The remaining structural question is whether every strong-to-strong mixed
blue edge admits a boundary hub, a square, an endpoint repair, or some
additional black compiler. No genuine mixed-to-mixed blue bridge is produced
here, and all-order black-core connectivity remains **OPEN**.

## 4. Reproduction

```
python3 05-knowledge/results/strong_blue_games_20261020f/boundary_transport_probe.py
python3 -O 05-knowledge/results/strong_blue_games_20261020f/boundary_transport_probe.py
python3 05-knowledge/results/strong_blue_games_20261020f/transport_identity.py
python3 -O 05-knowledge/results/strong_blue_games_20261020f/transport_identity.py
```

The identity checker independently visits every tournament with standard
directed HP on orders five and six and every legal relocation with the
stated nonempty blocks. It checks all six changed edges, directedness of
both triangles, the entire score vector, and the triangle count by direct
triple enumeration. It includes multiple orientations, without a blue,
strong, or isomorphism filter; its scope is the identity, not black detours.
