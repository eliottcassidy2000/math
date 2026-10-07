# Finite algebraic witnesses and exact triangle-packing certificates

**Status:** the three named OpenAI paper theorems are **accepted premises for this session**, following the owner's instruction; their proofs are not audited here. The finite-algebraic-witness consequence is a **classical deduction**, already present in the repository and in Madore §5.4. The geometric triangle-packing reduction and certificate checks below are **PROVED elementary**; no priority claim is made. The computation is **FINITE-EXACT**, and produces a four-chromatic Moser spindle control, **not an explicit six-chromatic graph**.

## 1. Inheritance and the constructive target

[THM-4558 — residue colourings of field planes](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md) and the [second OpenAI reading, section 4](oai2_openai_math_second_reading_20261006.md) already record the real-algebraic transfer. They also give the decisive field-search hostile: the entire old Heegner compositum is four-colourable, and the named Polymath field is five-colourable. Increasing the number of candidate points in one of those fields cannot produce the requested six-chromatic witness. MISTAKE-576/579 repaired the earlier class-number roadmap; no such roadmap is reused.

The live concept board is: finite non-colourability; real algebraic coordinates with a chosen embedding; independent proof checking; equilateral-triangle conflicts; random versus optimal removal; affine volume invariants. The least-used sidecar is the exact real-root isolation interval: a formal polynomial quotient alone does not identify the intended real coordinates.

The useful output object would be a finite packet

```
one algebraic real root + coordinate polynomials + exact edge list
    + a five-colouring refutation + a six-colouring witness.
```

Such a packet can be checked without trusting the search that produced it. This session implements the same interface for three-colour refutation and four-colour witness, then adds an independently certified triangle-packing operation.

## 2. What follows from the accepted plane theorem

[The Euclidean plane is not five-colorable, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/The-Euclidean-plane-is-not-five-colorable-September-23-2026/paper.pdf) states non-five-colourability with arbitrary colour classes, in ZFC. This is the premise used here; the measurable/unrestricted transfer is not re-proved.

First apply graph compactness: if every finite subgraph of the unit-distance graph were five-colourable, propositional compactness would give a five-colouring of the entire graph. Hence some finite point set induces a non-five-colourable graph. Repeatedly delete a vertex whenever the remaining induced graph is still non-five-colourable. The final graph `G` is vertex-critical: `G-v` is five-colourable for every vertex `v`. Giving `v` a new colour shows `chi(G)=6`, rather than merely `chi(G)>=6`.

For this fixed graph on `N` labelled vertices introduce real variables `(x_i,y_i)`. Require distinct points and, for every pair,

```
(x_i-x_j)^2+(y_i-y_j)^2 = 1       if ij is an edge,
(x_i-x_j)^2+(y_i-y_j)^2 != 1      if ij is a nonedge,
(x_i-x_j)^2+(y_i-y_j)^2 > 0       always.
```

This is one finite existential formula over the rationals. It has a solution in `R`; completeness/quantifier elimination for real closed fields gives a solution in the real algebraic numbers. All coordinates lie in one finite real number field. The exact graph, its six-colouring and its non-five-colourability survive unchanged.

This is the standard compactness/real-closed-field argument. [Madore, *The Hadwiger–Nelson problem over certain fields*, §5.4](https://arxiv.org/abs/1509.07023) explicitly gives both the field transfer and the enumeration/decision consequence. ZFC supplies the classical compactness framework; it is not, by itself, a numerical bound on a witness.

### Terminating search, with no supplied practical deadline

Enumerate finite simple graphs. Test five-colourability by a terminating finite procedure. For every non-five-colourable candidate, decide the displayed real-closed-field formula. A positive decision supplies algebraic sample coordinates. Under the accepted plane theorem this search eventually succeeds. Equivalently, enumerate finite tuples of real algebraic points and check their exact unit graphs. This is a computable search with a conditional halting proof, not a displayed witness or a useful runtime estimate.

The restriction to connected critical graphs gives a bounded spatial normalization: choose an edge, place its endpoints at `(0,0),(1,0)`, and use paths to bound every coordinate in `[-(N-1),N-1]`. It does not bound `N`, coordinate degree, or the bit size of a convenient algebraic description. Restricting to a single number field is not a complete enumeration.

Two cheap necessary filters are available. A six-critical graph has minimum degree at least five. Distinct plane points have at most two common unit-distance neighbours, since two unit circles have at most two intersections. Thus

```
sum_v binom(deg(v),2) <= 2 binom(N,2).
```

Minimum degree five gives `N>=11`. Equality at `N=11` would force every degree to be five, contradicting the handshake identity on eleven vertices. Therefore any such critical witness has at least twelve vertices. This is only a weak search filter, not an estimate of the true minimum.

## 3. A small independent algebraic and colouring verifier

A finite real-algebraic tuple can be written using one primitive element `theta`: a rational monic polynomial `f`, an isolating rational interval, and rational coordinate polynomials in `theta`. The checker uses rational polynomial arithmetic and Sturm sequences. Unit edges are exact zero tests, not distance tolerances.

The implementation also permits squarefree reducible `f`. In that case it tests whether `gcd(f,q)` has the selected root, rather than treating a nonzero remainder `q mod f` as a nonzero real number. The hostile is

```
f=(t^2-1)(t^2-2),       theta in (5/4,3/2),
q=t^2-2.
```

Here `q mod f` is nonzero but `q(theta)=0`. The selected real embedding is indispensable.

For the known Moser spindle choose

```
theta=sqrt(3)+sqrt(11),   f=t^4-28t^2+64,   5<theta<6,
s=(theta^3-20theta)/16=sqrt(3),
t=(36theta-theta^3)/16=sqrt(11).
```

With `O=(0,0)`, `A=(s,0)`, `B=(s/2,1/2)`, `C=(s/2,-1/2)`, rotate `A,B,C` about `O` by cosine `5/6` and sine `sqrt(11)/6`. The seven vertices have exactly eleven unit edges: two diamonds plus the edge joining their far tips. Each diamond forces its two nonadjacent tips to have the same colour in a three-colouring; the extra edge contradicts that assignment. This is the classical spindle, not a new chromatic construction.

The search returns a finite branching refutation: a node chooses an unassigned vertex and includes a child for every colour; a leaf names a monochromatic edge whose endpoints are already assigned. A separate recursive verifier checks every branch. The spindle refutation has 139 nodes, and an explicit four-colouring is also checked. Each one-vertex deletion is independently three-colourable. Coincident vertices, incomplete branches, false conflict leaves and false colourings are rejected.

The same proof format applies to five colours, but the checker has not been supplied with any six-chromatic algebraic graph. A practical large search could replace the tree by a standard independently checked SAT proof; that changes the compression of the refutation, not the geometric premise.

## 4. An exact geometric transfer from the triangle-removal paper

[The Sharp Terminal Leave in Random Triangle Removal, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/The-Sharp-Terminal-Leave-in-Random-Triangle-Removal-September-25-2026/The-Sharp-Terminal-Leave-in-Random-Triangle-Removal-September-25-2026.pdf) gives the accepted limit `F_n/n^(3/2) -> 1/(2 sqrt(2))` in `L_2` for uniform triangle removal from `K_n`. Its initial graph and random law are part of the statement. It is not an estimate for every graph or for an optimal packing.

The following elementary structure is specific to geometric unit-distance graphs. Let `G` have distinct vertices in the real plane. Form a graph `C(G)` whose vertices are the equilateral triangles in `G`, with two adjacent exactly when they share an edge. Sharing just one vertex creates no conflict.

**Triangle-conflict lemma.** `C(G)` is bipartite and has maximum degree at most three.

For each triangle `T`, regard its vertices as complex numbers and set

```
c_T=(z_1+z_2+z_3)/3,
eta(T)=(z_1-c_T)(z_2-c_T)(z_3-c_T).
```

This is nonzero and independent of vertex ordering. Move a shared unit edge to endpoints `0,1`. Its only possible equilateral third points are `1/2 +/- i sqrt(3)/2`; the corresponding cubics are respectively `-i sqrt(3)/9` and `+i sqrt(3)/9`. Restoring the common rotation multiplies both cubics by the same nonzero cube. Thus adjacent triangles have opposite `eta` values. Within each connected component all values are `+eta_0` or `-eta_0`, giving a bipartition. No absolute real orientation is assumed. Each of the three edges can have at most one other triangle, giving the degree bound.

### Exact optimal-packing certificates

An edge-disjoint triangle packing in `G` is exactly an independent set in `C(G)`. By bipartite matching/vertex-cover duality, a maximum packing has size

```
#triangles(G) - maximum_matching_size(C(G)).               (1)
```

A short certificate consists of the complete triangle incidence list, an independent set `S`, a matching `M`, and the complementary vertex cover `V(C) minus S`, satisfying

```
|S|+|M|=|V(C)|.
```

Every packing can take at most one endpoint of each matched edge, so the matching bounds its size above by `|V(C)|-|M|`. The supplied `S` attains that bound. The verifier reconstructs **all** triangles from the source graph, checks incidences, disjointness, the matching, cover, and equality. Omitting a triangle from the packet is explicitly rejected. The exact minimum number of edges left after a maximal removal sequence is consequently

```
|E(G)| - 3|S|.                                          (2)
```

A maximum packing is maximal; its removal leaves no triangle, so it is a legitimate terminal sequence. Conversely every terminal sequence gives a maximal packing. This proves (2), with no random asymptotic premise. Polynomial-time bipartite matching supplies the certificate once the triangle list is built. This session claims a useful derived instrument, not priority for its ingredients.

The sidecar needed to return to the original graph is the triangle-to-edge incidence list. The conflict graph alone forgets metric coordinates, shared vertices that are not edges, and the rest of `G`. Removing packed triangles can lower chromatic number; it is not a non-five-colourability-preserving simplification unless a separate retained refutation proves that fact.

### Minimal greedy control

Three consecutive equilateral triangles can have conflict graph a three-vertex path. Choosing the middle triangle first is a terminal packing of size one; choosing both ends gives size two. The source graph has five vertices and seven edges, so the two terminal leaves have four and one edges respectively. Therefore a random or merely maximal packing must not be reported as optimal. The checker realizes this chain exactly in the triangular lattice.

## 5. Mahler's affine invariant does not retain unit-distance incidence

[The Mahler Conjecture for General Convex Bodies, Theorem 1.1](https://github.com/openai/math/blob/main/preprints/The-Mahler-Conjecture-for-General-Convex-Bodies-September-22-2026/paper.pdf) gives the accepted bound

```
|K| |(K-s(K))^polar| >= (d+1)^(d+1)/(d!)^2,
```

with equality exactly for simplices. Its invariant survives arbitrary invertible affine maps; unit distance does not. This already blocks identifying the inequality with a graph-colouring certificate.

A stronger finite hostile keeps even the convex hull fixed. The four corners of `[-10,10]^2` have no unit edges. Adding all seven spindle points in its interior leaves that hull unchanged and changes the induced graph's chromatic number from one to four. Both configurations have the same centered volume product `400*(1/50)=8`, above the planar Mahler constant `27/4`. The checker verifies all unit edges exactly. Hence a hull, polar, or volume-product readout needs an additional incidence coordinate before it can certify the desired graph obstruction. The theorem can control a specified convex body; it neither decides the nonconvex unit-equation system nor recovers points discarded by taking a convex hull.

## 6. Reproduction and scope

[Checker](../../04-computation/experiments/finite_algebraic_witness_20261007.py) and [saved output](finite_algebraic_witness_20261007.out):

```
python -B 04-computation/experiments/finite_algebraic_witness_20261007.py
python -B -O 04-computation/experiments/finite_algebraic_witness_20261007.py
```

The explicit universe is the seven-vertex spindle and its seven vertex deletions; all 512 subsets of one nine-point triangular-lattice patch; independent exhaustive triangle-packing checks on those subsets; all triangle-order permutations for the spindle and full patch; the three-triangle greedy hostile; a fixed-hull control; a selected-root zero-test hostile; and fifteen malformed geometry/proof/packing inputs. Normal and optimized runs use exact rational arithmetic and the same proof verifiers.

Independent audit found and repaired two receipt-type boundaries before freeze: a point with a third coordinate was silently truncated, and Python boolean/float vertex aliases could pass equality comparisons. The verifier now checks two-coordinate shape and exact integer incidence types before comparing packets. Those demonstrated inputs are retained as rejection controls; the geometric theorems were unaffected.

The positive result is a checked finite-witness format and an exact geometric packing subroutine. The accepted plane theorem ensures that exhaustive complete search eventually finds a six-chromatic algebraic witness; neither our controls nor the two additional paper statements provide its coordinates, size, or a practical deadline.
