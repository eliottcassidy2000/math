# Graph/Pythagorean inheritance and a boundary-state right triangle

**Status: PROVED (elementary propositions below); CITED (the explicitly
identified Bessy--Havet counting bound); FINITE-EXACT (five small controls).
No Kuratowski, snark-minor, or general Tutte-flow theorem is proved here.**

Session: graph-pythagorean-20261018. This lane recovered existing mechanisms
before testing the proposed graph/primitive-triple dictionary. No theorem ID
is reserved and no novelty claim is made.

## 1. Inheritance, correction, and least-used sidecar

- **Closest proved mechanism:**
  [THM-261, Petersen root orthogonality](../../01-canon/theorems/THM-261-petersen-root-orthogonality.md).
  Map each edge `{i,j}` of `K_n` to the positive root `e_i-e_j` of `A_(n-1)`.
  For distinct edges, orthogonality is exactly disjointness of supports. Its
  graph is `KG(n,2)`; `n=5` gives Petersen. This is an actual representation
  map, although it does not assign a Pythagorean triple to a graph.
- **Canonical hostile / correction:**
  [MISTAKE-507](../../01-canon/MISTAKES.md) separates shared pair indexing
  from identification of `so(n)` with `sl(n)`, and cubic three-edge-coloring
  from a `Z/3Z` flow. The correct color-flow group is `F_2^2`; cubic
  nowhere-zero `Z/3Z` flows exist exactly for bipartite graphs. `K4` is
  three-edge-colorable and has no such three-flow. Historical arithmetic
  prime/snark analogies did not supply a factorization map.
- **Corrected near miss:** the
  [2026-09-25 Kuratowski/Tait note](procgen_kuratowski_20260925_tait_kempe_triples.md)
  explicitly distinguishes genuine `Delta-Y` moves from numeral triples.
  `Delta-Y(K5)=K3,3+e`, whereas the Petersen family is the `Delta-Y/Y-Delta`
  class of `K6`. Its source pins for several classical theorems were
  secondary, so those pins are not promoted here.
- **Least-used sidecar:**
  [THM-4113, half-Kempe atlas](../../01-canon/theorems/THM-4113-maximal-noncrossing-half-kempe-atlas.md)
  restores cyclic boundary order to root disjointness. On five boundary
  positions it removes a crossing `C5` from the 15 Petersen edges, leaving
  10 planar half-Kempe states. Orthogonality alone forgets planar order.
- **Exact gluing mechanism:**
  [THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md)
  counts three-edge-colorings as an inner product of ordered boundary
  extension vectors. Petersen's adjacent-vertex cut has disjoint supports
  of sizes 12 and 9, with multiplicities 1 and 2. The analogous `K4` cut
  has inner product 6. This is the natural place to seek a right triangle.
- **Another operation carrier:**
  [THM-4529, Petersen/Heawood Paley coordinates](../../01-canon/theorems/THM-4529-petersen-and-heawood-families-in-paley-coordinates-anti-circulant-all-odd-tournaments.md)
  realizes Petersen by four edge-disjoint `Delta-Y` moves from `K6`, and
  Heawood by seven from `K7`. It explicitly refutes transport of
  Hamiltonian-path parity along the family. The Fano incidence labels,
  rather than a shared numeral seven, are the surviving sidecar.

The concept board for this lane was: root support orthogonality; planar
cyclic order; ordered color-extension gluing; maximum-cut frustration;
Pell/Euclid arithmetic. The first three connect by explicit incidence maps.
The fourth produces the invariant below; the fifth remains a separate
arithmetic chart, with a precise incompatibility at `3:4:5`.

## 2. A choice-free right triangle for every finite simple cubic graph

Let `G` be a nonempty finite simple cubic graph, possibly disconnected. Put

```text
n = |V(G)|,
s = |E(G)| - MaxCut(G),
C = number of proper edge-colorings with three labelled colors,
M = 6^n / 3^s,
q(G) = C^2 / M = C^2 3^s / 6^n.
```

Here `s` is also the minimum number of edge deletions needed to make `G`
bipartite. Indeed a cut gives such a deletion set; conversely a bipartition
after deletions extends to a cut whose uncut edges are among those deleted.

**Proposition 1 (maximum-cut boundary norm).** Every maximum cut `X|Y`
has nonzero boundary extension vectors `f,g`, their entries are 0 or 1,
and

```text
<f,g> = C,                   ||f||^2 ||g||^2 = M.
```

The cut edges are ordered and use one common three-color gauge.

*Proof.* If a cubic vertex had at least two neighbors on its own side,
moving it would enlarge the cut. Thus each induced side is a disjoint union
of isolated vertices and single edges. Suppose side `X` has `i_X` isolated
vertices and `j_X` edges. An isolated vertex has six proper assignments to
its three boundary edges. A single internal edge has twelve: choose its
color in three ways, then order the two other colors at each endpoint in
two ways. The boundary determines every internal-edge color uniquely.
Consequently

```text
||f||^2 = 6^i_X 12^j_X,      ||g||^2 = 6^i_Y 12^j_Y.
```

Since `j_X+j_Y=s` and `i_X+i_Y=n-2s`, their product is
`6^(n-2s)12^s=6^n/3^s`. In particular both vectors are nonzero and the product
does not depend on the chosen maximum cut. The gluing bijection of
THM-4116 gives the inner product. QED.

**Proposition 2 (native boundary-state triangle).** In the tensor square of
the boundary-word space set

```text
U=f tensor g,                V=g tensor f,
P=U+V,                      Q=U-V.
```

Then

```text
<P,Q>=0,
||P||^2 = 2(M+C^2),
||Q||^2 = 2(M-C^2),
||P+Q||^2 = 4M.
```

Thus the normalized squared lengths of the two perpendicular legs are
`(1+q)/2` and `(1-q)/2`. They agree exactly when `C=0`.

*Proof.* `||U||^2=||V||^2=M` and `<U,V>=<f,g>^2=C^2`.
Expand the three norms. Cauchy--Schwarz gives `0<=q<=1`; no square root
approximation is used. QED.

This realizes the user's isosceles intuition in a genuine inner-product
geometry: a cubic graph with no three-edge-coloring gives the isosceles
endpoint. For a graph satisfying the usual additional snark hypotheses,
this is exactly its snark obstruction. For unrestricted cubic graphs it
also detects bridge obstructions and all other failures of colorability.

**Connection contract.** The source is the graph with an ordered maximum
cut and common color labels; the target is the two-dimensional span of
`P,Q`. The zero/nonzero coloring predicate survives. The final scalar
forgets boundary words, which colorings extend, component structure, and
planar cyclic order. Those must remain sidecars for any operation that
tries to prove a minor theorem. The tensor step squares the original
overlap, but loses no sign information here because coloring counts are
nonnegative. This is not an easier coloring algorithm: computing the
overlap still contains the original coloring problem.

The construction uses real lengths with rational squared lengths. A
nondegenerate isosceles right triangle has irrational hypotenuse/leg ratio,
so it cannot itself be a primitive integral Pythagorean triple.

## 3. Sharp range and a stronger asymptotic failure

**CITED input.** Bessy--Havet,
[*Enumerating the edge-colourings and total colourings of a regular graph*](https://www.lirmm.fr/~bessy/publis/3aretecol.pdf),
author-hosted manuscript dated 2011-11-05, checked 2026-10-18. Theorem 5
gives `C<=3*2^(n/2)` for a connected cubic loopless multigraph. Lemma 8
improves this to `C<=(9/4)*2^(n/2)` for connected simple cubic graphs.
Their proof uses an `s,t` ordering, grouping the first encountered cycle
before continuing the branching count. Parallel edges are allowed in the
general theorem but excluded in the sharpened lemma.

The color-count convention was checked directly: section 2, printed page 2,
defines `c_k(G)`; Theorem 4's proof on page 3 gives six choices at the first
cubic vertex, so colors are labelled. Lemma 8 on page 5 states:
“If G is a connected cubic simple graph, then c_3(G) <= (9/4) 2^(n/2).”
No additional 2-connectivity hypothesis occurs in that statement. Theorem
5's proof on page 4 explains why a colorable connected regular graph is
2-connected; the remaining case has count zero.

**Proposition 3 (sharp simple-graph range, with the cited input).**

```text
0 <= q(G) <= 1/4,
q(G)=1/4 iff G is K4.
```

*Proof.* The maximum-cut argument gives `s<=n/2`. For a connected graph
with `n>=6`, the cited sharper bound implies

```text
q <= (81/16) 3^(-n/2) <= 3/16 < 1/4.
```

The only nonempty simple cubic graph on four vertices is `K4`; it has
`C=6,s=2`, giving `q=1/4`. Counts multiply over disjoint union, while `n,s`
add, so `q` multiplies. Every component has at least four vertices, proving
the disconnected assertion. QED.

This supplies an exact incompatibility: the squared-leg asymmetry of
`(3,4,5)` is `(16-9)/25=7/25>1/4`. Therefore no nonempty simple cubic graph
has a boundary-tensor triangle similar to `3:4:5`.

Even more strongly, the cited bound shows `q(G_n)->0` for **every** sequence
of connected simple cubic graphs with order tending to infinity, including
colorable graphs. Isosceles convergence alone cannot classify snarks. It is
an entropy dilution effect of this normalization; only the exact endpoint
is a colorability certificate.

**A size-normalized repair (PROVED, no new geometry).** For a nonempty
finite simple cubic graph retain its order and define

```text
r(G)=q(G)^(1/n)=C(G)^(2/n) 3^(s/n) / 6,
```

with `r=0` when `C=0`. If the graph is colorable, the free action of global
permutations of its three colors gives `C>=6`. Since `s>=0`,

```text
r(G) >= 36^(1/n)/6 > 1/6.
```

Thus zero and colorability are separated uniformly in this derived rate;
growing colorable graphs cannot dilute toward zero after the root is
taken. Over disjoint union, `q` multiplies and `r` is the corresponding
vertex-weighted geometric mean. The rate is not the original boundary
angle or overlap: it uses the retained order to undo exponential scaling.
It still requires the same coloring count and therefore supplies no
shortcut to deciding colorability or proving a Petersen-minor theorem.

### An independent planar connected hostile

Let `R_k=C_k square K2` be the prism graph, `k>=3`. Cut it into vertical
columns and record the ordered colors `(a,b)` of the two horizontal edges
entering a column. If `a=b`, the outgoing pair must be `(c,c)` with `c!=a`.
If `a!=b`, the outgoing pair must be `(b,a)`. The transfer matrix therefore
splits into the adjacency matrix of `K3` and three copies of `K2`. Taking
its `k`th trace gives the elementary exact formula

```text
C(R_k) = 2^k + 3 + 5(-1)^k.
```

For even `k`, `R_k` is connected, simple, cubic, planar, and bipartite, with
`s=0`, hence

```text
q(R_k) = (2^k+8)^2 / 6^(2k) -> 0.
```

This hostile requires neither the cited counting bound nor any snark
classification.

Another elementary hostile starts with `K4` and repeatedly replaces one
vertex by a triangle, attaching the three former incident edges to its
three vertices. In a coloring of the new triangle, each color occurs an
odd number of times on its three external edges (count color incidences
inside the odd vertex set), hence each occurs once. Conversely three
distinct external colors force a unique internal coloring. Thus every
replacement preserves `C=6`, simplicity, cubicity, planarity and absence
of bridges. Its order grows while `q<=36/(2*sqrt(3))^n` tends to zero.

## 4. Exact small controls and a separate arithmetic chart

The old THM-4116 `boundary_color_counts` and `full_coloring_count` engines
were used independently: a direct maximum-cut scan supplied the boundary,
then both the inner-product and full-edge counts were evaluated. Every
listed side vector had only unit nonzero entries and its norm product
agreed with Proposition 1.

| Graph | n | MaxCut | s | C | q | normalized squared legs |
|---|---:|---:|---:|---:|---:|---|
| K4 | 4 | 4 | 2 | 6 | 1/4 | 3/8, 5/8 |
| K3,3 | 6 | 9 | 0 | 12 | 1/324 | 323/648, 325/648 |
| triangular prism | 6 | 7 | 2 | 6 | 1/144 | 143/288, 145/288 |
| cube | 8 | 12 | 0 | 24 | 1/2916 | 2915/5832, 2917/5832 |
| Petersen | 10 | 12 | 3 | 0 | 0 | 1/2, 1/2 |

Reproduce the broader session census with
`python 04-computation/graph_pythagorean_20261018_graph_probe.py`.
See the [graph probe](../../04-computation/graph_pythagorean_20261018_graph_probe.py)
and its [saved output](../../04-computation/graph_pythagorean_20261018_graph_probe.out).
The elementary prism formula above is an independent path.

There is also a literal complete-graph realization of the difference-one
Pythagorean family, with a different map. For `K_(2r+1)`,

```text
MaxCut=r(r+1),       s=r^2,
(a,b,c)=(2r+1, 2r(r+1), 2r^2+2r+1)
       =(|V|, 2 MaxCut, 2 MaxCut+1).
```

This is Euclid's formula at `(m,n)=(r+1,r)`, so `c-b=1` and `c-a=2s`.
It puts `K3` at `(3,4,5)` and `K5` at `(5,12,13)`. Separately `K_(m,n)`
has `mn` edges, giving Euclid's doubled leg `2mn=2|E|`. Coprimality,
opposite parity, and unequal positive parameters are still required for a
primitive nondegenerate triple. Balanced `K3,3` gives the degenerate
parameter pair `(3,3)`, not a primitive triangle.

Thus maximum-cut data make two parts of the analogy literal, but the
boundary-state isosceles construction is a different chart. A claimed
single isomorphism must still reconcile these maps and their lost data.

## 5. Tutte scope and source freshness

The word “Tutte” needs a specified theorem. In the snark context, distinguish
the elementary cubic three-edge-coloring / nowhere-zero four-flow
equivalence, the cubic Petersen-minor statement, and the general four-flow
conjecture. Tutte's excluded-minor theorem for regular matroids instead
uses `U_(2,4),F7,F7*`; it is a different carrier from this graph triple.

The repo's primary-source intake,
[apex cubic audit](../reference/CORE-PAPERS-APEX-CUBIC-2026-08-25.md),
records [arXiv:2608.22870](https://arxiv.org/abs/2608.22870) as an external
computer-assisted preprint under audit. Its Theorem 1.4 claims that every
2-connected apex cubic graph is three-edge-colorable and, together with
the prior reductions, claims the cubic Petersen-minor conclusion. This
lane did not replay that paper's finite obligations and does not import
them as internally proved. The paper does not settle the general
noncubic Tutte four-flow conjecture.

## 6. Independent reading of the other two lanes

The [arithmetic audit](graph_pythagorean_20261018_arithmetic_audit.md) was
read directly, including the square-shore integrality proof, the two
different unit-gap branches, the explicit Petersen switching map, the
Seidel spectral bounds and the full finite-field construction. The
character-sum proof, norm/discriminant calculation, and both switched
degree formulas check algebraically. The graph `H_t` should be described
as a representative of a symmetric-conference switching class, rather
than silently identified with the usual odd-order conference graph.

The [graph audit](graph_pythagorean_20261018_graph_audit.md) was also read
directly. Its five-pole support proof is exact: `5=2m_c+k_c` forces boundary
multiplicities `(3,1,1)`; the unique internal majority-colored edge joins
the two singleton positions; the remaining four-edge path has one forced
alternating coloring. Thus each adjacent singleton pair contributes six
words and two glued pentagons have coloring count `6|E(C) intersect E(D)|`.
This makes the complementary decomposition of `E(K5)` a literal Petersen
coloring obstruction. The flower maximum-cut proof also checks: three
edge-disjoint odd horizontal paths enforce three uncut edges when some
module's leaves are equal, while otherwise every module contributes one
uncut star edge. Neither argument relies on a numerical graph invariant
standing in for the missing extension addresses.
