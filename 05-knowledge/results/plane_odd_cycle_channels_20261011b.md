# Six colours as three labels with binary cycle obstructions

2026-10-11. **PROVED:** the elementary equivalences and finite search termination
below. **FINITE-EXACT:** the specified 1,099-graph census and geometric controls.
**OPEN:** whether the real plane admits six colours. The accepted lower-bound
premise H6 is not proved in this packet. Nothing here proves universal Collatz
coverage or a new colouring complexity bound.

## Inheritance and live board

The closest proved mechanism is the named-boundary relation in
[boundary palette compression](boundary_palette_compression_20261011.md).
The canonical hostile is an odd cycle: checking only triangles loses C5.
The corrected near miss is equating a count of possible states with their
complete, authenticated relation. The least-used sidecar is a witnessed
cycle-parity obstruction, which can be learned only when a candidate needs it.
The board is **pair bits / ternary labels / binary lifts / cycle witnesses /
boundary seams / geometric scale**.

The starting decimal claim needs a correction:

\[
 10^3+10^2-10=1090,\qquad
 1099=10^3+10^2-1=\sum_{n=1}^5 2^{\binom n2}.
\]

For a fixed labelled vertex set [n], a simple graph chooses one bit for each
unordered pair. Thus the census is 1,2,8,64,1024. A three-colour threshold is
a test on these graphs, not an extra condition in that count. This is the
first precise triangular-number connection: the dimension of the pair-bit
cube is the triangular number binom(n,2).

## The exact alternative statement

For a finite simple graph G=(V,E), and a map t:V->[r], retain just edges whose
endpoint labels agree:

\[
 E_t=\{uv\in E:t(u)=t(v)\},\qquad H_t=(V,E_t).
\]

**Proposition.** G has a proper 2r-colouring if and only if some t makes H_t
bipartite. Indeed, write each colour uniquely as 2t+b, with 0<=t<r and
b in {0,1}. Edges between different labels already have different colours;
within each label the bit must alternate. Conversely, choose a bipartition
bit on each component of H_t and obtain the full colouring by the same formula.
Isolated vertices count as components. The empty graph is permitted.

Equivalently, form the hypergraph on V whose hyperedges are vertex sets of
odd cycles in G. G is 2r-colourable if and only if this hypergraph admits
r labels with no monochromatic hyperedge. It suffices to retain chordless
odd cycles: a shortest monochromatic odd cycle has no chord. A chord would
split it into two cycles, exactly one odd, contradicting minimality.

Consequently, under H6 and the usual seven-colour upper bound:

* The plane has chromatic number six if and only if every finite algebraic
  unit-distance graph has a weakly three-colourable odd-cycle hypergraph.
* It has chromatic number seven if and only if one finite algebraic
  unit-distance graph has an odd-cycle hypergraph requiring at least four
  labels in this weak sense.

The reduction to finite algebraic coordinates is the already established
[finite witness mechanism](finite_algebraic_witness_20261007.md) and
[critical witness packet](plane_algebraic_witness_20261011.md). The word
"algebraic" does not supply a coordinate-degree or search-size bound.

## Counting, Fourier characters, and what parity remembers

Let c(H) count components, including isolated vertices. The same proof gives

\[
 P_G(2r)=\sum_{t:V\to[r]}2^{c(H_t)}[H_t\text{ is bipartite}].             \tag{1}
\]

For a graph H on n vertices and m edges, expand
\(\prod_{uv\in E(H)}(1-(-1)^{b_u+b_v})/2\) and sum over all bits b.
Character orthogonality cancels precisely those edge subsets with a vertex
of odd degree. Therefore

\[
 P_H(2)=2^{n-m}\!\sum_{F\subseteq E(H):\partial F=0}(-1)^{|F|}.          \tag{2}
\]

The inner sum is either 2^(m-n+c(H)) or zero. The sign character is trivial
on the binary cycle space exactly when every cycle is even. If it is
nontrivial, translation by one odd cycle pairs opposite signs. Combining
(1) and (2) is an exact nonnegative outer sum of individually nonnegative
binary-colouring counts, despite the signed inner expansions.

This is a finite Fourier/Potts identity, not a new general theory. It is
also not the automorphism-reversal notion of "even graph" occurring in
tournament equinumerosity. A triangle has even vertex degrees and odd edge
count. See the separate
[typed tournament and flow comparison](tournament_even_transfer_20261011b.md).

## A complete finite learner, with a checkable stopping condition

Start with no odd-cycle constraints. Seek an r-label assignment avoiding
the current monochromatic cycle constraints. Use breadth-first binary
propagation in H_t:

1. If it succeeds, retain all labels and bits: this is an actual full colouring.
2. If it fails, recover a simple odd cycle from the conflicting paths,
   shorten it through chords, and retain its actual vertices and edges.
3. If no label assignment survives the stored constraints, certify that
   the constraints cover every assignment, modulo permutations of the r labels.

Every added constraint excludes the current assignment, which satisfied all
previous constraints. Thus the process terminates on each finite graph.
The prototype uses exhaustive canonical label assignments to check negative
coverage. It is a complete finite proof-of-concept, not a scalable SAT proof
format or a speed improvement claim. The certificate is not complete just
because every cycle it lists is genuine: its coverage must also be checked.

An original-graph cycle basis is not a complete shortcut for these label
constraints. In K4 use the three star-basis triangles 012,013,023 and labels
(1,0,0,0). Every basis triangle is nonmonochromatic, but the remaining triangle
123 is monochromatic. A full binary parity basis is complete **after** forming
the residual graph H_t; the nonlinear operation of deleting differently
labelled edges does not preserve the old basis. This hostile separates the
linear cycle-space codec from the label-selection problem.

The [signed boundary compiler](plane_six_reframes_parity_20261011b.md) is the
compositional version for a fixed retained label assignment. It remembers
connectivity and XOR path parity at named ports and lifts compatible bits.
Gluing needs these path data; local bipartiteness alone is insufficient.

## Exact experiment and hostile controls

Run from the repository root:

```text
python -B -X utf8 04-computation/experiments/plane_odd_cycle_channels_20261011b.py
python -O -B -X utf8 04-computation/experiments/plane_odd_cycle_channels_20261011b.py
```

The universe is every labelled graph of orders 1 through 5, with r=1,2,3:
3,297 channel tables. An independent full-colour equality-partition count
checks (1); full edge-subset sums, literal bit assignments, and BFS counts
independently check (2). All learned negative certificates are checked over
their entire finite label universe. Total: **8,382 exact checks**.

Additional controls are an exact algebraic induced unit C5, the exact Moser
spindle inherited from the palette packet, abstract K6 and abstract K7.
C5 defeats the triangle-only test for a fixed label assignment. K7 learns
all 35 triples and fails if one is deleted from its negative certificate.
K7 is only an abstract graph control; no plane unit realization is asserted.
Malformed types, a false cycle, and incomplete negative coverage are rejected.
Checks use explicit exceptions and remain active under Python optimization.

## Connections with their losses and next target

The source is a 2r-colouring, target an r-label assignment plus binary lift,
map c->(floor(c/2),c mod2), preserved predicate proper colouring. Forgetting
the bit destroys within-label adjacency information; retained parity paths
and full witnesses restore it. The cheapest hostile is C5.

This makes the next geometric target concrete: build an algebraically
realized unit graph whose odd-cycle constraints cannot be weakly labelled
with three labels, or prove a compatible choice exists for every such graph.
Neither arbitrary abstract hypergraphs nor disjoint local colouring tables
solve that realization/gluing problem. The
[high-girth and crowding packet](plane_girth_crowding_20261011b.md) sharpens
the target: constraints must operate at unbounded cycle lengths and geometric
resolutions. No finite catalogue of short-cycle types is assumed complete.
