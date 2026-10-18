# Marked snark boundaries, counted six-colour defects, and the geometric gap

Status: **PROVED** elementary transfer and obstruction statements below;
**FINITE-EXACT** controls on the seven explicitly labelled graphs in §7.
No unit-distance six-equality fragment, seven-chromatic plane graph, or
geometry-preserving Petersen reduction is claimed.

## 1. Inheritance and the actual types

The starting mechanism is
[THM-4116, boundary-state gluing](../../01-canon/theorems/THM-4116-boundary-state-gluing-and-ap-odd-shell-tree-synchronizers.md):
ordered four-pole extension counts, with actual edge-colour witnesses, glue by
contraction of their full boundary vectors. Its Petersen/Blanuša/J5 atlas is
inherited here, not rediscovered as a new snark classification. The new use is
a counted lift into a *six-colour defect problem*, retaining pair labels,
binary shades, and the original graph.

The corrected flow group is \(\mathbb F_2^2\), not \(\mathbb Z/3\): see
MISTAKE-507. A cubic graph has a proper three-edge-colouring precisely when it
has a nowhere-zero \(\mathbb F_2^2\)-flow. These are edge colours. They are
vertex colours only after the explicit line-graph map.

The current geometric target is
[THM-4661, equality testers and the twenty-vertex exterior bound](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md).
The proved exterior bound is **20**, under its specified odd-degree carrier
and triangle-free base-field hypotheses; the finite SAT exclusion up to 21
has a different verification status. A field containing \(\sqrt{-3}\) does
not satisfy that base-field hypothesis. Neither bound applies merely because
a graph is cubic or has a Petersen minor. The target is a colourable unit
graph with distinct marked ports that are equal in *every* six-colouring.
The same theorem explains why an actual such fragment would settle the
seven-colour side of the plane problem.

[THM-4701, Petersen and finite unit circles](../../01-canon/theorems/THM-4701-kuratowski-graphs-the-petersen-graph-and-the-finite-unit-circle-graphs.md)
keeps unit embeddings, finite-field reductions, and minors separate. In
particular, contraction is not an operation on a fixed Euclidean unit drawing.
Its corrected pentagon statement concerns the side vectors, not an equality
between the pentagon vertex set and that side-vector set.

## 2. Pair defects: a useful six-colour invariant with its pairing retained

Let \(G\) be a finite simple cubic graph with \(m\) edges. Fix six colours
\(c\in\{0,\ldots,5\}\), paired by
\[
\pi(c)=1+\lfloor c/2\rfloor\in\mathbb F_2^2\setminus\{0\},
\qquad s(c)=c\bmod2.
\]
Here the nonzero vectors are represented by integers 1,2,3 and addition is
bitwise XOR. For a proper six-edge-colouring define
\[
d_v=\bigoplus_{e\ni v}\pi(c_e).
\tag{1}
\]
At a cubic star its three actual colours are distinct. Its pair labels are
either all distinct, giving \(d_v=0\), or have pattern \((a,a,b)\), giving
\(d_v=b\ne0\). Three copies of one pair are impossible because a pair has
only two shades. Every edge is counted twice, so
\[
\bigoplus_{v\in V(G)}d_v=0.
\tag{2}
\]
If \(G\) is a snark, meaning here a connected bridgeless cubic graph without
a proper three-edge-colouring, defects cannot all vanish. Equation (2)
excludes exactly one nonzero defect. Thus every proper six-edge-colouring
has at least two defective vertices. This statement depends on the fixed
pairing; it is not an obstruction to six-edge-colourability.

**The tetrahedral coordinate.** Identify the six colours with the six edges
of \(K_4\), pairing opposite edges. Among its twenty three-edge subsets,
four stars and four triangles have zero pair defect; the other twelve are
\(P_4\)'s and have nonzero defect. Across the fifteen perfect matchings of
the six palette colours, every fixed triple is defect-free six times and
defective nine times. Consequently summing the *number* of defective vertices
over all fifteen pairings gives \(9|V(G)|\) for every proper six-edge-colouring.
Averaging the pairings erases this diagnostic; a marked pairing is necessary.
Complete sums over all colourings are conjugate under palette relabelling;
the pairing sidecar matters for an actual witness and its boundary gluing.

## 3. Exact counted lift from an adjacent four-pole

Choose \(e=uv\), delete \(u,v\), and retain its four incident semiedges in
the fixed order \((uu_1,uu_2,vv_1,vv_2)\). Let \(m_e\) be the number of
three-edge-colour extensions of boundary word \((a,a,a,a)\), for fixed \(a\).
For a snark, THM-4116 gives the complete profile
\[
Z_e(a,a,b,b)=m_e\quad(a,b\in\{0,1,2\}),\qquad
Z_e(w)=0\text{ otherwise};\qquad m_e\text{ is even}.
\tag{3}
\]
Briefly, boundary parity leaves the four even-multiplicity types. Either
crossed two-colour type would restore a Tait colouring across \(uv\), so it
is absent. A two-colour Kempe-path swap bijects the all-equal and paired
types; colour permutation supplies the remaining equalities. Swapping the
two nonboundary colours is fixed-point-free, giving evenness. The code
independently checks the full 81-entry profile using matchings.

**Counted transfer theorem.** If \(m_e>0\), the minimum defect count is two,
and the exact number of proper six-edge-colourings whose defect set is
\(\{u,v\}\) is
\[
\boxed{3m_e2^m.}\tag{4}
Each of the three common nonzero defect values occurs \(m_e2^m\) times.

To prove it, take a full four-pole colouring with boundary \((a,a,b,b)\).
Give the restored edge \(uv\) a pair label \(c\notin\{a,b\}\). There are
two choices if \(a=b\), one otherwise, hence
\((3\cdot2+6\cdot1)m_e=12m_e\) pair-labelled objects. Choose one binary
shade per original edge. The only restrictions are opposite shades on the
two \(a\)-edges at \(u\) and opposite shades on the two \(b\)-edges at \(v\).
They concern four distinct edges and give two independent equations, so
there are exactly \(2^{m-2}\) choices. The resulting defects are precisely
\(d_u=d_v=c\). Conversely every colouring with this defect set restricts to
(3); the restored pair cannot equal either repeated pair, so the construction
is bijective. For a fixed \(c\), there are four allowed \((a,b)\)'s, proving
the refined count. Positive \(m_e\) is a real premise, not a property assumed
for every edge of every snark.

`lift_two_defects` is a source-bound implementation: it reconstructs the
ordered four-pole, validates the supplied full three-colouring and all shade
bits, restores the edge, and checks every actual vertex defect.

## 4. What transfers to vertex-six equality, and what does not

The line graph \(L(G)\) has maximum degree four. In any graph of maximum
degree \(\Delta\), a proper partial colouring extends greedily when
\(k\ge\Delta+1\): an uncoloured vertex sees at most \(\Delta\) forbidden
colours. Thus **every nonadjacent marked pair of \(L(G)\) admits all 36
six-colour assignments**; an adjacent pair admits precisely the 30 unequal
ones. This is true for every cubic source, not only the tested snarks.
No bare cubic line graph is a six-equality tester. Equivalently, a colourable
six-equality tester must give each marked port at least five neighbours;
otherwise one port can be recoloured away from the other.

There is an exact *abstract* palette lift. Let \(X_e\) be the line graph of
the four-pole, with adjacency only at the retained cubic vertices. Form
\(H_e=K_3\vee X_e\), retaining the three labelled palette vertices. Once
their distinct colours are chosen, \(X_e\) uses the other three colours.
Equation (3) therefore gives a six-equality tester at either paired port.
Its total six-colouring count is
\[
(6)_3\,9m_e=1080m_e.
\]
The **full** four-port weights, for fixed actual six-colours, are
\[
Z_{H_e}(a,a,a,a)=60m_e,
\qquad Z_{H_e}(a,a,b,b)=24m_e\ (a\ne b),
\tag{5}
\]
with all other words absent. The differing factors count the available
ordered palettes, \((5)_3\) and \((4)_3\). Thus even an exact support
transfer need not preserve uniform multiplicities.

For Petersen this is a 17-vertex, 69-edge abstract graph with 2,160 proper
six-colourings, and 360 for each fixed common colour at a marked port pair.
But the palette plus any vertex of \(X_e\) is a \(K_4\), which cannot be
unit-distance in the Euclidean plane. This positive algebraic transfer is
not a candidate drawing of the requested fragment.

## 5. Suppression must retain the shade holonomy

The relevant local reduction has an exact, limited survivor. On a degree-two
vertex impose both proper six-edge colouring and zero pair defect. Its two
incident colours must be paired but unequal, hence
\[
c_{\rm out}=\tau(c_{\rm in}),\qquad \tau(c)=c\mathbin{\mathrm{xor}}1.
\]
A path with \(r\) internal degree-two vertices has full boundary relation
\(c_{\rm out}=\tau^r(c_{\rm in})\), with exactly one interior extension
for every allowed boundary assignment. Relation matrices compose exactly:
\(R_rR_s=R_{r+s}\). Suppression is lossless when this twist is retained;
an ordinary untwisted wire is correct only for even \(r\). A closed
zero-defect wire has no extension when the total twist is odd.

The smallest shade-loss witness is one degree-two vertex: ports \((0,1)\)
extend and \((0,0)\) do not, although their pair labels agree. This statement
concerns a marked zero-defect constraint, not an automatic law of every
six-colouring. For ordinary *vertex* colouring, the two-edge path has five
extensions at equal fixed endpoint colours and four at different colours;
suppressing its middle vertex into an edge changes these numbers to zero
and one. Even boundary support is then lost.

**Full shade-lift theorem.** More generally, fix any pair assignment
\(\ell:E(G)\to\{0,1,2\}\) on a simple cubic source. If a star has three
edges of one pair, no lift exists. Otherwise build \(C_\ell\) on the
\(m\) *original edge labels*, joining two when they meet at a vertex and
have the same pair label. A star of pair multiplicity \(111\) gives no
constraint; a star of multiplicity \(210\) gives exactly one. Each original
edge meets at most one such constraint at each endpoint, so
\(\Delta(C_\ell)\le2\). The assignment lifts to a proper six-edge-colouring
if and only if \(C_\ell\) is bipartite. In that case the exact number of
lifts is \(2^{\kappa(C_\ell)}\), counting isolated vertices as components.
This follows directly: shades are precisely proper binary colourings of
\(C_\ell\), independently reversible on each component. The compiler
`shade_fibre` returns either a monochromatic star, a literal odd cycle of
original edge labels, or all components and a binary witness; `shade_lift`
then applies supplied component bits and rechecks the full source colouring.

In particular the exact defect partition function is
\[
\sum_{\text{proper six-colourings }c} z^{|\{v:d_v(c)\ne0\}|}
=\sum_{\substack{\ell:\text{no monochromatic star}\\C_\ell\text{ bipartite}}}
 2^{\kappa(C_\ell)}z^{|E(C_\ell)|}.
\tag{6}
\]
This is a finite identity, not a concentration or positivity assumption.
For the adjacent two-defect construction, \(C_\ell\) is exactly two
disjoint constraint edges plus isolates, recovering \(2^{m-2}\).

The triangular prism gives a cubic hostile to charge-only reconstruction.
Label all three top-triangle edges by pair 0, all three bottom-triangle
edges by pair 2, and the three spokes by pair 1. Every star has multiplicity
\(210\), every defect is the vector 2, and their total XOR is zero.
Nevertheless \(C_\ell\) consists of two odd triangles plus three isolates,
so no shade lift exists. The isolated spokes must still be retained in the
fibre even though they cannot repair either odd cycle. These graphs are
constraints on colours; no geometric picture is used to certify them.

**Operational next target.** A proposed geometric defect coupler must retain
its actual unit coordinates and all induced edges, expose the marked six
colours (pair and shade), and supply its complete boundary extension table.
Contracting that table with the counted packet (4), or a shade-twisted wire,
is a valid finite test: the desired result has a nonzero equal-port entry and
no separated-port entry. Retaining only (2), a Petersen minor, or the fact
that both tables are nonempty is insufficient. This specifies an executable
search interface; this packet supplies no geometric coupler that succeeds.

## 6. Petersen minors: one-way information, with an actual cubic hostile

The theorem that a 2-connected cubic graph without a Petersen minor is
three-edge-colourable is **CITED** from
[Robin Thomas's author statement](https://thomas.math.gatech.edu/FC/generalize.html).
The published [Edwards--Sanders--Seymour--Thomas doublecross paper](https://doi.org/10.1016/j.jctb.2015.12.006)
describes its role within that proof programme. We do not treat that one
paper as a standalone proof of the entire theorem. Its contrapositive
supplies a necessary minor for bridgeless cubic snarks, not a converse or a
vertex-six boundary equivalence.

The JSON contains independently checked, explicitly marked Petersen
subdivisions in Petersen, both Blanuša graphs, J5, J7, **and Desargues**.
Desargues is a connected bridgeless cubic graph with 192 proper
three-edge-colourings, including the displayed complete witness. Thus even
within that graph class a Petersen subdivision does not imply a snark.
The subdivision reader checks every source edge, ten distinct branch
vertices, and disjoint path interiors; it never contracts a drawing.

## 7. Finite universe, controls, and reproduction

All vertex and edge labels are fixed by `examples()`. One marked edge (the
first canonical edge) is tested per source. This is not a complete snark
census. Whole-graph Tait enumeration is exhaustive; each four-pole profile
is checked at all 81 words by a second, perfect-matching engine.

| Source | vertices / edges | \(m_e\) | six-colourings with defects exactly at the marked endpoints |
|---|---:|---:|---:|
| Petersen | 10 / 15 | 2 | 196,608 |
| Blanuša I | 18 / 27 | 4 | 1,610,612,736 |
| Blanuša II | 18 / 27 | 4 | 1,610,612,736 |
| Flower J5 | 20 / 30 | 12 | 38,654,705,664 |
| Flower J7 | 28 / 42 | 44 | 580,542,139,465,728 |

K4 and Desargues are positive Tait controls (6 and 192 colourings); all five
snark controls have zero. The earlier THM-4116 atlas supplies the first four
snark graph definitions and profiles; J7 tests the same transfer on a larger
source. All compatible marked pairs in all seven line graphs are extended
and replayed. One Petersen shade fibre is independently counted over its
entire \(2^{15}\) binary cube. Wire relations are composed for lengths
0 through 12 and literally counted through length 3. Twenty tetrahedral
triples and all fifteen palette pairings are checked. Malformed exact types,
port order, missing subdivision edges, improper shades and palettes reject.
The full K4 defect polynomial is independently computed by (6), over all
\(3^6\) pair assignments, and by exhaustive proper six-colouring; the
triangular-prism odd-cycle receipt is replayed separately.
The two engines give exactly \(384+1728z^2+1536z^3+432z^4\).

The implementation needs only Python's standard library. Production lifts
authenticate supplied witnesses; finite exhaustive discovery remains in
`main()` and the explicitly named enumerators. Reproduce with:

```text
python -B 04-computation/experiments/snark_boundary_transfer_20261018a.py
python -B -O 04-computation/experiments/snark_boundary_transfer_20261018a.py
```

The companion JSON is the deterministic normal output. Normal and optimized
runs agree at **91,470 checks**. LF SHA256:
`57e702ef97992a9d43dd30c908dffba0ab90116fc5d6682afc6d2978f83c40c2`.
Independent full proof/API review and optimized-versus-saved replay: **PASS**,
including the shade-fibre theorem, source reconstruction, exact multiplicities,
odd-cycle receipts and geometric scope. No correction was requested.
