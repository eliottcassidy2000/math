# Variable pair labels: exact boundary cells and the remaining equality test

Status: **PROVED**, the finite-interface statements below; **FINITE-EXACT**,
the declared witness-bank compression; **OPEN**, a unit-distance fragment
that keeps equal-endpoint six-colourings and eliminates every separating one.
This packet supplies neither such a fragment nor a value of the chromatic
number of the plane.

## 1. Inheritance and the current object

The closest proved mechanism is the pinned bipartite-fibre criterion in
[reanchoring kernel, section 7](reanchoring_kernel_20261018c.md).
The canonical hostile is its four-vertex graph whose chosen cycle basis
misses a monochromatic triangle. The corrected near miss is an isolated
vertex with a pin of the wrong ternary label: parity tests alone accept it.
The least-used sidecar is the complete affine response at the attachment
boundary, including component phases and colour names.

The anchor is six-colour boundary support. The niche is lossless finite
interface compression; the wildcard is odd-degree descent, in
[odd primes, functions and cycles](odd_prime_function_cycle_20261018d.md).
The live concept board is: variable pair labels; actual odd circuits; pinned
path parity; marked boundary; global colour permutation; field-coset contact.
This keeps a concrete map and its lost information attached to each analogy.

The source graph and 570 full witnesses are frozen in
[the preceding joint branch packet](coupled_support_branches_20261018c.md).
The new [variable-pair compiler](variable_pair_fibres_20261018d.md) searches
the label coordinate and retains learned actual-cycle/path guards.
The [cubic geometry packet](surviving_repair_geometry_20261018d.md) determines
which attachment edges these certificates must actually satisfy.

## 2. Fixed-label response is an affine binary code

Let \(G=(V,E)\) be a finite simple graph, \(B=(b_1,\ldots,b_k)\) a list of
distinct marked vertices, and \(\ell:V\to\{0,1,2\}\) a proposed pair label.
For the arithmetic palette convention, write

\[
 c(v)=3s(v)+4\ell(v)\pmod6,\qquad s(v)\in\mathbb F_2.
\]

This uses the particular pairs \(\{0,3\},\{1,4\},\{2,5\}\). Any other pairing
requires a retained palette relabelling.
Let \(H_\ell\) retain precisely those actual edges \(uv\in E\) with
\(\ell(u)=\ell(v)\). Edges of different labels are automatically proper.
Edges of \(H_\ell\) require \(s(u)+s(v)=1\).

**Boundary-cell theorem.** A fixed label word has a six-colour lift iff
\(H_\ell\) is bipartite. Choose a reference parity \(\pi_v\) on each component.
Every lift, uniquely, is

\[
 s(v)=\pi_v+\epsilon_{C(v)},\qquad
 \epsilon_C\in\mathbb F_2.
\tag{1}
\]

For named boundary colours \(a_i\), the exact additional conditions are

\[
 a_i\bmod3=\ell(b_i),\qquad
 (a_i\bmod2)+(a_j\bmod2)=\pi_{b_i}+\pi_{b_j}
 \quad\hbox{whenever } C(b_i)=C(b_j).
\tag{2}
\]

If \(h\) components avoid \(B\), each compatible named boundary has exactly
\(2^h\) full lifts. If \(v\) components meet \(B\), the named boundary response
has exactly \(2^v\) elements.

**Proof.** Propagate a binary sign across each same-label edge. An odd cycle
is exactly the obstruction to consistent propagation. In a bipartite component,
the only freedom is its initial sign, giving (1). A boundary pin fixes that
component's initial sign; two such fixes agree exactly when their path parity
agrees. Label agreement is separately necessary. Components untouched by the
boundary retain independent signs. This proves both directions and both counts.

The lossless interface is therefore: boundary labels, a partition of boundary
indices by component, parity offsets within each block, and \(h\) if counts are
needed. The full interior label word is also retained when an actual interior
colouring must be reconstructed.

## 3. What can be glued after compression

If a future graph meets the old graph **only through \(B\)**, the above response
is sufficient for testing its compatibility with this fixed label word. Match
the literal named boundary colours; the old interior can then be restored by
(1). If a new edge touches an old vertex outside \(B\), that vertex must first
be added to the boundary. An eight-port summary is not a certificate for an
attachment to arbitrary interior vertices.

This distinction matters in the cubic experiment, whose contacts include old
interior vertices. Its compiler pins all 714 old colours for its stronger
extension test. The eight-port atlas below must not replace that test.

For colour orbits, a permutation of all six colour names transports a full
proper colouring to another one. It can change the ternary label word.
We retain this permutation explicitly. No claim is made that a transported
word lies in the original fixed arithmetic-label fibre.

Different cells can overlap. Summing their \(2^h\) multiplicities would
double-count colourings. The atlas is an existence certificate, not a formula
for the total number of colourings of the 714-vertex graph.

## 4. Exact criterion for the unpinned equality target

Let \(s,t\) be distinct vertices of an **unpinned** finite graph \(G\).
There is a proper six-colouring and every such colouring has \(c(s)=c(t)\)
if and only if:

1. some label word has bipartite \(H_\ell\); and
2. for every such word, \(s,t\) lie in the same component of \(H_\ell\) at
   even path parity.

Indeed, different labels force different colours. Equal labels in distinct
components allow their two phases to be chosen independently. Equal labels
in the same component force equality precisely at even parity. Every
six-colouring determines such a label word, so this covers the whole palette
space, not a selected family of charts.

Adding \(st\) to an even same-label path closes an odd cycle. This gives an
exact connection between the terminal equality condition and the odd-cycle
cuts used by the compiler. It is not enough to find one even path or one
successful label word: the condition quantifies over every feasible word.
With additional pinned vertices, separately grounded components can have
forced relative phases, so use (2), not this unpinned connected-component
criterion.

Positive logical control: \(K_7-st\) forces equality at the missing edge.
Hostiles: \(K_7\) has no six-colouring, and a two-edge path has separating
six-colourings despite its even length. These are abstract graph controls,
not unit-distance realizations. The current geometric source terminals
\((2,4)\) are not unit-separated; an eventual equality certificate there
would still need the appropriate geometric connector.

## 5. The declared 570-witness bank

The inherited graph has 714 vertices, 2,914 actual edges, and eight marked
source vertices. Its source is the one-point union of a 5-cycle and a 4-cycle.
There are exactly 570 canonical proper six-colour source patterns, reproduced
independently by restricted-growth enumeration. Of these, 123 have equal
colours at source terminals \((2,4)\), and 447 separate them.

For each saved full colouring and each of the 15 perfect matchings of its
six colour names, assign one ternary label per matched pair. This supplies
a bank of 8,550 valid fixed-label cells. Each cell's component signs are
then varied, and its boundary words are canonicalized under global palette
permutations.

| Exact finite statistic | Value |
|---|---:|
| Candidate cells, before deduplication | 8,550 |
| Distinct named component/parity boundary diagrams | 3,427 |
| Distinct canonical source support sets | 3,176 |
| Selected cells covering all 570 source orbits | 49 |
| Certified packing lower bound for this bank | 24 |

The canonical support sizes 1, 2, 4, 8, 16, 32 occur respectively
329, 1,505, 2,398, 2,317, 1,399, 602 times among the candidate cells.
Greedy set cover followed by redundant-cell deletion supplies the 49-cell
upper bound. A saved set of 24 source patterns has at most one member in
every candidate cell, proving the lower bound **within this bank**.
Neither 49 nor 24 is asserted to be optimal over arbitrary ternary words.

Each selected cell retains its full 714-trit label word, its originating
colouring/pairing, boundary diagram, hidden-component count, and exact support.
Each of the 570 receipts retains a cell index, visible component phases,
a six-colour permutation, and a digest of the reconstructed full colouring.
The consumer recomputes components from the actual graph; it does not trust
the stored diagram. Thus the compression preserves a verified witness for
every boundary orbit without pretending to represent every possible interior
colouring.

## 6. Reproduction and audit boundary

Files:

- [Generator and receipt consumer](../../04-computation/experiments/boundary_fibre_atlas_20261018d.py).
- [Exact atlas and 570 reconstruction receipts](boundary_fibre_atlas_20261018d.json).

Run from the repository root:

~~~text
python -B -X utf8 04-computation/experiments/boundary_fibre_atlas_20261018d.py --out 05-knowledge/results/boundary_fibre_atlas_20261018d.json
python -B -O -X utf8 04-computation/experiments/boundary_fibre_atlas_20261018d.py --out 05-knowledge/results/boundary_fibre_atlas_20261018d.json --verify-only
~~~

The consumer binds the inherited input by its LF-normalized hash and its
coordinate/graph hashes. All checks use explicit exceptions and remain active
under optimization. Controls enumerate all 75 simple graphs on one through
four labelled vertices, every ternary labelling, and every binary shade word,
independently comparing existence, named boundary support, and \(2^h\)
multiplicity. Pinned-path, isolated wrong-label, empty-boundary and terminal
controls retain the failure boundaries above.

Final normal generation and the optimized saved-receipt consumer both pass.
An independent union-find parity consumer recomputed all 8,550 candidate
supports, the 3,427 diagrams and 3,176 support sets, the bank-specific
24-pattern packing, the 49-cell cover, and all 570 full colouring receipts
with every source pin and digest. Its audit also rejected malformed literal
label words and Boolean index aliases. No independent-audit correction to
the mathematical result remained.

## 7. The next construction gate

A retraction to the old graph preserves every old colouring by pullback,
regardless of the number of new vertices or retained core edges. The new
cubic packet authenticates just such a repair: each of its six fractional
layers projects back by forgetting its shift.

The next useful configuration must therefore couple exterior layers so that
this map fails on an actual edge or identification. Such a failure only
defeats that particular repair. The decisive consumer remains the complete
variable-label criterion: keep at least one equal-terminal lift and exclude
every separating one with actual odd-cycle/pinned-path or complete search
certificates. More edges, a failed fixed chart, a nontrivial quotient cycle,
and a field extension by themselves do not imply this conclusion.

## 8. Incoming Hamming relation: a sharper joint target

Incoming commit 785a59499e supplies
[the ten-label Hamming packet, section 7](ten_label_matching_hamming_20261018n.md).
Its four-chart result has an exact consumer in the present boundary-cell
language. On ten ordered binary coordinates let

\[
 C=\operatorname{span}_{\mathbb F_2}(170,204,240,255,768),\qquad
 D=\operatorname{span}_{\mathbb F_2}(15,240,768),
\]

where integers are bit masks, with coordinate zero the least significant bit.
The 32-word code \(C=H_8\oplus I_2\) is the disjoint union of
\(r+D\) for \(r=0,170,204,102\).

Each chart can be constructed as an abstract fixed-label fibre of (1).
Use the boundary blocks \(\{0,1,2,3\}\), \(\{4,5,6,7\}\), \(\{8,9\}\).
Join each block's first vertex to each other vertex by one edge if their
required parity difference in \(r\) is one; otherwise use a two-edge path
with a fresh hidden vertex. All ternary labels are zero. The resulting named
boundary colours are \(c_i=3b_i\), with shade masks exactly \(r+D\).
The first graph has 17 vertices and 14 edges; the other three each have
13 vertices and 10 edges. Each boundary mask has one hidden extension.

The lower bound of four charts is structural. Any free visible-component
indicator of an affine component chart contained in \(C\) belongs to \(C\).
There is only one weight-two nonzero word, 768; all other nonzero words have
weight at least four. Disjoint nonempty supports therefore permit at most
three free components on ten coordinates. Each chart has at most eight
words, so at least four are needed to cover \(C\). The four constructions
achieve the bound.

The script reproduces these four **named** responses with the same fibre
engine as the 714-vertex atlas. An independent exhaustive shade enumeration
confirms the graphs, all 32 boundary masks, and unique hidden extensions.
Canonicalizing under arbitrary palette permutations would lose this marked
binary relation and is not performed in this control.

This supplies a joint-support test object, not a geometric selector. The
four graphs are separate fixed-label contexts; their union of response sets
does not specify a single graph with unrestricted ternary labels, much less
a unit-distance realization. The next transfer must realize the selection
between charts while preserving the ordered ports and excluding unwanted
boundary words. Covering every individual source word, as the 49-cell atlas
does, is a different predicate.

There is a small decisive hostile: simply gluing all four chart graphs on
their ten ports imposes intersection of their responses, not union. The
zero chart contains an even path from port 0 to port 1; the 170 chart adds
the edge between them. Together they give a same-label triangle and empty
fixed-label support. The script authenticates this triangle. A selector
must implement the disjunction without imposing all contexts simultaneously.

The companion [cube clock packet](cube_unitary_clock_packets_20261018n.md)
also explains a useful difference between two cycle tests. Integer clock curl
is linear and must vanish in \(\mathbb Z\), so a suitable full cycle basis
can certify it. A clock change by two can pass its parity test and still fail
the integer test. Monochromatic-cycle membership in the present label search
is not linear in the original cycle space: the four-vertex hostile requires
recomputing the same-label graph before reading its circuits. Both methods
retain the value group and the actual source, but their guards are not
interchangeable. No Collatz energy bound follows from this comparison.
