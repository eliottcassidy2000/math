# Source-fixed difference tables, palette gauges, and shade repair

Status: **PROVED** finite-interface theorems below; **FINITE-EXACT** the declared
controls. A failed gauge search is not a six-colour obstruction. No new
six-equality fragment or seven-chromatic unit-distance graph is claimed.

## 1. Inherited mechanism and retained source

[The marked snark packet](snark_boundary_transfer_20261018a.md) separates pair
labels from binary shades and proves the exact bipartite shade-fibre count.
[The signed-boundary compiler](plane_six_reframes_parity_20261011b.md) retains
actual paths, component connectivity, internal odd cycles, and pinned source
bits. These mechanisms are reused here on **vertices of the actual graph**;
there is no implicit edge-colouring-to-unit-geometry transfer.

The target remains the marked six-equality problem of
[THM-4661, equality testers](../../01-canon/theorems/THM-4661-odd-degree-carriers-a-failed-six-colour-exterior-has-at-least-twenty-vertices-and-equality-testers-are-seven-colour-certificates.md).
Its field and exterior hypotheses are not inferred from an abstract graph.
The current counterpressure is more specific: at least thirteen additional
projected pairs and absence of a source-fixing retraction do not by themselves
remove any separating source colouring.

Fix a finite labelled graph \(B\) on \(n\geq1\) vertices, a proper colouring
\(\phi:V(B)\to\{0,\ldots,5\}\), and \(L\) named copies. A `Layers` record
retains all \(Ln\) vertices, the exact induced graph of every copy, and every
cross edge. The first copy and its colouring remain literal, not merely
isomorphic. The abstract constructor makes no claim of metric realizability.
The rational-coordinate constructor checks all pairs exactly, rejects
coincident vertices, and requires zero translation of the original copy.

Live interfaces are: geometric difference; palette permutation; vertex pair
label; pinned shade; complete boundary support. Passing through any one of
these interfaces requires its stated sidecar.

## 2. Difference tables and their exact palette quotient

For marked source coordinates \(C_a\), retain

\[
D_{xy}=\{(C_b-C_a,a,b):\phi(a)=x,\ \phi(b)=y\}.
\]

The ordered pair \((a,b)\) is part of each entry. For layers translated by
\(t_i,t_j\), an entry is an actual unit edge exactly when
\(\|C_b-C_a+t_j-t_i\|^2=1\). Thus the complete colour-difference matrix is

\[
M_{ij}[x,y]=1\quad\Longleftrightarrow\quad
\text{at least one such authenticated edge has source colours }(x,y).
\tag{1}
\]

For the restricted colouring ansatz
\(c(i,a)=\sigma_i(\phi(a))\), with \(\sigma_i\in S_6\) and
\(\sigma_0=\mathrm{id}\), (1) is a lossless quotient:

\[
M_{ij}[x,y]=1\ \Longrightarrow\ \sigma_i(x)\ne\sigma_j(y).
\tag{2}
\]

Internal edges remain proper because permutations are injective. Conversely,
(2) checks every cross edge. `find_gauges` is complete for this ansatz only.
Unused colours are retained but permutations agreeing on the used source
colours are identified when counting actual colourings.

The subgroup preserving the three pairs \((01),(23),(45)\) has order 48.
It is strictly weaker than \(S_6\). Six forbidden assignment cells already
separate them: force \(\sigma(0)=0\) by forbidding images 1 through 5, and
forbid \(\sigma(1)=1\). Every paired gauge fails; 96 full permutations pass.
This six-cell bound is sharp for covering the 48 gauges: each single
assignment prohibition hits exactly eight of them, so five cannot cover all.
The witness is an abstract graph with six isolated source anchors, not an
unverified unit-distance configuration.

Even full \(S_6\) gauges can fail while actual six-colourings exist. Take seven
copies of six isolated vertices, source colour all zero, and connect each
pair of copies by one edge using a different vertex at each incidence.
The actual graph is a 42-vertex matching. Monochromatic-copy gauges demand a
seven-colouring of the copy quotient \(K_7\), while source-pinned binary shade
repair has exactly \(2^{15}\) extensions. A quotient obstruction is not an
obstruction in the original graph.

For two copies there is a simple complete alternative. Enumerate every
canonical proper partition \(\kappa\) of the child copy, and every injection
of its used classes into six colour names. Retain exactly the injections
satisfying the cross edges. This is a disjoint parametrization of **all**
source-fixed six-colourings: each colouring determines its unique canonical
partition and class injection. `complete_two_layer_count` implements this
complete finite sum; it does not restrict the child to the source partition.

## 3. Exact pinned shade fibres and all fifteen pairings

First fix the pairing \((01),(23),(45)\). At vertex \((i,a)\), retain the
pair label \(\ell(i,a)=\lfloor\phi(a)/2\rfloor\). Let \(H\) contain exactly
the actual edges whose endpoints have equal pair labels. A shade is a bit
\(b(i,a)\), with the original source bits prescribed. A proper colouring
with these retained labels exists **if and only if**:

1. \(H\) is bipartite; and
2. on every path joining two source vertices, path-length parity equals the
   XOR of their prescribed source bits.

Necessity is immediate along edges. For sufficiency choose the binary phase
of each connected component. A source-containing component has exactly one
allowed phase if its pins agree; a component without source has two. Hence
the exact number of colourings in this fibre is

\[
2^{\#\{\text{components of }H\text{ containing no source vertex}\}}.
\tag{3}
\]

Isolated vertices count as components. The compiler returns either an actual
colouring, an actual odd cycle, or an actual source-to-source path with
incompatible parity. All these refer to the original graph; a pinned path
need not have odd length. `audit_repair` reconstructs the whole instance.

Now run the same construction for each of the 15 perfect matchings of the
six colour names. Each pair is given a canonical orientation only to encode
bits. The witness is decoded back into the original names, and the entire
source colouring is checked unchanged. One successful chart certifies a
repair; fifteen failures reject only this fifteen-chart cover. Counts from
different charts overlap and cannot be added. A one-anchor source and one
isolated child illustrate the loss: the 15 fibres have total size 30, but
their union has only six child colours.

There are two explicit negative-scope controls. In an equilateral unit
triangle regarded as three copies of a single source anchor, every chart
retains one common pair label and finds an odd cycle, whereas a proper
three-colouring exists. Separately, an abstract three-copy six-anchor graph
forces one copy of source colour 0 to become 1 and another to become 2.
No single perfect matching can pair 0 with both 1 and 2. All fifteen charts
return pinned-path obstructions, while independent full \(S_6\) gauges pass.

The exact rational translate with source points \(0,2\), source colours
\((0,1)\), and offsets \(0,1\) is a unit path on positions \(0,1,2,3\).
The source pins at positions 0 and 2 disagree along an even path. Exactly
three matching charts fail and twelve succeed. Thus even within exact
geometry, a failed chosen pairing cannot be treated as a terminal result.

## 4. Copy-level holonomy and a support-preserving forest theorem

A still smaller repair uses one bit \(t_i\) per copy and
\(c(i,a)=\phi(a)\mathbin{\mathrm{xor}}t_i\), with \(t_0=0\).
For an actual cross edge with equal pair labels, the exact equation is

\[
t_i\oplus t_j=1\oplus(\phi(a)\bmod2)\oplus(\phi(b)\bmod2).
\tag{4}
\]

Cross edges between different pairs impose no condition. A cycle of these
equations is consistent precisely when its signed XOR is zero. The compiler
retains actual cross-edge indices for every equation. Its inconsistent
cycle is a **copy-gauge** obstruction, not necessarily a cycle of \(H\).
Three disjoint edges on three two-vertex copies give an explicit odd copy
cycle but a perfectly valid vertex-shade repair.

**Forest theorem.** Suppose copies are joined according to a forest, with
exactly one actual cross edge per adjacent copy pair and no other cross
edges. Every source colouring extends, for every chosen palette pairing.
Indeed the constraint subgraph of (4) is a forest, so propagate its bits.
For a connected tree of \(L\) copies, if \(s\) bridge edges have matching pair
labels, the exact number of uniform-copy flip repairs is \(2^{L-1-s}\).
The first copy fixes one binary phase; every other constraint component is
free. Vertex-level repair may have additional freedom. No source retraction
is required by this proof.

Apply this to the source \(B\), the one-point union of \(C_5\) and \(C_4\) at vertex 0,
with edge set \(01,12,23,34,04,05,56,67,07\). Add thirteen copies and one
bridge to each. Projected pairs complete the graph on anchors
\(\{0,2,3,4,5,6,7\}\) except \(24\), adding exactly thirteen missing pairs.
Orient the pair \(26\) as source 6 to copied 2. The physical graph has
112 vertices and 139 edges.

It has no source-fixing homomorphism to \(B\): every homomorphism from the
copied \(C_5\) into \(B\) maps it onto the original \(C_5\), while vertex 6
has no neighbour on that cycle. The finite check enumerates its ten closed
five-walks. Nevertheless the forest theorem repairs every one of the 570
canonical source patterns, including all 447 separating patterns at ports
\((2,4)\). There are 123 equal patterns. Every one of the 15 matching charts
succeeds: 8,550 witness obligations in total, 6,705 of them separating.

[The exact metric control](shared_phase_threshold_hostile_20261018b.json)
supplies the actual coordinates and complete edge census for this same
graph. Its translations are engineered mixed radii; the result is not a
claim about the restricted shared \(b/c\) phase search bank. Our abstract
support proof consumes the authenticated edge bank, not a projected graph.
The literal pentagon ports \((2,4)\) are nonunit-separated. Even hypothetical
equality forcing at that pair would still require a geometrically legal
terminal connector before giving the deleted-unit-edge target.

## 5. Operational consequence and verification

A candidate six-equality fragment must retain an equal source witness and
eliminate every separating source pattern. For a separating source pattern,
any successful gauge or shade chart supplies an explicit separating witness
and rejects that candidate. For an equal pattern it supplies the desired
positive witness instead.
For a pattern failing all fifteen charts, the returned odd cycles and pinned
paths are useful retained obligations, but a full boundary-colouring solver
is still needed. There is no completeness claim for bounded chart families,
bounded-width quotients, or snark minors. In particular, coupled cross edges
must defeat the forest repair; merely breaking a preferred projection or
every source retraction is insufficient.

The exact test universe is every simple source graph on one through three
vertices, every cross-edge subset between two named copies, and every
canonical proper source colouring using at most six names. It compares the
complete partition/injection sum to direct literal six-colour enumeration,
and compares (3) to the independently counted fixed-pair subset. Separate
controls check rational geometric difference provenance, unused colours,
the sharp six-cell subgroup obstruction, all fifteen chart outcomes,
copy-holonomy false negatives, source reauthentication, and bool/float
rejections. The thirteen-bridge control verifies every source pattern.

Reproduce from the repository root:

```text
python -B 04-computation/experiments/shade_difference_support_20261018b.py
python -B -O 04-computation/experiments/shade_difference_support_20261018b.py
```

The companion JSON is the exact standard output, including source patterns,
repair counts, copy bits, and every failed chart in the two-star hostile.
Normal and optimized output agree: **31,940 checks**, including **10,802**
complete small source instances. The LF-normalized JSON SHA256 is
`f0f931d747061254db32a49af78372a554ddc17508011a32a355605691606557`.
An independent proof/API audit and optimized-versus-saved replay passed,
including all fifteen charts and the source-pattern qualifier in §5.

The companion metric control also received an independent coordinate and
support audit: every coordinate was reconstructed from its marked rational
unit phase, scalar modular checks excluded all 6,077 nonedges, exact field
norms authenticated all 139 edges, and all 570 stored colourings were replayed.
This audit preserves the inherited marked field and its good-place premises.
