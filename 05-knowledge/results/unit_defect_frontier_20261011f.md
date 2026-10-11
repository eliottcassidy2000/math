# Prime-screen escape, attainable displacements, and the colouring-defect frontier

**PROVED:** the elementary deletion-context antichain identities, finite-prime
escape supplier, necessary 2-adic norm guard and Mycielski geometry obstruction.
**INHERITED:** the conjugation-fixed-prime colouring mechanism of THM-4558.
**FINITE-EXACT:** declared parameter searches, finite profiles and three exact
point configurations. **OPEN:** a six-colour forced-equality unit fragment.
All three newly constructed configurations below have explicit four-colourings.

## Inheritance and the actual construction objective

The immediate suppliers are [conditional counts](boundary_context_weights_20261011d.md),
[the fixed-universe zero threshold](unit_equality_quantized_guards_20261011e.md),
and [the Moser-field obstruction](unit_equality_candidates_20261011e.md).
The current field theorem is
[THM-4558, residue colourings](../../01-canon/theorems/THM-4558-residue-colourings-of-field-planes-the-heegner-tower-is-at-most-four-chromatic.md),
including its repaired rotation/Heegner claims (MISTAKE-576/579).
The new incoming [golden/taxicab boundary note](golden_taxicab_boundary_20261011.md)
constructs nontrivial precoloured interfaces; its private leaves do not force
their own prescribed colours in an unprecoloured graph.

Closest mechanism: preserve the complete boundary support before composition.
Hostile: larger field degree and more points can leave the chromatic number
unchanged. Corrected near miss: a norm-one rotation need not close a displacement
attainable in its proposed seed field. Least-used sidecar: the actual norm
valuation at the prime behind a field-colouring certificate.

Portfolio: six-colour equality / finite-prime escape / defect-antichain recursion.
Board: **unit incidence, marked boundary, defect support, prime decomposition,
attainable displacement, actual colour witness**.

## 1. Retain the obstruction support, not just its count

Let G be a fixed finite graph, k a palette and c any assignment of at most k
colour names. Its defect set is

\[
 M_G(c)=\{xy\in E(G):c(x)=c(y)\}.
\]

The largest spanning subgraph preserving c as a proper colouring is exactly
F_c=G-M_G(c). Consequently, when searching for an equality fragment on a
fixed point set with a chosen colouring, it is sufficient to test F_c:
omitting additional compatible unit edges cannot eliminate more separating
colourings. This does not supply a suitable c or G.

Let A_k(G) be the inclusion-minimal members of the family of all M_G(c).
For any deleted edge set D,

\[
 G-D\text{ is k-colourable}\quad\Longleftrightarrow\quad
       \exists A\in A_k(G),\ A\subseteq D.                 \tag{1}
\]

Every defect contains an inclusion-minimal one, proving both directions.
For any nonempty A in this antichain, G-A has a colouring, and every one
of its k-colourings is monochromatic on **every** edge of A: otherwise
its defect set would be a proper subset of A. Thus each deleted edge is
an equality port pair. If G is a unit-distance graph, these are genuine
unit-separated equality fragments.

Conversely, a forced-equal omitted unit edge lies in a minimal defect
contained in the deleted set of its full unit graph. This is a reformulation
of critical obstruction extraction, not a new existence proof. When G already
has a proper k-colouring, the antichain is just {empty}, which excludes
forced equality on every edge removed from E(G). If G is the **full induced
unit graph**, this excludes every unit-separated equality fragment on that
point set. Full incidence is necessary for the latter inference: the two-edge
path 0,a,1 is two-colourable and forces its unit-separated ports equal, but
their unit edge is absent from that source graph. This is the same positive
control retained in the companion compiler, not a counterexample to (1).

The antichain is lossless for all spanning edge-deletion colourability queries;
it discards extension multiplicities and the actual colour assignments. Keep
the graph, palette, port labels and a literal witness as sidecars. Its size
can be exponential. A short description here is not a general efficient
colouring algorithm.

### A recursion that preserves this support exactly

With a boundary B, keep an antichain A_G(pi) separately for every equality
partition pi of its labelled ports. For two graphs with exactly B in common,
the antichain for their union is

\[
 A_{G\cup H}(\pi)=\min_{\subseteq}
        \{D\cup E:D\in A_G(\pi),E\in A_H(\pi)\}.          \tag{2}
\]

Align colour names on B before gluing; private interiors are otherwise
independent. Taking minimal unions then proves (2). Forgetting a port takes
the minimal union over its permitted extension partitions. The operations
are set union and removal of dominated supports, rather than scalar addition
or multiplication of counts. This gives a concrete recursive representation
of the proof obstruction.

The smallest hostile is a two-edge path and one edge sharing their endpoints
at palette two. Individually they both have empty minimal defect. Gluing while
forgetting the marked boundary would incorrectly retain that empty defect.
Keeping the two port patterns instead recovers the triangle obstruction.
The Moser spindle has eleven singleton minimal defects at palette three,
but only the empty defect at palettes four and six. Abstract K7 has21
singleton minimal defects at palette six and no injective plane unit embedding.

## 2. One rotation can evade any prescribed finite prime bank

For a finite set S of rational primes put

\[
 P=\prod_{p\in S,\ p\ne2}p,\quad t=2P,\quad D=4t-1=8P-1,
 \qquad \omega_t=\frac{2t-1+\sqrt{-D}}{2t}.                \tag{3}
\]

The empty product is one. The norm of omega_t is exactly one. The imaginary
quadratic field K=Q(sqrt(-D)) splits at every prime in S. At odd p in S,
-D is1 modulo p and p does not divide D; the two roots lift by Hensel's
lemma. At2, -D is1 modulo8 and is a square in Q_2. Square factors in D do
not change this conclusion or its quadratic field.

Let L be any conjugation-stable number field containing K. There is no
complex-conjugation-fixed prime of L over any p in S. Such a prime would
restrict to a conjugation-fixed prime of K, whereas the two primes of the
split quadratic extension are exchanged by conjugation.

Thus (3) removes **all instances of the inherited fixed-prime residue screen
over S**, including higher residue-degree versions. It says nothing about
primes outside S or other colouring methods. It supplies no chromatic lower
bound. It does give an explicit algebraic rotation beyond any preselected
finite bank of these certificates.

For S consisting of all primes at most19, t=9,699,690 and D=38,798,759.
This universal construction is intentionally simple, not optimized.

## 3. Escape angles must also fit a real geometric displacement

For the original Moser field M=Q(sqrt(-3),sqrt(-11)), any nonzero z with
z conjugate(z) rational satisfies

\[
                         v_2(z\bar z)\equiv0\pmod2.     \tag{4}
\]

Indeed primes over2 are unramified and fixed by complex conjugation in M,
as retained in THM-4558 and the explicit two-F4-frame construction. Hence
v_P(z conjugate(z))=2v_P(z), while for a rational number its valuation is
its ordinary v_2. This is a necessary norm condition, not a complete norm
classification.

The script checks even integer t through120. It searches for a **displayed
split imaginary quadratic subfield** among Q(sqrt(-3)), Q(sqrt(-11)) and
Q(sqrt(-(4t-1))) at every prime at most19. The passing values are

\[
                    22,32,38,46,60,66,68,82,90.
\]

This particular sufficient split-certificate test is not advertised as a
complete local decomposition algorithm. Its first entry22 violates (4),
as do32,38 and46. Its first entry surviving the norm guard is60, and this
one has an actual displacement witness:

\[
 z=\sqrt{33}+3\sqrt{-3},\qquad z\bar z=60,
 \qquad\omega_{60}=\frac{119+\sqrt{-239}}{120},
 \qquad|z-\omega_{60}z|=1.                               \tag{5}
\]

With a=(1+sqrt(-3))/2 and b=(5+sqrt(-11))/6, the same endpoint is
z=16a-8-12ab+6b. More economically, writing s=sqrt(-3),

\[
                  z=3s(1-b+\bar b).
\]

It can therefore be reached by nine consecutive unit rhombi with opposite
port displacements s,s,s,-sb,-sb,-sb,s conjugate(b),s conjugate(b),s conjugate(b).
This certifies geometric reachability, not a six-colour equality relation.

## 4. Three actual configurations, with every unit edge retained

The [candidate script](../../04-computation/experiments/unit_prime_escape_candidate_20261011f.py)
constructs the bent chain above in the exact eight-dimensional rational basis
of Q(sqrt(-3),sqrt(-11),sqrt(-239)). It then joins its omega_60-rotated copy
and finally takes the Minkowski sum with the unit triangle {0,1,a}.
Coincident points are identified exactly before finding **all** unit pairs.
There is no approximate distance tolerance or assumed Cartesian-product
incidence.

| Actual point set | Vertices | Unit edges | Certified chromatic number |
|---|---:|---:|---:|
| Bent nine-rhombus chain |28|46|4|
| Chain and its omega_60 copy |55|93|4|
| Triangular thickening of the joined set |153|405|4|

The base chain has an extra unit edge between corners6 and12, which every
three-colouring of its constituent rhombi would identify. It is therefore
not three-colourable. The joined and thickened configurations retain this
obstruction. For each set the script supplies an independently checked
four-colouring of the complete exact unit graph. These are successful
geometric screen escapes and unsuccessful six-equality candidates.
Every omitted unit pair is separated by that same four-colouring.

The solver is bounded DSATUR used only to find witnesses. Its lack of a
witness would be labelled unresolved, not an uncolourability proof. In these
three cases it finds a four-colouring after28,55 and211 search nodes.
The saved output records every rational coordinate and every colouring.

## 5. Two quantified recursion obstructions

First, the [diamond-chain packet](unit_field_escape_recursion_20261011f.md)
proves an all-m family with6m+1 vertices,10m+1 unit edges and chromatic number4.
Take m=2P for the prime product in (3). Its new rotation discriminant is
D=12m^2-1=48P^2-1, so its quadratic field splits at every prime in S by
the same proof. Thus **for every finite prime bank there is an actual
four-chromatic unit configuration whose field evades all those fixed-prime
screens**. This simultaneously gives a norm-compatible escape supplier and
an all-size hostile to inferring chromatic growth from screen escape.

Second, ordinary Mycielski lifting of any graph with a vertex v of degree
at least three cannot have an injective planar unit-distance realization.
The old vertex v and its clone share all original neighbours of v, giving
three distinct common unit neighbours. Two distinct unit circles meet in
at most two points. In particular the familiar abstract chromatic-raising
recursion cannot simply be iterated from a four-chromatic unit seed.
Changing the recursion must also change its local incidence, not merely
its vertex coordinates.

## Connection ledger and remaining construction requirement

| Source -> target | Map / preserved predicate | Lost information and sidecar | Cheapest decisive test |
|---|---|---|---|
| Colourings -> defect antichain | Minimal monochromatic edge sets preserve all deletion feasibility queries | Counts, geometry and literal witnesses | Compare every small deletion context |
| Two interfaces -> union | Marked antichain product(2) | Port forgetting loses compatibility | Path plus edge at palette2 |
| Prime bank -> rotation | Formula(3) removes its conjugation-fixed primes | No unit incidence or chromatic lower bound | The all-m four-chromatic family |
| Rotation -> spindle endpoint | Authenticate norm and closing distance | Norm-one alone need not supply a seed displacement | t22 fails(4); t60 has(5) |
| Abstract chromatic recursion -> geometry | Mycielski preserves the abstract construction | Three-common-neighbour obstruction | One degree-three source vertex |

The missing object is a rigid six-colour relation in an actual geometric
source, or an assembly whose retained full boundary relations eliminate
every separating assignment. The all-distance compiler in
[the recursion atlas](recursion_escape_atlas_20261011f.md) removes the need
to prescribe its nonzero port distance in advance. It does not supply that
rigidity. The new candidates demonstrate why increasing field size, vertex
count, or path length does not itself create the missing relation.

## Reproduction and scope

Run both named scripts in `04-computation/experiments/`, with ordinary Python
and `python -O`. The defect script checks all labelled graphs on0..4 vertices
at palettes2..4 and every edge-deletion context:2,283 independent literal
comparisons, plus320 marked gluing cases, Moser/K7 profiles, exact prime
suppliers, norm identities and typed hostiles. Total2,901 checks.
An independent audit caught Python equality aliases in the initial profile
validator: Boolean pattern labels, float edge labels and mutable defect sets
could compare equal to an authentic packet. The repaired validator checks
every exact entry type and source edge before reauthentication; nine new
hostiles retain this repair. The mathematical profile identities were unchanged.
The geometric script independently authenticates its basis arithmetic,
rhombi, exact full incidence, closing edge, colourings and malformed inputs:
109 checks. Neither finite census is a universal graph search.

LF-normalized output SHA256:

- `unit_defect_frontier_20261011f.out`:
  `1731a3f774f854458f70189bb3ef5c16e6103e8d22213c928521c3ed7738a264`.
- `unit_prime_escape_candidate_20261011f.out`:
  `3f664ada1cfa2544b64fb16f840429c76bcc9f72ee0c45c36504d454fdc066e8`.
