# Counting exact port-observer fibres by colored modular trees

**Status: PROVED using the cited cyclic-triple realization theorem; independent
review pending.** This counts **labelled** relative fibres, not unlabelled
isomorphism classes, new A000568 values, or attainable odd counts.

The new observer rigidity theorem in [PROOF.md](PROOF.md), Theorems9–11,
classifies equality of complete eliminated boundary matrices, with named
interior vertices and fixed arcs to/from the retained ports. Here we count
those equivalence fibres exactly from the old modular decomposition.

## 1. The native colors

Each interior vertex v has color sigma(v), its full binary vector of arc
orientations to the retained ports. No arbitrary orientation of ties is used.
Fix the color at every labelled vertex. A pair u,v is **cross-color** if
sigma(u) differs from sigma(v). A port distinguishes u,v iff this holds.

Under fixed port incidences, equality of the relative augmented cyclic-triple
hypergraph means exactly

1. the interior tournaments have the same labelled cyclic triples;
2. every cross-color pair has the same orientation.

Indeed if a port distinguishes u,v, reversing uv toggles that triple between
cyclic and transitive. If it does not distinguish them, the port is a source
or sink on the triple, so flipping uv does not change its cyclic status.
Triples with two or three ports are already fixed. By the new observer
rigidity theorem, these are exactly the full boundary-equivalent interiors.

## 2. Inherited modular fibre theorem

We use the published cyclic-triple realization theorem as routed through
[THM-4803, corrected scope](../../../../01-canon/theorems/THM-4803-intransitivity-is-the-relation-modulo-its-modular-tree-and-a000568-from-the-prime-cycle-index-census.md)
and [the corrected path-cover composition proof](../../tournament_composition_classification_20261026/fibre/RESULT.md).
Its realization fibre is the common canonical strong-module tree with:

- each prime quotient independently retained or dualized;
- each linear quotient's children independently permuted.

These labelled operations are independent and faithful. This published
classification, and not a new proof of its structural direction, is the
input here. The old assertions about automorphism orders and ordinary SCC
isotypes were refuted; neither is used.

Let a tree node be **monochromatic** if all its labelled leaves have one
color. At a linear node, read its child modules in their original total
order. Call a child pure if it is monochromatic. Split this list into maximal
consecutive runs of pure children of one color. Mixed children act as barriers.
Let R range over all such pure runs at all linear nodes. Let p_mono be the
number of monochromatic prime nodes.

**Theorem (relative fibre count).** The exact number of labelled interiors
boundary-equivalent to T with those port incidences is

    F_boundary(T,sigma) = 2^p_mono · product_R |R|!.

This is also an explicit recursive grammar for the entire fibre: dualize
only monochromatic prime nodes and permute only the children inside the
specified pure runs. It is not merely a numerical count.

**Proof.** At a prime node, dualization reverses every arc between distinct
children. Such a reversal preserves all cross-color pairs iff every pair
of vertices in different children has equal color. With at least two
nonempty children this is equivalent to the entire node being monochromatic:
a mixed child cannot agree with every vertex of another child. Thus exactly
the p_mono prime choices survive.

At a linear node, any inversion of two children reverses all arcs between
them. It is allowed iff the two children are both pure and of the same color.
A mixed child therefore cannot cross any other child, and children of different
colors cannot cross. Pure children of the same color separated by such a
barrier cannot exchange positions without crossing the barrier. Consequently
the permissible permutations are exactly the independent full permutations
inside maximal consecutive pure runs, giving |R|! each.

Every vertex pair's orientation is controlled at its least common ancestor.
Changes strictly below that ancestor cannot cancel its reversal. Hence the
nodewise restrictions are both necessary and sufficient, remain independent,
and cannot give duplicate labelled outputs. Multiplying their sizes proves
the formula and the generating grammar. □

## 3. Controls and what marking does

With no ports there is one empty-vector color. Every node is monochromatic;
the formula reduces to the inherited unmarked cyclic-triple fibre count
2^(#prime nodes) product_(linear nodes) arity!.
If every vertex has a different color, the relative fibre is a singleton.
For a transitive interior, the tree has one linear node with singleton
children: its fibre size is product of factorials of the consecutive color
run lengths. Thus the placement of equal colors, not just their multiplicity,
matters. The words0011 and0101 give sizes4 and1, respectively.

Adding more retained ports refines the vertex colors, so it can only remove
allowed node operations. This is a literal recursive symmetry-breaking law.
If a newly added port is constant on every monochromatic prime node and on
every pure linear run of length at least two in the old fibre, its bit does not split any allowable
operation and the fibre size is unchanged. Otherwise at least one old
operation is killed and the fibre shrinks. The test is exact and finite on
the marked tree; it is not an asymptotic estimate of automorphism order.

These boundary fibres are not automorphism groups: their members may have
different automorphism orders (the corrected THM-4803 order12 hostile).
They are classes of replacement interiors with identical response to all
continuations attaching only at the retained ports. Scalar H alone identifies
many more objects and is not the equivalence used here.

## 4. One designed port, or a sharp logarithmic random budget

**Theorem (one engineered discriminator).** For every interior tournament T
there is a choice of one new retained port's incidence bits for which
F_boundary(T,sigma)=1. For n>=2, one port is optimal for this objective:
without ports the inherited cyclic-triple fibre has at least two members.
The one-vertex case already has a singleton fibre without a port.

**Construction and proof.** Color the modular tree recursively with two colors.
Require every nontrivial child module to be mixed, using this same construction
inside it. At a linear node, alternate the colors on each consecutive run of
singleton children; mixed children separate those runs. At a prime node with
only singleton children, use both colors; if it has a nontrivial child it is
already mixed, and remaining singleton colors can be arbitrary. This makes
every nontrivial node mixed. No prime node is monochromatic, and the only pure
children at any linear node are singletons, whose adjacent colors differ.
The fibre formula is therefore1. The bit assignment is an actual orientation
of the new port's edge to each interior vertex, hence a physical tournament
extension. No probability or unspecified external gadget is needed. □

This **does not** make the extension prime as an unmarked tournament. New
modules containing the new port can occur. It makes the *relative fibre*
unique when its prescribed port incidences and augmented cyclic triples are
observed. It also does not preserve H when the port is added or select an
arbitrary target count. This is a discriminator construction, not a +2 gadget.

**Theorem (sharp random marking scale).** Add k independently uniform random
ports, meaning each of the kn incidence bits is independent and fair. Arcs
among the new retained ports can be fixed arbitrarily. For n>=1,

    Pr(F_boundary>1) <= (n−1) 2^(−k).

For transitive interiors the exact uniqueness probability is

    Pr(F_boundary=1) = (1−2^(−k))^(n−1).

Thus log2(n)+O(1) ports is the sharp worst-case scale for a fixed nontrivial
success probability; uniform success tending to1 requires an unbounded
additive margin in the worst-case transitive family. One designed port and
randomly chosen ports therefore have very different information efficiency.

**Proof.** Write q=2^k. The vertex colors are independent uniform values in
[q]. A nontrivial fibre occurs exactly when some prime node is monochromatic,
or some adjacent pair of children of a linear node has a monochromatic union.
Every such support has at least two leaves, so its probability is at most1/q.
The number of events is at most

    #prime nodes + sum_(linear nodes)(arity−1)
      <= sum_(all internal nodes)(arity−1)=n−1.

The union bound proves the first statement, without pretending these events
independent. For a transitive interior the only requirement is that consecutive
singleton colors differ. There are q(q−1)^(n−1) such color words among q^n,
proving the exact formula. □

A more precise additive quantity also has an exact expression. Write log=log2.
For each prime node C contribute q^(1−|C|). For every consecutive interval J
of ell>=2 children at a linear node contribute

    log(ell/(ell−1)) q^(1−sum_(child in J)|child|).

The sum of these contributions is exactly E[log F_boundary]. Indeed an
interval contributes precisely when its leaves have one color, and
log(r!)=sum_(ell=2)^r (r−ell+1)log(ell/(ell−1)). This is an expectation of
relative-fibre information, not a formula for the tournament's automorphisms,
Hamiltonian count, or A000568. It gives a structural rather than size-only
measurement of which symmetries persist under marking.

## 5. Reproduction and continuation

`colored_fibers.py` independently builds all labelled interior tournaments
through five vertices, computes every strong module by literal subset tests,
and compares the displayed product against direct grouping by cyclic-triple
bits and fixed cross-color arcs. It checks all2- and3-colorings through order4,
and all2-colorings at order5:39,287 labelled color/tournament cases.
Each of the1,099 labelled tournaments through order5 also gets its explicitly
constructed discriminator; the event criterion, random union bound and exact
transitive probability are checked over every tested coloring.
Colors are abstract distinct port signatures; any finite palette can be
encoded by sufficiently many binary port bits. Port–port arcs are irrelevant
to this count once fixed.

The constructive lesson is exact. Refining marks predictably shrinks the
replacement grammar, and full relative-C3 conservation freezes all endpoint
responses, not merely H. To gain two paths, a proposed move must leave this
relative fibre by changing a cyclic triple visible to the retained boundary.
The five-term elimination atlas and signed reciprocal response law in
PROOF.md can track the resulting change, but do not yet construct the needed
renewable target landing. Chan–Pak remains OPEN.

## 6. Deliberately relaxing guards produces a comparison-poset grammar

**Theorem (sparse relative guards).** Keep the full labelled interior
cyclic-triple hypergraph and fixed port incidences, but require equality only
of a chosen collection of quadratic port-start coefficients. Discard chosen
queries whose port does not split the queried interior pair, since these are
vacuous. Each remaining query fixes the orientation of its interior pair.
Assign that comparison to the least common ancestor of the pair's leaves
in the inherited common modular tree.

At a prime node, retain its two orientation choices iff no assigned comparison
crosses two of its children; otherwise the original quotient orientation is
forced. At a linear node, put the assigned child comparisons into a partial
order P_v, consistent with the original total order. Its allowed quotient
orders are exactly the linear extensions of P_v. Therefore the exact labelled
fibre size is

    2^(#unconstrained prime nodes) · product_(linear nodes v) e(P_v),

where e(P_v) counts linear extensions. This also explicitly generates the
physical interior tournaments satisfying these sparse guards.

**Proof.** A prime quotient's only choice is complete dualization, reversing
every child-to-child arc, so one comparison rules it out. For a linear
quotient, satisfying the retained comparisons means exactly choosing a total
order extending their acyclic transitive closure. Each leaf pair's orientation
is controlled at its LCA, so the constraints separate by node; lower changes
cannot undo an ancestor comparison. The inherited faithful fibre
parametrization gives the product and no duplicate labelled outputs. □

Full quadratic port data recovers the preceding colored-run formula. Sparse
guards preserve interior H but need not preserve the full boundary or the
Hamilton count of an attached system: that is precisely
why this is a lawful space of candidate compensating moves rather than the
refuted frozen-observer seam. The missing responses must be carried and
controlled, not assumed constant.

**At most n−1 targeted quadratic queries suffice with the engineered port.**
For every prime node, select one cross-color pair in distinct children; such
a pair exists because that node is mixed. At every linear node, select one
cross-color pair from each adjacent pair of children; their union is not
monochromatic by construction. Querying these pairs locks every prime
orientation and every adjacent comparison in every linear order. The
number of queries is

    #prime nodes + sum_(linear nodes)(arity−1) <= n−1.

With the internal cyclic-triple inventory and fixed port incidences given,
these particular quadratic observations therefore determine the entire
labelled interior, not only its boundary response. Selection uses the known
modular tree, so this is a designed certificate, not a claim that arbitrary
n−1 queries suffice. For a transitive interior all n−1 adjacent comparisons
are necessary: leaving one unqueried permits swapping that adjacent pair
without changing any other pair comparison. This contrasts a deliberately
chosen certificate with indiscriminate scalar sampling.

The next constructive experiment can now retain a chosen subset of this
certificate and enumerate the remaining linear extensions/prime reversals.
Every candidate is a genuine tournament, so positivity and missing-edge
problems of formal network synthesis are avoided at the generation stage.
Exact count/endpoint gains and renewable target-decreasing readiness remain
the hard obligations. No closure of the473 family repairs is asserted.
