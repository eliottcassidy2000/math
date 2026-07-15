---
id: HYP-6825
title: Canonical merged-metagraph addresses and exact tiling fibers
status: PROVED general Hamiltonian-path fibre inverse (THM-781), cyclic-triangle flow laws (THM-785), and three-sorted tiling/line/node recursion (THM-796) + VERIFIED FINITE objective atlas n=3..7; general structural-address completeness conjecture open
source: codex-2026-07-14-S4
artifacts:
  - 04-computation/tournament_tiling_metagraph_address_codex_S4.py
  - 04-computation/tournament_tiling_metagraph_address_n7_codex_S4.py
  - 05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.out
  - 05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.json
  - 05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.out
  - 05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.json
  - 03-artifacts/visualizations/tournament-tiling-explorer.html
  - 07-reflections/the-metagraph-needs-a-stalk-canonical-addresses-tiling-fibres-and-the-4d-lrc-object-codex-S4.md
  - 00-navigation/METAGRAPH-PRESERVATION-AVENUES-2026-07-14.md
  - 01-canon/theorems/THM-773-prime-seven-sheet-monodromy-and-tournament-fibre.md
  - 04-computation/lrc14_prime7_sheet_monodromy_metagraph_codex_S6.py
  - 05-knowledge/results/lrc14_prime7_sheet_monodromy_metagraph_codex_S6.json
  - 07-reflections/the-heptagon-node-is-a-basepoint-not-the-sheet-state-codex-S6.md
  - 01-canon/theorems/THM-778-centered-christoffel-endpoint-skew-product.md
  - 04-computation/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.py
  - 05-knowledge/results/lrc14_centered_christoffel_endpoint_skew_product_codex_S7.json
  - 00-navigation/LRC14-CONTINUED-FRACTION-FRONTIER-2026-07-14.md
  - 01-canon/theorems/THM-781-hamiltonian-path-inverse-metagraph-fibre.md
  - 04-computation/tournament_tiling_metagraph_fibre_inverse_codex_S8.py
  - 05-knowledge/results/tournament_tiling_metagraph_fibre_inverse_codex_S8.json
  - 01-canon/theorems/THM-785-cyclic-triangle-flow-on-blue-black-metagraph.md
  - 04-computation/merged_metagraph_transitivity_flow_codex_S9.py
  - 05-knowledge/results/merged_metagraph_transitivity_flow_codex_S9.json
  - 07-reflections/symmetry-lives-before-the-quotient-black-drift-as-disintegration-codex-S9.md
  - 01-canon/theorems/THM-796-three-sorted-recursive-tiling-line-node-incidence.md
  - 04-computation/three_sorted_metagraph_recursion_codex_S9.py
  - 05-knowledge/results/three_sorted_metagraph_recursion_codex_S9.json
  - 04-computation/test_tournament_tiling_explorer_line_api_codex_S9.js
related:
  - HYP-2245
  - HYP-2989
  - HYP-3808
  - HYP-3809
  - HYP-3813
  - HYP-6815
  - THM-550
  - THM-646
  - THM-773
  - THM-778
  - THM-781
  - THM-785
  - THM-787
  - THM-796
---

# HYP-6825 — Canonical metagraph addresses and tiling fibers

Every tournament isomorphism-class node in the explorer convention's
converse-merged metagraph at `n=3..7` now has an executable objective address
rooted at the transitive class, and every fixed-path tiling maps exactly to and
from the corresponding node fibre.  The browser explorer itself stops at
`n=6`; an exact refinement-based continuation supplies `n=7` offline.

The address is not raw graph distance.  It is the tuple

~~~text
(one-tile quotient depth,
 lexicographically least spine/rib/sea path word,
 rooted blue/black complement-line distance and word,
 stable weighted blue/black colour-refinement class,
 stable combined-channel colour-refinement class,
 recursive all-vertex-deletion parent address,
 canonical converse-orbit tournament code).
~~~

The canonical code is an admitted final tie-breaker.  Every preceding
controlled-forgetting stage is audited for collisions, so a profile is never
silently promoted to a complete invariant.

The exact relation under study is two-sorted:

~~~text
labelled tiling/tournament realization
  -> canonical tournament isomorphism class
  -> merged complement-orbit node,
~~~

with inverse maps returning fibers, not a fictitious unique tiling.

## The distinction that makes the address work

Two graphs had previously been described with overlapping colour language.
They must remain separate.

1. A **local flip** changes one staircase tile.  Its quotient gives the radial
   coordinate from the transitive class.  Its channels are named `F_S`, `F_R`,
   `F_G` for SC--SC spine, SC--NS rib, and NS--NS sea moves.
2. A **complement line** is the antipodal pair `{t,t xor (2^m-1)}`.  It is
   `L_B` (blue) when the tiling is anti-diagonally grid-symmetric and `L_K`
   (black) otherwise.  These are the blue/black lines drawn by the
   tournament-tiling explorer.

Blue/black complement lines do not express one-step distance from the
transitive tiling.  Conversely, a local-flip distance forgets the antipodal
line parity/folding structure.  The address is therefore bi-axial.

## Exact finite result

The explorer convention fixes the Hamiltonian path `n -> ... -> 1` and orders
the `m=C(n-1,2)` non-path tiles row by row.  Exhaustive enumeration gives:

| n | tilings | classes | merged nodes | blue/black line components | rooted line-WL cells |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 1 | 2 |
| 4 | 8 | 4 | 3 | 1 | 3 |
| 5 | 64 | 12 | 10 | 1 | 10 |
| 6 | 1024 | 56 | 34 | 1 | 34 |
| 7 | 32768 | 456 | 272 | 1 | 272 |

Thus, through `n=7`:

- the projected weighted blue/black complement-line graph is connected;
- raw blue/black incidence is incomplete (`n=7`: 159 cells for 272 nodes,
  largest collision cell 7);
- rooted stable colour refinement of that same weighted graph separates all
  272 merged nodes in three rounds at `n=7`;
- the local quotient distance agrees with minimum labelled tiling Hamming
  weight at every node; and
- self-converse nodes have odd tiling fibre, non-self-converse merged nodes
  have even fibre, supplying an independent checksum.

At `n=7`, the recursive local-depth/path plus all-deletion parent atlas gives
270 cells and leaves exactly two twin pairs.  Within each pair the nodes also
have the same all-black line incidence (`62` for one pair and `86` for the
other), but the next blue/black neighbourhood refinement separates them.  This
is a concrete witness that recursive line placement adds information not
present in deletion ancestry or raw degree.

The line-WL completeness is a finite observation, not a claim for all `n`.
The canonical orbit code makes the address exact even if refinement develops
twins later.

## Exact forward and inverse indices

For a tiling mask `t`, the JSON stores

~~~text
t -> canonical class code -> merged address n-aXX
  -> (global tiling index, fibre-local index).
~~~

Global tilings are ordered by

~~~text
(Hamming distance from all-off transitive tiling,
 merged-node address rank,
 canonical unmerged class code,
 explorer mask).
~~~

The fibre-local order uses `(Hamming weight, class code, mask)`.  Conversely,
each node record lists every tiling mask and global index in its fibre.  The
committed JSON atlases contain all `33866` tilings at `n=3..7`.  The browser
explorer loads the `n<=6` atlas to display the address and both indices; the
compact `n=7` atlas stores mask-indexed class/node/global/fibre arrays and
complete inverse node fibres.

## Direct inverse theorem: paths modulo automorphisms

THM-781 removes the exhaustive-cube scan from the mathematical definition of
the inverse.  If a merged node carries class representatives `T` (one for a
self-converse class, two for a converse pair), enumerate their directed
Hamiltonian paths, relabel each path to the explorer's fixed path, and read the
remaining arcs as tile bits.  The exact functions are

```text
tiling t -> A_n({[tau(t)],[tau(t)^op]}),

node u -> union_[T] HP(T)/Aut(T).
```

Two paths produce the same tiling exactly when an automorphism relates them.
Thus an unmerged class has `H(T)/|Aut(T)|` tilings; a non-self-converse merged
pair has twice this number.  This is a theorem for every size once its class
and node atlas is defined, not merely an observation through `n=7`.

The executable API `MetagraphFibreAtlas.tiling_to_node` computes the forward
quotient and `node_to_tilings` constructs the inverse directly from paths.
All `33,866` stored tilings, `530` classes, and `321` nodes through `n=7` pass
both directions with zero mismatches.  The browser exposes the analogous
`tilingToMergedNode` and `mergedNodeToTilings` functions through `n=6`.

For the prime-seven pullback this explains the previously mysterious fibre
size without scanning masks:

```text
n7-a267: H=175, |Aut|=7, H/|Aut|=25.
```

The 25 masks are precisely automorphism orbits of Hamiltonian observer cuts on
the cyclic heptagon.  THM-773's 5,040 owner-sheet assignments select these cut
orbits with nonuniform multiplicity; the node alone still forgets the owner
assignment and transport future.

## Prime-seven sheet pullback: the node is a base, not the state

THM-773 gives the first exact functor from an LRC sheet tiling into this
atlas.  At `c=7`, seven unramified owner tokens form an exact cover precisely
when they are a permutation of `F_7`.  The marked-cut tournament sends all
`7!=5040` assignments to transitive node `n7-a000`; circular three-step
precedence sends all of them to heptagon node `n7-a267`.  Choosing the least
Hamiltonian path maps the 5040 assignments onto exactly the 25 masks in
`a267`'s stored inverse fibre.  The two tournament gauges differ on the six
long chords, matching `a267`'s local depth six.

This exact success is also a sharp controlled-forgetting failure.  Two speed
rows can have the same node, same mask, same six finite-field moments, and the
same owner-to-sheet assignment but different next event owner and free sheet.
Thus the inverse fibre is not an LRC continuation state.  Inverse windings,
endpoint order/phase, metric base, and the global `x -> x+1` sheet carry remain
mandatory transport fields.

## Continued-fraction path stalk

THM-778 now gives an exact address for the missing endpoint-order field.  The
merge of every two owner midpoint clocks is a centered rational mechanical
word; a one-bit parity cocycle survives its Euclidean shears.  For owner-local
event `(a,i)`, the centered Beatty rank

```text
i + sum_(b != a) ceil((w_b(2i+1)-w_a)/(2w_a))
```

is its exact global simultaneous-wall index.  Thus the full wall schedule can
be mapped forward and backward without discarding owner labels or tie blocks.

The natural next-event tournament on owner clocks is transitive in every
chamber.  On HYP-6835's eight-owner example it therefore stays at one ordinary
isomorphism node while its labelled Hamiltonian path takes `948` values and
flips `6,620` pair edges.  This is a particularly sharp answer to the original
mapping question: the node is objective but static; continued-fraction
substitutions act on the labelled path/tiling stalk above it.  The promising
ordering is consequently lexicographic in

```text
(merged-node address, Euclidean endpoint-block address,
 labelled path/mask index, owner-token/redundancy state),
```

not in merged-node rank alone.

The ten covered walls of the named eight-owner movie have mask word
`(25773,32153,31115,14635,615,30093,31115,615,14233,6035)`.  Its owner word is
palindromic but its mask word is not.  Exhausting all 5,040 assignments proves
that the reflected-mask relation is multivalued: only 9/25 masks have a unique
image and the largest image fibre has size 7.  Therefore even fibre-local mask
index is not a reflection-equivariant state unless the owner-labelled lift is
retained.

THM-779 supplies the predicate on this path stalk: continued blocking requires
the chronological next-wall owner to agree with the current collision-pair hop
target.  In the ten-wall movie every covered-wall gap contains at least three
other walls, so these are isolated first returns rather than a consecutive
blocking chain.  A continuation-complete address must therefore distinguish
the current `SURJ`/hop state from the later return mask.

## Cyclic-triangle flow order

THM-785 adds the missing mathematical horizontal coordinate to the exact
node/fibre dictionary.  If `d_i` are the scores, then

```text
sum_i (d_i-(n-1)/2)^2 = n(n^2-1)/12 - 2 C3.
```

Thus `C3` orders the merged nodes from the unique transitive root to exactly
the regular/near-regular score locus.  On an explorer complement line,

```text
Delta C3 = d_0-d_(n-1)-1.
```

Blue symmetry forces `d_i+d_(n-1-i)=n-1` and reduces this to `2d_0-n`;
black lines retain the transverse defect
`epsilon=d_0+d_(n-1)-(n-1)`.  THM-785 proves closed binomial laws for the blue
step distribution, the full step distribution, and the black absolute-defect
distribution for every `n`.

The categorical line topology is exactly

```text
pure_blue --blue-- mixed --black-- pure_black,
```

but increasing `C3` need not follow that categorical arrow.  At `n=7`, all
six blue boundary supports point from pure-blue to mixed, whereas black
boundary flow has `2798` pure-black-to-mixed instances against `1254` in the
opposite direction (`1075` versus `522` projected supports).  The labelled
line ensemble remains centrally symmetric; the asymmetry appears after
iso-class projection and `C3` orientation.

Every node now has a stored `flow_rank` ordered first by `C3`, then phase and
rooted blue/black position, score/Landau shape, tiling-fibre mass, weighted
line refinement, and finally the exact HYP-6825 address.  All balanced nodes
through `n=7` have a
nondecreasing `C3` path whose colour word is `B* K*`.  Extending that reach law
beyond `n=7` remains open.

THM-787 independently recovers the same axis as the integral energy
`E4=E4(transitive)-8C3` and supplies the complementary unmerged-class
histograms and depth-pair current table.  Its `368` pure-black classes at
`n=7` are the `184` converse-merged pure-black nodes here.  THM-785's closed
blue law promotes THM-787's formerly finite parity/max pattern to all sizes.

## Preservation statement

This solves the combinatorial address problem, not the LRC quotient problem.
The raw tournament node still destroys the observer cut, gap metric, exact
threshold side, endpoint owner, scale/residue, wall chronology, and line-orbit
identity.  The tiling fibre repairs the fixed-cut address but not those LRC
sidecars.  HYP-6815's slope suspension should therefore carry the metagraph
address as a constructible chamber label with metric/owner/monodromy stalks,
not identify the 4-coordinate object with a bare node.

## Three-sorted recursive incidence

THM-796 separates the atlas into tiling endpoints `T_n`, complement-line
instances `E_n`, and converse-merged nodes `N_n`, with exact inverse fibres at
every arrow.  Its tiling recursion is the two-face pullback

```text
T_n = (T_(n-1) x_[T_(n-2)] T_(n-1)) x {0,1}_apex.
```

Choosing the apex-zero endpoint identifies `E_n` with the compatible tiling-
face pairs.  Passing both faces to bare lower lines loses one coherent endpoint
phase: for `n>=5`, `E_n` is a uniform two-sheeted torsor over
`E_(n-1) x_[E_(n-2)] E_(n-1)`.

At node level the face is a weighted correspondence `D_n`, not a function.
Its primitive row separates all nodes through `n=7` (support alone leaves eight
twin pairs there), but strong lumpability fails on `1206/1312` nonzero `n=7`
node blocks.  Thus a complete static node fingerprint is still not a complete
continuation state: the lower tiling/Hamiltonian-path orbit must survive.

The blue/black face recursion is closed for all `n`.  Upper-, low-face-, and
high-face-blue are pairwise independent, although upper blue forces the two
face colours equal; the remaining information is a pure three-way interaction.
This exact pre-quotient symmetry locates THM-785's black left/right imbalance
in the disintegration over unequal node fibres rather than in the raw line
tower.

HYP-6865's concurrent Smith-diagram audit gives the horizontal flow order a
second interpretation: harmonic voltage from the transitive class to the
distributed rail is perfectly concordant with score variance through `n=6`
and `99.9%` concordant at `n=7`.  It lives on the unmerged local-flip graph, not
the complement-line graph; pairing that voltage with THM-796's primitive face
vector is a promising two-axis address rather than a proved graph identity.

## Tournament Analysis

Vertices are candidate information carriers rather than runners:

~~~text
raw SC type, radial depth, blue/black incidence, rooted blue/black WL,
local parent atlas, combined coloured WL, structural address, exact address.
~~~

The pairwise observable is how many unordered node pairs each carrier
separates.  The retention switch and retention-per-description-bit switch both
produce transitive tournaments at `n=3..7`, but their edge order differs (22
edge flips at `n=7`).  At `n=7`, separated-pair counts advance

~~~text
16192 (SC type), 27974 (depth), 36636 (line incidence),
36856 (rooted line-WL, the information ceiling).
~~~

The tie Hamiltonian path is the carrier list above.  At `n=7` each gauge has
score histogram `{0:1,...,7:1}`, zero directed 3-cycles, eight singleton SCCs,
and one Hamiltonian path.  The challenged assumption is that tournament
vertices must be runners or arcs; here they are quotient carriers/proof
obligations.

## Open boundary

1. Test whether rooted weighted blue/black 1-WL remains complete at `n=8`;
   if not, classify the first twins by line-class, spectrum, or 2-WL.  Existing
   nauty class machinery may avoid a full `8!` per-tiling canonicalization.
2. Prove or refute connectivity of the projected blue/black line graph for all
   `n`.
3. Prove that local quotient depth always equals minimum tiling weight/MFAS;
   the finite equality alone does not justify the general statement.
4. Replace raw all-deletion parent codes by a genuine functor on the whole
   recursive metagraph tower, retaining deletion multiplicities and orbits.
5. Pull the exact node/fibre atlas back over the LRC slope suspension and test
   fibre-purity of threshold nonemptiness after adding metric gap, owner,
   carry, scale, and wall-monodromy sidecars.
