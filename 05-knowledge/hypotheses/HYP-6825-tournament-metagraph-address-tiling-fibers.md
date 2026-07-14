---
id: HYP-6825
title: Canonical merged-metagraph addresses and exact tiling fibers
status: VERIFIED FINITE ATLAS n=3..7; general completeness conjecture open
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

## Preservation statement

This solves the combinatorial address problem, not the LRC quotient problem.
The raw tournament node still destroys the observer cut, gap metric, exact
threshold side, endpoint owner, scale/residue, wall chronology, and line-orbit
identity.  The tiling fibre repairs the fixed-cut address but not those LRC
sidecars.  HYP-6815's slope suspension should therefore carry the metagraph
address as a constructible chamber label with metric/owner/monodromy stalks,
not identify the 4-coordinate object with a bare node.

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
