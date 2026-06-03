# Unit Distance n=22, Tournaments, LRC, and the Grid Disproof, S614

The current exact small-`n` state is crisp.  Alexeev and Tikhonov determine
`u(n)` through `n=21`, with `u(21)=57`, and put

```text
60 <= u(22) <= 61.
```

They also identify why `n=22` is the first awkward case: the graph-only
`F`-free bound gives `overline u(22)=62`, but the two `62`-edge candidates
contain totally unfaithful subgraphs, so geometry pushes the upper bound down
to `61`.  The lower bound `60` comes from the Moser-ring search/database of
Engel, Hammond-Lee, Su, Varga, and Zsamboki.  Thus the open problem is one bit:
find an embeddable `61`-edge graph, or kill all `61`-edge candidates.

Rebase integration: incoming HYP-2170/S599 should be kept as the
triangular/Harborth lattice baseline.  Its `49` is not the planar optimum claim
for `n=22`; it is the best value in a much thinner carrier, and the contrast
with the published `60`-edge Moser-ring lower bound is exactly the point.

S614 records the first useful reduction from the exact `n=21` frontier.  Any
`61`-edge UDG on `22` vertices has a vertex of degree at most `5`; deleting it
leaves at least `56` edges on `21` vertices.  Since `u(21)=57`, such a graph
has minimum degree `4` or `5`.  So a `61` proof/hunt should not search the plane
raw.  It should search `56/57`-edge `21`-cores plus one degree-`5/4` embeddable
extension vertex.

This is the same shape as the recent LRC `n=14` progress.  A coarse quotient is
good, but only if it is a sufficient quotient.  LRC's `Res_27` shadow must carry
owner labels, carry cocycles, and endpoint/pinch channels.  Unit-distance
graphs need graph coimages plus Moser coordinates, embeddability data, dense
deletion decks, and totally-unfaithful obstruction labels.

The "grid conjecture disproof" should also be read in this language.  The
recent OpenAI-originated, human-verified disproof of Erdos's asymptotic
`n^{1+o(1)}` expectation does not beat the grid by finding a clever low-degree
planar doodle.  The expository note by Alon, Bloom, Gowers, Litt, Sawin,
Shankar, Tsimerman, Wang, and Wood explains the carrier: CM fields, many
norm-one elements, a pigeonhole argument over split prime ideals, and
Golod-Shafarevich class-field towers with bounded root discriminant.  Sawin's
explicit refinement gives more than `n^1.014` unit-distance pairs infinitely
often.  The plane sees one projection coordinate; the counting happens in a
high-degree arithmetic lattice.

That is the shared lesson:

- The visible grid is not the master object.
- The graph-only quotient is not the master object.
- The tournament on points/runners is often not the master object.
- The master object is the retained carrier plus the probes that detect what
  the quotient forgot.

For `u(22)`, that carrier is likely a Moser-coordinate extension/deletion
system with unfaithful-subgraph probes.  For LRC `n=14`, it is the
`Res_27`/owner/carry/pinch statistic.  Tournaments still matter, but the
vertices should be proof obligations, dense cores, obstruction filters, or
side-channel routes.

Sources used this session:

- Alexeev-Tikhonov, "Erdos' Unit Distance Problem for Small Point Sets",
  https://arxiv.org/abs/2412.11914.
- Engel-Hammond-Lee-Su-Varga-Zsamboki, "Diverse beam search to find
  densest-known planar unit distance graphs", https://arxiv.org/abs/2406.15317.
- Alon-Bloom-Gowers-Litt-Sawin-Shankar-Tsimerman-Wang-Wood, "Remarks on the
  disproof of the unit distance conjecture", https://arxiv.org/abs/2605.20695.
- Sawin, "An explicit lower bound for the unit distance problem",
  https://arxiv.org/abs/2605.20579.
