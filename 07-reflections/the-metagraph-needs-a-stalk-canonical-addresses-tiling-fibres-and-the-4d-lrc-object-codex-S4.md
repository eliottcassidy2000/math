# The metagraph needs a stalk: canonical addresses, exact tiling fibres, and the four-coordinate LRC object

**codex-2026-07-14-S4. Status:** exact finite tournament atlas at `n=3..7`,
information-preservation synthesis, and a research program; not a proof of
LRC(14).

## The new finite fact

The tournament-tiling explorer has always contained an exact map

```text
fixed-path tiling -> labelled tournament -> isomorphism class.
```

What it lacked was a structural address on the converse-merged class nodes and
an explicit inverse-fibre index.  HYP-6825 now supplies both.

At `n=7`, the fixed-path cube has `32768` tilings, `456` tournament classes,
and `272` converse-merged nodes.  The weighted blue/black complement-line graph
on those 272 nodes is connected.  Raw blue/black incidence produces only 159
profiles, with collision cells as large as 7.  Root the graph at the transitive
node and repeatedly refine a node by the multiset of colours and multiplicities
of its neighbours.  Three rounds of rooted weighted 1-WL refinement produce
**272 different colours**.  They separate every node.

The same is true at every explorer size:

| n | merged nodes | raw line-incidence cells | rooted line-WL cells |
|---:|---:|---:|---:|
| 3 | 2 | 1 | 2 |
| 4 | 3 | 3 | 3 |
| 5 | 10 | 9 | 10 |
| 6 | 34 | 26 | 34 |
| 7 | 272 | 159 | 272 |

This is not a theorem for arbitrary `n`.  It is nevertheless meaningful: the
blue/black structure is more than a decorative parity overlay.  Once rooted
and read recursively, it is a complete fingerprint of every merged node the
explorer currently knows.

## The colour collision that had to be resolved first

Two different graphs in the repository have both been narrated as a blue/black
geometry.

1. The **local flip graph** changes one non-path staircase tile.  It describes
   how classes spread from the all-off transitive tiling.  After converse
   merging, its three honest channels are
   `F_S=SC--SC`, `F_R=SC--NS`, and `F_G=NS--NS`: spine, rib, and sea.
2. The explorer's **complement-line graph** pairs a tiling with the tiling
   obtained by flipping every tile.  Its line is blue when the tiling is fixed
   by anti-diagonal grid reflection, black otherwise.

The second graph is antipodal, not local.  It cannot by itself mean “one step
farther from transitive.”  The first graph does not retain the antipodal
line/fold parity.  Treating them as one graph erases exactly the distinction the
user wanted to use.

The correct coordinate system is therefore two-dimensional before it is
anything else:

```text
radial coordinate = local one-tile-flip quotient depth;
fold coordinate   = rooted weighted blue/black complement-line structure.
```

At `n<=7`, the local quotient is connected and its distance from the transitive
node equals the minimum Hamming weight of any fixed-path tiling in the node
fibre.  This is the computed finite bridge to MFAS/flip depth.  It still needs
a general proof.

## The objective address

For each merged node the executable address is

```text
(d,
 w_local,
 d_line and w_line,
 colour_WL(line),
 colour_WL(all channels),
 recursive deletion-parent address,
 canonical converse-orbit code).
```

Here:

- `d` is local quotient distance from transitive;
- `w_local` is the least shortest path word with the convention
  spine `<` rib `<` sea;
- `(d_line,w_line)` is the least rooted blue `<` black complement-line word;
- both WL colours retain edge multiplicity and count a self-loop twice in the
  incidence degree;
- the parent coordinate uses the multiset of **all vertex deletions**, not the
  old explorer convention of deleting the top vertex of one representative;
- the final canonical tournament code is an explicit exact tie-breaker.

That final coordinate matters philosophically.  A structural profile should
not be declared canonical merely because no collision has yet appeared.  The
right approach is to report collision counts at every controlled-forgetting
stage, then retain a mathematically exact code at the end.  If `n=8` produces
WL twins, the address remains valid and the twins become a new object to
study.

The objective node order is lexicographic in this tuple.  It begins at the
transitive class and spreads through flip-depth shells; inside a shell it uses
the blue/black rooted structure, recursive parents, and finally the exact
code.  It is objective relative to the explorer's declared Hamiltonian-path
and tile-order gauge.  It does not pretend that erasing that gauge is free.

## Tiling-to-node and node-to-tiling are a relation, not a bijection

The exact forward map is

```text
tiling mask
  -> canonical unmerged class code
  -> converse orbit
  -> merged node address n-aXX.
```

The inverse of a node is its whole fibre.  There is no canonical unique tiling
unless one deliberately chooses another gauge.

Two indices are useful:

```text
global index = sort by
  (Hamming weight, node address, unmerged class code, explorer mask);

fibre index = within one node, sort by
  (Hamming weight, unmerged class code, explorer mask).
```

The first makes the transitive tiling index zero and displays the outward
spread.  The second identifies a tiling inside its exact node fibre.  The JSON
atlases store both directions for all `33866` tilings at `n=3..7`; the browser
explorer displays them through `n=6`, and the compact offline continuation
stores mask-indexed maps plus complete inverse fibres at `n=7`.

THM-781 subsequently replaced the extensional inverse list by an intrinsic
one.  A fixed-path tiling presentation of a class `[T]` is exactly a directed
Hamiltonian path of `T`, modulo the free action of `Aut(T)`.  Therefore

```text
class tilings = HP(T)/Aut(T),
node tilings  = union over the one or two converse classes,
|class tilings| = H(T)/|Aut(T)|.
```

This matters conceptually: the observer cut is not an awkward label attached
after the tournament is built.  It is a Hamiltonian-path orbit.  For the
heptagon node `n7-a267`, `H=175` and `|Aut|=7`, so the explorer's 25 masks are
forced as `175/7`.  The inverse fibre is now generated from the object itself,
not recovered by searching all `2^15` tilings.

The `n=7` comparison isolates why the blue/black recursion is useful.  Local
depth/path plus the all-deletion parent address gives 270 cells and leaves two
twin pairs.  The nodes in each pair even have identical raw all-black incidence
(`62` or `86`).  Their rooted blue/black neighbourhoods differ, and the next WL
round separates both pairs.  Line placement therefore contains exact
information that neither deletion ancestry nor raw degree sees.

There is one further refinement the node fibre does not replace.  The line
metagraph of klein-S163 remembers the simultaneous-isomorphism orbit of the
ordered pair of endpoint tournaments.  Several line classes can lie over the
same unordered class pair.  Thus the hierarchy is genuinely

```text
tiling -> line orbit -> endpoint pair -> merged-node incidence.
```

Aggregated blue/black multiplicity is enough for the finite node fingerprint,
but not enough to reconstruct a particular line orbit.

## Precisely what information must survive for LRC

There is no context-free answer such as “preserve the tournament.”  The right
answer is predicate-relative and recursive.

For a representation map `q:S->Q` and a theorem-facing predicate `P`, forgetting
is sound only if at least one of the following is true:

1. `P` is constant on every fibre of `q`;
2. the discarded coordinate is reconstructed from retained coordinates;
3. a dual certificate annihilates every possible effect of the coordinate;
4. transport across all future wall/lift/deletion operations is fibrewise
   compatible; or
5. the ambiguous fibre is emitted as a named residual obligation.

The fourth clause is the one static quotient audits usually miss.  A state can
be pointwise pure and still have different futures.  This suggests a precise
Nerode-style definition of the minimal proof state:

> Two states are equivalent only if every legal continuation—wall word,
> residue lift, scale sheet, peel, deletion, or threshold change—has the same
> terminal witness outcome.

The minimal theorem-safe quotient is the quotient by this continuation
equivalence.  The raw tournament class is much coarser.  The exact metagraph
address plus tiling fibre repairs the combinatorial cut, but it still omits
several proven-to-matter LRC coordinates.

### A theorem-safe sufficient packet currently needs

- the normalized speed/core shape and the scale/residue fibre;
- the chosen observer/cut or fixed Hamiltonian-path gauge;
- exact circular gap widths, or equivalently sector potential plus microphase
  order at the required resolution;
- the closed-threshold wall flag distinguishing equality from dangerous
  interior;
- blocker/endpoint ownership, including tie ownership;
- wall chronology or monodromy over the whole arithmetic clock period;
- the slope-sheet or lift index when scale is being quotiented;
- deletion/peel transport data;
- the complement/time-reversal involution and which side of it is occupied;
- the metagraph node, tiling fibre address, and line-orbit identity when those
  combinatorial carriers are used; and
- a named certificate or residual proof obligation.

Not every proof branch must carry every field.  A branch may prove that a field
is irrelevant.  What is forbidden is dropping it without a fibre-purity,
reconstruction, annihilation, or residual argument.

## The four-coordinate object is a suspension with a constructible atlas

For an affine family `V(c)=cP+R`, HYP-6815 defines

```text
Phi_{P,R}(u,t) = min_i ||p_i u+r_i t||,
X_{P,R} = {(u,t,c,lambda): u=ct and Phi(u,t)>=lambda}.
```

This is a four-coordinate incidence object inside
`(R/Z)^2 x N x [0,1/2]`.  It is not a smooth four-manifold, and it is not the
full 13-dimensional LRC configuration space.  It is the exact stratified
suspension of one affine family.  LRC(14) asks whether every integer-slope
section at `lambda=1/14` is nonempty.

The metagraph should not replace this object.  It should be its **constructible
atlas**:

```text
base chamber label: merged tournament node address;
stalk over a node: exact tilings / observer cuts / line orbits;
metric stalk: gaps, sector potential, microphase order;
threshold stalk: owner-labelled clearance filtration in lambda;
transport: wall-crossing and lift/scale monodromy;
section: the arithmetic slope u=ct.
```

Then the lonely-runner question is not “which tournament class occurs?”  It is

```text
does every arithmetic section of this labelled constructible stack meet the
closed safe substack at lambda=1/14?
```

This makes the role of the exact address clear.  An address tells us which
base chart we are in.  It does not replace the stalk.

## The new recursive sheet signal

THM-761 arrived while this atlas was being built.  For a scaled core plus
exceptions, its sheets

```text
t_k=(t_0+k)/c, k in Z_c
```

preserve the core margin exactly.  Each exception removes an arc-preimage from
the sheet cycle.  With at most six exceptions the bad-sheet budget leaves a
free sheet under the theorem's exact gcd inequality.  At seven exceptions the
plain union bound reaches density one, and its tight method-wall is a cyclic
tiling of `Z_c`.

This is the strongest recursive bridge found in this session:

```text
continuous affine suspension
  -> integer slope c
  -> finite sheet fibre Z_c
  -> equal-length cyclic cover/tiling
  -> a new quotient-and-fibre problem one level down.
```

The sheet tiling is not the staircase tournament tiling.  But both demand the
same preservation discipline: a quotient node plus a fibre address, owner
labels, and transport.  The `r=7` wall should therefore be attacked as a
fibrewise cyclic-tiling classification, not as another unlabelled runner
tournament statistic.

## Preservation ledger across the main quotients

| Carrier | Preserves | Destroys | Repair |
|---|---|---|---|
| labelled fixed-path tiling | every non-path arc and cut gauge | alternate Hamiltonian gauges, metric phases | retain orbit/fibre maps |
| tournament class | full combinatorial tournament up to relabelling | which tiling/cut, observer, metric, owners | exact tiling and line stalk |
| converse-merged node | complement-invariant class content | chirality/which converse side | side bit or full two-class orbit |
| blue/black incidence profile | antipodal line degrees and loops | recursive neighbourhood placement | rooted weighted WL |
| rooted weighted line-WL | all finite 1-WL data | possible higher-order twins/line orbits | canonical code, 2-WL, line metagraph |
| sector-carry state | 14-sector potential and microphase order | global clock history if static | wall word/monodromy |
| raw round tournament | semicircle order | six-sevenths of threshold walls, exact clearance | sector/gap/threshold stalk |
| prime-grid cover | base lift blocking and private owners | higher-lift carry unless recorded | persistent owner/carry tree |
| affine slope section | exact family maximin | other affine families | normalize and atlas family space |
| sheet bad-set quotient | shared scale clock and exception cover | endpoint owners if only union size kept | labelled bad arcs and gcd strata |

## Tournament Analysis

The vertices in the finite audit are information carriers, not runners:

```text
SC type,
radial depth,
blue/black incidence,
rooted blue/black WL,
local deletion-parent atlas,
combined coloured WL,
structural address,
exact address.
```

The pairwise observable is the number of unordered merged-node pairs separated
by the induced partition.  The switch is either total retention or retention
per `log2(number of cells)` description bit.  At `n=7`, separation rises

```text
16192 -> 27974 -> 36636 -> 36856
```

from SC type through depth and line incidence to rooted line-WL.  Each gauge is
a transitive eight-vertex tournament with score histogram `{0:1,...,7:1}`, no
directed 3-cycle, singleton SCCs, and one Hamiltonian path; changing gauge flips
18 edges.  Compactness changes the order of carriers but not the fact that the
exact recursive line structure reaches the information ceiling.

## What is genuinely new, and what remains open

New and exact:

- an objective finite address for every merged explorer node;
- exact tiling-to-node and node-to-tiling indices;
- rooted weighted blue/black 1-WL completeness at `n<=7`;
- a clean separation of local flip depth from antipodal complement lines;
- an all-deletion recursive parent coordinate;
- direct explorer integration; and
- the constructible-atlas interpretation of the affine slope suspension,
  sharpened by THM-761's finite sheet fibre.

Open:

- blue/black connectivity and 1-WL completeness at `n=8` and beyond;
- a proof that quotient depth equals minimum tiling weight/MFAS in general;
- reconstruction of line-orbit fibres from node/edge data;
- the exact Nerode quotient for LRC wall/lift continuations;
- attaching metric/owner/carry/monodromy stalks to the metagraph address;
- classifying THM-761's seven-exception cyclic tiling wall; and
- using any of this to close the remaining LRC(14) scale-normal residual.

The structural lesson is simple: the underlying object is not one graph.  It
is a graph of charts with fibres and transport.  The merged node gives a place;
the tiling gives an address inside the place; the metric and owner labels say
what the place means; the wall/lift monodromy says what it can become; and the
arithmetic slope selects the section the proof must follow.

## Artifacts and principal sources

- `04-computation/tournament_tiling_metagraph_address_codex_S4.py`
- `04-computation/tournament_tiling_metagraph_address_n7_codex_S4.py`
- `05-knowledge/results/tournament_tiling_metagraph_address_codex_S4.{out,json}`
- `05-knowledge/results/tournament_tiling_metagraph_address_n7_codex_S4.{out,json}`
- `03-artifacts/visualizations/tournament-tiling-explorer.html`
- HYP-6825, HYP-6815, HYP-2245, HYP-2989, HYP-3106, HYP-3513
- THM-760, THM-761, HYP-6785
- `three-observer-categories-tiling-is-relative-tournament-is-blind.md`
- `the-line-metagraph-is-rigid-klein-S163.md`
- `the-blue-subgraph-IS-the-half-tiling-metagraph-...md`
- `two-axes-of-the-tournament-metagraph-...md`
- `lrc-restricted-tournament-quotient-stack-s535.md`
- `lrc14-prime-grid-cover-sector-carry-and-threshold-tournament-...md`
