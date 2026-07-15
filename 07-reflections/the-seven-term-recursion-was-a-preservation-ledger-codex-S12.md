# The seven-term recursion was a preservation ledger

**Session:** codex-2026-07-15-S12  
**Companion theorem:** THM-801  
**Exact audit:** `mobius_cech_metagraph_codec_codex_S12.py/.out/.json`

The old formulas looked like three ways to count the same sort of triangle:

```text
A+B+C-D-E-F+G,
A+B-C,
A+B-C+D-E-F+G.
```

They are not interchangeable formulas.  They are declarations of what a
recursive chart has kept, what it has double-counted, and which deepest slot
can still witness phase.  The important advance is to stop reading the minus
signs as arithmetic and start reading every letter as an owned overlap.

## How the view sharpened

THM-442 correctly identified the full word as the third finite difference of
the quadratic staircase cell count.  Its warning was equally important:
Hamiltonian-path count is not cell-affine.  THM-549/550 then exposed the
mirror-folded square/pronic parity charts.  The odd chart taught the first
preservation lesson: a corner and an overlap can have equal cardinality and
cancel as scalars while remaining different geometric slots.  THM-553 made
the slots local through the birth/crossing address `(beta,tau)`.

HYP-3234 supplied the right language—local signed charts with sidecars—but it
was still mostly a legality principle.  THM-796 finally separated the three
actual object sorts:

```text
tiling endpoint -> complement line -> converse-merged node.
```

It showed exactly where recursion changes character: a pullback on tilings, a
phase torsor on bare lines, and a non-lumpable weighted correspondence on
nodes.  THM-801 now closes the loop back to the seven old letters.  The three
full subtriangles are literal restriction maps; the formerly neglected middle
one contracts the interval-root gap.  Their triple overlap is the phase
witness that makes the full complement-line descent exact.

The progression was therefore:

```text
scalar identity
 -> local tile slots
 -> typed tiling/line/node incidence
 -> overlap phases
 -> exact Cech gluing
 -> transported interaction stalk.
```

## What the object actually has to remember

There is no single best metagraph.  There is a family of observations of the
same marked-path object:

| observation | exact information | first thing destroyed |
|---|---|---|
| literal tiling | every interval-root orientation | isomorphism invariance |
| complement line | unordered antipodal endpoints | endpoint phase |
| high/low face pair | endpoint-deletion restrictions | apex bit |
| full high/gap/low triple | every tile and overlap owner | nothing at tiling level |
| `Xi` | upper and two face node pairs, colours | gap face and literal lower lines |
| `Omega` | upper and all three face node pairs, colours | literal lower lines |
| `S2` crossing layers | mirror-pair phase counts | positions within a layer |
| merged node | converse class | marked Hamiltonian path orbit |
| primitive face row | aggregate branching ratios | individual child lift law |
| `K_u(z)` | boundary cyclicity over the whole node fibre | full lower-node destination |
| corewise Möbius stalk | exact leg/apex membership over a literal core | compactness |

This table explains why several apparently successful coordinates do not
compose.  A statistic can identify every current node and still fail to say
how one particular child continues.  Static separation and recursive
sufficiency are different predicates.

## Two Möbius operations, two different discoveries

The word “Möbius” was doing double duty in past threads.

First, Boolean subset inversion turns intersection moments of colour events
into exact colour atoms.  Upper, low-face, and high-face blue have all product
marginals of order at most two, yet their law differs from the product law by

```text
ab(1-b)(x-1)(y-1)(z-1).
```

The entire missing colour datum is one top Boolean coefficient.  In exact-set
coordinates its signs are `++-+--+`, the old odd Legendre chart realized with
constant amplitude.  Marginalize any role and the interaction disappears.
Condition on upper colour and it reappears as equal positive and negative
`B2` face charges that cancel globally.

Second, the legal marked-Hamiltonian-path faces form an interval poset, not the
full vertex Boolean lattice.  Its whole top Möbius function is the four-term
endpoint square

```text
Omega f=f-full_low-full_high+common_core.
```

Applied to `C3`, it counts cyclic triangles meeting both path endpoints.
Applied to the primal and dual Smith currents, its antisymmetric part is the
black endpoint defect `epsilon` and its symmetric part is the longitudinal leg
current after an apex correction.  The two-dimensional transitivity-flow
cloud is literally an electrical face-curvature cloud.

These are not two metaphors for one calculation.  One is the incidence
algebra of colour subsets; the other is the incidence algebra of legal path
faces.  Their meeting is useful precisely because they preserve different
things.

## The “middle face” changes the tournament question

The high and low faces delete path endpoints, so they are induced
subtournaments with inherited marked paths.  The gap face sends
`(a,b)->(a-1,b)` whenever the interval has length at least three.  It shortens
every retained interval root.  The result is again a valid tournament tiling,
but not a fixed induced subtournament of the original.

This challenges the default choice of tournament vertices.  A recursive
vertex can be an interval root, a gap layer, an overlap obligation, or a
marked-path minor operation.  The immediate mathematical task is to make the
gap face intrinsic: identify the category in which it is a natural morphism
and compute how path automorphisms transport it.  Until that is done, its
weighted node row is exact but presentation-dependent.

## Symmetry cannot itself explain the black drift

The line-level colour masses are left/right symmetric.  The new Smith identity
is sharper: reflection sends `epsilon` to `-epsilon` inside every fixed black
projected node-pair fibre, and fixes the boundary-curvature pair `(q0,q1)`.
Every signed histogram is symmetric; every fixed coefficient is even.

So the observed left/right imbalance of black flow cannot be a broken line-
level symmetry.  It is a disintegration effect.  The only candidates are
unsigned strata such as `|epsilon|`, `(q0,q1)`, stabilizer size, path-orbit
multiplicity, and unequal source mass.  This removes a large class of false
explanations and makes MPA-36 a sharply finite conditional-current problem.

## Continued fractions fit as transport, not as another rank

The continued-fraction archaeology reaches the same conclusion from the LRC
side.  A digit or convergent is only a value address.  THM-778 became exact
when it carried the centered phase, tie blocks, metric wall ranks, owner
labels, inverse steps, and the induced substitution on the token fibre.

The metagraph analogue is now explicit.  A node or face signature is a base
address.  Its continuation lives in the corewise Boolean function

```text
I_(u,c)(p_L,p_H,a).
```

The next useful continued-fraction experiment is not to correlate partial
quotients with node ranks.  It is to let centered Christoffel substitutions
act on the low/high leg variables and then on this function's Möbius
coefficients.  A finite invariant coefficient sector would be a genuine
transported address.  Failure would say exactly which interaction degree the
continued-fraction action creates.

## Recursive structural picture

The metagraph can now be read in four coupled directions:

```text
horizontal: C3 / score variance / Smith potential / leg current,
vertical:   endpoint and gap face correspondences,
transverse: epsilon and boundary C3 curvature,
stalk:      literal core + path-orbit leg/apex interaction.
```

Blue is the clean self-dual locus: zero antisymmetric Smith curvature and a
single rigid top colour interaction.  Black is not random residue.  It is the
union of signed curvature pairs whose signs cancel before node projection but
whose absolute strata distribute unequally over fibres.  Mixed nodes are the
interfaces where those distributions carry the largest electrical/current
effect, matching the incoming Smith “battery” observation.

This is closer to the four-dimensional object the LRC work keeps asking for,
but it is still not the LRC object.  Tournament coordinates are combinatorial
and observer-relative.  LRC needs metric scale, threshold, owner, carry, and
wall transport.  The right bridge is a bundle of transported stalks over this
four-direction tournament base, not a claim that one node coordinate already
contains loneliness.

The next pull arrived while this write-up was live.  THM-802 gives the
unmarked cell-network Tutte product
`prod_(k=1)^(n-2)(x-1+[k]_y)`, so MPA-39 no longer starts from zero: the open
step is to mark that product by `A/B/C` face ownership and ask whether a
deletion-contraction coefficient recovers `q` or the colour cubic.  At the
same time, the exact class-network Smith audit found that the 132 apparent
`n=7` potential/score discordances are real but confined to adjacent score
levels, and that the proposed primal/dual resistance reciprocity is false.
That is useful separation, not bad news.  The cell network, complement-line
network, and class-flip network are three electrical functors with different
edge sets; THM-801's local curvature identities survive without identifying
their global resistances.

## Pullable frontier

MPA-33 through MPA-39 now give bounded tasks for other agents:

1. generalize the line-gluing phase theorem to arbitrary covers;
2. decide `Omega+B2` injectivity at `n=8`;
3. minimize the corewise Möbius stalk by interaction bidegree;
4. disintegrate black current by unsigned Smith/Möbius curvature;
5. give the gap face an intrinsic path-minor semantics;
6. make Christoffel substitutions act on the stalk;
7. compute the marked Tutte/Smith specialization of the staircase cover.

The most informative failure would be as valuable as success: the first
`n=8` codec collision, saved with its literal core and lowest missing
coefficient, would tell us exactly what the current picture still cannot see.
