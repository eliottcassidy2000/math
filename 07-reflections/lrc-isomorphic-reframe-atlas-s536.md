---
source: codex-2026-06-01-S536
status: exploratory synthesis after repo and web search
tags:
  - lonely-runner
  - isomorphic-problems
  - tournament-analysis
  - endpoint-incidence
  - zonotopes
  - distance-graphs
  - finite-sieves
---

# LRC Isomorphic Reframe Atlas

This session searched the repo and current web literature for problems that may
be genuinely isomorphic, boundary-isomorphic, or at least profitably
tournamentizable versions of the Lonely Runner Conjecture.

The working answer is:

```text
LRC is a one-parameter subgroup asked to hit a threshold chamber.
The exact obstruction is a finite protected boundary of forbidden arcs.
Tournament Analysis becomes useful whenever a reframe exposes pairwise
protection, pairwise pressure, or pairwise survivor comparison.
```

This refines S430/S450/S512.  S430/S450 identified the protected anti-Bohr
boundary.  S512 gave the exact observer-source tournament.  The new S535
mainline work adds a compression ladder: distance-rank, layered, modular, and
apex-reduced tournament mappings.  So the right question is no longer "is LRC a
tournament?" but:

```text
Which isomorphic LRC language gives the smallest tournament state space while
still preserving the arithmetic constraints that make the walk nontrivial?
```

## Web Snapshot, 2026-06-01

Primary sources checked or reused:

- Perarnau-Serra survey, "The Lonely Runner Conjecture turns 60":
  https://arxiv.org/abs/2409.20160
- Barajas-Serra, seven runners and regular chromatic number:
  https://arxiv.org/abs/0710.4495
- Bienia-Goddyn-Gvozdjak-Sebo-Tarsi, flows and view obstructions:
  https://www.sciencedirect.com/science/article/pii/S0095895697917706
- Henze-Malikiosis, covering radii, view obstructions, billiards:
  https://arxiv.org/abs/1609.01939
- Beck-Hosten-Schymura, Lonely Runner polyhedra:
  https://arxiv.org/abs/1606.01783
- Giri-Kravitz, Lonely Runner spectra:
  https://arxiv.org/abs/2304.01462
- Fan-Sun, amended spectrum conjecture:
  https://arxiv.org/abs/2306.10417
- Malikiosis-Santos-Schymura, finite checking:
  https://arxiv.org/abs/2411.06903
- Rosenfeld, eight runners:
  https://arxiv.org/abs/2509.14111
- Trakulthongchai, nine and ten runners:
  https://arxiv.org/abs/2511.22427
- Rosenfeld, nine runners:
  https://arxiv.org/abs/2512.01912
- Sungkawichai-Trakulthongchai, eleven through thirteen total runners
  (`k in {10,11,12}` in stationary-runner notation):
  https://arxiv.org/abs/2604.23906
- Jensen, mixed thresholds:
  https://arxiv.org/abs/2605.27941
- Blanco-Criado-Santos, shifted LRC and zonotope counterexamples:
  https://arxiv.org/abs/2603.24784
- Goncalves-Ramos, linear programming:
  https://arxiv.org/abs/2010.02271
- Bedert, Riesz products:
  https://arxiv.org/abs/2511.16636
- Chow-Rimanic, function fields:
  https://arxiv.org/abs/1711.01207
- Czerwinski, random runners:
  https://arxiv.org/abs/1102.4464
- Perarnau-Serra, correlations and dynamic interval graphs:
  https://arxiv.org/abs/1407.3381

I treat the 2025-2026 runner-frontier papers as current preprint claims, not
settled canon.  They matter strategically because they make the repo's `n=14`
total-runner focus adjacent to the public frontier.

## Grading

```text
A = exact restatement for the same speed instance.
B = exact after adding a standard restriction, compactification, or variant.
C = boundary-preserving analogue; could transport proof ideas.
D = loose analogy; useful only if it produces a concrete pairwise observable.
```

## Atlas Of Candidate Isomorphisms

| # | Reframe | Grade | Tournament connection | Possible leverage |
|---|---|---:|---|---|
| 1 | Stationary-runner Diophantine norm: find `t` with all `||v_i t|| >= 1/n`. | A | Exact observer-marked tournament: `observer -> i` iff runner `i` is safe. LRC iff observer is a source. Tie path: speed order. | This is the base source-reachability problem from S511/S512. |
| 2 | One-parameter subgroup hits the safe cube in `(R/Z)^k`. | A | Vertices are coordinates; orient `i -> j` by which coordinate is deeper inside the safe cube or reaches the wall first. | Converts "line misses box" into a chamber walk with edge flips at coordinate walls. |
| 3 | Forbidden open-arc cover of the multiplier circle. | A | Vertices are forbidden arcs or endpoints; orient arc `i -> j` if `i` protects more of `j`'s boundary than conversely. | A counterexample should be a strongly connected all-protected core. If the protection tournament is acyclic, peel a leaf endpoint. |
| 4 | Endpoint-protection hypergraph. | A | Rows=endpoints, columns=speed arcs. Pairwise observable is private protection debt between rows or between protectors. | The cleanest exact finite object. Target: prove every all-covered core has a source/sink/private row. |
| 5 | Regular circular coloring of integer distance graphs `G(Z,D)`. | A | Edges of the distance graph become runner constraints; orient distances by which one forbids the other endpoint multiplier. | Imports graph-coloring language but keeps the multiplier restriction, which ordinary coloring would forget. |
| 6 | Anti-Bohr set avoidance. | A | Each speed defines a Bohr neighborhood of zero; orient two speeds by whose anti-Bohr complement survives more intersections. | Makes the boundary object a protected anti-Bohr frontier. Good for Fourier and additive-combinatorics tools. |
| 7 | View obstruction: a ray avoids periodic obstacles. | A/B | Obstacles become vertices; orient obstacle `i -> j` if it shadows `j` along more directions or first-contact times. | Replaces runner language by visibility. Tournament cycles are mutual shadowing cycles. |
| 8 | Billiard motion in a torus with forbidden walls. | B | Wall families are vertices; orient by first-hit/last-hit order along the billiard orbit. | May expose monotone return maps and interval-exchange forbidden words. |
| 9 | LR zonotope covering radius. | A/B | Facets or generators are vertices; orient `i -> j` by larger covering debt/facet slack after projecting to the opposite generator. | Counterexample becomes a covering-radius debt cycle. Uses convex geometry while preserving endpoint debt. |
| 10 | Lonely Runner polyhedron contains an integer point. | A/B | Inequality pairs `(i,j)` are vertices; orient by which inequality is tighter under a lattice candidate. | Turns integer-point search into a chamber tournament over active constraints. |
| 11 | Coefficient of asymmetry / deep lattice point in zonotopes. | B/C | Facets around an interior lattice point are compared by asymmetry ratio. | Good if endpoint debt can be read as an asymmetric facet load. |
| 12 | Finite checking with bounded speed boxes. | B | Candidate residue classes are vertices; orient `A -> B` if `A` has smaller survivor burden after a prime/projection step. | Creates a finite "sieve tournament" ranking survivors; SCCs are the only hard components. |
| 13 | Proper/improper tuple sieve modulo primes. | B | Vertices are tuple classes modulo `p`; edge compares lift-survivor count or earliest proper witness. | For `n=14`, use mod `2` and `7` channels as a tournament on CRT survivor classes. |
| 14 | Maximum loneliness spectrum `ML(v)` and `S(k)`. | B/C | Vertices are deletion coordinates or spectral accumulation strata; orient by which deletion preserves larger maximum loneliness. | The recursion `S(k)` accumulating on `S(k-1)` smells like tournament vertex deletion. |
| 15 | Mixed thresholds `d_i` instead of uniform `1/n`. | B | Vertex `i` has its own threshold. Orient `i -> j` by excess slack `||v_i t||-d_i` versus `||v_j t||-d_j`. | Natural induction route: known `13` total runners plus one mixed `1/14` constraint. |
| 16 | Shifted LRC with nonzero starting phases. | B/C | Vertices are shifted arcs; orient by whose shift exports boundary debt to the other. | Recent counterexamples warn which tournament decorations are essential. |
| 17 | Lonely Vector Property in dimension two. | C | Vertices are rational vectors; orient by which vector sum blocks a target direction. | If it fails, identify the missing LRC structure by the first directed cycle in the vector tournament. |
| 18 | Function-field lonely runners. | C | Speeds are polynomials; orient by valuation/degrees of boundary collisions. | Gives an algebraic toy model where p-adic/depth tournaments may be exact. |
| 19 | Invisible runners in finite fields. | C | Vertices are finite-field speeds; orient by character sum or visibility deficit. | Finite-field tournaments may reveal which characters create endpoint protection. |
| 20 | Random runners are very lonely. | C | Orient `i -> j` by larger random phase gap contribution. | Helps distinguish generic transitive/no-core behavior from structured near-counterexample SCCs. |
| 21 | Lacunary-sequence / Local Lovasz Lemma regime. | C | Orient pairs by dependency strength; tie path by speed order. | If dependency tournament has bounded out-neighborhoods, LLL-style proof may be expressible as acyclic peel. |
| 22 | Riesz products and improved universal gap bounds. | B/C | Fourier modes are vertices; orient mode `a -> b` if its Riesz weight suppresses more bad mass. | Pairwise Fourier interactions can be studied as an Ising/tournament quadratic form. |
| 23 | Linear programming bounds. | B/C | Test functions or constraints are vertices; edge compares dual slack. | Produces certificate tournaments: acyclic dual slack means a simple bound, SCC means real obstruction. |
| 24 | Correlation/dynamic interval graphs. | C | Vertices are runners; orient by positive/negative correlation of bad intervals. | Negative-correlation tournaments could certify overlap tension in forbidden arcs. |
| 25 | "Almost lonely" and weak variants. | C | Orient blockers by who is the unique allowed close neighbor. | May prove that any source-avoiding walk must carry a blocker cycle. |
| 26 | CRT residue-cover problem. | A/B for fixed denominator quotients | Vertices are residue channels; orient by which channel exports more endpoint debt. | For `n=14`, the `2 x 7` split suggests a small channel tournament, not a 13-runner raw state. |
| 27 | p-adic or Bruhat-Tits frontier. | B/C | Vertices are p-adic leaves; orient by depth of uncovered/debt-export boundary. | A counterexample must be recurrent in the depth tournament; transient leaves peel. |
| 28 | Ostrowski/continued-fraction normal form. | C | Digits/carries are vertices; orient carry `i -> j` if it forces/blocks the next carry. | Useful if endpoint debt has a unique normal form with no adjacent repair carries. |
| 29 | Sturmian/Beatty return-word scheduling. | C | Return words are vertices; orient by which word delays a safe visit longer. | Converts irrational rotation intuition into a finite word tournament at rational approximants. |
| 30 | Circular-arc graph / nerve cover. | A for cover structure | Arcs are vertices; orient by left endpoint order unless one arc strictly contains/protects the other. | Circular-arc cover theory gives a finite endpoint certificate; directed cycles represent mutual protection. |
| 31 | Hall/set-cover duality. | A for endpoint IP | Rows and columns can both be tournamentized by pairwise Hall deficit. | If every all-protected core violates a Hall-like inequality, LRC follows in endpoint language. |
| 32 | Nowhere-zero flows and matroids. | C | Flow edges or circuits are vertices; orient by modular flow pressure. | May import circuit/cocircuit duality; use tournaments to detect cyclic pressure cores. |
| 33 | Oriented matroid/topes of the wall arrangement. | B/C | Topes/chambers are vertices; orient by wall-crossing direction relative to a base Hamiltonian path. | Gives a natural chamber-walk tournament; source target is a tope sign pattern. |
| 34 | Toric arrangement chamber walk. | A/B | Chambers are vertices; edge direction is the time orientation of crossing a wall. | This is the exact finite dynamical system behind each speed instance. |
| 35 | Cyclic permutohedron forbidden-word problem. | B/C | Adjacent transpositions/chambers are vertices; orient by whether swap moves observer toward source. | Source-avoidance becomes a forbidden word in chamber adjacency. |
| 36 | Symbolic dynamics/subshift of safe-bad words. | C | Words or states are vertices; orient by extension dominance. | A proof could show no primitive speed subshift avoids the source symbol forever. |
| 37 | Tropical/max-plus minimax: maximize `min_i ||v_i t||`. | C | Vertices are active constraints; orient by which active constraint wins after a local perturbation. | Tropical active-set tournaments may isolate the small set of true bottlenecks. |
| 38 | Thermodynamic energy `E(t)=-sum log ||v_i t||`. | C | Runners are vertices; orient by marginal energy relief when deleting the other. | Energy minima should have no cyclic blocker pressure unless a true obstruction exists. |
| 39 | Rapidity/formal-group load. | C | Orient by rapidity load `atanh(1-2||v_i t||)` or relief under formal addition. | Links LRC to the repo's formal group, but must keep endpoint decoration to avoid becoming scalar-only. |
| 40 | Phase synchronization/desynchronization of oscillators. | C | Oscillators are vertices; orient by phase lead or by who opens/closes pairwise separation. | Natural source of opening/closing tournaments and cyclic spread fingerprints. |
| 41 | Coding/covering-code language in the torus. | C | Code coordinates are vertices; orient by larger Hamming/torus margin from the bad codeword. | Could import covering radius and perfect-code intuition, especially for Paley-like arithmetic sets. |
| 42 | Scheduling/quorum anti-rendezvous. | D | Agents/tasks are vertices; orient by whose safe window dominates. | Mostly applied metaphor unless it yields a finite circular-arc protection matrix. |

## The Compression Ladder From S535

The most useful new mainline signal is that different tournament maps collapse
different amounts of structure:

```text
distance-rank mapping:
  pointed class count = n
  runner sub-tournament is transitive
  LRC = observer reaches bottom/top of a total order
  problem: arithmetic constraints are mostly lost

layered mapping:
  O(n) coarse shell states
  LRC = blocking shell empty
  keeps threshold geometry, loses fine runner order

modular/CRT mapping:
  tournament on residue channels
  preserves arithmetic gates
  promising for n=14 because 14 = 2*7

apex-reduced mapping:
  condition on an apex/source-sink alignment, then study n-2 interior runners
  preserves interior geometry with one less runner

standard half-turn / observer-source mapping:
  largest state space, but exact and richly constrained
```

This suggests a proof architecture with two projections:

```text
raw LRC walk -> compressed observer-rank/layer walk
             -> arithmetic channel tournament
```

The distance-rank projection supplies the tiny target.  The channel tournament
supplies the arithmetic constraints lost by total ordering.  A proof would show
that the compressed walk cannot avoid the source layer once the residue-channel
walk obeys the CRT/unit constraints.

The S534o rebase signal adds an important caveat.  Unit-channel tests are not
automatically strong at large composite denominators: for `n=18`, `n*=9`, the
single-level parity/unit obstruction is vacuous for primitive sets because the
unit group and the `3`-adic sublevel make zero representable almost always.
Thus the channel tournament should be treated as a compression layer, not as a
complete proof invariant.  It must retain coupled higher-order debt, or it will
repeat the same "parity fires at n=4, fades by n=8, vanishes at n=18" failure.

## Why Tournaments Are The Natural Sink

Every serious reframe produces one of five pairwise observables:

```text
1. protection:      which constraint covers whose endpoint?
2. pressure:        which runner is an irreplaceable blocker?
3. slack:           which coordinate/facet has larger safety margin?
4. survivor burden: which residue/chamber survives more sieve lifts?
5. opening motion:  which pair is creating separation rather than consuming it?
```

With a switch and a tie Hamiltonian path, each becomes a tournament.  The
fingerprints then have direct LRC meaning:

```text
source/sink      = peelable row, lonely observer, or decisive certificate
3-cycle          = mutual protection/pressure, first true obstruction atom
SCC              = non-peelable all-protected core
score histogram  = distribution of debt or survivor burden
edge flips       = wall-crossing events in the LRC clock
Hamiltonian paths = number of compatible global protection orders
```

This is the strongest conceptual bridge to the repo's tournament program:

```text
LRC counterexample should look like a strongly connected tournament core in
at least one faithful protection/pressure/channel language.
If every faithful tournamentization is peelable, LRC follows.
```

## Best Leverage Leads

1. Endpoint-protection tournament.
   Build the protection tournament on endpoint rows for hard `n=14` candidates.
   If it is acyclic, turn a source/sink into an unprotected boundary witness.
   If it has SCCs, classify the smallest labelled SCCs and ask which are
   realizable by primitive speeds.

2. Distance-rank plus channel constraints.
   S535 found the distance-rank map has only `n` pointed transitive classes.
   Alone it is too weak; paired with CRT channels it could become the exact
   "small state plus arithmetic constraint" proof system.

3. CRT/unit-channel tournaments.
   The pulled S533c result says inside debt is governed by unit-group balance
   modulo `n*`.  For `n=14`, the relevant split is parity plus mod `7`.
   Treat residue classes as vertices, and signs/units as edge labels.  The
   S534o n=18 attempt warns that this must be a coupled debt tournament, not
   only a single-character parity test: at `n*=9`, the unit test alone becomes
   vacuous on primitive sets.

4. Zonotope facet-debt tournament.
   Translate endpoint debt into covering-radius facet debt.  A counterexample
   should require a cyclic facet-debt core.  If LR zonotopes are always
   facet-peelable in the primitive/cosimple class, this becomes convex geometry.

5. Mixed-threshold induction.
   Jensen's mixed-threshold paper is almost tailor-made for using current
   total-13 progress against `n=14`: twelve speeds at threshold `1/13` plus one
   runner at threshold `1/14`.  Tournamentize omitted-runner choices by which
   omission leaves the largest safe interval.

6. Finite-sieve survivor tournament.
   The Rosenfeld/Trakulthongchai finite-checking method already has survivor
   states and projections.  Orient survivor classes by lift burden, properness
   time, or deletion effect.  Hard SCCs are then small explicit objects.

7. Riesz/LP Fourier interaction tournament.
   Fourier approaches give pair interactions between modes and runners.  The
   repo already has Ising/tournament machinery for quadratic interactions.
   Use it to find which Fourier pairs create the obstruction to improving the
   bound to `1/n`.

8. Function-field laboratory.
   Use valuations and degrees as exact depth labels.  If the function-field
   analogue has a clean source/peel theorem, pull that proof shape back to
   integer p-adic frontiers.

9. Counterexample-first generation.
   Instead of searching speeds first, enumerate small labelled all-protected
   tournament/SCC cores, then test speed realizability.  This reverses the
   search and could rapidly reveal forbidden protection shapes.

10. Chamber-walk forbidden words.
    In the toric arrangement or cyclic permutohedron, source avoidance is a
    word-avoidance problem.  Edge-flip fingerprints and Hamiltonian path counts
    may identify source-avoidance words that no primitive speed walk can spell.

## A Working Unification

The exact LRC formulations all preserve:

```text
one clock
many forbidden neighborhoods
finite signed boundary
strict protection relation
source target
```

The promising speculative isomorphisms are exactly those that preserve the
last three.  Any reframe that forgets the signed boundary becomes a scalar
loneliness meter; useful for intuition, but not enough to prove the conjecture.

The tournament slogan I would keep:

```text
LRC is not "a tournament" in one privileged way.
It is a source-reachability problem whose faithful shadows are tournaments.
The proof likely lives in showing that every source-avoiding shadow develops
a forbidden SCC, or else has a peelable source/sink that exposes a witness.
```
