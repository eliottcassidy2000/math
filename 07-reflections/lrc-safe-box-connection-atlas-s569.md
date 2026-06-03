---
source: codex-2026-06-03-S569
status: long repo archaeology / connection atlas
tags:
  - lonely-runner
  - safe-box
  - clock-triage
  - anti-bohr
  - endpoint-core
  - pinch
  - resonance
  - primal-dual
  - assumption-challenge
---

# LRC Safe-Box Connection Atlas

This pass searched broadly around the prompt:

```text
Does the orbit hit the safe box even once?
Which clocks matter?
Which speed families can be ignored, and which remain dangerous?
Can the model be changed so the check becomes efficient?
```

The repo already points to a clean architecture:

```text
closure/rank dismissal
  -> resonance/folding classification
  -> cheap pinch-witness search
  -> endpoint-cover peel
  -> labelled protected endpoint core, if no witness appears
```

So the efficient checker should not ask for "the whole orbit" first.  It should
return one of two certificates:

```text
primal:  a safe time, usually found by pinch/antipode/pair-sum clocks
dual:    a minimal labelled endpoint core proving why the current clock failed
```

The open proof work is then to show that the dual object cannot persist after
the allowed descents, projections, or refinements.

## Search Scope

The search swept across safe boxes, torus geodesics, forbidden covers, endpoint
protection, Bohr/anti-Bohr language, pinch and pair-sum sieves, resonance and
Fourier clocks, Lee/coding analogies, zonotopes, view obstruction, sandpile and
sector occupancy models, symbolic dynamics, braid persistence, and LP/Farkas
certificates.

No web search was needed for this pass.  The strongest anchors were internal:

- `07-reflections/lrc-safe-box-clock-triage-and-wild-remodels-s568.md`
- `07-reflections/lrc-reset-horizon-speed-family-gauge-s567.md`
- `07-reflections/lrc-time-as-orbit-the-reset-commensurability-categorizer-s563.md`
- `05-knowledge/hypotheses/HYP-2080-lrc-reset-commensurability-resonance-categorizer.md`
- `07-reflections/diophantine-approximation-lonely-runner-s361.md`
- `07-reflections/lonely-runner-tight-stratum-s357.md`
- `07-reflections/kernel-residue-trick-atlas-2026-05-30.md`
- `07-reflections/lrc-lens-atlas-s430.md`
- `05-knowledge/hypotheses/HYP-1900-lrc-tournament-incidence-core.md`
- `07-reflections/lrc-multi-sieve-recursive-pinch-moduli-no-apex-s562.md`
- `07-reflections/lrc-recursive-multisieve-lower-bound-ledger-s563.md`
- `07-reflections/lrc-zero-branch-star-formalization-s548.md`
- `07-reflections/lrc-tree-entropy-the-order-parameter-s543.md`
- `07-reflections/lrc-symbolic-dynamics-target-shift-s541.md`
- `07-reflections/lrc-spacetime-braids-and-persistence-vineyards-with-a-wild-multitude-s540o.md`
- `07-reflections/lrc-sectors-as-nodes-dual-mapping-occupancy-and-the-dft-duality-s536o.md`
- `07-reflections/lrc-zeckendorf-bridge-s451.md`

## Executive Synthesis

The "total reset" clock is a category label, not the main witness finder.
For integer-normalized LRC instances every rank-one primitive orbit physically
closes after period `1`, so raw reset length collapses.  What remains useful is
the amount of folding inside that period: short resonances, repeated critical
moments, endpoint coincidences, and how much of the safe box boundary is hit
only as a wall.

The repo's strongest split is:

```text
dense or high-rank closure with an interior box point  -> dismiss as easy
closed rank-one or low-rank orbit near boundary        -> keep
pinch or antipode gives a time                         -> done
all open gaps covered                                  -> peel endpoint core
endpoint core has no labelled descendants              -> false obstruction
endpoint core survives all descents                    -> real proof target
```

This makes the question "does the orbit hit the safe box once?" feel more like
a finite primal-dual cover problem than a continuous-time search problem.

## Clocks That Matter

1. **Closure/rank clock.**  This decides whether the orbit is dense in a torus,
   closed rank-one, or trapped in a proper low-rank subtorus.  It can dismiss
   many cases if the closure has an interior safe-box point.  It does not by
   itself locate the witness in the hard integer case.

2. **Resonance/folding clock.**  This is the real replacement for raw reset
   length after integer normalization.  Count short relation vectors, repeated
   endpoint moments, distinct critical moments, and safe-measure collapse.
   HYP-2080 says the arithmetic progression is the maximally folded boundary
   model: many resonances, few distinct moments, zero open safe measure.

3. **Endpoint cover clock.**  The S357/S361 anti-Bohr frame says the decisive
   trichotomy is positive open gap, boundary witness, or full open cover with
   every endpoint protected.  This is the clock that can produce a dual
   certificate.

4. **Pinch/pair-sum clock.**  HYP-2075/S562 say pair sums `v_i+v_j` are natural
   finite coordinates.  Pinch candidates are cheap primal witnesses because
   they test places where two runners are forced against opposite danger walls.

5. **Sector/discrepancy clock.**  S536 dualizes runners into fixed circle
   sectors.  Loneliness becomes "the two observer-adjacent cells are empty."
   The DFT of the occupancy vector is the resonance/character-sum picture.
   This clock is valuable if it keeps observer-cell labels.

6. **Symbolic/wall clock.**  S541 codes the orbit as a compactified word with
   open target `G` and wall target `W`.  A counterexample would be an arithmetic
   periodic word avoiding both.  Coarse bad-symbol cycles are too lossy unless
   they retain owners, carries, and gate labels.

7. **Braid/vineyard clock.**  S540o reframes runners as spacetime braids and
   loneliness as a fat observer tube or persistence-vine threshold.  This may
   be best for perturbation stability: if a witness bar survives a bottleneck
   move, nearby speed sets can be dismissed in families.

8. **Half-turn clock.**  Useful as a diagnostic for wall alignment and tension,
   but not decisive unless translated into endpoint-labelled data.  It is easy
   to overvalue because it is visually simple.

## Cases To Ignore First

These cases should usually not consume proof-search time until they re-enter
through a labelled endpoint core:

- **Dense full-torus closures with interior safe-box intersection.**  Kronecker
  style closure already supplies a hit after approximation.
- **Any closure with certified interior box point.**  The search should project
  to the orbit closure first and stop if the safe box has positive relative
  interior there.
- **Low-tier pinch witnesses.**  If the pair-sum/antipode clock finds a witness,
  do not keep classifying the whole speed family.
- **Endpoint covers that peel to empty.**  THM-391's zero-branch star theorem is
  the warning: many scary p-adic covers have no strict endpoint core.
- **Bare half-turn anomalies.**  Without endpoint ownership they are signals,
  not certificates.
- **Bare entropy or p-adic occupancy signals.**  S543/S548 show they are good
  order parameters, but not proof objects unless the descendants are labelled.

## Cases To Keep

These are the families that deserve retained state in a checker or proof
ledger:

- **Rank-one closed orbits with high resonance folding.**  Raw reset is
  normalized away; folding is not.
- **Boundary-only tight families.**  Arithmetic progressions are the model:
  zero open safe measure but wall hits.
- **Full covers with nonpeeling endpoint cores.**  This is the anti-Bohr dual
  obstruction.  It is the most important object to export.
- **Low-rank subtori whose closure misses the safe-box interior.**  These need
  a lower-dimensional endpoint/core check, not generic equidistribution.
- **High-denominator residual packets after pair-sum pinch.**  S562/S563 show
  residuals can move to generated moduli, not only prime factors of the base
  denominator.
- **Endpoint debt with owners.**  Keep endpoint owner, event owner, compactified
  wall, cross-prime coordinate, Gabor label, and pair-tension label.

## Strong Connections Found

**HYP-2080 / S563: reset becomes resonance.**  The user's "total reset moment"
is a strong category idea, but in normalized integer LRC the reset period is too
coarse.  The useful invariant is resonance-folding inside the period.  Count
short relations and critical-moment collisions.

**S357/S361: LRC is finite anti-Bohr, not classical approximation.**  The usual
Diophantine-approximation story points in the wrong direction.  We want all
runners away from zero at the same time.  The finite endpoint problem is:
uncovered endpoint, boundary witness, or protected full cover.

**S430/HYP-1900: the central object is labelled incidence.**  Rows are
constraints/endpoints/cells/cuts.  Columns are speeds/protectors/arcs/gates.
The bad object is not "many close runners"; it is a leafless labelled
protection hypergraph.

**S526 permutohedron frame: chambers linearize loneliness.**  Inside each braid
or permutohedral chamber, the loneliness condition is a linear interval
condition.  This suggests storing chamber handoff and endpoint debt instead of
sampling times.

**HYP-2075/S562: pair sums are witness coordinates.**  Pair sums `v_i+v_j`
define natural pinch clocks.  The one-modulus apex obstruction dissolves under
multi-sieve recursion, so the checker should use generated residual packets.

**HYP-2076/S563: lower-bound proof as obligation descent.**  The proof should
reduce hard cases by a trichotomy: local tier witness, exported frontier mass,
or conserved descent to child obligations.  The existing proof-obligation
tournament is too transitive unless it remembers endpoint owners.

**S548: zero-branch p-adic covers are not enough.**  Covered zero branches peel
to empty endpoint core.  This blocks a tempting false path: p-adic coverage
alone is not the obstruction.

**S543: entropy is an order parameter.**  Tree entropy distinguishes generic
safe behavior from tight AP walls, but it is not the dual certificate.  Use it
to sort regimes, then hand off to endpoint-core checks.

**S541: symbolic dynamics needs a richer alphabet.**  The compactified word
model is promising, but coarse bad-symbol cycles are spurious.  Add owner,
gate, carry, gap-race, and pair-tension labels.

**S536: sector occupancy is the real-space dual of resonance.**  Empty observer
cells and Fourier character debt are the same information in different bases.
This is a strong candidate for a remodel because it replaces moving runners by
a fixed cell cycle plus boundary-crossing events.

**S540o: persistence can make family certificates.**  If safe-box hits are bars
in an observer-vine, bottleneck stability may prove that whole speed families
remain safe under perturbation.

**S451: endpoint debt may want a normal form.**  Zeckendorf/Ostrowski language
is not literally LRC, but it suggests a canonical no-adjacent-carry normal form
for endpoint debt transfer along recursive denominators.

**New coding-theory tangent: rank-one LRC as a Lee code.**  For integer speeds,
the orbit `t -> t v mod 1` is a one-generator cyclic code in the torus, and the
safe-box question asks whether some codeword has every coordinate at Lee
distance at least `1/n` from zero.  Short resonance vectors are dual codewords.
This suggests computing dual Lee distance / short-relation spectra against
safe-measure collapse.

## Wild But Testable Remodels

**1. Witness-or-core certificate pipeline.**

Implement a checker that never returns "unknown" without structure:

```text
input speed set
  normalize / compute closure rank
  if closure interior intersects safe box: return closure certificate
  run pair-sum and antipode pinch candidates
  if safe time found: return primal certificate
  build forbidden endpoint cover
  peel endpoints with private gaps
  return minimal labelled core
```

The returned core should include endpoint owner, protector interval, wall side,
denominator, residual packet, and any cross-prime or pair-tension labels.

**2. Resonance-dual-code profiler.**

For a speed set, compute:

```text
short integer relations r with r dot v = 0
relation support histogram
distinct endpoint/critical moments
safe measure
dual Lee-distance proxy
```

Prediction: tight or near-tight families have unusually many short relations,
few distinct critical moments, and low safe measure.

**3. Sector sandpile / exclusion process.**

Replace runners by chips on the `n` sector cycle.  Boundary crossings move one
chip to an adjacent sector.  LRC asks whether the two observer cells are empty.
Try recurrent-state or sandpile-style invariants, but keep observer labels and
crossing owners.  If the quotient forgets which chip moved, it likely destroys
the endpoint-core predicate.

**4. Vineyard stability prover.**

Build the observer tube persistence vineyard for a speed set and perturb the
speed vector.  A long safe bar should certify nearby families; a bar ending at
a compact wall should export endpoint debt.

**5. Endpoint debt normal form.**

Construct the recursive denominator/product-depth transfer graph for endpoint
debt.  Test whether surviving debts admit a path-like, Fibonacci-cube, or
Zeckendorf-style no-adjacent-carry representation.  If yes, proof search gets a
canonical simplifier.

**6. Dual Farkas/IP core search.**

Translate finite endpoint covers into an incidence matrix.  Either solve for a
safe time / uncovered row, or return a dual row-weight certificate.  Report
rank, Smith factors, private pivots, and whether the core is leafless after
label-preserving peels.

## Tournament Analysis / Assumption Challenge

Do not default to runners as tournament vertices.  The better vertex sets for
this problem are:

```text
pinch candidates
endpoint owners
protected endpoints
short resonance vectors
sector cells
boundary-crossing events
proof obligations
braid crossings
vineyard events
dual-code checks
```

Candidate tournament for the next script:

- **Vertices:** labelled endpoint obligations after pinch search.
- **Pairwise observable:** which obligation has the earlier/private removable
  endpoint under the current peel order, with labels for shared protector and
  wall side.
- **Switch/gauge:** orient `A -> B` when peeling `A` first strictly reduces the
  protected interval count or residual modulus depth for `B`; reverse when `B`
  does.
- **Tie Hamiltonian path:** denominator depth, then endpoint coordinate, then
  owner speed index.
- **Fingerprints:** score histogram, directed cycles, SCCs, edge flips after
  adding pair-sum labels, Hamiltonian-path count, and sink/private-endpoint
  concentration.

Predicate preserved: whether the endpoint cover has a nonempty labelled core
after allowed peels.

Information destroyed: exact continuous time order inside chambers, unless wall
events are retained as labels.

Challenged assumption: the proof object does not need to be a runner graph.  It
is probably an endpoint/protector incidence core with optional runner labels,
not a graph whose primary vertices are runners.

## Practical Next Session

The highest-value script is `find_witness_or_core`:

```text
04-computation/lrc_witness_or_core_s570.py
05-knowledge/results/lrc_witness_or_core_s570.out
```

Start with small bounded integer speed sets and known hard rows.  The first
success criterion is not proving LRC; it is returning interpretable certificates
that classify every tested set into:

```text
closure-easy
pinch-witness
boundary-wall
empty-core-after-peel
nonempty-labelled-core
```

If nonempty labelled cores appear, those are the next real objects to study.
If none appear outside known boundary families, the repo's many clocks have
collapsed into a much smaller proof problem.
