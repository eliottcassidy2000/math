---
source: codex-2026-06-01-S537
status: computational probe plus creative hypothesis zoo
tags:
  - lonely-runner
  - tournament-analysis
  - restricted-quotients
  - isomorphism-classes
  - source-reachability
  - crt-channels
---

# Restricted Tournament Functors for LRC

The user's request asks for mappings of LRC into tournament space that still
make the conjecture a statement about which isomorphism classes an arithmetic
walk can exhibit, but where the exhibited set is not the bloated raw
observer-marked quotient.  This session turns that into a functor problem.

For a fixed denominator `n` and primitive speed set `V={v_1,...,v_{n-1}}`, a
restricted tournament functor is:

```text
F_n(V,t) = a marked tournament or stratified marked tournament class
```

constant on the LRC wall cells, with a target family `Target(F,n)` such that:

```text
V satisfies LRC(n)  iff  the F_n(V, -) walk exhibits some class in Target(F,n).
```

The ranking criteria are not just "few classes."  A quotient is useful only if
it has:

```text
fidelity:     LRC target visible in the class;
purity:       target classes are not mixed safe/unsafe fibers;
compression:  class count grows much slower than A000568;
memory:       enough arithmetic/geometric constraint survives to prove hitting.
```

Tiny but impure is bad.  Exact but huge is bad.  The sweet spot is the smallest
pure quotient that still remembers why the arithmetic walk is constrained.

## Computation

Artifact:

```text
04-computation/lrc_restricted_tournament_mappings_s537.py
05-knowledge/results/lrc_restricted_tournament_mappings_s537.out
```

Tournament Analysis declaration:

```text
pairwise observables:
  half-turn phase, safety danger shell, safe-arc sentinel separation,
  CRT channel majority, safe-only subconfiguration

switch/gauge:
  observer-source threshold, lexicographic danger override,
  fixed boundary sentinels, residue-channel majority

tie Hamiltonian path:
  fixed vertex order, with observer/sentinels marked first

fingerprints:
  total pointed classes, good-only/bad-only/mixed fiber counts
```

Open wall-cell probe, no boundary compactification:

```text
n=4
  standard_source          classes=8    good_only=1  mixed=0
  binary_danger_lex        classes=4    good_only=1  mixed=0
  ternary_danger_lex       classes=4    good_only=1  mixed=0
  two_sentinel_safe_arc    classes=8    good_only=1  mixed=0
  safe_only_subtournament  classes=4    good_only=1  mixed=0
  crt_channel              classes=6    good_only=2  mixed=0

n=5
  standard_source          classes=22   good_only=2  mixed=0
  binary_danger_lex        classes=7    good_only=2  mixed=0
  ternary_danger_lex       classes=8    good_only=3  mixed=0
  two_sentinel_safe_arc    classes=72   good_only=23 mixed=0
  safe_only_subtournament  classes=7    good_only=2  mixed=0

n=6
  standard_source          classes=50   good_only=4  mixed=0
  binary_danger_lex        classes=11   good_only=4  mixed=0
  ternary_danger_lex       classes=11   good_only=3  mixed=0
  two_sentinel_safe_arc    classes=135  good_only=17 mixed=0
  safe_only_subtournament  classes=11   good_only=4  mixed=0
  crt_channel              classes=16   good_only=3  mixed=0

n=7
  standard_source          classes=96   good_only=6  mixed=0
  binary_danger_lex        classes=17   good_only=6  mixed=0
  ternary_danger_lex       classes=12   good_only=3  mixed=0
  two_sentinel_safe_arc    classes=283  good_only=40 mixed=0
  safe_only_subtournament  classes=17   good_only=6  mixed=0
```

The first conclusion is sharp:

```text
binary danger-lex and safe-only collapse the standard source quotient by
roughly 4x-6x in these probes, while preserving pure good target classes.
```

The second conclusion is cautionary:

```text
two-sentinel safe-arc is geometrically faithful but expands class count.
Its value is not compression.  Its value is making the safe arc a marked
interval object whose target can be attacked by circular-order geometry.
```

## External Hooks Checked

A short web pass did not change the computation, but it sharpened which
decorations are natural rather than arbitrary.

1. Renault's six-runner proof emphasizes congruence classes modulo `6`; that
   backs the idea that a useful restricted class must remember residue-channel
   data, not just a raw tournament class.
   Source: https://doi.org/10.1016/j.disc.2004.06.008

2. Henze--Malikiosis and the later Alcantara--Criado--Santos shifted-LRC work
   push LRC toward covering radii of lattice zonotopes and dyadic fundamental
   domains.  That makes the zonotope-facet and permutohedral-chamber functors
   more than metaphors: they are candidate discretizations of the same
   covering problem.
   Sources: https://doi.org/10.48550/arXiv.1609.01939 and
   https://arxiv.org/abs/2506.13379

3. The Perarnau--Serra survey frames LRC as a problem with many equivalent
   languages.  This is exactly the right setting for a "restricted functor"
   search: the goal is not a pretty encoding, but an encoding whose target
   classes are small, pure, and proof-bearing.
   Source: https://arxiv.org/abs/2409.20160

4. Tournament Hamiltonian-path machinery, including fixed-endpoint variants,
   justifies keeping a canonical tie Hamiltonian path in the functor definition:
   it is a structural tournament primitive rather than bookkeeping noise.
   Source: https://www2.eecs.berkeley.edu/Pubs/TechRpts/1987/6111.html

5. The zonotopal proof of tournament score-sequence theorems gives a bridge in
   the other direction: tournament score data already has a convex-geometric
   zonotope language, so LRC-zonotope data and tournament-quotient data should
   be forced to talk to one another.
   Source: https://arxiv.org/abs/2009.09322

Concurrent upstream S535o arrived while this session was closing and claimed
HYP-2018/HYP-2019 for the metric-retention spectrum, especially the near-graph
mapping whose realizable classes are circular-indifference graphs.  A later
upstream S535 commit claimed HYP-2020 for the broader restricted quotient
stack.  The BLEX stack here is therefore renumbered as HYP-2021 and should be
read as the oriented tournament refinement of those symmetric threshold/stack
pictures: near-graph says "who is dangerously close," quotient-stack says
"which pure source layer," while BLEX says "how close-danger orients around
the observer-source target."

## The Six Tested Functors

### 1. Standard Observer-Source Functor

Vertices are observer plus runners.  Observer arcs are threshold arcs:

```text
0 -> i  iff  ||v_i t|| >= 1/n.
```

Runner-runner arcs are half-turn phase arcs.  LRC is source exhibition.

This is exact and pure, but it is still close to the raw marked quotient.
It is the baseline, not the final compression.

### 2. Binary Danger-Lex Functor

Assign danger levels:

```text
blocker = 2, observer = 1, safe = 0.
```

Higher danger beats lower danger; within equal runner danger, use half-turn
phase.  Thus every blocker is forced above every safe runner.  The observer is
a source exactly when all runners are safe.

This is the best computational hit of the session.  It keeps the LRC target
pure and collapses many unsafe circular classes into a small ordered
composition:

```text
blocker round tournament  >  observer  >  safe round tournament.
```

The target is the top safe layer:

```text
empty blocker block + safe round tournament in the fixed safe arc.
```

Hypothesis: this quotient is the minimal fixed-vertex pure quotient that still
remembers both the observer-source target and the runner half-turn geometry.

### 3. Safe-Only Subtournament Functor

Discard blockers and keep the half-turn tournament induced by currently safe
runners.  The class key is:

```text
(safe_count, safe_runner_tournament_class).
```

LRC is exhibition of the top stratum `safe_count=n-1`.

This has the same sampled class counts as binary danger-lex.  It is harsher:
it forgets the blocker arrangement entirely.  Its virtue is that it turns LRC
into a pure "reach the top stratum" problem.  Its weakness is that the proof
needs an external invariant to explain why blockers must disappear.

### 4. Ternary Danger-Lex Functor

Assign:

```text
blocker > observer > near-safe > deep-safe.
```

Within equal runner shell, use half-turn phase.  This is a source-neighborhood
functor rather than a pure compression functor.  It sees whether the witness is
barely lonely or robustly lonely.

At `n=7` it had fewer classes than binary danger-lex in the sample, but only
three good-only target classes.  That means it merges some source geometry.  It
may be excellent for proving wide-gap cases and poor for wall-only cases.

### 5. Two-Sentinel Safe-Arc Functor

Add two fixed vertices:

```text
L = 1/n,    R = 1 - 1/n.
```

Build the half-turn tournament on `{L,R}` plus runners, with `L` and `R`
fixed in the isomorphism.  LRC is target membership:

```text
all runners lie in the marked arc from L to R.
```

This expands the class count in the probe.  It is still important because it
turns the target into a two-marked interval tournament class.  The safe target
is no longer an observer score; it is a separation statement by fixed sentinels.

This functor is the geometric adjoint of the Ferrers interval menu S525.

### 6. CRT Channel Functor

For composite `n`, group speeds by the largest proper divisor channel.  The
observer sees a channel as safe only if all runners in that channel are safe;
channel-channel arcs are majority half-turn arcs.

This is the arithmetic quotient.  It is small at `n=6`, and for `n=14` should
be a seven-channel object.  But S534o warns that single-level unit/parity tests
become vacuous at `n=18`.  Therefore the CRT functor must be decorated by
coupled debt tensors; raw channel majority is a shadow, not a proof invariant.

## Untested Functors That Now Look Worth Building

### 7. Sentinel Danger-Lex Hybrid

Use `L,R` sentinels and binary danger levels.  The class remembers both:

```text
which side of the safe-arc bracket a runner violates,
and the round tournament structure inside each side.
```

Target is all runners between sentinels.  This should be less bloated than the
pure two-sentinel half-turn functor and more geometric than binary danger-lex.

### 8. Left-Right Defect Functor

Split blockers into left-defect and right-defect:

```text
left blocker:  0 < v_i t < 1/n,
right blocker: 1-1/n < v_i t < 1.
```

Functor class:

```text
left-defect block > observer > safe block > right-defect block
```

or the circular analogue with left/right sentinels.  This should interact
directly with THM-384/387's two adjacent observer gaps.

### 9. Endpoint-Owner Protection Tournament

Static functor on a speed set, not a time-cell functor.  Vertices are speeds.
Orient:

```text
i -> j  iff endpoints owned by i are protected by j more than
          endpoints owned by j are protected by i.
```

A counterexample requires every endpoint protected.  If the owner-protection
tournament is acyclic, a source/sink should expose an unprotected endpoint or a
private protection leaf.  A non-peelable counterexample must carry an SCC in
this static tournament.

This is the dual of the time-cell source functors.

### 10. Endpoint-Row Tournament

Vertices are endpoint rows rather than speed owners.  Orient endpoint `e -> f`
if every protector of `f` also protects or shadows `e`, or if `e` has greater
private slack after deleting common protectors.

Target is not a source class but a peelable row class.  This maps LRC to:

```text
every endpoint-protection class exhibited by a primitive speed set is
either not full-cover or has a peelable tournament row.
```

### 11. Omission-Witness Tournament

For each runner `i`, delete `i` and use the known `n-1`-runner witness region
for the remaining set.  Orient omissions:

```text
i -> j iff omitting i leaves a larger or more robust witness region than
          omitting j.
```

LRC at `n` becomes a class-exhibition statement about whether one omission
source survives intersection with the omitted runner's safe set.  This is the
tournament version of mixed-threshold induction.

### 12. Permutohedral Chamber Functor

Vertices are local chambers or adjacent swaps in the circular order of runners
and sentinels.  Orient a chamber edge by whether it moves the observer toward
source status.  Classes are not tournaments on runners but tournaments on
chamber moves.

This could turn source avoidance into a forbidden word problem in the cyclic
permutohedron.

### 13. Apex-Bracket Functor

Use the source-sink/apex arc from S530.  Condition on the apex arc state, then
map the residual interior to a smaller danger-lex or safe-only class.  S535's
apex-reduced mapping was extremely small; the next version should retain the
apex bracket and the residual blocker count instead of discarding too much.

### 14. Coupled CRT Debt Tournament

Vertices are residue channels; edge labels are not just majority directions
but resonance-debt tensors:

```text
within-channel debt,
cross-channel pair debt,
triple and higher coupled debt.
```

The tournament edge is chosen by the first nonzero debt order.  This is the
version that could survive the S534o warning: parity alone dies, coupled debt
does not.

### 15. Fourier Mode Tournament

Vertices are Fourier modes or Riesz-product factors.  Orient mode `a -> b` if
mode `a` suppresses more forbidden mass or contributes more negative dual
slack than `b`.  The target is a dual certificate class.

This is not a runner tournament.  It is a proof-certificate tournament, useful
if the class-exhibition problem moves to LP/Riesz dual space.

### 16. Zonotope Facet Tournament

Vertices are LR zonotope generators or facets.  Orient by larger covering debt
after projection away from the other generator.  Target class is a facet-peel
or deep lattice point certificate.

This may be the cleanest bridge from endpoint debt to convex geometry.

## Abstract Hypotheses

### H1. Binary Danger-Lex Is The Minimal Pure Fixed-Vertex Quotient

Among fixed `n`-vertex tournament functors whose observer score alone detects
LRC and whose runner-runner relation is a deterministic coarsening of half-turn
phase plus safety, binary danger-lex is minimal up to reversal:

```text
blocker > observer > safe.
```

Any further collapse either mixes safe and unsafe classes or forgets all
runner geometry.

### H2. Binary Danger-Lex Class Counts Are Ferrers-Convolution Counts

The full binary danger-lex menu should be a convolution:

```text
C_BLEX(n) = sum_{b+s=n-1} R_b * R_s * glue(b,s),
```

where `R_k` is the round tournament class count on `k` vertices and
`glue(b,s)` records how the blocker and safe arcs interleave around the
observer wall.  The target count is exactly the fixed-safe-arc Ferrers menu:

```text
Target_BLEX(n) = open_source_menu(n).
```

### H3. Safe-Only Is The Quotient Of Binary Danger-Lex By The Blocker Ideal

Safe-only and binary danger-lex had the same sampled class counts because
blockers collapse to a single upper ideal whose internal round class is forced
or irrelevant in the sampled open cells.  If this equality fails at larger
boxes, the failure measures the first nontrivial blocker-geometry obstruction.

### H4. Ternary Shell Is A Robust-Loneliness Quotient

Ternary danger-lex should prove wide-gap/lacunary/random cases more naturally
than tight cases.  It separates:

```text
source by a wall witness,
source by a shallow open interval,
source by a deep robust interval.
```

The missing good target classes at `n=7` are not a bug; they are wall geometry
being intentionally crushed.

### H5. Two-Sentinel Is The Geometric Adjoint Of The Ferrers Menu

Two-sentinel classes expand because they mark the proof boundary.  The
trade-off is worth it if the target admits a clean theorem:

```text
all runners in marked arc L->R  iff  class lies in a two-sentinel Ferrers ideal.
```

This would make S525's Ferrers interval menu a theorem about marked circular
tournaments rather than a computed target list.

### H6. Sentinel Danger-Lex Is The Best Next Computational Functor

Hybridizing sentinels with danger shells should dominate both parents:

```text
fewer classes than two-sentinel,
more boundary memory than binary danger-lex.
```

Prediction: it has pure good fibers and class counts between BLEX and
two-sentinel, with left/right defect classes matching the THM-384 two-gap
states.

### H7. CRT Channels Need Coupled Debt Labels

Raw CRT channel tournaments are meaningful only when decorated by resonance
debt.  S534o says single-character parity becomes vacuous at `n=18`; therefore
the useful channel class is:

```text
(majority tournament, unit-balance tensor, higher-order debt profile).
```

This is the arithmetic partner of binary danger-lex.  BLEX supplies purity;
CRT-debt supplies memory.

### H8. Endpoint-Protection SCCs Are Dual To Source-Avoiding Walks

A time-cell source-avoiding walk should project to a static endpoint-owner SCC.
If the endpoint-owner tournament is acyclic, source avoidance is impossible.
Thus a counterexample must exhibit two classes at once:

```text
no source in BLEX time quotient,
nontrivial SCC in endpoint-owner protection quotient.
```

### H9. Omission-Witness Tournaments Are The Right Induction Interface

Given progress through total `n-1`, the `n`-runner problem should be read as a
tournament on omissions.  A source omission is the runner whose deletion leaves
the largest safe region.  Mixed-threshold induction is the claim that some
source omission survives its own runner.

### H10. Meaningful Restriction Is A Galois Connection

There is a contravariant relation:

```text
coarser class quotient  <->  stronger decoration needed for purity.
```

Distance-rank is maximally coarse and needs almost all arithmetic externally.
Standard observer-source is maximally faithful and needs little decoration.
The proof should sit at the fixed point of this connection:

```text
binary/sentinel danger quotient + coupled CRT/endpoint debt decoration.
```

## The Conclusion

The best new mapping is not the distance-rank `n`-state walk and not the
two-sentinel geometry by itself.  The best mapping is a stack:

```text
Layer 1: binary danger-lex tournament
  blocker > observer > safe
  small pure class set; LRC = exhibit observer-source class.

Layer 2: sentinel or left/right defect marks
  remembers which side of the safe arc failed.

Layer 3: CRT/coupled-debt channel labels
  remembers why the arithmetic walk is constrained.

Layer 4: endpoint-owner protection tournament
  dual static certificate; counterexample requires an SCC.
```

So the next proof-shaped statement is:

```text
For every primitive speed set, the decorated binary danger-lex walk exhibits
a good-only source class before the endpoint-owner protection tournament can
close into a realizable coupled-debt SCC.
```

That is complicated, but it is the right kind of complicated: a restricted
isomorphism-class target, a small pure quotient, and the arithmetic memory
needed to make "must exhibit" plausible rather than wishful.
