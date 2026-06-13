---
source: codex-2026-06-03-S568
status: speculative synthesis and triage rules
tags:
  - lonely-runner
  - safe-box
  - clocks
  - resonance
  - endpoint-cover
  - tournament-analysis
  - assumption-challenge
---

# Safe-Box Clock Triage and Wild Remodels

The question is:

```text
does gamma(t) = (v_1 t, ..., v_d t) mod 1 hit
B_n = {x : ||x_i|| >= 1/n for every i} even once?
```

That is cleaner than asking for "loneliness over time."  It is a hitting
problem for one orbit and one box.  The right response is to demote most clocks
to diagnostics and keep only the clocks that can change the yes/no answer.

## The Clock Hierarchy

### 1. Closure clock: only a regime splitter

The first object is the orbit closure

```text
K(v) = closure {t v mod 1 : t in R}.
```

If `K(v)` intersects the interior of `B_n`, Kronecker gives a hit.  If `K(v)`
misses `B_n`, the speed family is a genuine obstruction inside its subtorus.
If `K(v)` only touches the boundary of `B_n`, it is a compactified wall case.

So the closure clock matters only to split:

```text
full/high-rank closure with an interior box point: safe, ignore
boundary-only closure: compactification case, keep
rank-one closed orbit: hard LRC case, keep
low-rank subtorus disjoint from the box: possible structured obstruction, keep
```

For standard primitive integer LRC, this clock stops distinguishing after the
reduction: every orbit is rank-one and resets at `t=1`.  HYP-2080 is right:
the useful invariant is not reset length, but how the closed orbit folds.

### 2. Endpoint clock: decisive for certification

For rank-one/integer systems, the direct safe-box complement is an interval
cover on the time circle:

```text
D_i = {t : ||v_i t|| < 1/n}
B-hit exists iff union_i D_i does not cover the time circle,
with closed-wall repair handled by endpoint compactification.
```

This is the clock that cannot be ignored.  Its wall endpoints are

```text
t = (k +/- 1/n) / v_i.
```

Every exact certificate should ultimately reduce to one of:

```text
open gap:       a component of time not covered by danger intervals
wall witness:   a boundary point covered only weakly, compactified safe
bad core:       full open cover and every endpoint strictly protected
```

The first two cases are safe.  The only real counterexample-shaped object is
the third: a full open cover with a nonpeeling endpoint-protection core.

### 3. Pinch clock: the elegant witness extractor

The endpoint clock certifies, but it can be too large.  The witness clock should
be the lower-envelope maximum of

```text
F(t) = min_i ||v_i t||.
```

A maximum of `F` occurs at a runner antipode or at a balanced pinch where two
distance functions meet.  Those times live in

```text
t = (2k+1)/(2v_i)              single-runner antipode
t = m/(v_i + v_j)              opposite-side pair pinch
t = m/|v_i - v_j|              same-side pair tie
```

The pair-sum branch is the natural LRC witness branch already visible in
HYP-2075/S557.  So an efficient checker should not ask "scan the whole period."
It should ask:

```text
Does a low-complexity pinch candidate already give F(t) >= 1/n?
If not, which endpoint-cover core explains the failure?
```

This suggests a two-output oracle:

```text
find_witness_or_core(S):
  1. test pair-sum/antipode pinch candidates by increasing denominator tier
  2. if found, return witness t and margin
  3. else build only the protected endpoint core needed to explain failure
  4. return residual core packets, not the full endpoint universe
```

### 4. Half-turn phase clock: mostly ignorable for yes/no

The half-turn tournament clock is excellent for measuring spread, entropy,
walk shape, and quotient class.  But the safe-box predicate is anchored at
`1/n`, not `1/2`.

For the yes/no question, the half-turn clock can be ignored unless it is
decorated with endpoint data.  It becomes useful only as:

```text
phase diagnostic:       how folded/spread the orbit looks
corridor diagnostic:    how many phase cells a safe interval crosses
quotient compressor:    which marked class/fiber the endpoint state lives in
```

Bare half-turn tournaments are therefore not proof objects.  Endpoint-labelled
half-turn corridors might be.

### 5. Resonance clock: the hard-family sorter

For primitive integer sets, all runners reset at `t=1`, so the reset length is
constant.  The sorter is the resonance lattice:

```text
Lambda(v) = {m in Z^d : m dot v = 0}.
```

Short resonances make critical times coincide.  Coincidence compresses the
movie.  HYP-2080's punchline is the right one:

```text
generic low resonance: many moments, positive safe measure, easy
high resonance: few distinct moments, small or zero safe measure, hard
AP-like maximal folding: boundary-tight, proof-critical
```

So resonance matters after the closure clock, not as a separate wall clock but
as a compression meter for the endpoint and pinch clocks.

## Cases We Can Safely Ignore

1. Full-dimensional incommensurate orbits.

If the speeds are Q-independent enough that `K(v)=T^d`, then the orbit hits the
positive-volume safe box.  Ignore.

2. Any orbit closure with an explicit interior safe point.

This includes many low-rank noninteger families.  Once a point of the closure
has every coordinate distance `>1/n`, density inside `K(v)` finishes it.
Ignore after the closure feasibility certificate.

3. Integer/rank-one systems with a low-denominator pinch witness.

If pair-sum or antipode candidates hit `F(t) >= 1/n`, the endpoint cover does
not need to be built.  Ignore the rest of the period.

4. Full danger covers with empty endpoint-protection core.

These are wall-only or peelable cases.  They are safe after compactification.
Do not confuse "open cover" with "counterexample."

5. Half-turn-only anomalies.

A strange half-turn tournament, high `H`, low `H`, or many phase crossings is
not dangerous unless it also carries endpoint-cover debt.

## Cases We Need To Worry About

1. Rank-one closed orbits with high short-resonance count.

These are the true LRC hard cases.  They fold many endpoint/pinch moments onto
few actual times.

2. Boundary-only tight families.

AP and V*-type rows can have zero open safe measure but still have closed wall
witnesses.  These are proof-critical because any open-only method falsely sees
danger.

3. Full open covers whose endpoint core does not peel.

This is the counterexample-shaped object.  The core must be owner-compatible,
pressure-cyclic, and resistant to the known zero-branch/star peel mechanisms.

4. Low-rank irrational subtori with no obvious interior box point.

Full irrational is easy, but a dense orbit in a proper subtorus is only as good
as that subtorus's intersection with the box.  These should be handled by
subtorus-box feasibility before being dismissed.

5. High-denominator pair-sum residuals.

If low-tier pinch candidates fail, HYP-2075/HYP-2076 say the remaining object
is not "all time"; it is a residual packet ledger.  These packets need owner,
residue, and pressure labels.

## A More Elegant Decision Model

Replace "simulate the runners" by a three-stage certificate pipeline:

```text
Stage A: closure feasibility
  Input: real/projective speed vector.
  Compute integer relations m dot v = 0.
  Decide whether K(v) intersects int(B_n).
  If yes, done.

Stage B: pinch maximum
  Input: rank-one/integer residual.
  Search the antipode and pair-pinch clocks by denominator tier.
  If max F(t) >= 1/n, return t.

Stage C: endpoint-core dual certificate
  Input: no found pinch witness.
  Build danger interval cover only as much as needed.
  Peel exposed endpoints.
  If core empty, return wall/gap certificate.
  If core nonempty, return labelled owner-pressure core.
```

This is primal-dual:

```text
primal certificate:   a time t in the safe box
dual certificate:     a full danger cover with a nonpeeling protected core
```

The conjectural theorem is that the dual certificate cannot exist for
primitive LRC speed sets.  Computation should stop producing scalar "near
misses" and instead output which side of this primal-dual fork it lands on.

## Wild Remodels

### Tropical lower-envelope model

View `F(t)=min_i ||v_i t||` as a tropical object.  The safe-box question is

```text
is max_t F(t) >= 1/n?
```

The active set at an optimum is usually two runners, sometimes more in a
resonant case.  Hard families are exactly those where many tropical vertices
coalesce.  This turns "resonance folding" into degeneracy of a tropical
arrangement.

### Coding-theory model

The rank-one integer orbit is the cyclic code

```text
C_v = {t v mod 1}.
```

LRC asks whether the code has a word whose Lee distance from zero in every
coordinate is at least `1/n`.  The obstruction is a covering-radius statement
for a one-generator code in the Lee metric.  Resonance is low dual distance:
short `m dot v=0` relations are dual codewords causing bad covering.

### Matroid/flow model

Endpoint protection is an incidence matroid on interval owners and boundary
points.  A counterexample is a leafless protection matroid with an integer
realization.  The efficient checker should compute rank, pivots, private
endpoints, and owner cycles instead of sampling time.

### Packet automaton model

Use pair-sum pinch packets, not time cells, as states.  A transition is legal
only if endpoint owners and CRT residues agree.  LRC becomes: every legal
packet automaton either emits a witness packet or descends to a smaller packet.
This is the recursive multisieve picture in finite-state form.

### Tournament Analysis remodel

Do not put runners as vertices by default.  Better vertices here are:

```text
pinch candidates
endpoint owners
protected endpoint packets
short resonance vectors
subtorus relations
proof obligations
```

Pairwise observable:

```text
which object removes more uncovered safe margin or protects more endpoint debt?
```

Switch:

```text
orient toward the object that explains the obstruction at lower denominator.
```

Tie Hamiltonian path:

```text
increasing denominator, then owner speed, then index.
```

Fingerprints:

```text
score histogram, SCCs, directed cycles, Hamiltonian paths,
edge flips under dyadic/CRT lift, surviving endpoint-core rank.
```

This preserves the predicate better than runner tournaments because the target
is not "which runner beats which runner."  The target is "which obstruction
packet prevents the orbit from entering the box."

## Bottom Line

The clocks that matter for the safe-box question are:

```text
closure/rank clock:  dismiss dense/interior cases
pinch clock:         find witnesses efficiently
endpoint clock:      certify gaps, walls, or protected cores
resonance clock:     sort hard closed orbits by folding
```

The half-turn phase clock matters only after endpoint decoration.

The whole scenario should be remodeled as a primal-dual hitting problem:
either produce a safe-box time from the pinch/lower-envelope clock, or produce
a minimal protected endpoint core.  Everything else is instrumentation.
