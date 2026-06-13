# Kernel-Residue Trick Atlas

**Instance:** codex-2026-05-30-kernel-residue-atlas
**Date:** 2026-05-30
**Status:** creative synthesis; not canon

## Premise

The Fejer session produced a useful correction:

```text
Fejer/IPR concentration
  = additive energy
  = variance of pinned triangle localization.
```

So the first bridge from Fourier mass to geometry is not mystical.  It is a
kernel identity in disguise.  The hard part begins after the kernel has become
local geometry:

```text
kernel / phase shadow
  -> pinned local statistic
  -> exclusion-corrected packets
  -> incidence / packing / rank
  -> scalar invariant or obstruction.
```

This note asks how that move can be reused across the repo and in famous
problems.

## The Meta-Trick

The reusable trick is:

```text
Do not try to prove a scalar invariant directly.
First find the kernel whose second moment is already the desired local
obstruction.
```

Then apply a four-step lift:

```text
1. Shadow: find the compressed quantity people already know how to bound.
2. Kernel: rewrite it as a convolution, autocorrelation, transfer kernel, or
   finite endpoint cover.
3. Pinning: evaluate the kernel at a marked vertex, endpoint, root, residue,
   color, or deletion.
4. Exclusion: subtract the diagrams that violate the actual combinatorial
   condition.
```

The Fejer example is:

```text
Shadow:    IPR / additive energy.
Kernel:    difference autocorrelation of S.
Pinning:   J_3(0,v), the through-0-and-v triangle profile.
Exclusion: none for triangles; starts at J_5.
```

The new HYP-1804 says the first three lines are exactly one axis for circulant
tournaments.  This tells us where the real proof labor begins: not at
spectral-to-spatial, but at spatial-to-exclusion-packing.

## What The Main Objects Represent

### Tournament

A tournament is a complete binary measurement field.  It is a sign assignment
on all positive roots of type A, but it is also an experimental design in which
every pair has been forced to choose a side.

The danger: many applications do not have this completeness.  Ties, missing
comparisons, weighted probabilities, and noisy preferences first need a
completeness-defect ledger before hard tournament invariants can be trusted.

### Score

Score is the Cartan shadow: the first-order divergence of the sign field.  It
sees source/sink pressure and coarse ranking, but it throws away curl.  Every
time a score-based heuristic fails, look for a residue in the two-row
representation or in the odd-cycle packet algebra.

### 3-Cycle

A 3-cycle is the first curl.  Under Fejer, pinned triangle variance is the
spatial incarnation of additive energy.  Under representation theory, it is
the first root relation.  Under applications, it is the smallest ambiguity
packet: A beats B, B beats C, C beats A.

The trick: use 3-cycles as a cheap diagnostic for where a kernel has already
localized, but do not confuse them with the final obstruction.  OCF lives in
all odd-cycle packets.

### Odd-Cycle Conflict Graph

`Omega(T)` is an exclusion graph of relation packets.  Its independent sets are
compatible packets of nontransitivity.  `H(T)=I(Omega,2)` says Hamiltonian path
count is a high-fugacity hard-core gas on ambiguity packets.

The lesson from Paley/Interval is brutal and useful:

```text
more cycles can lose to more packable cycles.
```

So famous-problem analogues should look for a transition from abundance to
packability, not just from small count to large count.

### Deletion Residue

Deleting a vertex, endpoint, color, speed, or residue class is a projection.
What survives is the residue.  Exact kills are easy to classify and often
misleading.  Dangerous examples are near-kills with just enough rank,
independence, or phase imbalance left to survive cancellation.

### Phase

Phase is what remains when the obstruction is not localized.  Paley/Interval is
mostly phase: the difference sits in Fourier modes, trace signs, and fiber
multiplicities, not in a tiny deletion residue.  Whenever deletion features are
broad, switch to phase channels.

### Incidence

Incidence is the proof layer.  If scalar features and support shadows know the
right answer but cannot prove rank, build the matrix or hypergraph:

```text
private pivots, endpoint collisions, Smith factors, protection graphs,
colored transversals, repeated-vertex correction diagrams.
```

The repo repeatedly learns that support matching is not rank.

## Trick 1: Kernel-To-Pinned-Variance

General form:

```text
sum_x |K(x)|^2  <->  variance of pinned local count.
```

Repo instance:

```text
additive energy E(S)
  <-> Var_{v != 0} J_3(0,v)
```

Potential reuses:

- For path homology, try to express boundary-rank variance as the second moment
  of a pinned path kernel.
- For real-root failures, express Newton-defect terms as second moments of
  pinned Omega residues.
- For active ranking, use the variance of pinned ambiguity packets to choose
  the next comparison.

Famous-problem analogue:

- Union-closed sets: coordinate frequency is the first shadow; the kernel might
  be pairwise union co-occurrence.  A second moment of "sets containing i and
  j after union" may reveal the residue when every coordinate frequency is
  below `1/2`.
- Graph reconstruction: deck statistics are shadows; pinned subgraph
  co-occurrence across deleted vertices may be the kernel whose rank determines
  reconstruction.

## Trick 2: Boundary Protection Graph

Lonely Runner has a clean endpoint language.  The forbidden set is a finite
union of open intervals.  A counterexample is an open cover in which every
endpoint is strictly protected by another interval.

That is exactly an incidence object:

```text
vertices: forbidden endpoints
edge/hyperedge: interval j strictly protects endpoint e
counterexample: every endpoint has protection indegree > 0
```

This is the same shape as endpoint-transfer collisions in the repo:

```text
private endpoint     -> witness survives
collision endpoint   -> needs incidence rank / peeling
all endpoints covered -> possible obstruction
```

Current external context: Sungkawichai-Trakulthongchai's 2026 paper gives a
computer-assisted proof of Lonely Runner for reduced `k in {10,11,12}` speed
sets, after Rosenfeld/Trakulthongchai checked lower cases.  Jensen's 2026
mixed-threshold work uses safe/unsafe Fourier indicators and arithmetic
progressions.  Both fit the same finite-cover plus Fourier-kernel split.

Repo trick to try:

```text
Build the endpoint-protection graph for every tight scanned set.
Find the local motifs that occur in boundary-only examples.
Then prove a full open cover would need an impossible motif.
```

## Trick 3: Return Residue

Caccetta-Haggkvist is a return problem:

```text
large outdegree gives expansion budget;
no short directed cycle means return residue is killed;
eventually expansion must collide with the root.
```

For Cayley digraphs, this becomes zero-sum additive growth.  The repo's cyclic
probe already found the interval model:

```text
A={1,...,d}, n=d(g-1)+1
|j(A union {0})| = jd+1 until first return.
```

The Aharoni-Berger-Chudnovsky-Guo-Zerbib nonuniform/rainbow theorem gives the
parallel harmonic checksum:

```text
directed girth <= 2 * sum_v 1/(outdeg(v)+1)
```

That suggests a weighted return-residue feature for general digraphs:

```text
R(v) = (
  layer sizes,
  second-neighborhood surplus,
  backward collisions,
  return walks,
  harmonic outdegree weight
).
```

Repo trick to try:

```text
Endpoint-transfer private pivots are to rank
as first-return pivots are to directed girth.
```

Can leaf-peeling of collision hypergraphs become a proof model for rainbow CH?

## Trick 4: Exclusion-Diagram Expansion

Fejer gives pinned walk/localization counts cheaply.  Simple cycles require
subtracting repeated-vertex diagrams.

This is a general technique:

```text
easy kernel count
  - illegal coincidence diagrams
  = actual object
```

Examples:

- `J_5`: pinned closed-walk count minus the first repeated-vertex diagrams.
- Path homology: all paths minus non-allowed face diagrams.
- Graph reconstruction: all deck co-occurrences minus coincident-deletion
  diagrams.
- LRC: total forbidden endpoint cover minus protected-endpoint overlaps.
- CH: all rooted walks minus back-return and layer-collision diagrams.

This is probably the most transferable proof technology from the Fejer note.
It gives a way to use spectral/analytic kernels without pretending they already
count simple combinatorial objects.

## Trick 5: Phase-Residue Split Before Optimization

Before attacking any maximization or famous conjecture, ask:

```text
Is the obstruction localized?
```

If yes, use residue:

- deletion residue,
- endpoint gap,
- exact kill,
- near-kill rank,
- private pivot.

If no, use phase:

- Fourier character,
- trace sign,
- Walsh degree,
- Krawtchouk layer,
- eigenvalue argument,
- root phase / Newton margin.

This prevents mixing proof methods.  H=63 and THM-025 are residue examples.
Paley/Interval is phase/fiber.  H=37 is landscape.  Lonely Runner has both
finite endpoint residue and Fourier safe/unsafe phase.  CH has return residue
and harmonic phase weights.

## Trick 6: Packability Beats Abundance

OCF teaches:

```text
scalar count is not enough;
the compatibility graph of packets matters.
```

This should be imported into famous problems whenever a classical approach
counts many local objects but cannot assemble them.

Possible translations:

- Erdős matching / disjoint cycles: local cycles are packets; the real object is
  their conflict graph.
- Union-closed: sets are packets; closure operations create compatibility
  constraints; high coordinate frequency may be a packing statement.
- Graph reconstruction: cards are packets; Kocay-style equations are incidence
  compatibility, not just counts.
- Stanley-Stembridge/Hikita-style positivity: monomial/elementary terms are
  packets; positivity may be hidden in a compatible-packet basis.

The repo-native move is to build an `Omega` analogue before optimizing.

## Trick 7: Local Finite-Difference Theorems

Global extremal statements are hard.  The interval Fejer data suggests a more
local route:

```text
every single signed-pair flip away from interval lowers a pinned kernel
in majorization order.
```

If true, this gives a proof skeleton:

1. Prove local monotonicity of the kernel.
2. Prove exclusion corrections preserve enough monotonicity.
3. Integrate local moves along a path in the orientation cube.

Analogues:

- LRC: one speed replacement changes endpoint protection by a controlled local
  move.
- CH Cayley: replacing one generator changes sumset-growth slack.
- Single-core signatures: appending or flipping one bit changes reachable
  `r_core` targets through a finite-state transition.
- Active ranking: one comparison flip changes Omega packets locally; use that
  derivative as acquisition.

## High-Value Concrete Experiments

### A. `J_5` Correction Atlas

Compute pinned spectral walk profiles and exact `J_5` profiles for circulants,
then enumerate the correction diagrams.  If the correction variance is also a
low-degree expression in autocorrelation data, the Fejer proof gap shrinks
again.

### B. `phase_profile(S)` For Circulants

Create a reusable profile:

```text
top_q, ipr, entropy, additive_energy,
trace_signs, eigen_argument_moments,
J3var, J5var, alpha=(alpha1,alpha2,...)
```

Use it on Paley/Interval/circulant maximizers and on path-homology symbol
matrices.  This would join HYP-1801 and HYP-1804.

### C. Endpoint Protection Hypergraph

Extend the Lonely Runner scan with endpoint-protection indegrees, peelability,
and unprotected endpoint motifs.  Compare tight examples with endpoint-transfer
collision hypergraphs.

### D. Return-Residue Extractor

For arbitrary digraphs, compute rooted layer profiles, backward collisions, and
first-return data.  Benchmark CH-tight Cayley examples against random and
tournament-derived digraphs.

### E. Weighted-Language Automaton For Single-Core Gaps

Build a finite-state or dynamic-program target search for the `r_core`
signature map.  Treat `r=3,10` as language non-membership, not as tournament
folklore.

### F. Packet-Compatibility Feature Block

For any problem object, build:

```text
packets, conflict graph, independence polynomial, deletion residues.
```

This generalizes `Omega(T)` and might make the repo's methods portable to
union-closed families, reconstruction decks, rainbow cycles, and active ranking.

## A Famous-Problem Betting Table

| Problem | Repo translation | Best trick |
|---|---|---|
| Lonely Runner | finite forbidden endpoint cover | boundary protection graph |
| Caccetta-Haggkvist | delayed return residue | layer growth + return pivots |
| Union-closed | coordinate projection residue under union | kernel-to-pinned frequency variance |
| Graph reconstruction | deletion deck projection | incidence rank / Kocay matrix as residue |
| Erdős matching / disjoint cycles | packet conflict graph | OCF-style packability over abundance |
| Stanley-Stembridge positivity | compatible packet basis | phase-to-incidence lift |
| Circulant H maximization | Fourier phase/fiber transition | Fejer-to-exclusion expansion |
| Active ranking / PRP | incomplete tournament completion | expected H-drop + ambiguity packets |

The common advice is the same: do not stare at the scalar.  Build the packet
space, find the kernel, pin it, and only then count.

## External Anchors

- Sungkawichai-Trakulthongchai, "Eleven, twelve, and thirteen lonely runners",
  https://arxiv.org/abs/2604.23906.
- Jensen, "Lonely runner mixed thresholds" paper,
  https://arxiv.org/abs/2605.27941.
- Aharoni-Berger-Chudnovsky-Guo-Zerbib, "Nonuniform Degrees and Rainbow
  Versions of the Caccetta-Haggkvist Conjecture", SIAM J. Discrete Math.,
  https://epubs.siam.org/doi/10.1137/22M1529658.
