# Applied Residue/Phase/Incidence Programs

Source: `opus-2026-05-30-S5`

This is a second pass over the residue/phase/incidence synthesis, but with the
direction reversed.  Instead of asking "what do the abstractions mean?", ask:

```text
If the abstraction is right, what should we build or prove next?
```

The useful object is not just a tournament, or `H`, or a homology group.  The
useful object is a complete pairwise world together with three failure
detectors:

```text
completeness: is the complete binary pairwise model still valid?
residue:      after a projection, what structured survivor remains?
phase:        which orthogonal/eigen/character channel carries the anomaly?
incidence:    what matrix/hypergraph/transport interface proves the rank fact?
```

The order matters.  Completeness is the admission test.  Residue and phase are
the two main explanatory axes.  Incidence is the proof layer when scalar
features and adjacency shadows have become too compressed.

## A Small Contrast Experiment

I added `04-computation/residue_phase_incidence_contrast_s5.py`, which reuses
the existing `TournamentFeatures` extractor and prints a compact table for:

- transitive `T7`;
- Paley `T7`;
- interval `T7`;
- the two `H=63` single-core tournaments;
- the `n=6`, `H=37` arc-flip local trap;
- THM-025, the `n=9` non-real-rooted Omega counterexample.

Stored output: `05-knowledge/results/residue_phase_incidence_contrast_s5.out`.

The table is intentionally not a classifier.  It is a feedback loop: pick
standard examples, inspect which coordinates move, then revise the next
hypotheses.

Key rows:

```text
Paley T7       H=189 alpha=(80,7,0)  residue=(20,1,0) rank=2 rho=.750
Interval T7    H=175 alpha=(59,14,0) residue=(16,2,0) rank=2 rho=.729
H=63 cores     H=63  alpha=(31,0,0)  residue=(0,0,0)  rank=0 rho=1.000
H=37 n6 trap   H=37  alpha=(14,2,0)  residue=(2,0,0)  rank=1 rho=.857
THM-025 n9     H=237 alpha=(94,10,1) residue=(2,1,0)  rank=2 rho=.979
```

This makes four different mechanisms visible.

The `H=63` examples are exact projection kills.  One vertex supports every odd
cycle; deleting it leaves no Omega residue.  The anomaly is not that the
residue is complicated.  The anomaly is that the complete-core signature can
realize `r=31` while apparently avoiding `r=3` and `r=10`.

THM-025 is not an exact kill.  Almost all odd cycles pass through one vertex,
but two cycles survive deletion, and those two are disjoint.  That tiny
survivor has rank 2 and `I(R,2)=9`.  This is exactly the kind of residue that
can defeat a too-simple real-rootedness story: small enough to be missed by
bulk statistics, structured enough to preserve an obstruction.

Paley and interval `T7` are not near-kills.  Both are regular, both have the
same three-cycle count, and both have broad deletion residues.  The difference
is fiber geometry: Paley has more odd cycles and fewer disjoint pairs, interval
has fewer cycles and more disjoint pairs.  This is the phase/fiber regime, not
the localized-residue regime.

The `H=37` trap is different again.  It is neither exact kill nor dangerous
near-kill.  It has a moderate residue and nonzero score variance, and its
interesting property is dynamical: greedy arc flips can stall on a plateau.
That makes it a phase-landscape example, closer to an energy barrier than to a
projection defect.

## The Abstract Picture

The repo has been repeatedly compressing tournaments in ways that are useful
but lossy:

```text
T -> H(T)
T -> score sequence
T -> Omega alpha-vector
T -> support shadow
T -> even graph
T -> quotient bucket
T -> endpoint child class
T -> path-homology ranks
```

Every compression has two outputs.  The first is the visible quotient.  The
second is the residue that the quotient did not preserve.  The mistake is to
record the first output and forget the second.

The newer lesson is that residue alone is still not enough.  Some phenomena do
not look like a small survivor after projection.  They look like phase:

- an eigenvalue sign pattern changing with `k mod 4`;
- a Walsh/Krawtchouk level vanishing;
- a circulant connection set moving mass between Fourier modes;
- a root channel becoming badly imbalanced;
- a local-search landscape changing without a one-vertex kill.

And some phenomena require a third category.  Endpoint transfer is the model:
support matching, adjacency, and scalar counts all see shadows, but rank lives
in an incidence matrix.  The proof object is a private child, a nonzero minor,
a leaf-peelable collision hypergraph, or an odd Smith profile.

So the practical decomposition is:

```text
complete binary pairwise data
  -> residue_vector(T)
  -> phase_vector(T)
  -> incidence_vector(operation on T)
```

## Problem Program 1: n=9 Real-Root Failure Classification

Current state: THM-025 refutes universal real-rootedness of `I(Omega(T),x)`.
The counterexample has:

```text
I(Omega,x) = 1 + 94x + 10x^2 + x^3
max-loss deletion residue = (alpha1, alpha2, alpha3) = (2,1,0)
residue rank = 2
near_kill_rank2_vertices = 1
```

This suggests a concrete classification program:

1. Enumerate or sample `n=9` tournaments.
2. Compute `omega_near_kill_vertices`, `omega_near_kill_rank2_vertices`,
   max-loss residue alpha-vector, and the root/Newton defect of `I(Omega,x)`.
3. Split failures into:
   - exact-kill or complete-core cases;
   - near-kill rank-2 cases;
   - broad-residue phase/root cases.
4. Check whether every first failure has a near-kill rank-2 certificate.

The sharper hypothesis is not merely "real-root failures have small residue."
The sharper hypothesis is:

```text
failure = small deletion residue with rank at least 2
          plus an imbalanced root/phase channel.
```

This is immediately testable using the existing residue features.  The missing
piece is a root-defect helper that records Newton inequalities and root
statistics next to the residue table.

## Problem Program 2: Single-Core Signature Gaps

The `H=63` examples say that complete Omega is not forbidden.  They say
something more precise:

```text
single-core complete Omega gives H = 1 + 2*r_core(signature)
```

The observed gaps `r_core != 3,10` are now more important than the old false
"63 is forbidden" story.  If `r=3` and `r=10` are permanently absent in the
single-core signature language, then the complete-core routes to `H=7` and
`H=21` are blocked for a concrete automaton reason.

This should be treated as a finite-state arithmetic problem.  The signature
formula weights each `1...0` inversion by a power of two depending on the gap.
That is close to a transducer:

```text
binary word -> weighted inversion sum -> reachable target set
```

Concrete next experiment:

1. Build a dynamic program over suffix state, not over full signatures.
2. Track reachable residues and exact small targets.
3. Prove that targets `3` and `10` are impossible by invariant, or find the
   first counterexample.
4. Then separate non-single-core complete Omega as a second family.

This is a place where the abstract residue story becomes a proof obligation:
the exact kill has reduced the tournament problem to the image of a weighted
binary language.

## Problem Program 3: Circulant Maximizer Phase Transition

Paley and interval `T7` are a warning against overusing deletion residue.
Their deletion residues are broad, not small, and both have rank 2.  The
interesting contrast is:

```text
Paley T7:    alpha=(80,7,0),  support_excess=44
Interval T7: alpha=(59,14,0), support_excess=23
```

Paley wins at `p=7` by producing many odd cycles while keeping disjointness
low.  Interval eventually wins in larger regimes by packing more ordered
structure.  That crossover should be modeled as a phase-channel problem:
Fourier characters of the circulant connection set, trace signs, additive
energy, and Fejer-kernel-like concentration.

Concrete next experiment:

1. For circulant tournaments `C_p(S)`, compute a `phase_profile`:
   eigenvalue magnitudes/arguments, low trace signs, additive energy of `S`,
   and cycle-count deviations.
2. Join it with the existing residue vector.
3. Test whether the Paley/interval crossover is predicted by phase features
   before residue features.
4. Push the same phase profile into path-homology symbol matrices for Paley
   and interval connection sets.

The working prediction: the circulant maximizer transition is not a localized
residue threshold.  It is a phase-dominance transition.

## Problem Program 4: Endpoint Transfer As Incidence, Not Adjacency

Endpoint transfer now has several small-n facts that look separate but are
probably one theorem:

- private child witnesses imply full `F_2` row rank;
- support matching is not enough;
- merged SC collision columns have support 3;
- those support-3 columns are incidence hyperedges, not necessarily triangles
  in the parent metagraph;
- the collision hypergraph leaf-peels completely through the tested range;
- even-graph transfer exposes 2-primary Smith factors.

The abstraction says: stop looking for an adjacency explanation.  The object is
an incidence system.  The theorem should be stated in that category:

```text
endpoint transfer rank = private pivots + peelable collision block
                         + absence of 2-primary torsion.
```

Concrete next experiment:

1. Compute the next endpoint level if feasible.
2. Track the collision hypergraph's 2-core and leaf-peeling order.
3. Compute Smith factors for the same matrix.
4. Try to prove peelability by an endpoint interval or SCC-defect ordering.

This also suggests engineering features for `tournament_tda.py`: when an
operation is part of the data, emit incidence summaries such as row rank,
private-column count, collision support spectrum, peel-core size, and small
Smith factors.

## Problem Program 5: OCF-Guided Active Ranking

The practical application that survived the application probe is active
pairwise ranking.  The residue/phase/incidence synthesis makes it more
concrete.

For a complete tournament, `H(T)` measures how many total orders remain
consistent with the pairwise world.  For an incomplete or soft pairwise system,
we should first measure completeness defect.  Then query comparisons that are
expected to reduce the remaining OCF ambiguity.

A prototype can be simple:

```text
state:
  observed pairwise graph, optional weights, current tournament completion(s)

features:
  H(T)
  Omega cycle packets
  residue vector by vertex deletion
  completeness defect if data is partial/weighted/tied

acquisition:
  for candidate pair (i,j), estimate expected drop in H
  break ties by expected drop in near-kill residue rank
  expose named Omega packets resolved by the comparison

output:
  ranked list plus ambiguity certificate
```

The point is not just to sort faster.  The point is to return a certificate:
"these three cyclic packets are the reason the ranking is unstable."  That is
useful for PRP/RAG reranking, LLM evals, A/B testing, sports schedules, network
meta-analysis, and any domain where pairwise comparisons are expensive or
noisy.

The first prototype should not wait for all theory to be perfect.  It should
use small `n` exact `H`, Monte Carlo or local approximations for larger `n`,
and a clear JSON report.

## Problem Program 6: Completeness Defect For Ties And Partial Data

The trienerment lesson is that some tournament obstructions vanish when ties or
missing arcs are allowed.  So before using `H`, OCF parity, or residue rank on
real data, ask whether the data is close enough to a tournament.

The parabolic law

```text
E[Delta c3 | score=s] = s(n-1-s)/2
```

is the most concrete checksum currently available.  It says what frustration
injection should look like in a complete binary pairwise world.  A tied,
missing, weighted, or thresholded comparison system should carry a
`completeness_defect` vector measuring deviation from this forced parabola.

Concrete next experiment:

1. Define the defect for partial/tied data.
2. Test it on trienerment examples realizing tournament-forbidden `H=7` and
   `H=21`.
3. Perturb weighted comparison matrices, threshold to tournaments, and compare
   variance of `H`, residue rank, and phase features against completeness
   defect.

If this works, it gives the applied ranking prototype a principled warning
label: when completeness defect is high, hard tournament invariants should be
reported as unstable.

## Problem Program 7: Paley/Circulant Path Homology

The GLMY/Tang-Yau-style path-homology thread needs both residue and phase.
Residue tells us which supports survive projection.  Phase tells us how
circulant characters make the boundary matrices split.

Concrete next experiment:

1. Implement symbol matrices for circulant digraph path homology.
2. Specialize to Paley connection sets and interval connection sets.
3. Compare Betti ranks, torsion/rank drops, and phase profiles.
4. Ask whether Paley flatness forces closed formulas for low-degree Betti
   numbers.

This is also a bridge back to coding theory: Paley support shadows, QR/NQR
orbits, and circulant symbol matrices are all different faces of the same
finite character decomposition.

## Revised Feedback Loop

Future sessions should use this loop:

```text
1. Pick a standard contrast set.
2. Compute residue, phase, and incidence features separately.
3. Ask which lens moves and which lens stays flat.
4. Generate a hypothesis only in the lens that moved.
5. Add one script or table that can refute that hypothesis cheaply.
```

The current standard contrast set should be:

```text
transitive;
Paley and interval circulants;
H=63 single-core complete-Omega classes;
THM-025 real-root failure;
n=6 H=37 local trap;
endpoint-transfer even-graph and merged-tournament matrices;
trienerment/tied examples realizing forbidden tournament values.
```

The main danger is overfitting one lens.  Residue explains many recent wins,
but Paley/interval and local-search traps are not primarily exact-kill
phenomena.  Phase explains broad spectral behavior, but it misses private
endpoint pivots.  Incidence explains rank proofs, but it is the wrong first
language for application data until completeness is checked.

The abstraction is useful only if it keeps producing smaller, more falsifiable
programs.  The table from this session does that: it turns a loose slogan into
four separable mechanisms and a concrete queue of experiments.
