# Breakthrough Candidate Map

**Session:** codex-2026-05-30-breakthrough-map
**Mode:** speculative research triage, not canon
**Prompt:** look outward for famous or lesser-known open problems that this
workspace might unexpectedly touch.

This is a map of possible breakthrough targets, from conservative publishable
extensions to genuinely sideways bridges.  The point is not that the repo is
already close to solving all of these.  The point is that several mature local
objects now have enough structure to act as probes of external problems:

```text
Hamiltonian path count        H(T)
Odd-cycle conflict graph      Omega(T)
OCF                           H(T)=I(Omega(T),2)
Residue-rank/deletion profile projection -> residue
Tiling/Krawtchouk cube        Boolean Fourier/compression
GLMY path homology            beta constraints
Paley/interval circulants     arithmetic vs additive structure
```

## Round 1: Nearest Targets

### 1. Which Odd Values Occur As H(T)?

External problem: characterize the odd integers that occur as the number of
Hamiltonian paths in a tournament.  This appears explicitly in community
questions around Redei parity and the missing values 7 and 21.

Local lever:

- OCF turns the value problem into an inverse independence-polynomial problem:
  `N = 1 + 2 alpha_1 + 4 alpha_2 + ...`.
- THM-343/THM-079 prove `H != 7,21`; THM-344 refutes the old hope that `63`
  is forbidden.
- The single-core complete-Omega signature formula explains exactly how `63`
  unlocks and why the same mechanism cannot immediately produce `7` or `21`.
- The arithmetic-candidate count `f(N)` already separates scalar arithmetic
  load from structural kill load.

Possible breakthrough:

Prove a first general "forbidden value mechanism" theorem: every permanently
forbidden value is explained by a small finite set of impossible alpha/residue
shapes.  The first realistic target is narrower: prove the single-core image
gap `r_core(s) notin {3,10}` for all signatures, closing the complete-core
parts of `H=7` and `H=21` in one stroke.

Next experiment:

Search the weighted-signature image modulo powers of 2 and 3, then try to
prove absence of 3 and 10 by a finite automaton or recurrence.  In parallel,
enumerate non-core complete-Omega candidates through n=9 and classify their
core/intersection patterns.

### 2. Maximal Hamiltonian Paths: Beyond Szele-Alon

External problem: determine or sharpen the maximum `P(n)` of Hamiltonian paths
in tournaments.  Szele's lower bound was the first probabilistic-method result;
Alon proved the maximum is within a polynomial factor of the random mean.

Local lever:

- The repo extends OEIS A038375 and has exact/new lower-bound data for `n=12,13`.
- OCF decomposes `H` into total odd cycles versus disjoint odd-cycle packings.
- Paley wins small `p == 3 mod 4`, but interval/cyclic tournaments overtake
  among circulants at larger primes.
- THM-138 gives a mechanism: Paley wins `alpha_1`; interval wins `alpha_2+`.
- Krawtchouk coordinates and Walsh degree bounds give a compressed search
  basis for extremizers.

Possible breakthrough:

Show that Alon's cyclic/interval construction is asymptotically optimal among
circulant tournaments, with an explicit constant or explicit crossover
criterion.  A more ambitious version would give a nontrivial improved upper
bound on `P(n)` for all tournaments by bounding high-order OCF packing terms.

Next experiment:

For prime circulants, express the first few `alpha_k` in additive-energy or
Fejer-kernel terms and prove the interval advantage for `alpha_2+` eventually
dominates Paley's `alpha_1` advantage.  Use p=19/23 computations as the
calibration points.

### 3. Classify Real-Rootedness Failure, Not Real-Rootedness

External problem family: real-rootedness, ultra-log-concavity, Lorentzian
polynomials, and independence-polynomial zero locations.

Local lever:

- Universal real-rootedness is false: THM-025 gives an n=9 counterexample.
- The counterexample is not generic.  It is a near-kill residue object:
  deleting one vertex leaves only two old cycles, but those two cycles have
  the exact disjoint residue needed to preserve the bad independent triple.
- ULC at `k=1` has an unconditional Turan-style proof.  ULC at `k=2` looks
  plausible for tournament co-conflict graphs even when real-rootedness fails.

Possible breakthrough:

Classify the first failure family at n=9 and prove an unconditional `k=2`
ultra-log-concavity theorem for tournament conflict graphs.  This would be a
clean result even though universal real-rootedness is dead.

Next experiment:

Generate all n=9 non-real-root examples, not just targeted samples.  Record
`max deletion loss`, residue alpha-vector, co-conflict graph shape, and
three-disjoint-triangle partition imbalance.  Try to prove that every failure
contains a near-kill residue of rank at least 2.

### 4. Paley/Circulant GLMY Path Homology

External problem family: GLMY path homology and newer Fourier/symbol-matrix
methods for circulant digraphs.

Local lever:

- THM-108 proves `beta_2(T)=0` for every tournament.
- The remaining even-Betti and seesaw questions are sharply constrained but
  not fully closed.
- Paley/circulant tournaments have automorphism enough for Fourier reduction.
- The repo already identified Tang-Yau symbol-matrix methods as the right
  external language for circulant path homology.

Possible breakthrough:

Prove a closed Betti formula for Paley tournaments `T_p`, or prove a general
even-Betti vanishing range for all tournaments.  A Paley formula would be a
strong bridge between arithmetic tournaments and directed algebraic topology.

Next experiment:

Implement Tang-Yau symbol matrices for Paley connection sets up to p=19/23 and
compare the ranks with the existing `circulant_homology.py` data.  Search for a
constant-symbol or empty-bad-locus theorem analogous to the repo's T_7/T_11
observations.

## Round 2: Erdős Doors

### 5. Erdős-Moser Transitive Subtournaments, But Soft

External problem: Erdős-Moser/Stearns-type bounds ask how large a transitive
subtournament every tournament must contain.  Erdős Problem #1216 asks a sharp
form and is now marked solved in the negative, but the directed Ramsey problem
behind it remains a central tournament theme.

Local lever:

- A transitive subtournament is a local acyclic island; `H(T)` is global mass
  over all transitive total orders compatible with arcs.
- Paley/random tournaments are extremal for avoiding large transitive sets, but
  they can still have large or structured `H`.
- Krawtchouk low-degree features may estimate "soft transitivity" without
  finding the largest transitive subset.

Possible breakthrough:

Define and prove an entropy/Renyi-style soft Erdős-Moser theorem: every
tournament has either a large transitive subtournament or a large low-order OCF
packing signature.  This would not replace classical directed Ramsey bounds,
but it would give a new quantitative invariant between largest-subset and
total-order count.

Next experiment:

For known extremal Erdős-Moser constructions, compute `H`, root spectrum, and
`alpha_1,alpha_2` against Paley, interval, and random tournaments.  Look for a
separation not visible from largest transitive subset alone.

### 6. Erdős-Schütte Domination and Cycle-Domination Residues

External problem: Erdős Problem #902 asks for the least size of a tournament in
which every n-set is dominated by some outside vertex.  The known bounds are
still far apart.

Local lever:

- The beta_2 proof architecture already uses dominated/free 3-cycles and good
  vertices.
- The SCC/good-cut theorem reduces path-cut obstruction to component
  condensation.
- Projection-residue features distinguish localized cycle mass from flat
  cycle mass.

Possible breakthrough:

Not solve #902 outright, but prove a cycle-domination analogue: bounds on
families of 3-cycles or odd cycles in which every small support set is
dominated.  This could feed back into beta_1/bad-vertex bounds and perhaps
give new tournament-domination lemmas.

Next experiment:

Define `S_k^cycle`: every k-set of vertices, or every k-set of cycle supports,
is dominated by an outside vertex.  Test random, Paley, interval, and H=63
families.  Compare to `bad vertex <= 3`, beta_1, and deletion-loss residue.

### 7. Erdős-Hajnal Via H-Mass

External problem: tournament Erdős-Hajnal asks for polynomial-size transitive
subtournaments in every fixed-H-free tournament.

Local lever:

- `H(T)` counts all Hamiltonian transitive orders, not just the largest
  transitive subtournament.
- OCF converts H-mass into a hard-core partition function of odd-cycle
  packings.
- Forbidden-subtournament families may impose root-spectrum or alpha-vector
  constraints.

Possible breakthrough:

An "H-Erdős-Hajnal" theorem: for a fixed forbidden tournament `F`, every
`F`-free tournament has anomalously high or low normalized H-mass, or a
controlled OCF alpha spectrum.  This would be a soft version of EH and may be
more tractable computationally.

Next experiment:

For all forbidden `F` on 4,5,6 vertices, sample `F`-free tournaments and
measure normalized `H/(n!/2^(n-1))`, alpha-vector depth, and root/residue
features.  Look for monotone families.

## Round 3: Sideways Bridges

### 8. Directed Odd-Cycle Erdős-Pósa and Bermond-Thomassen

External problem family: packing vertex-disjoint directed cycles versus hitting
all directed cycles.

Local lever:

- OCF alpha_k literally counts k-packings of vertex-disjoint odd directed
  cycles.
- The conflict graph `Omega(T)` packages the entire packing problem.
- The H=21 proof already uses 3-cycle matching bounds and poisoning graphs.

Possible breakthrough:

Prove tournament-specific odd-cycle packing/hitting inequalities sharper than
general directed graph theorems.  Even a 3-cycle-only result could strengthen
the forbidden-H machinery.

Next experiment:

For tournaments with bounded `alpha_k` or bounded 3-cycle matching number,
compute minimum odd-cycle transversals.  Compare with Frankl-style matching
bounds and Lichiardopol/Bang-Jensen-Bessy-Thomasse disjoint-cycle results.

### 9. Rédei-Berge Positivity and Noncommutative Deletion-Contraction

External problem family: symmetric functions, chromatic symmetric functions,
Stanley-Stembridge-style positivity, and noncommutative deletion-contraction.

Local lever:

- Grinberg-Stanley express tournament Rédei-Berge functions using odd
  power-sum variables and prove the parity/refinement backbone behind OCF.
- Mitrovic's noncommuting Rédei-Berge function has deletion-contraction, the
  exact inductive property missing in many local proof attempts.
- Hikita's Stanley-Stembridge proof shows positivity can move through
  probability distributions rather than geometry alone.

Possible breakthrough:

Build a tournament positivity/failure atlas: identify which local OCF,
Krawtchouk, or path-homology invariants are shadows of coefficients in the
noncommutative Rédei-Berge expansion.  A clean deletion-contraction proof of
one repo theorem would be publishable.

Next experiment:

Compute the noncommutative Rédei-Berge expansion for the n=8 H=63 classes and
the THM-025 n=9 real-root failure.  Check whether exact-kill and near-kill
residues are visible as small coefficient patterns.

### 10. Quadratic-Residue Codes, Delsarte Bounds, and Circulant Extremizers

External problem family: QR codes, association schemes, circulant Hadamard
shadows, and Delsarte linear programming.

Local lever:

- Paley tournaments are QR objects; their spectra are Gauss sums.
- T_7 and T_23 connect to Hamming/Golay-side code data.
- Krawtchouk coordinates already compress tournament tilings.
- Interval tournaments beat Paley in larger circulant H counts, producing an
  "anti-Ramanujan" extremal story.

Possible breakthrough:

Use coding-theoretic LP bounds to certify circulant H-maximizers or prove a
Paley-to-interval phase transition.  Conversely, tournament H may define a new
weight enumerator for QR/circulant codes.

Next experiment:

Turn circulant connection sets into Krawtchouk weight distributions and fit
`H` as a low-degree functional.  Compare Delsarte feasible regions against
actual p=7,11,13,17,19,23 extrema.

### 11. Condorcet Thermodynamics and Kemeny Algorithms

External problem family: social choice, rank aggregation, Kemeny/feedback arc
set, and property testing of preference cycles.

Local lever:

- Majority preferences form tournaments.
- `I(Omega,2)` is a partition function over independent packets of Condorcet
  cycles.
- `Delta H` and deletion-loss features locate high-leverage comparisons.

Possible breakthrough:

Make a practical theorem/algorithm: small OCF features approximate distance to
transitivity or identify comparisons whose re-evaluation most reduces cycle
temperature.  This is likely more engineering breakthrough than pure theorem,
but it would give the mathematical framework a visible application.

Next experiment:

Implement a weighted Omega approximation for real preference datasets and
compare OCF features to Kemeny local-search improvements and FAST heuristics.

### 12. Quiver/Cluster Mutation as a Toy Exchange Graph

External problem family: quiver mutation, cluster exchange graphs, maximal
green sequences, and scattering diagrams.

Local lever:

- The tiling metagraph is an arc-flip graph with H as an energy.
- Quiver mutation preserves H in the restricted source/sink tournament cases.
- Bucket transport already gives quotient conservation laws for move systems.

Possible breakthrough:

Show that a restricted tournament mutation category has a cluster-like
potential whose specialization is H or OCF.  This is a moonshot, but the
source/sink mutation H-invariance is too clean to ignore.

Next experiment:

Enumerate all vertex mutations that preserve tournament status through n=8.
Record their action on H, Omega alpha-vector, SCC/good-cut height, and
Krawtchouk coordinates.  Look for exchange relations.

## What I Would Bet On

If the goal is a real mathematical result in the next few sessions, the best
targets are:

1. **n=9 real-root failure classification** via deletion residue rank.
2. **Single-core signature image gap** `r_core notin {3,10}`.
3. **Circulant H-maximizer phase transition** with additive-energy/Fejer tools.
4. **Paley path-homology formula** using symbol matrices.

If the goal is a surprise bridge to a famous problem, the best bets are:

1. **Erdős-Schütte domination** through cycle-domination and beta proof lemmas.
2. **Erdős-Hajnal soft transitivity** through H-mass and OCF spectra.
3. **Rédei-Berge deletion-contraction** through noncommutative symmetric
   functions.
4. **Delsarte/QR code bounds** through Krawtchouk tournament coordinates.

The shared meta-pattern is not "tournaments solve everything."  It is more
specific: this repo has become unusually good at finding the residue left by a
projection.  Famous problems often become hard exactly at the place where a
coarse invariant forgets too much.  Our strongest chance is to build residue
invariants that remember just enough.

## Source Trail

- Erdős Problem #1216, transitive subtournament/direct Ramsey question:
  <https://www.erdosproblems.com/1216>
- Erdős Problem #902, Schütte domination problem:
  <https://www.erdosproblems.com/902>
- Noga Alon, "The Maximum Number of Hamiltonian Paths in Tournaments":
  <https://web.math.princeton.edu/~nalon/PDFS/hamilton.pdf>
- Grinberg-Stanley, "The Redei-Berge symmetric function of a directed graph":
  <https://arxiv.org/abs/2307.05569>
- GLMY, "Homologies of path complexes and digraphs":
  <https://arxiv.org/abs/1207.2834>
- Pure pairs / tournament strong Erdős-Hajnal:
  <https://arxiv.org/abs/2202.13977>
- Mitrovic, noncommuting Rédei-Berge deletion-contraction:
  <https://arxiv.org/abs/2504.20968>
- Hikita, Stanley-Stembridge proof:
  <https://arxiv.org/abs/2410.12758>
- Tang-Yau circulant path homology:
  <https://arxiv.org/abs/2602.04140>
