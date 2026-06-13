# Formalization Root Signs And Packet Modules

Session date: 2026-05-30  
Prompt: spend a long formalization session, use it for inspiration, search around, think freely, and pursue tangents.

## Executive Thesis

The useful formalization move was not to encode a large theorem from the
outside of the project.  It was to turn the smallest representation-theoretic
atom from the previous session into Lean:

```lean
root i j = e_i - e_j
```

An oriented tournament arc is a signed type-A root.  The root is not yet the
tournament, and that is the point.  It is the linear atom from which many of
the repo's compressed invariants can be rebuilt: reversal is negation,
consecutive arcs telescope, closed walks have zero total root, and the directed
triangle is the first nontrivial closed-walk relation.

This changes the formalization target.  The primitive should not be "prove a
fact about a 3-cycle."  The primitive should be:

> Paths carry endpoint charge.  Closed paths carry no endpoint charge.

Once that is in the system, a 3-cycle is just the first visible relation in the
kernel of the endpoint map.

## What Was Formalized

Added `TournamentH7.RootSigns` with namespace `TypeA`.

The ambient lattice is:

```lean
abbrev RootSpace (n : Nat) := Fin n -> Int
```

The signed type-A root is:

```lean
def root {n : Nat} (i j : Fin n) : RootSpace n :=
  fun k => (if k = i then 1 else 0) - (if k = j then 1 else 0)
```

The module proves, without project axioms:

- `root_self`: the degenerate root `e_i - e_i` is zero.
- `root_swap`: reversing orientation negates the root.
- `root_backtrack_sum`: a two-edge backtrack has zero total root.
- `root_add_root`: consecutive roots telescope.
- `root_eq_zero_iff`: the only zero signed root is degenerate.
- `root_cycle_sum`: the directed triangle relation telescopes to zero.
- `walkRootSum_append_single`: any finite vertex walk telescopes to the root
  from first vertex to last.
- `walkRootSum_closed`: any closed vertex walk has zero total root.

The important upgrade is `walkRootSum_append_single`.  The triangle theorem is
now a corollary of a path theorem rather than a bespoke arithmetic exercise.

## Formalization Lesson

The repo has many scalar invariants: `H`, `alphaCount`, bucket counts, SCC
defects, residue ranks, phase profiles, and transport row sums.  Formalizing
the scalar identities first is tempting, but it risks freezing late-stage
compression as primitive.

The better order is:

1. Formalize atoms.
2. Formalize maps.
3. Formalize kernels/residues.
4. Only then formalize scalar evaluations.

`RootSigns.lean` is an atom-and-map layer:

- atom: `root i j`,
- map: finite walk to root sum,
- kernel: closed walk has zero total root,
- scalar evaluations: none yet.

This is why the module feels smaller than the idea it unlocks.  It is not an
OCF proof.  It is the beginning of the linear language in which OCF may become
less mysterious.

## Representation-Theory Reading

Type A is the right first formal layer because the tournament vertex set
already carries the Weyl action of relabelling.  A tournament is a sign choice
on the positive roots:

```text
for each unordered pair {i,j}, choose either e_i - e_j or e_j - e_i.
```

Relabelling is the Weyl group action.  Reversal is global sign negation on all
chosen roots.  Score is a divergence-like statistic.  Directed cycles are curl
or kernel witnesses.  Hamiltonian paths are ordered chains of roots whose sum
collapses to a single endpoint root.

That last sentence is the new formal door:

```text
(e_v0 - e_v1) + (e_v1 - e_v2) + ... + (e_v{k-1} - e_vk)
= e_v0 - e_vk
```

So a Hamiltonian path has a huge internal ordering but a tiny endpoint shadow.
All of the path multiplicity counted by `H(T)` lives in the fiber over endpoint
roots plus the compatibility conditions imposed by the tournament chamber.

This reframes Hamiltonian-path counting:

- endpoint root is the coarse quotient,
- path order is the fiber,
- tournament signs are chamber inequalities,
- odd-cycle packets are obstructions or residues inside those fibers.

## The Three-Cycle Reinterpreted

The directed triangle is usually introduced as the smallest odd directed cycle.
The root formalization says something more structural:

```text
(e_i - e_j) + (e_j - e_k) + (e_k - e_i) = 0
```

So the 3-cycle is the minimal exactness relation, not merely the smallest
cycle.  It is the first place where local tournament orientation becomes a
closed root packet.

That matters for OCF.  OCF says:

```text
H(T) = I(Omega(T), 2)
```

where `Omega(T)` is built from odd directed cycles.  If cycles are kernel
relations in the type-A root lattice, then `Omega` is not only a conflict graph
on cycle supports.  It is a conflict graph on root-kernel packets.  The scalar
evaluation at `2` may be shadowing a module or chain-complex identity.

The wild target becomes:

> Replace scalar OCF by a packet-module OCF: Hamiltonian path fibers decompose
> into independent odd-cycle kernel packets, and evaluating a graded character
> at fugacity `2` recovers `H`.

This is not proved.  But it is now a clearer formalization target because the
kernel relation has a Lean seed.

## Packet Module Hypothesis

The phrase "packet module" should mean something precise enough to formalize:

- A packet is a finite closed-walk support with zero total root.
- An odd packet is a closed-walk packet whose combinatorial support is an odd
  directed cycle or a generated odd-cycle relation.
- A packet module is the free abelian group on such packets modulo boundary
  and support-identification relations.
- `Omega` becomes an incompatibility graph on packet basis elements, not merely
  on cycles as opaque objects.
- `I(Omega, 2)` becomes a character/evaluation of a packet exterior or
  independence algebra.

This suggests a stronger version of the representation-lens target from the
previous session:

```text
scalar OCF:
    H(T) = I(Omega(T), 2)

packet OCF:
    HP(T) as a chamber fiber has a filtration whose associated graded object
    is assembled from independent odd-cycle packets, and the dimension/trace
    evaluation at fugacity 2 gives the scalar OCF.
```

The Lean path to this should not start with `H`.  It should start with closed
walk packets and their supports.

## Formal Roadmap

The next Lean layers should be incremental.

1. Root support.

Define the support of `root i j` as `{i,j}` when `i != j`, and empty or
singleton-degenerate when `i = j`.  Prove `root_eq_zero_iff` is compatible with
support emptiness.

2. Directed edge as a root sign.

Bridge `Tournament.arc i j = true` to choosing `root i j`.  Keep this separate
from `RootSigns.lean`, which should remain type-A infrastructure independent
of tournaments.

3. Walk compatibility.

Define a tournament-compatible vertex walk: every consecutive pair follows an
arc.  Then reuse `walkRootSum_append_single` without reproving algebra.

4. Directed cycle as closed root packet.

For `DirectedCycle T k`, prove the associated `walkRootSum` is zero.  This
should be a structural theorem, not arithmetic.

5. Odd-cycle packet.

Add a small record containing a directed cycle, oddness, support, and zero-root
certificate.

6. Packet conflict.

Relate cycle-support disjointness or intersection to packet support.  This is
where the current `Omega` object can start being upgraded.

7. Packet algebra.

Define independence monomials over packets.  The scalar independence
polynomial is then a specialization of a packet algebra, not the primitive.

8. Character/evaluation.

Ask what the evaluation `x = 2` actually counts.  It may be a two-state choice
per independent packet, a sign representation fiber, or a boundary-orientation
choice.

## Tangents From The Formal Seed

### Incidence Matrices

The root map is the incidence boundary map of the complete directed graph:

```text
edge (i,j) -> e_i - e_j
```

Closed walks lie in the kernel.  This is ordinary graph homology, but in this
repo the tournament orientation selects one signed representative per
undirected edge.  That means many "tournament" questions are really questions
about chambers in a fixed incidence complex.

This gives a practical target: compute Smith normal forms of packet incidence
matrices.  Endpoint-transfer sessions already warned that support matching is
not enough; rank and torsion matter.  Root signs make that warning canonical.

### Score As Divergence

For a tournament, sum all chosen roots:

```text
R(T) = sum_{chosen arcs i->j} (e_i - e_j)
```

The coordinate at vertex `v` is:

```text
outDegree(v) - inDegree(v) = 2 outDegree(v) - (n-1)
```

So the score sequence is the divergence of the global root flow.  Regular
tournaments are divergence-free.  Paley tournaments are not just score-flat;
they are root-flow divergence zero with extra character symmetry.

This suggests a phase/residue split:

- divergence: score channel,
- curl/kernel: cycle packet channel,
- harmonic/character part: Paley/circulant phase channel.

The repo's residue-phase-incidence split may literally be a Hodge-style split
on the complete graph/tournament chamber.

### Hamiltonian Paths As Endpoint Fibers

Every Hamiltonian path from `a` to `b` has root sum `e_a - e_b`.  Thus `H(T)`
splits into endpoint-root fibers:

```text
H(T) = sum_{a != b} H_{a,b}(T)
```

This is a natural representation-theoretic refinement because the Weyl group
acts on endpoint roots.  It suggests replacing raw `H` by an endpoint-root
distribution.  For self-complementary or regular tournaments, this distribution
may have symmetries invisible in the scalar.

Practical computation: add `endpoint_hamiltonian_path_counts(T)` and compare
Paley, interval, transitive, H=63 single-core, and THM-025 examples.  The
single-core cases should have concentrated packet residues but maybe still
structured endpoint fibers.

### OCF As A Boundary Formula

OCF may be a disguised boundary formula:

```text
path fiber count = boundary-compatible independent packet count
```

Odd cycles are not arbitrary cycles.  They are the packets that change parity
or orientation choice in a way that affects Hamiltonian path counting.  The
evaluation at `2` says each independent packet contributes a binary degree of
freedom.  In formal terms, this smells like an exterior algebra or Boolean
algebra of packet toggles.

### Applications

The applied ranking idea gets a cleaner explanation.  Pairwise comparison data
is a partially observed chamber in the type-A root arrangement.  Completing a
comparison chooses a root sign.  Inconsistent preference cycles are closed
root packets.  Active ranking should query the comparison that most reduces
packet uncertainty, not merely the pair with the highest local entropy.

This gives a feature language for product systems:

- endpoint-root fiber uncertainty,
- divergence/score imbalance,
- closed-packet residue,
- packet conflict independence,
- expected packet collapse after a query.

That is a more algebraic version of the earlier OCF-guided active ranking
proposal.

### Speedups

A formal root layer also hints at computational speedups:

- Cache Hamiltonian path counts by endpoint root.
- Use divergence channels to prune impossible regular/SC candidates.
- Use closed-packet support to update `Omega` incrementally after edge flips.
- For circulants, diagonalize the root-flow and packet channels by characters.
- For active ranking, maintain packet residues rather than recomputing all
  odd cycles after each queried pair.

The main computational bet is incremental packet maintenance.  Edge flips
alter only packets crossing that root.  If `Omega` can be maintained as a
packet conflict graph with local updates, `H`-related acquisition functions may
become usable beyond tiny tournaments.

## Formalization Principles Going Forward

The session suggests a rule for the whole repo:

> Before formalizing a theorem, classify every object in the theorem as an
> atom, map, kernel, quotient, residue, phase channel, incidence matrix, or
> scalar evaluation.

For example:

- `root i j`: atom.
- `walkRootSum`: map.
- closed walks: kernel.
- endpoint-root fibers: quotient.
- odd-cycle packets: kernel basis candidates.
- `Omega`: incidence/conflict graph on packets.
- `I(Omega,2)`: scalar evaluation.
- score sequence: divergence/projection.
- good-cut defect: quotient residue.
- endpoint transfer: incidence matrix/rank.

This classification prevents premature scalar formalization.

## Verification Notes

`lake build TournamentH7.RootSigns` completed successfully.  A direct
`#print axioms` check for the main root theorems reports only Lean/Mathlib
foundation axioms (`propext`, `Quot.sound`) and no project axioms.

An all-up `TournamentH7.Verify` build was attempted because the audit file was
updated with `#print axioms` wrappers for the new root facts.  That full audit
could not be completed in this session: the local Mathlib/build cache was
cold, and overlapping Lake processes initially fought over the same `.olean`
outputs.  After killing the stale process, a direct single-process Lean check
of `Verify.lean` still stopped before the new root audits because an existing
dependency artifact, `TournamentH7.GoodCuts.olean`, was absent.  The root module
itself is verified independently.

This is a process lesson too: future Lean sessions should avoid starting more
than one Lake build in this workspace, and should build the narrow module first
before attempting the all-up audit.

## Concrete Next Session

The next formalization session should add a separate module, perhaps
`TournamentH7.RootPackets`, with no new grand theorem:

1. Import `TournamentH7.RootSigns` and `TournamentH7.Cycles`.
2. Define the vertex list associated to a `DirectedCycle`.
3. Prove its `walkRootSum` is zero by `walkRootSum_closed`.
4. Define `CycleRootPacket`.
5. Add support and disjointness lemmas.
6. Only then touch `Omega`.

If that works, the repo will have crossed from "root signs are a metaphor" to
"odd cycles are certified kernel packets."  That is the first real foundation
for representation-refined OCF.
