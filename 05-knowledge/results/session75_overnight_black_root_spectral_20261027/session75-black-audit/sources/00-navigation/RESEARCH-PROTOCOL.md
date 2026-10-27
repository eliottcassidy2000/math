# Research Protocol

**Status:** CURRENT
**Scope:** mathematical research sessions in this repository
**Purpose:** make exploration self-directed, connective, rigorous, and cumulative without forcing every idea through the current flagship problem.

The main prize is **LRC(14)**. It is not the only worthwhile mathematics in
the repository. Direct progress, a decisive refutation, a reusable instrument,
an exact reformulation, a recovered niche thread, and a genuine bridge between
subjects are all valid session outcomes.

## 1. Truth and freshness

Fetch before relying on the frontier and again before filing a result. If the
tree is clean, rebase/pull before building on incoming work. If it is dirty,
inspect `HEAD..origin/main` without rebasing, checkpoint coherent owned work,
then synchronize; never hide or sweep another session's changes. Treat incoming
work as mathematical signal: ask whether it changes the current object,
invariant, proof route, computation, or vocabulary.

During long computations, fetch and reassess at natural checkpoints; rebase
only from a clean tree. Compare new work before choosing the next batch; do not
silently change the inputs of an already running reproducibility job.

Use this precedence order:

1. explicit correction or retraction together with its repaired theorem;
2. current canon with a proof and scoped statement;
3. reproducible computation with its matching source, output, and controls;
4. hypothesis detail and current navigation ledgers;
5. session logs, reflections, drafts, and historical outputs.

`01-canon/MISTAKES.md` overrides older claims. Session logs and reflections
are idea provenance, not truth authorities. Before inheriting a target, search
canon for its **statement, constants, quantifiers, theorem identifiers, and
synonyms**. Searching only for the proposed method is insufficient.
Frontmatter status controls even inside canon: a `RESERVED` proof candidate
stays outside the proof graph until adversarial audit and explicit promotion.

## 2. The session portfolio: Anchor / Niche / Wildcard

Every open-ended mathematical session maintains three lanes:

- **Anchor:** the principal objective, normally a live LRC(14) obligation or
  another explicitly assigned problem.
- **Niche:** an underexplored or orthogonal thread with reusable artifacts,
  a corrected-but-close conjecture, an empty operation/dual cell, or a result
  that has been proved but not connected forward.
- **Wildcard:** a freely generated mathematical compulsion, analogy, or
  unusual object with a cheap decisive probe.

The lanes are not quotas. A niche or wildcard may become the anchor when it
produces an exact object, a strong counterexample, a transferable lemma, or a
new instrument. Do not manufacture relevance to LRC(14); do ask whether a real
map, obstruction, extremal, or proof mechanism transfers.

When a promising compulsion appears, record it immediately with its cheapest
falsification test. Run that first probe when affordable. Positive signals
deserve pursuit even when the route is indirect.

## 3. Keep an active concept board

Maintain a small board of roughly three to seven live concepts. For each one,
record:

```text
object / representation
predicate or question
known invariant or extremal
operation currently being applied
lost coordinate or obstruction
cheapest next test
```

Whenever a result, pull, paper, or new object arrives, compare it against every
board item:

- Is there a common invariant, extremal family, recurrence, or obstruction?
- Is one a quotient, lift, dual, limit, localization, or composition of the
  other?
- Does a proof or counterexample transport?
- Does one object retain information the other discarded?
- Does an operation law explain both?

This comparison is recursive: every real connection may create a new board
item, operation, or sidecar.

## 4. Generate perspectives procedurally

Use the research surface

```text
objects
  x representations/lenses
  x invariants
  x operations
  x symmetries/quotients
  x scales/regimes.
```

Typical representations include sets, indexed sequences, supports with
multiplicity, graphs, tournaments, hypergraphs, lattices, automata, codes,
measures, generating/Dirichlet functions, event words, and proof obligations.
Typical operations include deletion, adjoining an observer, complement,
duality, composition, product, quotient/lift, localization, reduction modulo
a prime, dilation, degeneration, differentiation, Fourier/Mellin transform,
and passage to a limit.

Blank cells are prompts. Generate several candidate cells, rank them by
novelty, structural fit, and cost of falsification, and test at least one when
the session is exploratory. Do not assume the natural vertices, variables, or
quotient are the faithful ones.

## 5. The connection contract

A connection is mathematical only when it specifies:

```text
source object -> target object
the map or correspondence
the predicate or theorem preserved
the information destroyed
the sidecar needed to restore that information
one concrete transfer, prediction, or falsification test
```

Shared vocabulary, a matching decimal, or a familiar shape is an analogy, not
a bridge. Keep analogies as tangents until the contract is filled.

For wall/arrangement models, the contract must also say whether the target is
the bare wall set, ordinary complement, thickened complement, or selected
inequality cells, and whether orientation, owner, height, or deletion labels
are part of the state.

## 6. Explain why, not only whether

For a **true** statement, record:

- the mechanism that makes it true;
- equality and boundary cases;
- minimal dependencies and exact scope;
- the coordinate in which it becomes natural;
- the strongest evident generalization and its likely failure boundary.

For a **false** statement, record:

- the smallest or cleanest witness;
- the first invalid implication, type change, or lost coordinate;
- the strongest part that survives;
- a repaired statement or replacement object;
- the error genus and where else it may recur;
- the new question exposed by the failure.

A negative verdict without this anatomy is unfinished research.

## 7. Validation gate

Before promoting a claim, complete the applicable checks.

### Mathematical typing

- Name the objects and domains precisely.
- Separate indexed sequence from support, scalar equality from polynomial
  identity, labels from isomorphism classes, node values from exponents, and
  affine dimension from Poisson or Weyl rank.
- State whether each implication is necessary, sufficient, or an equivalence.
- Audit quantifiers and exhibit a non-vacuity witness.

### Symmetry and quotient discipline

- List translation, dilation, relabeling, complement, reflection, and other
  relevant actions.
- Search orbit representatives rather than an arbitrary coordinate slice.
- State which predicate a quotient preserves, what it forgets, and which
  sidecar carries the forgotten coordinate.

### Hostile probe

- Attack the claim for at least a few minutes before extending or naming it.
- Test canonical constructions, equality cases, boundaries, degenerate cases,
  structured adversaries, and scale continuation.
- For any proposed uniform bound, first try to construct an indefinitely worse
  family.
- Cross-check against already proved bounds and known extremals.

### Computation

- State the exact universe, filters, and exclusions; negative results inherit
  every filter.
- Include positive controls and, when feasible, a second implementation,
  pruning configuration, relabeling, or optimized/unoptimized replay.
- Make the program print the theorem's **consequence object**, not merely an
  intermediate ratio or the assumptions injected into the model.
- Save the source and matching output together with reproduction commands and
  hashes for load-bearing results.

### Formalization

- A successful build proves only the statement encoded.
- Check that hypotheses are satisfiable and match the paper theorem.
- Check the intended module is in the root import closure.
- Record the build command, axiom audit, and whether any finite decision is
  kernel, native, or externally certified.

### External literature

- Read the relevant primary-source theorem and its hypotheses, not only an
  abstract or secondary summary.
- Check the source version and equality-case content.
- Treat very recent provenance and priority as uncertain until a stable source
  exists.

## 8. LRC-specific route triage

Before investing in a new near-floor LRC route, state whether it:

1. breaks translation invariance when the target does;
2. respects dilation invariance by working on primitive orbits;
3. controls a maximum or tail rather than only a mean;
4. retains phase/sign information or is exact enough to survive cancellation;
5. adapts across moduli or performs an honest complete finite enumeration.

Failure of an axis does not make the mathematics useless, but it sharply limits
what the route can prove. Record that limit before computation.

## 9. Tournament Analysis without forcing a tournament

Tournament Analysis is a preferred lens when an intrinsic pairwise relation
preserves a named target predicate. It is not a requirement to manufacture a
tournament.

- Declare the vertices, pairwise observable, switch/gauge, and tie behavior.
- Exact ties remain ties, a preorder, or confluent data unless a mathematically
  meaningful rule resolves them. Do not add a cosmetic tie-breaker merely to
  obtain a tournament.
- Explain why the chosen vertices are more faithful than alternatives such as
  gaps, events, residues, blocks, modes, circuits, or proof obligations.
- Report what the tournament quotient preserves and destroys.
- If the quotient yields only a transitive ranking with no target content, say
  so and use a better object.

## 10. Make the process cumulative

Log hypotheses whether confirmed or refuted. Preserve scripts and outputs.
Write a reflection when a mechanism connects domains or changes how future
work should be done. Add a method card to `META-PATTERNS.md` only when it has a
clear trigger, action, counterindication, and evidence from at least two
distinct threads, or from one severe failure followed by a demonstrated
repair.

At close-out, report separately:

```text
direct progress
indirect/niche progress
connections created
claims refuted or narrowed
new reusable methods or artifacts
honest remaining frontier
```
