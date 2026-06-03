---
source: codex-2026-06-03-S574
status: synthesis + finite prototype
tags:
  - lonely-runner
  - clock-blockers
  - obligation-hypergraph
  - private-pivots
  - open-gap
  - burnside-labels
  - matroid-shadow
---

# Rearranging the Clock Ledger as an Obligation Hypergraph

S573 gave the right necessary conditions for a row to sit below

```text
E_n = 2/(2n-1).
```

But the `D/U/N` list is still too flat.  The more useful object is a labelled
hypergraph.

```text
vertices  = proof obligations
edges     = speeds
incidence = speed v blocks obligation o
```

The obligation vertices are:

```text
D_q  : some speed is divisible by q, for 2 <= q <= n-1
U_a  : unit antipodal shell {a,2n-1-a} is hit
N_j  : j/n has a runner at distance 0 or 1/n
```

So a possible strict sub-edge row is not primarily a speed set.  It is a full
cover of this obligation hypergraph.

## Why This Rearrangement Helps

Each uncovered obligation is already a proof:

```text
uncovered D_q -> t=1/q clears at or above E_n
uncovered U_a -> S553 unit-inverse clock clears at E_n
uncovered N_j -> j/n clears above E_n
```

That means the only rows worth worrying about are full covers.  The proof
problem becomes:

```text
full D/U/N cover -> M(S) >= 1/n.
```

This is weaker than the false spectral-gap claim, and stronger than simply
enumerating speed sets.

## What S574 Measured

`04-computation/lrc_obligation_hypergraph_s574.py` audits representative floor
rows, open-gap lifts, and nonunit-hole rows.

It prints:

```text
coverage counts by D/U/N layer,
critical obligations with exactly one blocker,
private obligations owned by each speed,
exact maximin witness,
active runners at that witness,
obligation-tournament fingerprints.
```

The useful surprise is that open-gap rows are not weak covers.  They are full
covers with their own private pivots.

For example:

```text
AP n=7 floor:              17 obligations,  9 critical
lifted flip {2} n=7:       17 obligations, 10 critical, M=5/33
AP n=8 floor:              17 obligations, 10 critical
n=8 nonunit floor:         17 obligations,  7 critical
n=8 nonunit open rows:     17 obligations, 10 critical, M=3/23
AP n=14 floor:             34 obligations, 22 critical
V* n=14 floor:             34 obligations, 21 critical
```

The open-gap rows survive the whole D/U/N ledger.  Their exact witnesses are
pair-sum refinements inside the full-cover cell.

## Private Obligations Are the Descent Handles

An obligation is private if exactly one speed covers it.  If a speed owns a
private obligation, any attempted exchange has to preserve that label or the
row immediately exposes a cheap clock witness.

This suggests an exchange/descent proof:

```text
1. Start with a full D/U/N cover.
2. If it has removable nonprivate lift mass, lower or exchange it.
3. If every move hits private labels, record a cover-circuit.
4. Prove every cover-circuit has an n-clock floor witness or a pair-sum witness
   above 1/n.
```

This is close in spirit to endpoint-private-pivot work: a counterexample would
need a leafless, private-pivot-free incidence core.  The S574 rows are not
counterexample-like; they have many private labels.

## Rearrangement 1: Cover-Circuit Normal Form

Treat a row as a cover of obligations.  A minimal full cover is a circuit-like
object: every speed has at least one private obligation.

The AP is the clean lower circuit.  The n=7 lifted flip row

```text
(1,5,6,11,16,17)
```

is a different circuit preserving all unit shells and small denominators while
moving the exact witness to `10/33`.

The next test is to enumerate minimal full covers, not all speed sets.  This
should be vastly smaller and should sort floor circuits from open-gap circuits.

## Rearrangement 2: Labelled Burnside Event Words

HYP-2085/HYP-2087 say the binary lonely time word is too coarse.  It remembers
which grid slots are lonely but forgets ownership.

The obligation hypergraph is the compressed labelled event word:

```text
D_q owner labels,
U_a owner labels,
N_j owner labels.
```

Instead of asking whether `1_X(t)` is invariant, ask whether the labelled
obligation cover is invariant under reflection, lift, or unit action.  The
binary word is the shadow; the D/U/N cover is the proof object.

## Rearrangement 2b: Strong-Lens Ownership

Incoming HYP-2089 says the hard LRC picture is a regular strong encirclement:
the observer is trying to escape to sole-source status, while runners surround
it in a nearly rotational strong block.

The D/U/N private labels say who prevents which escape clock.  That connects
directly to THM-397's endpoint-blocker lemma: a collective cover is not just
"covered"; some runner owns an endpoint or clock obstruction.  In the
hypergraph language:

```text
observer escape clock = obligation vertex
runner preventing escape = hyperedge owner
private blocker = proof debt / descent pivot
```

So the strong-lens proof can be made incidence-theoretic: a regular
encirclement can only remain hard if its escape obligations form a private,
label-stable cover-circuit.

## Rearrangement 3: Three-Gear Automaton

The three layers are three different gears:

```text
D layer: divisibility gear
U layer: odd unit-shell gear
N layer: even floor-clock gear
```

A strict sub-edge row has to synchronize blockers across all three gears.
Open-gap lifts are synchronization states that do not return to the floor, but
still cannot descend below it.

This suggests a finite automaton whose state is:

```text
(covered D labels, covered U labels, covered N labels, private labels)
```

and whose transitions are speed insertions, lift exchanges, and gcd descents.

## Rearrangement 4: Tropical Face Complex

The original LRC maximin is the lower envelope of runner tents.  The D/U/N
ledger is a coarse face complex cut out by clocks that would otherwise clear
the edge.

Inside a full-cover cell, the exact witness is a refined face, often a pair-sum
pinch:

```text
n=7 open lift: t=10/33
n=8 open lift: t=4/23
```

So the proof can split:

```text
coarse cover cell -> no cheap clock
refined tropical face -> pair-sum or floor witness
```

## Tournament Analysis / Assumption Challenge

The S574 tournament is intentionally almost transitive.

- **Vertices:** obligations `D_q`, `U_a`, `N_j`.
- **Observable:** coverage count, then layer, then label.
- **Switch:** orient toward smaller coverage count, i.e. more fragile
  obligation.
- **Tie Hamiltonian path:** `D`, then `U`, then `N`, increasing label.
- **Fingerprint:** transitive score histograms, no directed 3-cycles,
  singleton SCCs, one displayed Hamiltonian path.

This is not meant to discover chaos.  It is meant to expose private pivots.
The challenged assumption is that a useful tournament must have runner
vertices.  Here runner vertices destroy the predicate; obligation vertices
preserve it.

## The New Proof Target

The S573 target was:

```text
D/U/N ledger -> M(S) >= 1/n.
```

S574 makes it more operational:

```text
full D/U/N cover
-> private-pivot exchange or cover-circuit normal form
-> floor or pair-sum witness
-> M(S) >= 1/n.
```

That is the creative rearrangement I trust most right now.
