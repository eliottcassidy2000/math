---
id: HYP-2701
title: LRC14 sparse residual tails should be absorbed by a generated-context deficit automaton
status: OPEN; namespace reserved for exact scout
source: codex-2026-06-20-S64
tangent: T936
depends_on:
  - HYP-2697
  - HYP-2698
  - HYP-2700
related:
  - HYP-2699
  - HYP-2694
  - THM-557
  - THM-558
  - THM-555
---

# HYP-2701 - Sparse Tail Deficit Automaton

## Claim

After HYP-2700, consecutive blocks already dominate arbitrary bounded cluster
shapes on the top residual layers `|R| >= 5` and on the true full-cover
functional `p_0`.  The remaining HYP-2697/HYP-2698 obstruction is the sparse
tail `|R| <= 3`, where nonconsecutive shapes can beat the consecutive block on
named residual packets.

The expected repair is that those sparse-tail wins are absorbed by the
generated context word:

```text
x -> w_x(R),  R subset {1,...,6},
```

where `w_x` comes from the finite sector-mask OR/deletion automaton and is best
read in miss-zeta product coordinates.  The finite proof object should be a
**deficit automaton** whose states record

```text
(|R|, sector geometry of R, hit-count deficit, miss-zeta context weight).
```

The automaton vertices are not runners or sectors alone.  They are proof-state
deficits: a residual packet that wins for an arbitrary shape, together with the
minimum generated context pressure needed to neutralize that win.  Edges compare
state transitions under adding singleton context carriers, merging singleton
context carriers into coherent blocks, or lifting `|R|` toward the already-safe
large-residual layers.

## What Is Claimed Now

Only the namespace and route are claimed in this stub.  The immediate scout
will:

1. enumerate the exact `|R|=3,4` violators from HYP-2700;
2. decompose each violator by residual-sector geometry and hit-count gaps;
3. test singleton-product context kernels `g_r` as a neutralizing automaton
   step;
4. compare singleton kernels against coherent generated contexts from
   HYP-2698;
5. report tournament fingerprints for the proof-state automaton: score
   histogram, directed cycles, SCCs, edge flips, and Hamiltonian-path count.

## Assumption Challenge

The session explicitly rejects the easy vertex choices:

- **not runners:** raw runner vertices forget the generated residual language;
- **not residual coordinates:** arbitrary positive weights make the theorem
  false;
- **not binary tournament arcs:** the `1/2`-scale tournament is only a correlate
  of the exact `Z/7` coloring.

The quotient preserves the LRC predicate "context-generated residual pressure
times cluster capacity" and destroys the individual carrier phases.  The
challenged assumption is that sparse residual packets must be handled by a
separate ad hoc finite check; HYP-2701 tries to make them a finite automaton
frontier attached to HYP-2698's miss-zeta product word.

## Status

No LRC14 proof is claimed.  This is a live stub reserved before running the
exact sparse-tail scout.
