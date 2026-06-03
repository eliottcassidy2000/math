---
source: codex-2026-06-03-S582
status: exploratory bridge from Tournament Analysis to finite automata
tags: [Tournament-Analysis, automata, LRC, endpoint-cover, middle-state, HYP-2109]
---

# Tournament automata with left, middle, and right

The clean connection is this:

```
tournament = the binary shadow of many pair automata.
```

A normal Tournament Analysis switch says "compare two things and orient the
edge."  A three-state automaton says "watch the comparison evolve, keep the
middle state visible, and orient the edge only at the end."  The middle state is
where the proof lives.

## The three states

- `L`: the left endpoint/owner/chart currently has the obstruction.
- `M`: the pair is on the wall, tie seam, midpoint corridor, or live residue.
- `R`: the right endpoint/owner/chart currently has the obstruction.

The smallest useful transition rule is the wall rule:

```
M --L--> L      M --R--> R
L --R--> M      R --L--> M
L --M--> M      R --M--> M
```

So an edge cannot flip from `L` to `R` directly.  It must pass through `M`.
That is exactly the behavior we want for LRC: a proof obstruction should not
teleport from one owner to the other; it should cross a visible wall where
endpoint-cover, residue, or tie-wall data can be inspected.

## Why HYP-2108 makes this sharper

HYP-2108 says the endpoint-cover residual has one scalar obstruction

```
P(S) = max_i ( ||v m_i|| + (v/2) l_i - 1/n ).
```

Tightness would mean every component remains in the middle corridor:

```
||v m_i|| <= 1/n - (v/2) l_i.
```

In automaton language, a counterexample is trying to keep all component cells in
terminal `M` around a closed circuit.  The summed corollary says that would
force the average midpoint phase below `1/n`, while experiments sit near the
generic `1/4`.  So the hard target becomes a no-closed-middle theorem:

```
all local automata terminal M around the circuit is impossible.
```

Then all non-middle exits become owner obstructions, where Lemma C, B', private
pivots, or CRT/dominance can act.

## Addition versus multiplication

Addition creates the summand graph: shells `{a, C-a}`, pair sums, endpoint
pinches, and cheap-pair witnesses.

Multiplication by units acts on that additive graph by permuting shells and
turning local witnesses into clocks.  In the automaton frame, addition supplies
event words and multiplication supplies time changes of those words.

Odd versus even explains why `M` keeps reappearing.  Odd `C=2n-1` has no
self-antipodal midpoint shell, so the summand graph tends to force a side.
Even denominators or composite odd `C` create apex/nonunit holes where a pair can
sit in `M` longer.  The n=14 story has both: clean odd geometry at `C=27`, but
composite unit/nonunit arithmetic that preserves middle states until a
certificate kills them.

## Algorithmic consequence

For a proposed safe-box orbit hit, do not sample all times first.  Compile each
local proof obligation into a three-state cell:

1. event alphabet from endpoint owners, residues, arc centres, or gluing charts;
2. transition rule that forces side changes through `M`;
3. terminal projection to a tournament with a declared tie path;
4. fingerprints of the binary shadow plus a separate audit of terminal and
   transient `M` states.

This is cheaper than exact time enumeration because most pairs become stable
`L/R` early.  Only cells that enter or remain in `M` require arithmetic detail.
For Cprime, those cells are exactly the all-short large-owner residue/circuit
residuals of HYP-2107 and HYP-2108.

## Assumption challenge

The vertices should not default to runners.  For this problem the stronger
vertices are:

- components of `G(S')`;
- endpoint-owner clauses;
- cover arc centres;
- residue automaton states;
- certificate germs in the apex/sheaf picture;
- proof obligations in the certificate calculus.

The quotient preserves owner side and middle survival.  It destroys magnitude:
lengths, residue moduli, and phase distances must stay as labels if they are
used by B', Lemma C, dominance, or `P(S)`.

## Wild extension

Build a "middle graph" in parallel with the tournament.  Vertices are local
cells currently in `M`; edges record shared owners, shared residues, or shared
arc centres.  A tight counterexample would need the tournament shadow to avoid a
private-pivot exit while the middle graph contains a closed circuit.  HYP-2108
then attacks the circuit by midpoint average; HYP-2107 attacks it by bounded CRT
emptiness.

That would turn the current frontier into a two-layer proof object:

```
binary tournament for exits
middle automaton graph for the impossible residual
```

This feels like the right simple machine: left, right, and the one place where
the proof can still hide.
