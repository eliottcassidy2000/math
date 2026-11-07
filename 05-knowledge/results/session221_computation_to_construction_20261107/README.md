# Session221: from machine semantics to tournament constructions

**Checkpoint 1: audited single-arc difference construction, bounded-target compiler and recovered provenance. Research continues. Chan–Pak completeness remains OPEN.** No literature-priority claim or unrestricted spectrum extension is made.

The session starts at main commit `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`, with Session220's exact binary envelopes, selector obstruction and bounded ribbon screen. The aim is a literal construction, not an inference from computational universality.

## Every even gap from one physical arc flip

For every integer N>=1, the construction gives **two strongly connected tournaments on 5+ceil(log2 N) vertices**, differing in exactly one arc, whose Hamiltonian-path counts differ by exactly 2N. This order is optimal within the adjacent-chain, one-apex flip class, not among all tournaments.

For N>=2 the actual counts are 12N-9-f(N-1) and 14N-9-f(N-1), where f(0)=0, f(2u)=f(u), and f(2u+1)=f(u)+u+(u mod2). The N=1 case has counts9 and11 on five vertices. The explicit adjacency constructor is logarithmic in the target gap. The baseline is controlled, not arbitrarily selectable: this is **difference-spectrum coverage**, not absolute spectrum coverage.

The underlying weighted tape formula and adjacent response were already proved in the October20 feedback-word packet. The contribution is their explicit all-even inverse, single-chain-edge realization, sharp size within this class, strong lift and parallel interaction law. A radius-one exclusion cellular-automaton step corresponds to simultaneous disjoint chain-edge flips plus relabeling. Its exact count change includes positive cross terms; adding isolated increments would be false.

[Physical construction, complete marked gap intervals and proofs](physical/PROOF.md) · [Independent audit](audit/PHYSICAL_AUDIT.md).

## A precise machine-to-tournament transfer

The inherited nonnegative-polynomial saturation lemma for register machines and the complete capped inverse method from Sessions203/216 apply to the inherited positive ribbon profile automaton. For each fixed span and target cap B, coordinatewise truncation at B+1 gives a finite, exact automaton for **all word lengths**. Breadth-first search either returns a shortest ribbon word for each count at most B or certifies absence within that ribbon family.

This is a target-dependent compiler, not one fixed finite-state selector for all targets. Its worst-case state bound is enormous. It does not show that any given unbounded sequence of targets is attained or that witnesses have logarithmic order.

The new implementation and proof application are separated from the inherited inputs and prior fixed-point computation. Span-two cap99 was already completely known; reproducing it is a control, not a coverage advance. Signed reduced transfer matrices cannot be capped: their cancellations destroy the quotient law. Likewise saturation does not commute with arbitrary machine decrement microsteps.

- [Exact compiler theorem and limits](compiler/THEOREM.md).
- [Recovered source map and transfer contracts](compiler/PROVENANCE.md).
- [Independent audit](audit/CAP_AUDIT.md).

## Inheritance and live concept board

- **Anchor:** physical target selection for Hamiltonian-path counts.
- **Niche:** exact local cellular-automaton moves on a marked apex tape.
- **Wildcard:** variable matching seams and positive continuants.

The five active comparisons are machine macrostep versus decrement microstep; tape state versus scalar path count; positive saturation versus signed quotient; final-state compression versus continuation; and marked seam connectivity versus scalar gate arithmetic.

Closest inherited mechanisms are Session133's 32-fraction one-apex counter, Session135's macrostep saturation lemma, Sessions203/216's capped inverse, THM-4880's positive ribbon transfer and THM-4881's cap99 fixed point. The canonical hostile is equal H with different next-step H. Session220's finite-state witness obstruction is a corrected near miss, not a prohibition on target-dependent algorithms. The least-used sidecars are the actual tape origin/marked boundary, positive path-cover coordinates and transient endpoint roles.

The older TT3/C3 block encoding faithfully carries tape Hamming distance but has H=3^(number of ones), so it forgets word order. The TM-to-CA/FRACTRAN simulations preserve machine execution, not tournament-count surjectivity. The present work keeps those different predicates explicit.

## Open boundary

A bounded inverse is not an all-target proof. A machine simulation is not a size-efficient tournament constructor. The physical difference theorem has passed independent proof and literal-graph audits; fixed-core inverse compression and a controlled seam sublanguage are still being finalized. The missing global step remains a construction selecting a valid tournament for every required odd target.

## Checkpoint reproduction

Run `python3 verify_all.py` from this directory. Python standard library only. It checks three completed author fixed points and three independent BFS/graph audit runs, with all generated receipts confined to a temporary directory.

The author runs contain 5,102 states at span2/B99, 14,618 at span2/B199, and 7,025 at span3/B25. The independent span3/B47 run closes at 22,399 nonempty states and realizes47. These ranges and exclusions are controls against inherited results, not new coverage records. The auditor's BFS includes one empty initialization state; state-count conventions are stated in its report.
