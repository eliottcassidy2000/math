# Session221: from machine semantics to tournament constructions

**COMPLETE: audited single-arc difference construction, fixed-apex inverse, controlled seam language and recovered machine/CA provenance. Chan–Pak completeness remains OPEN.** No literature-priority claim or unrestricted spectrum extension is made.

The session starts at main commit `38fdf9f234159aaedb2d0ba1e0f7bcdefb49e817`, with Session220's exact binary envelopes, selector obstruction and bounded ribbon screen. The aim is a literal construction, not an inference from computational universality.

## Every even gap from one physical arc flip

For every integer N>=1, the construction gives **two strongly connected tournaments on 5+ceil(log2 N) vertices**, differing in exactly one arc, whose Hamiltonian-path counts differ by exactly 2N. This order is optimal within the adjacent-chain, one-apex flip class, not among all tournaments.

For N>=2 the actual counts are 12N-9-f(N-1) and 14N-9-f(N-1), where f(0)=0, f(2u)=f(u), and f(2u+1)=f(u)+u+(u mod2). The N=1 case has counts 9 and 11 on five vertices. The explicit adjacency constructor is logarithmic in the target gap. The baseline is controlled, not arbitrarily selectable: this is **difference-spectrum coverage**, not absolute spectrum coverage.

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

A bounded inverse is not an all-target proof. A machine simulation is not a size-efficient tournament constructor. All four packets have passed independent mathematical and executable audits. The missing global step remains a construction selecting a valid tournament for every required odd target.

## Fixed-apex inverse with a logarithmic physical horizon

For a fixed labeled apex tournament A on k vertices, every attained count m has a same-apex tape witness with at most

    L(m)=max{2a+3b : 3^a 5^b<=m} <= floor(3 log_5 m)

chain vertices. Delete singleton chain strong components. Retaining at most one chain source and one chain sink additionally preserves all apex start/end and joint endpoint counts, at cost at most two extra vertices.

The classical Busch bound and global SCC compression are inherited. The specific construction retains the labeled apex core, supplies an explicit mask deletion test and endpoint-aware normal form, and combines them with the inherited positive transfer to give an exact inverse in O_k(m^(3k log_5 2) log m) bit operations. This is polynomial in the numerical target, not its binary length. It decides membership in this family; it does not prove that the family accepts every allowed target. Compression is a final-word operation, not a valid prefix-state quotient.

[Proof and compiler](automata/PROOF.md) · [Independent endpoint-DP audit](audit/FIXED_CORE_AUDIT.md).

## A programmable seam macro with an honest inverse

The recovered O^m J^2 operation, m>=1, has the exact two-register action

    (h,R) -> (h+qR,2h+qR),  q=2^(m+1),

where R counts root-starting two-path covers. A chain of these gates is a literal crossing matching of backward arcs. Its counts obey a variable continuant, each reachable pair has a unique reversible program, and output witnesses have logarithmic order.

This is a reduced sublanguage of existing O/Z/J operations, not a new primitive. Its scalar spectrum has density zero: at most O(X^(log_2 phi)) values up to X, and only residues 1/5 modulo 8. Scalar H is insufficient even here: two programs give H=89 with R=98 or 106, and the same next gate produces 481 or 513. The exact pair decoder therefore does not supply an inverse for all scalar targets.

[Physical macro, decoder, collision and sparsity proof](gates/REPORT.md) · [Independent audit](audit/GATES_AUDIT.md).

## What the merge accomplishes

The paired single-arc constructor cancels the same nonlinear f(M) term in two actual path counts. This gives a simple target-gap selector even though each absolute count retains that curvature. The direct Hamilton-path injection identifies the surplus paths, so the result is not only a subtraction of evaluator outputs.

It cannot be converted into absolute-count coverage by simply adding one to the surplus: N=3 gives strong counts 27 and 33, a gap of 6, while 7 is globally forbidden (THM-1370, `h-spectrum-omits-7-21-all-n`). An exception-aware physical assembly would need a new proof. There is also no global order-optimality claim: the five-vertex outgoing masks [20,1,2,7,14] and [20,9,2,5,14] are strong, differ only on arc1–3, and have counts 9 and 13. Our adjacent-chain construction uses six vertices for that gap 4. The portable replay checks this explicit boundary hostile.

The concrete next target is a physical way to use the shared surplus while controlling the baseline and all continuation coordinates. The fixed-apex inverse and capped edit observer can test such proposals exactly. The seam macro demonstrates why a programmable, uniquely decodable pair language is still insufficient when its scalar support is sparse. Universal simulation, finite computation and all-target arithmetic coverage remain distinct.

## Verified publication and final replay

[Checkpoint 1](https://github.com/eliottcassidy2000/math/commit/1546ced59d4482bd491fe575906a8032e1a03817) published the single-arc theorem, positive ribbon specialization and provenance. All 32 files were read back exactly at that commit; the full packet and documentation checks passed.

Checkpoint 2 adds the audited fixed-apex and seam results, the independently checked path injection, complete portable replay and this synthesis. It preserves the previously read root-change/Carmichael packet and all other concurrent main history. The current file's immutable commit URL is the final index.

Run `python3 verify_all.py`. It verifies hashes, three author and three independent bounded compiler fixed points, all three physical/inverse/macro verifiers in normal and optimized modes, and five independent audit scripts. Only Python's standard library is required. All generated audit receipts use a temporary directory. Every local theorem, cited dependency, bounded computation and unresolved coverage claim is separately scoped in its source report.
