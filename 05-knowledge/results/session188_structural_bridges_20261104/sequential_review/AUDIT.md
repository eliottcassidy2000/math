# Independent review of the sequential fanout budget

**Verdict: PASS. No discrepancy found.**

I reviewed `circuits/SEQUENTIAL_BUDGET.md` and `verify_sequential.py`. Producer normal and optimized outputs agree and match `sequential_receipt.json`.

## Potential monotonicity

For a gate-parent alias g->p with positive gate-child degrees D=d(p), d=d(g), rerouting gives d'(p)=D-1+d-t, where t counts already-shared children. Since the new degree is at least one, p's potential change is d-1-t. Removing g contributes -(d-1), giving the claimed combined change -t. Every other gate parent loses a child; the drop is one exactly when its old distinct-child degree was at least two. Repeated pins must be collapsed, and the proof does so.

If the predecessor is external, its added fanout is excluded from the gate potential. If g has no gate child, deleting it can only decrease parent degrees. Constant deletion and dead cleanup only remove gate edges. These cases establish the stated monotonicity under the explicitly bounded, no-new-gate operation set.

## Sequential constant-region charge

Deleting the certified constant region in topological order preserves every such gate's original outgoing child set until its own deletion: its C-children occur later, and boundary reroutes are deferred. Each nonoutput live C gate therefore consumes exactly its previously measured d-1 contribution. Other changes cannot increase the potential. This proves the stepwise E_j bound and telescoping inequality.

For a live single-output circuit, the output has zero gate children and every other gate has at least one. Hence Phi=e_gate_gate-(s-1), and the optional Phi<=s-n_used+1 bound follows by counting distinct input incidences. This assumes the ordinary gate/input coordinates explicitly stated in the note; it does not make affine coordinate changes free.

## Counterexample scope

In the descending alias ladder, G_i genuinely has two current children L,R at each step. Its surplus transfers to G_(i-1), so the same unit of potential can be charged k times if one mistakenly charges unary-boundary fanout. Exactly two gates are deleted at each of those steps. The final constant deletion consumes the single unit and removes three gates, leaving y AND z.

This is a correct counterexample to an unrestricted boundary-charge rule for these literal simplification operations. It does not refute a claim restricted to globally minimum circuits or supply a new lower bound. The packet states this limitation accurately.

## Independent controls

`independent.py` does not import producer functions. Across 19,866 DAG presentations with zero, one, or two pins through four gates, including repeated pins and dead vertices, it checks 97,098 literal reroutes, 78,426 constant deletions, and 7,680 exact gate-parent transfer identities. It separately constructs and simplifies the sharp ladder for every k through 128. All results match the proof, and normal/optimized independent outputs agree.

Artifacts: `normal.json`, `optimized.json`, `independent.py`, `independent.json`, and `independent.optimized.json`.
