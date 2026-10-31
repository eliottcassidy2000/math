# Scope correction: cyclic products do not exhaust the sewn graph

**Status: CORRECTION / PROVED explicit counterexample / FINITE-EXACT exhaustive check.** This is the latest overnight packet and overrides any graph-complete reading of packet j's cyclic product formula. Packet j's tree theorem and packet k's finite105-state matrix certificate are unchanged. Their frozen bytes are retained for provenance.

In packet j section2, W_v(P_v) sums over local **Hamiltonian cycles** containing the removed ports. The connected-matching formula counts exactly the sewn products of those local cycles. It does **not** in general count every Hamiltonian cycle of the sewn union graph. Restoring a cyclic assembly's cut ports may produce several local cycles rather than one. On a tree the two-edge-cut induction excludes this problem; on a circuit it does not.

Here is a small explicit graph hostile. Source A is the triangular prism on0..5:

    edges01,12,02,34,45,35,03,14,25; ports01,34.

Source B is the four-cycle on6..9, with edges67,78,89,69 and ports67,89. Remove all four ports and insert seams06,17,38,49. There are exactly three Hamiltonian cycles in the resulting ten-vertex graph. Source A has two Hamiltonian cycles containing both ports; source B has one. Their two sewn products are Hamiltonian, but a third global cycle is

    0,2,1,7,8,3,5,4,9,6,0.

Restoring A's ports in this cycle produces the two triangles012 and345. Thus its local matching is excluded by the condition that P_v union M_0 be one cycle. The final graph is simple; its seams are distinct; the obstruction is not duplicate labels or edge reuse.

**Repaired complete representation.** Use all spanning linear-forest boundary states, including those whose reclosure by the original ports has several components. This is precisely the stronger named-port transfer already proved in Session32 `ports/PORT_TRANSFER.md`, recovered and cited in packet j. The local cyclic-matching domain remains useful for an explicitly selected product language, but cannot replace that full transfer without an additional reclosure theorem.

The example is graph-theoretic. No integer square-sum realization of these particular labels is claimed. It suffices to refute an unconditional graph-level completeness inference. The48-of105 map and pointwise finite kernel concern the selected single-cycle boundary domain and stay valid in that scope.

Reproduce with `python3 05-knowledge/results/square_sum_spectrum_overnight_20261031l/audit_scope.py`; optimized Python gives identical output. The checker enumerates every Hamiltonian cycle of all three explicit graphs and distinguishes the two products from the extra cycle.

## Final handoff

Use this correction together with [packet k's handoff](../square_sum_spectrum_overnight_20261031k/HANDOFF.md) and [packet j's broader handoff](../square_sum_spectrum_overnight_20261031j/HANDOFF.md). The decisive open cyclic target now has two separate obligations: realize the arithmetic seams and retain all required path-cover states, rather than only states obtained by cutting one local Hamiltonian cycle. Interval absorption and universal raw tournament-H completeness remain OPEN.

Cutoff:2026-10-31 13:00UTC =07:00 America/Denver. At/after it, stop research, finish Git close-out and pause `overnight-square-sum-and-tournament-spectrum-research`. No new chats or subagents are authorized.
