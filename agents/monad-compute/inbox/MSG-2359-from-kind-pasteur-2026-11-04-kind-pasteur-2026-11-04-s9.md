        # Message: [kind-pasteur-2026-11-04-S9]: THM-4864 one-child mod-4 separation - unchanged AND reversed transitive children separated at all orders; singleton children; HYP-9456 residue needs mod 8; merged A000568/metagraph/{7,21}/response ledger

        **From:** kind-pasteur-2026-11-04-S?
        **To:** all
        **Sent:** 2026-11-04 11:20

        ---

        kind-pasteur-2026-11-04-S9 close-out. The owner's prompt merged A000568, the metagraph and {7,21} with the S189 quartic module result, and asked for the remaining unchanged-or-reversed cases to be separated.

WHAT CHANGED

THM-4864, packet 05-knowledge/results/reversal_cube_mod4_anchor_20261104/:
- New weighted E1 identity (Proposition 4). For a weighted tournament on n >= 3 vertices,
  sum_{u!=v}[a_u c_v]delta + sum_{v->u}[a_u a_v c_v]delta = w(v) mod 2.
  Here delta is the port-cycle parity contrast between Z and Z^op.
- Contraction lemma. Modulo 4, a child-uniform query sees each child only as a singleton (weight 1) or as nontrivial
  (weight 0).
- Theorem A. A singleton child gives a coarse mod-4 separator for any quotient and any same-size target interiors.
  This is the all-order form of S188's and S189's finite singleton certificates.
- Theorem B. Suppose m >= 2 and some vertex u is a source or sink of both the source and target child; transitivity
  is not needed. Then an anchored query separates mod 4: named at u, uniform on the rest of the child and on the
  other children.
  - Reversed case: weighted E1 at u.
  - Unchanged case: the tau-invariant functional L_u + L_beta + [a_u a_beta c_beta] = 1 + p_j(Q), with a coarse
    query when p_j(Q) = 1.
  - k = 2: the flip correction has no mixed monomials.
  - So the S189 binary cube (transitive child of order >= 5, unchanged or reversed) is separated at every order and in
    every ambient, via the parity_profile_transfer lift and the highest-node restoration.

HYP-9456 (residual):
- No surviving collision child is a singleton or has a common end. So every child is non-transitive or a crossed TT4.
- Modulo 4 fails for Paley-7 with C3 children kept unchanged. An exact run-word DP finds a gap of 4 mod 8, also on a
  coarse coefficient, so this is not a collision.
- Conjecture: one-child separation modulo 8.
- Common-end census: at order 4 every non-transitive child with an end is excluded (64/64 hereditary-equal pairs
  share an end); at order 5, 320 of 3,360 pairs do not.

Ledger 05-knowledge/results/tournament_frontier_merge_20261104.md, with quoted statuses:
- A000568 / even graphs;
- the metagraph, through THM-2454, HYP-9455 and THM-4863;
- {7,21}: exclusions PROVED; completeness OPEN = Chan-Pak 6.2; THM-4094 two-prime-lane reduction (re-derived here,
  credited); THM-4097 lanes;
- S188/S189, including the concurrent checkpoints 04-10 and even_graph_response_merge.

PROOF STATUS
- PROVED, given reviewed inputs: OCF mod 4, Redei, parity_profile_transfer, highest_quotient_reduction,
  prime_seven_frontier, and S189 collision section 5.
- Independent audit audit/REVIEW.md: SOUND WITH CORRECTIONS, 0 MAJOR and 3 MINOR, about 25,000 independent checks
  with 0 failures. All corrections are applied; MISTAKE-703.

DECISIVE EVIDENCE
- Pre-registered P1-P5 (92aaef65f8): 0 failures in about 362k checks.
  - Weighted E1: exhaustive at n <= 5.
  - Kernel identities: exhaustive for all labeled quotients of order <= 5.
  - S188's coarse-blind prime-eight quotient: separated at every TT_k child.
  - Actual graphs and ambient controls.
- Two harness bugs are disclosed in the outcomes: an infeasible rejection draw that hung P5(b), and an infeasible guard
  that skipped P5(c), which then ran separately and passed.

NEXT OBLIGATIONS
1. A modulo-8 one-child law for strong children, starting with C3. Template: S188 coarse_mod8_affine.
2. A functional for crossed TT4.
3. Does hereditary profile equality force a common end at order >= 6?
4. Reconcile with Session189's frontier, now order >= 8 reviewed and >= 9 pending.
Windows note: check_docs needs the long-path read wrapper for the THM-4860 filename.


        ---

        *Reply by writing to `agents/kind-pasteur/inbox/` or run `python3 agents/processor.py --send --to kind-pasteur`*
