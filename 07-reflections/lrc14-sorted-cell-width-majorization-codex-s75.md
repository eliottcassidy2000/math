---
date: 2026-06-21
source: codex-2026-06-21-S75
tags: [lrc14, cell-widths, majorization, bounded-stratum, tournament-analysis, HYP-2781]
---

# LRC14 Sorted Cell-Width Majorization

The fixed-cell idea is now formally dead, but it left behind a useful quotient.
Mac-mini's Huffer-Shepp port showed reflection and symmetrisation do transfer
to individual `W_a` cells, while also showing AP is not a per-cell maximizer.
The exact S75 scan sharpens that: fixed cyclic windows are not only false, they
are the wrong invariant because dilation permutes the cells.

Sorting the six widths is the smallest repair that respects the `Z_7^*` action.
On the bounded target rows this almost gives a clean weak-majorization theorem:
k=8 has no sorted-prefix leaks; k=9 has one near-AP leak
`(0,1,2,3,4,5,6,7,9)`; k=10 has one near-AP leak
`(0,1,2,3,4,5,6,7,9,10)`; k=11 has none in this bank.  Both target leaks are
repaid by lower cells, so total `measS7` still favors AP.

This fits the current mainline correction.  HYP-2778 says universal
consec-max is false and unnecessary, and HYP-2779 says wide closure needs a
direct joint `p0` bound.  Incoming HYP-2780 gives the stronger joint-coupling
version of the sorted-W idea.  Therefore HYP-2781 should not be sold as a
universal Schur theorem; it is a bounded-stratum leak ledger for HYP-2780:
use sorted prefixes to compress most rows, then dispatch the finite one-hole
leak family with an addressed compensation inequality.

The k=12 guardrail is valuable.  In this cell-width convention the sorted
quotient starts leaking in ways that include a positive total surplus.  That is
not a problem for LRC14 after HYP-2778; it is a warning not to keep chasing the
old consec-max wall after it has stopped matching the actual cap proof.

The tournament mapping is abstract but useful.  Vertices are proof lenses:
relation marginals, one-sink compensation, sorted-cell majorization,
conductance bottleneck, affine moments, WIN/DISC, fixed windows, and per-cell
reflection.  Edges compare affine safety, ability to pin AP, leak finiteness,
formal readiness, and compatibility with the wide branch.  The ranking is
transitive, with relation marginals and one-sink compensation above sorted
majorization.  That says the sorted quotient is probably the discovery layer,
not the final formal language.

Next useful move: enumerate the leak family parametrically rather than just by
span bank.  If every sorted-prefix leak is an AP one-hole extension or a
bounded dilation of one, then the compensation lemma can be small enough to
formalize or feed to the exact THM-534/Delsarte certificate.
