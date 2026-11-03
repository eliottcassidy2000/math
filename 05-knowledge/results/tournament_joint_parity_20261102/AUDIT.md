# Audit and reproduction record

**PASS for the written scoped results.** All-order statements are supported
by proofs or explicitly cited theorems; the controls below are FINITE-EXACT.
No global spectrum completeness or unproved dependency is promoted.

## Independent audits

- [Class-pairing audit](audit_parity.md): independent adjacency-matrix
  canonicalization, literal path orders, and component-BFS reconstruction
  of the alternating trails. It checks the general odd-fibre lemma on all
  odd-fibre partitions of sets of sizes2,4,6,8,10. Closed components and
  parallel edges are retained; no representative lies in a cycle.
- [Growth audit](audit_growth.md): independent named joins, all literal
  six-vertex orders for64 children, explicit cover cuts and vertex transport.
  It exhausts all labelled parents through order4 to establish sharpness.
- [Cross-lane audit](audit_cross.md): physical toothpick ports, exact clocks,
  reflection fixed sets, insertion/root weights, and self-converse two-tail
  formula. The parent additionally checked the latter directly using unique
  ordered SCC decomposition and symmetric singleton padding.

## Universes and replays

| Program | Exact universe / positive and hostile controls |
|---|---|
| parity_probe.py | Every ordinary class n0..6 via fixed-spine words; union of labelled orbits covers every labelled mask. Literal and DP H, odd fibres, both trail pairings, converse transport. Rational Burnside n0..40, including n0,1,2 exceptions. |
| spectrum_probe.py | Complete unfiltered class augmentation n0..6; labelled orbit masses, literal/DP H, SCC star recurrence, square-central recurrence. Strong H<=21 globally complete using the cited f(7)=25 cutoff. |
| joint_audit.py | Every labelled tournament n0..6, independently canonicalized; all F_n and converse-fixed coefficients agree with augmentation. One-birth support hostile with literal children. |
| toothpick_probe.py | One classical seed through integer stage128; exact endpoints/midpoints, transported parents and birth times. Independent dyadic recurrence; stage3 blocked ports. Not an arbitrary-seed census. |
| audit_growth.py | Separate named-join and literal-order implementation; exhaustive H-only sharpness through order4. |
| audit_parity.py | Independent graph-component pairing and matrix canonicalization, plus exhaustive small odd-fibre partitions. |

The four producer programs were replayed by the parent under `python3 -I -B`
and `python3 -I -B -O` (joint_audit's already saved normal result was compared
with the optimized replay). JSON/text outputs match byte for byte. Both
independent audit programs were also run normally and optimized by their
auditors. Assertions use explicit exceptions and survive optimization.
The published receipts and source hashes are included in this packet.

`python3 agents/check_docs.py` passed the bounded-startup, warning-position,
frontmatter-status, correction-sentinel and router-smoke checks. Both PNG
diagrams were inspected visually by their producers and by the parent;
no clipped labels or missing fixed-point annotations were found.

## Corrections and scope checks

Historical TANGENTS T1544 had a missing cumulative `-1` and an unrelated
THM-893 citation. Both are repaired and recorded immediately below the
maintained banner in `01-canon/MISTAKES.md`. A draft replacement sentence
about the path cube was narrowed during review: the fixed-spine cube does
surject onto ordinary classes, with fibres H/Aut; it simply does not carry
the ambient full S_n action. No valid cube quotient was withdrawn.

Inherited parity/valuation, SCC series, the two gaps, prior finite coverage
and response-kernel corrections are credited separately. The new numerical
experiments do not upgrade them to larger universes. Source weights in all
attachment tests are the original labelled-cut weights. The reflection
pairing is explicit and choice-dependent; H is preserved on non-self-converse
converse pairs, but not required on the self-converse locus. The unique-H1
obstruction precludes any fully H-preserving free class pairing.

The incoming six-checkpoint Session149 packet was reconciled in
[INCOMING_RECONCILIATION.md](INCOMING_RECONCILIATION.md) before rebase/push.
No claim of literature priority is made.
