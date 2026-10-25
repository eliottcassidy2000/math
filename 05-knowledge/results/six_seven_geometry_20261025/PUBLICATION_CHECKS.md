# Proof review, reproduction, and publication boundary

2026-10-25. The mathematical packet was produced in an isolated current
worktree; unrelated shared-checkout files were not modified. The starting
revision was82691d391b, and concurrent mainline checkpoint f7041582a8 was
read as incoming mathematical evidence. Its overlapping results are cited.
The containing Git commit and final remote ancestry check establish publication;
this prose does not certify a push that has not happened.

## Independent proof gates

| Packet | Independent review | Result |
|---|---|---|
| Root deletion, branch components, cube and divisor relations | [Fano-lane review](audits/ROOT_DELETION_AUDIT.md) | PASS |
| Fano gauge/deck exact sequence, free quotients, singularities and marking census | [Root-agent review](audits/FANO_AUDIT.md) | PASS |
| Canonical characters, double-sextic quotients and arithmetic count identity | [Root-deletion-lane review](audits/CHARACTER_AUDIT.md) | PASS |
| Resolved support reflection, actual Cremona involutions and nonlift, auxiliary family | [Root-deletion-lane review](audits/SUPPORT_REFLECTION_AUDIT.md) | PASS |

The common free-action/regular-component-action lemma was separately checked
algebraically by the parent and both relevant lanes, and on all5040 markings.
No substantive mathematical defect was found; no inherited theorem or
correction lineage needed retraction. Candidate statuses were promoted only
after the corresponding reviews. No Lean formalization is claimed.

## Executed exact replays

Run from repository root:

    python3 -B 05-knowledge/results/six_seven_geometry_20261025/replay.py

The harness runs each of seven scripts normally and with `-O`, requires exit0
and byte-identical stdout, and writes [REPLAY_RECEIPT.json](REPLAY_RECEIPT.json)
with source and output SHA-256 hashes. All seven pairs passed.

The source universes include252 seven-to-six deletions; explicit eight-conic
intersection computations over F17 and independently F41; all5040 Fano
markings; all30 subgroup choices; all5596 partial Steiner systems generated
by independent paths; canonical character computations through n12;35
root-surface point counts and all36 seven-subsets at two primes; and126
support/Cremona controls. The independent character audit also enumerates
every P6 point over F3, F5 and F7 and tests the equations directly.

Hostile controls are structural, not just extra positive cases: wrong smooth
branch genus, no global-sign descent for odd characters, all nonfree code
markings, the changed support rank/order and F3 count, the noncanonical
quotient discrepancy, the new-axis odd valuation, and all three exceptional
parameters of the auxiliary sextic family.

The support producer uses SymPy; the other producers and independent
polynomial/finite-field audit paths use the Python standard library.
Independent checkers import no producer implementation. Finite-field counts
refer to specified equations, not good reduction of complex resolutions or
the existence of distinct positive integer-square points.

## File integrity and document gate

`MANIFEST.json` records every packet file except itself, with exact byte size
and SHA-256. Run

    python3 -B 05-knowledge/results/six_seven_geometry_20261025/verify_packet.py

to verify the complete file set and bytes. Replays should leave those bytes
unchanged on this checked environment. Different Python/platform stdout newline
conventions can change replay metadata without changing the mathematics.

The bounded navigation is updated in START-HERE and CURRENT-FRONTIER and the
result is indexed. `python3 agents/check_docs.py` and `git diff --check`
both passed. The first document pass exceeded the total startup-byte budget;
shortening redundant displayed link labels in START-HERE, retaining their
destinations and all mathematical claims, repaired that gate. Existing
freshness/warning banners remain at their original positions.

The final packet leaves quotient isomorphism classification, nonfree minimal
models, integral saturation, and arithmetic lifting open. It supplies no
resolution of LRC(14) or the integer-square existence questions.
