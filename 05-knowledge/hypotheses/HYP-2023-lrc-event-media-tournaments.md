---
id: HYP-2023
status: OPEN
source: codex-2026-06-01-S539
related:
  - HYP-2018
  - HYP-2019
  - HYP-2020
  - HYP-2021
  - HYP-2022
  - HYP-2017
  - HYP-1961
---

# HYP-2023: LRC can be encoded by event-media tournaments on holes, gates, certificates, and carries

**Claim.** The next useful tournament mappings should put vertices on the
medium of the LRC obstruction rather than on runners: empty sectors, boundary
gates, crossing events, danger-cover endpoints, carry digits, or rotor-like
vacancy states.  LRC becomes forced exhibition of an observer-clasp two-hole
class, or equivalently the impossibility of a closed kinetic/carry/certificate
orbit avoiding that class.

**Evidence.** `lrc_event_media_tournaments_s539.py` compared six exact open-cell
functors for `n=4..6`:

```text
hole_only_bare:            classes 3,4,6     mixed 2,3,3
hole_only_anchored:        classes 14,31,68  mixed 0,0,0
sector_survival_bare:      classes 2,3,10    mixed 1,1,2
sector_survival_anchored:  classes 24,61,118 mixed 0,0,0
gate_priority_bare:        classes 3,6,13    mixed 1,2,6
gate_priority_anchored:    classes 78,373,626 mixed 0,0,0
```

The law is:

```text
bare media tournaments are tiny and mixed;
observer-anchored media tournaments are pure and larger.
```

Thus target anchoring is not cosmetic.  It is the gauge that turns a small
event-media quotient into a proof object.

**Interpretation.** HYP-2022 makes sectors the nodes and identifies the DFT
duality with resonance.  This HYP goes one layer deeper: nodes can be the holes
that move through sectors, the gates whose certificates fail, the endpoints of
the danger-cover arcs, or the carry digits of the sector odometer.  The
tournament relation can encode survival, next-event priority, certificate
shadowing, or forced carry precedence.

**Predictions.**

1. Every rotation-invariant bare event-media quotient without observer anchoring
   has mixed good/bad fibers for sufficiently large `n`.
2. The anchored hole-only quotient is the smallest natural pure event-media
   target; it should grow much slower than full sector occupancy.
3. Gate-priority tournaments need present-state colors: future pressure alone
   cannot certify present loneliness.
4. A counterexample is a closed kinetic tournament orbit avoiding the
   observer-clasp two-hole class and maintaining a nonpeeling endpoint/carry SCC.
5. For `n=18`, the missing memory should live in a prime-power carry tournament,
   not in a single-level residue/parity quotient.

**Next tests.**

1. Identify the uncertified speed sets in S539 and separate boundary-only from
   genuine open-cell misses.
2. Extend hole-only and gate-priority counts through `n=7` with optimized
   canonicalization.
3. Build the event-word tournament on the first `n` boundary crossings after a
   cell and search for forbidden observer-empty factors.
4. Build an endpoint-shadow tournament for the danger-cover arcs `{B_i}` and
   test whether acyclicity peels a non-cover certificate.
5. Add a `3`-adic carry tournament for the `n=18` gate frontier.

**Files.** `04-computation/lrc_event_media_tournaments_s539.py`;
`05-knowledge/results/lrc_event_media_tournaments_s539.out`;
`07-reflections/lrc-event-media-tournaments-holes-gates-certificates-s539.md`.
