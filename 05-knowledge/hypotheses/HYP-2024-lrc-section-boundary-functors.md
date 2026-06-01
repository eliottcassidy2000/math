---
id: HYP-2024
status: OPEN
source: codex-2026-06-01-S539
related:
  - HYP-1981
  - HYP-1982
  - HYP-1986
  - HYP-1999
  - HYP-2001
  - HYP-2008
  - HYP-2009
  - HYP-2018
  - HYP-2019
  - HYP-2020
  - HYP-2021
  - HYP-2022
  - HYP-2023
  - THM-382
  - THM-383
  - THM-384
  - THM-387
---

# HYP-2024: Section-boundary tournament functors challenge the runner-vertex default

**Claim.** In LRC tournament analysis, the vertex set should be treated as a
choice of proof microscope, not as an inherited obligation.  Runners are only
one possible vertex set.  Fixed circle sections, fixed section boundaries,
empty sections, wall-crossing events, cover arcs, residue channels, and even
proof obligations can also be tournament vertices, provided the quotient keeps
the LRC target visible.

The section-boundary version says:

```text
At q=n, divide the time circle into n fixed sections.
The LRC target is exhibition of a state in which the two observer-adjacent
danger sections are empty.
```

This extends HYP-2022's sector-occupancy dual mapping and HYP-2023's
event-media tournament mapping.  HYP-2022 says fixed sectors are already a
highly restricted real-space dual of the resonance picture; HYP-2023 says
holes, gates, certificates, and carries can be vertices when the right anchor
is retained.  HYP-2024 says these were not just clever alternatives, but
evidence that the vertex set itself belongs inside the quotient-stack design
space.
It is also the section-node analogue of HYP-2021's BLEX source target and
HYP-2018/HYP-2020's threshold-decorated quotient stack.  Instead of asking
"which runner tournament class is exhibited?", ask:

```text
Which fixed-section or fixed-boundary tournament classes can the arithmetic
clock exhibit, and which of those classes certify empty danger sections?
```

**Evidence.** `lrc_section_boundary_functors_s539.py` audits four families of
non-runner tournament functors on bounded primitive speed sets:

```text
section_pressure:       vertices are fixed sections; edges compare occupancy
                        count, speed flux, and distance from the observer.

section_empty_colored:  same vertices, but danger and occupied status are
                        encoded as section colors.

boundary_flux:          vertices are fixed section boundaries; edges compare
                        incoming and outgoing flux across adjacent sections.

void_pressure:          vertices are fixed sections; empty sections beat
                        occupied sections.
```

The audit uses cyclic/dihedral canonicalization rather than arbitrary
tournament isomorphism, so the fixed circular scaffold is retained.  It reports
Tournament Analysis fingerprints: class counts, pure-good and mixed fibers,
certification rates, score histograms, directed triangles, SCCs, Hamiltonian
paths, and a meta-tournament over functors.

Bounded results:

```text
n=4, max speed 10:
  section_empty_colored_q4  classes=26  pure=4   mixed=0  cert=108/109
  boundary_flux_q8          classes=376 pure=42  mixed=0  cert=108/109
  void_pressure_q8          classes=111 pure=15  mixed=0  cert=108/109

n=5, max speed 8:
  section_empty_colored_q5  classes=133 pure=15  mixed=0  cert=67/69
  boundary_flux_q10         classes=3487 pure=360 mixed=0 cert=67/69
  void_pressure_q10         classes=589 pure=91  mixed=0  cert=67/69

n=6, max speed 7:
  section_empty_colored_q6  classes=456 pure=52  mixed=0  cert=20/21
  boundary_flux_q12         classes=2710 pure=258 mixed=0 cert=20/21
  void_pressure_q12         classes=1114 pure=141 mixed=0 cert=20/21
```

The missing speed sets in these open-cell certifications are the expected
wall-only extremal behavior: the AP/regular-polygon type is lonely only on
threshold walls, so open section states do not certify it.  This mirrors the
THM-382/THM-383 compactification issue and should be handled by a boundary
layer rather than interpreted as a failure of the functor.

**Interpretation.** HYP-2022 independently confirms that evenly spaced sector
vertices are not a decorative analogy: their DFT is the real-space dual of the
resonance/inside-debt picture.  HYP-2023 independently confirms that event and
media vertices can become pure after the correct observer anchor is attached.
HYP-2024 broadens those into a protocol for choosing and auditing non-runner
tournament vertices.

The pure section-empty quotient is a direct fixed-frame
form of the LRC statement.  It remembers exactly the two facts the raw
runner-tournament quotient hides:

```text
1. where the observer danger region sits in the circle;
2. whether that region is occupied.
```

Boundary-flux is not merely a static quotient.  Its edges change when runners
enter or leave fixed section boundaries, so it behaves like a derivative layer
for HYP-2021's BLEX and HYP-2018's near-graph: it records how the dangerous
sections become empty, not only that they are empty.

Void-pressure is the strangest but most suggestive map.  It treats absence as
the active object.  In this view the LRC witness is not a runner configuration
but a moving vacuum: a protected empty sector wins a tournament against
occupied sectors.

**Assumption-challenge protocol.** For future LRC/tournament sessions, before
settling on a map, explicitly answer:

1. What are the tournament vertices: runners, gaps, sections, boundaries,
   wall-crossing events, residue classes, cover arcs, Fourier modes, matroid
   circuits, proof obligations, or something else?
2. What pairwise observable is being compared?
3. What switch or gauge turns the observable into an edge?
4. What LRC predicate is visible in the quotient?
5. What information is destroyed, and is that loss deliberate?
6. Why is this vertex choice not just inherited habit?

**More radical vertex sets.**

1. **Event tournaments:** vertices are times when a runner crosses a section
   boundary or danger endpoint; edges compare which event relieves more
   blocker pressure.
2. **Cover tournaments:** vertices are danger arcs in the time circle; edges
   compare endpoint ownership, uncovered gap contribution, or Helly defect.
3. **Residue-channel tournaments:** vertices are CRT channels or p-adic
   depth layers; edges compare exported endpoint debt.
4. **Fourier-mode tournaments:** vertices are resonance characters; edges
   compare contribution to inside-debt cancellation.
5. **Proof-obligation tournaments:** vertices are local lemmas or obstruction
   clauses; edges point from the obligation that discharges another to the one
   it controls.  A proof becomes a forced source or sink in this meta-space.

**Predictions.**

1. `section_empty_colored_qn` is the smallest fixed-section pure quotient for
   open cells; compactified wall states add exactly the THM-383 boundary menu.
2. `boundary_flux_q2n` is a derivative refinement of BLEX/two-sentinel maps:
   it separates left/right entry events at the cost of a larger label space.
3. Void-pressure classes should correlate with HYP-2008 apex/largest-gap data:
   a lonely witness is a high-rank empty-sector/apex class.
4. Event and cover tournaments will connect HYP-2001 permutohedral handoff
   circulations with HYP-2021 endpoint-owner SCCs.
5. At hard denominators such as `n=14` and `n=18`, the useful finite target is
   likely an intersection of section-empty, boundary-flux, CRT-debt, and
   endpoint-owner proof-obligation classes.

**Next tests.**

1. Add compactified wall states so AP/regular-polygon witnesses certify in
   the section-empty quotient.
2. Compare section-empty classes directly against BLEX classes on identical
   open cells.
3. Build an event tournament whose vertices are boundary crossings and check
   whether source exhibition predicts THM-387 gap-race success.
4. Build a cover-arc tournament over danger intervals and compare SCCs with
   endpoint-pressure cores from THM-379/THM-380.
5. Test whether section-empty class counts have a closed necklace/Ferrers
   formula.

**Files.** `04-computation/lrc_section_boundary_functors_s539.py`;
`05-knowledge/results/lrc_section_boundary_functors_s539.out`;
`07-reflections/lrc-section-boundary-functors-and-assumption-challenge-s539.md`.
