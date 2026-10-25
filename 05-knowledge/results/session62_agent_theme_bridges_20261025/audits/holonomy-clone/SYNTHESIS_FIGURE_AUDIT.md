# Final synthesis and figure readback

2026-10-25. **PASS after summary-only qualification repairs.** This closeout
reviews `SYNTHESIS.md` sections 2, 4 and 5, together with
`figures/common-action-map.svg`, its PNG rendering and its generator.
Sections 1 and 3 belong to the other independent lane audits; hashing the
whole synthesis does not imply this auditor independently reviewed them.

## Mathematical readback

- Section 2 faithfully retains odd group order for the subset-root theorem,
  distinguishes exact fixed-XOR rank from vanishing order, states the
  simultaneous-element criterion, and accurately summarizes the 24-vertex,
  33-vertex and unbounded tournament constructions. The separate real
  exterior theorem does not acquire an unnecessary odd-order hypothesis.
  Negative coefficients after removing all three root factors correctly
  obstruct a positive degree-one exterior-factor interpretation, not
  unrelated nonlinear operations on an abstract set.
- Section 4 correctly states the equivariant parity splitting and the
  holonomy-dependent dimensions. The one-cycle/two-independent-cycle
  distinction, marked conic-label comparison, section versus trivialization
  boundary and nonidentification with the palette automaton are preserved.
- Section 5 now explicitly requires a **connected** graph with **holonomy
  image A**. Its regular nine-gap example explicitly uses **F2** and two
  cycle holonomies generating the whole group. Consequently its four
  quantities are exactly `r=1`, `d=3`, `dim H^1(V)=10`, and
  `dim H^1(ker parity)=8`. The graph remains additional data rather than
  an asserted intrinsic tournament-move graph.

The only changes requested during this summary audit were the connectedness,
full-holonomy and coefficient-field qualifications just recorded. They were
made in summary/figure files; no frozen source proof required correction.

## Figure and reproducibility

The SVG visibly states odd-order A, distinguishes individual and common
fixed spaces, labels the lower branch "Connected, holonomy A", and labels
the numerical example as connected gluing over F2. The formulas and the
four numerical values agree with the independently audited proof.
The text explicitly says the gluing graph is additional data and that no
plane-coloring bound follows.

An earlier PNG still had the pre-qualification lower-branch label. The
parent regenerated both formats. This auditor then copied the final
unchanged generator into its own scratch directory and ran it with `--png`.
The regenerated SVG and PNG were **byte-identical** to the parent's final
files. Both PNG versions were actually inspected; the final image has
legible equations and no detected overlap, clipping or title/body collision.
No image editing or theorem-source modification was performed by this audit.

Matplotlib reported unavailable default home cache directories and used its
normal temporary-cache fallback. Rendering exited successfully and the
byte comparisons passed; this was not a validation blocker or a permission
escalation.

Final reviewed SHA-256 pins are in `synthesis-figure-manifest.json`.
The prior theorem audits, source pins and 123,237-gate receipts remain
unchanged. No unresolved issue remains in the reviewed summary sections
or diagram.
