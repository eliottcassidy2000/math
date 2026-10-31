# Session 122: sharp periods, five arbitrary tails, and the E14 closure

2026-10-31. **Scoped proofs, exact finite certificates, and independent audits.**
Baseline: `d345f993475dd7b5272ac6c3b4ab693899515a38`.

## Results

1. **[E14 finite-content closure](collatz/E14_CLOSURE.md).** The remaining
   finite-depth Collatz price is paid by a complete near-tight clock analysis.
   [Independent audit](audit_collatz/E14_AUDIT.md),
   [reproduction](collatz/REPRODUCE.md),
   [final proof review](final_e14_review.md).
2. **[Five arbitrary added speeds](lrc/FIVE_TAILS.md).** All1287 eight-element
   cores in {1,...,13} admit any five positive integer attachments at radius1/14.
   This includes the new [715 four-tail families](lrc/LRC.md).
   Independent audits: [four tails](audit_lrc/AUDIT.md),
   [five tails](audit_lrc5/AUDIT.md).
3. **[Uniform exact dyadic periods](ramsey/DYADIC_PERIOD.md)** for the
   restricted C5/join/union graph class, with a recursive sharpness criterion.
   **[Exact rows5..8](ramsey/ROWS_5_TO_8.md)** hold from budget1.
   **[Finite generators](ramsey/FINITE_GENERATORS.md)** give a concrete
   all-row threshold algorithm. Audits: [periods/rows](audit_ramsey/AUDIT.md),
   [generators](audit_ramsey/FINITE_GENERATORS_AUDIT.md).
4. **[Improved interval conductors](conductor/INTERVAL_CONDUCTOR.md)** include
   a sharp interval-only constant in the small-radius range and sufficient
   conductors through six tails at radius1/14.
   [Independent audit](audit_conductor/AUDIT.md).

**Not claimed:** Collatz convergence or all-depth content control,
unrestricted LRC(14), unrestricted Ramsey improvements, a general optimized
capacity compiler, Lean formalization, human refereeing, or priority.

## Reproduce and inspect

Run `python3 -B verify.py` for both normal and optimized declared Python
pipelines, including complete regenerated four/five-tail trees and independent
reconstruction. The C++ large-potential and corridor commands are separate in
[REPRODUCE.md](collatz/REPRODUCE.md). The inherited rational Collatz certificate
and the preceding Session121 LRC certificate are stated dependencies rather
than silently counted as newly rerun results.

Compact JSON receipts are distributed as gzip/base64 text archives. Run
`python3 -B materialize_receipts.py` to recover their exact original JSON
bytes with length/hash validation. The master verifier decodes in memory.
Both the1.96GB potential and larger full LRC trees can be regenerated; none
is required as a giant repository blob.

[Session synthesis and remaining frontier](SYNTHESIS.md).

## Inheritance and concept board

* Anchor: complete closed core intervals, exact widths, original clocks and
  joint attachment. Hostile: {1,2} at1/4 admits3 and12 separately but not jointly.
* Niche: sharp C5 capacity and actual child realizability. Hostile: integral
  slope w(3)=4 does not permit sharpness outside multiples of4.
* Native-price lane: retain ordered words, first returns, translation and
  clock. The abstract score30 endpoint is not an actual identity word.
* Wildcard: exact residual-arc measure plus connectedness at equality.
  Phase-independent sufficiency is not full-core characterization.

Method cards: search the statement before the method; retain observer,
recurrence class and finite head; correct the object before sharpening the
technique. See the synthesis for the distinct maps and information losses.
