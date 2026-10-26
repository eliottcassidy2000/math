# Proof scope, independent reviews and incoming work

2026-10-26 UTC. This is the audit for [README.md](README.md), not a claim
that bounded testing settles any listed open problem. No global theorem
identifier was reserved or promoted during this session.

## Status and independent proof review

| Item | Status / evidence | Boundary |
|---|---|---|
| Root-27 to root-7 ballot equivalence | PROVED elementary arithmetic; [independent guard/source audit](anchor/INDEPENDENT_REVIEW.md) | Identification with the inherited native machine is CITED; macro counts are not native transition clocks |
| First-mismatch leaf formula | PROVED for every stated word, leaf and positive source; same audit | Only the canonical frontier guarantees that the first mismatch is a leaf |
| Fixed pumped-family criterion | PROVED for every exponent of every fixed admissible block | The selected bank of 344 blocks is FINITE-EXACT, not the universe of all blocks |
| Actual root-27 rank hostile | FINITE-EXACT literal reconstruction of odd return 19 | Refutes this charge, not all possible ranks |
| Normalized five-state macro forest | PROVED elementary normalization, unique inverse, cylinder and periodic-address facts; [independent review](atlas/INDEPENDENT_REVIEW.md) | All positive-integer survivors are not excluded by measure zero or lack of periodic addresses |
| Michel macro durations | CITED primary formula; 93 independent native branch/halting checks | The finite native controls do not replace the cited all-input configuration proof |
| Bouncer and two Bell block formulas | PROVED native tape inductions; [independent dense-tape review](machines/INDEPENDENT_REVIEW.md) | No unproved permanent chaotic classification is inferred for irregular tables |
| Goldbach25 equivalence | CITED external Lean theorem inspected in both directions; ten native checkpoints | Not locally built; no full axiom/import audit; no minimal-state claim |
| Program serialization and padding | PROVED elementary encoding/induction; [independent review](compiler/INDEPENDENT_REVIEW.md) | Marked path and decoder retained; no semantics assigned to arbitrary graph moves |
| Configuration tournament rewrite | PROVED exact canonical coloured-carrier conjugacy; machine review | Unbounded tape, origin and program remain in the data |
| Reduction hierarchy | Standard computability results and elementary deductions; machine review | Presentations are not minimality claims, and common targets do not imply mutual reductions |
| Legendre extension theorem | CITED known theorem, independently rederived CRT/Dirichlet proof; [sources and hostile](anchor/LEGENDRE_SIDECAR.md) | No assertion of increasing-prime prefix rates or a faster A000568 algorithm |

The primary author of each of the anchor, atlas and machine lanes was
different from its named independent reviewer. Reviews used independent
implementations for the decisive native/tape/guard controls. This is a
human-readable proof audit with exact experiments, not a new formal proof
in Lean or another proof assistant.

The root also checked the leaf-update algebra, both valuation cases in
the pump test, the Bell interval/count formulas, and the prime-extension
reciprocity sign. Two anchor wording repairs were made before acceptance:
pointwise `V(z)>z` is not order preservation, and z=1,2 are outside the
cited native parameter sector. Both are explicit in the review. The
external Bell correction and pointwise-growth mechanism are recorded in
[MISTAKES.md](../../../01-canon/MISTAKES.md) below its maintained header.

## Reproduction universes

Run from this directory:

```
python -B replay.py
python -B verify.py
python -O -B verify.py
```

`replay.py` runs the following nine jobs in both normal and optimized
Python. It compares normalized JSON to saved receipts and checks normal
equals optimized. Only its explicit `--refresh` mode changes receipts.
All mathematical gates use exceptions, not removable `assert` statements.

| Producer/reviewer | Explicit finite universe or decisive controls |
|---|---|
| anchor/check.py | 567 forest lifts, 510 cylinders, 3,928 leaf updates, 8,000 guarded returns, 8,600 pump controls, 32 frontier comparisons |
| anchor/review_check.py | 23,023 original-parameter macro comparisons; 2,304 guards; 62 section edges; 1,270 mismatch and 3,845 leaf cases; 27,264 independently generated pump cases |
| anchor/LEGENDRE_SIDECAR.md | All eight prescribed extensions of (3,7,11); sign and prime-1-mod-4 hostiles |
| atlas/check.py | Supplied orbit and full clock sum; complete survivor cylinders through depth 10; predecessor controls through 20,000; 1,022 periodic words; ternary digit-height and divisor-flow controls |
| atlas/INDEPENDENT_REVIEW.md | No-import five-state native simulator on 93 branches/halts; 6,144 divisor valuation comparisons |
| machines/check.py | 194 native blocks, 6,640 configuration rewrite cases and ten Goldbach checkpoints |
| machines/INDEPENDENT_REVIEW.md | Independent dense tape: 32 exponential blocks, 560 intermediate invariants and seven checkpoints through step 48,012 |
| compiler/check.py | All 20,800 complete one-/two-state tables; all 40,320 frame relabellings; opposite-outcome one-bit pair; 262 padding trace controls |
| compiler/INDEPENDENT_REVIEW.md | 4,040 field-edit commuting squares and 24 frame/complement controls |

`verify.py` pins packet files and six inherited mathematical sources with
SHA-256 after CRLF-to-LF normalization. It checks local Markdown links,
Python parsing and all nine saved receipt pairs. Its explicit `--refresh`
mode writes the manifest. A digest confirms identity of the audited
text or receipt; it is not a substitute for the mathematical proof.

## Incoming work read during this session

The live shared branch was fetched during the research. These are precise
connections to incoming work, not a claim that every new file was audited:

* Commit `fcc112a69d`,
  [twist-conductor INCOMING](../twist_conductor_reverse_20261025/INCOMING.md),
  gives a guarded affine-ray criterion. A fixed native word may be pumped
  indefinitely only if its ray stays in every guard half-space; finite
  observed repeats alone are insufficient. Our fixed-tree pump criterion
  is a different exact certification of a chosen repeating shape. Both
  preserve source and guards; neither implies coverage of every shape.
* Commit `a0645dc58c`,
  [quadratic-residue tournament](../analogy_reduction_atlas_20261025/QR_TOURNAMENT.md),
  motivated the Legendre sidecar. The known extension theorem proves
  structural universality, while the packet's numerical-prefix hypothesis
  remains a separate obligation. The shared typed-analogy ledger agrees
  with our distinction between encoding, operation and proof predicate.
* Commit `1b12f4e513`,
  [S56 preregistration](../bb6_cryptids_odd_orbits_20261025/PREREGISTRATION.md),
  records prospective deep-run and basin/hazard tests. A preregistration
  is not an experimental outcome. No unrun deep census is imported as
  support for the root-27 inequality or blank nonhalting here.

The old reflection routes in [HISTORY_ATLAS.md](atlas/HISTORY_ATLAS.md)
were used to find neglected objects. Their current repaired sources,
rather than their historical prose, control all claims. The supersingular
4+5+6 lists and Fibonacci/heap indexing were checked directly against
their maintained sources. All thirty atlas document links resolved in
the bounded source audit.

## What remains unresolved

The root-27 inequality is still equivalent to preventing a fixed source's
matching precision from reaching a finite tree's length. Aperiodic shapes
are not controlled by the fixed-pump certificates. The map here is not
ordinary Collatz, and no transfer between their barrier conditions has
been proved. The ternary survivor set's dimension, the Haar/Kraft value
phi^-6, and the five-state sequence are useful exact comparison objects,
not pointwise convergence arguments.

No generic undecidability theorem for generalized affine machines is
claimed to apply to the one classical Collatz map. No 25-state number,
Platonic count, moonshine prime list or matching tournament statistic
replaces an explicit reduction preserving the designated source and
remaining proof obligation.
