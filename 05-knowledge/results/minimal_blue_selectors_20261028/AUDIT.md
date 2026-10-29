# Proof and computation audit

2026-10-28. Local proved statements refer to M1--M5 in THEORY.md. These are
written proofs with exact computational checks, not formal verification or
an independent-agent audit. Published/inherited material is separately
credited in CONNECTIONS.md. No global theorem IDs are reserved.

## Validity gates

M1 distinguishes native path complement from converse, observed blue from
trivial interval symmetry, common-polarity runs from arbitrary constant runs,
and the Boolean verdict from a full reflection syndrome. Cross-pair equations
force the common polarity. Source/sink peeling proves uniqueness of the
offset. Explicit realizations and distinguishing contexts prove the sharp
state count; a finite census alone is not used to assert all-size minimality.

M2 tests structural recognition rather than a decoder that assumes original
labels. It includes the coincident-port singleton and nonstrong cores. The
inherited two-marker decoder is credited. Four-marker optimality is OPEN.

M3 proves both directions. The finite horizon bounds number of inserted
vertices, not number of complement steps. Fresh reset rows are class-constant;
an interior-conditioned reset is outside that old menu. Minimal quotient
does not mean every retained coordinate is independently necessary as an
encoding, or that the full relation is stored in O(log n) bits.

M4 specifies the selector cost and separates marked observations, conditional
reset laws and named guards. Different eligible counts alone would not prove
different normalized laws; the H=29 event checks that difference explicitly.
Guard support/span minimum is over the5040 permutation words at the original
order7, not over all possible later attachment gadgets. All7 reconstruction
with a successor mark is a finite statement, immediately challenged by M5.

M5's padding preserves both endpoint maps, fixes queried boundary vertices,
and leaves a nonempty old interval. At the first fresh reset the entire
literal law agrees. Strongness follows from an exhibited spanning cycle.
The theorem does not apply to deleting padding or retaining a moving cursor.

The81/27 pair tests retaining actual relation/symmetry information; its
ordinary classes already differ. The C15=32/31 mismatch is an inherited
counterexample, not recomputed here. The Borwein margins are re-certified,
but the tiny integral deficits and the all-integer crossing conjecture are
not recomputed or newly proved.

## Exact finite universes

* blue_minimum.py: all33868 presentations of orders1..7, all three updates
  (101604 total), and a separating pure-context profile for each realized B.
* selector_probes.py: all33868 presentations decoded from the four-marker
  context, also after cyclic relabeling; all5040 guard words on the7-pair;
  all seven potential extra marked positions; full selected-reset partner
  distributions before/after complement; every nonsingleton cell in the
  inherited exact7-bank, refined by the successor mark; buffered strong
  controls for widths1..6 and85 contexts each.
* bridges.py: certified interval sums through m=40249 with60-digit directed
  rounding, pi bounded by a rational Machin alternating series; all ordinary
  presentation fibers through7; all21 Paley7 automorphisms and seven explicit
  anti-involutions; inherited exact81/27 automorphism/attachment algorithms
  replayed with an independent constraint-search group check, cyclic-triple
  equality and subset-DP Hamilton counts.

Shared ancestry is disclosed: the state classifier uses the old boundary
primitive, and comparison uses literal relation construction; the native
class bank and isomorphism oracle are inherited. verify.py also rebuilds
the entire rooted7-bank rather than trusting only its saved cell count.
The all-size conclusions rest on the proofs in addition to these tests.

## Reproduction

From the repository root, standard-library Python3 suffices:

```text
python -B 05-knowledge/results/minimal_blue_selectors_20261028/verify.py
python -B -O 05-knowledge/results/minimal_blue_selectors_20261028/verify.py
```

The verifier checks source/artifact seals, compares fresh results with saved
receipts, checks the selected H=29 probabilities, rebuilds the old endpoint
bank and checks its equality, and rechecks seals after the runs. It does not
rerun the larger inherited A000568 enumeration or the full C15 census.
Text hashes normalize CRLF to LF; gzip hashes use original bytes.

Individual producers can be run directly. They replace their corresponding
receipt; they do not silently update any seals. The default verifier must
reject a changed source or receipt. Scripts use active RuntimeError checks
under both normal and optimized Python. Normal and optimized run reports
are saved outside the sealed manifest because their optimization flag differs.
