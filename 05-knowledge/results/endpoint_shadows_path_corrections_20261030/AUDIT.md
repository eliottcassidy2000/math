# Scope and reproducibility audit

2026-10-30. Same-agent proof review and independent computational paths;
**no independent-agent audit, formal proof, or Lean build is claimed**.

## Mathematical gates

**Collatz.** The prefix is an exact odd-valuation word, not an integral
affine endpoint alone. P>Q makes the anchor negative; b/D is reduced.
The clean h_j range requires M>b, and the mass/moment estimates require
M>=2b. Sources retain w(n)=4/((n+1)(n+3)); endpoint weights are transported,
not re-evaluated. The j-block integer clock is exactly dj. The proof of the
moment boundary uses the whole infinite source progression; the supplied
finite checks verify rational enclosures, not an infinite sum by truncation.
The phase-wise limits are proved by dominated convergence. No irrationality
claim for Hurwitz zeta is made. Fixed-word constants are not uniform over
all words, nor is a post-mismatch ROOT certificate constructed.
The source-dependent changing-word lower bound (C10) follows directly
from the already checked positive integer phase identity Dn+b=hM; no
unproved growth estimate for a particular infinite orbit is imported.

**Tournament.** Odd m ensures the common path edges agree in direction and
both ports coincide. Every unspecified-edge completion gives two legal
moves; blackness is not assumed. The rank is that of K2-K0, not K0, and
its grid principal submatrix suffices. Two computations (diagonal gauge
and finite-field row reduction) support the algebraic rank argument.
Random completion preserves the SAME free coordinates in all three states.
The module union bound does not assume independence of the three prime
events. For the general r-path compiler, repeated edge prescriptions are checked
for consistency before any completion is permitted, and N>=max(64,32r)
is retained. Its module failure sum is less than1/64; the norm failure is below
2exp(-8N). The existence theorem is therefore joint, rather than two
separate existence claims. The near-ramified overwrite variant is NOT
asserted prime or multiplier invariant.

**Quotient.** The edge-parity carrier is exact with a labelled base. An
Aut(base) orbit is not necessarily an ordinary target isomorphism class.
The construction excludes bounded-rank and bounded-treewidth correction
carriers, not all compact symbolic grammars. Short serpentine formulas
already encode these particular large grids. No new A000568 term follows.

## Finite universes and controls

`check.py` performs73,281 always-active gates:

-45 distinct expanding words: all words of length1..5 over valuations1..3
  with expanding full slope, plus the stated signed-cycle controls.
  At repetition depths1..8, four sources per cylinder are replayed literally:
  1,440 source instances with exact valuations and endpoint identities.
  Rational partial sums with an integral tail bound verify both the original
  mass estimate and transported weak-L1 bound at four tail cutoffs.
- Grids of sides3,5,7,9,11,17,25: legal P/Q, shared orientations, both endpoint
  ports, exact scores,2(m-1)² toggles, Eulerian support, diagonal gauge and
  rank over fields3 and5. The m=3 zero-rank interior is retained as a boundary
  case, not silently omitted.
- Every labelled tournament of orders3 and4, and all its legal two-move
  histories:792 histories check the parity-carrier decoder and converse.
- All858 proposed path programs: depths1..3 on3 vertices and1..2 on4
  vertices. Exhausting all labelled source tournaments independently matches
  the compiler's consistency decision and its exact2^(free edges) completion
  count. The immediate reversed-path conflict and the four-move return
  P,Q,Q,P are retained. Six size ranges check the general union-bound budget.
- Explicit prime triples at81 and121 vertices, seeded281 and321, both
  accepted at the first candidate. Pair-generated module closure certifies
  primality; that routine is independently compared with exhaustive proper
  subset enumeration on all labelled tournaments at orders3 and4. Exact
  Gram row-sum bounds770 and1328 certify the source norms. Deterministic
  source SHA256 values are in the receipt; no floating eigensolver is used.
- All8,192 binary colourings of13 vertices under the39 explicitly listed
  affine permutations yield7,800 with trivial stabilizer, hence200 orbits.
  An independent polynomial convolution gives the physical Gamma-bank
  coefficients, including200 twice and total4,319 nonempty states.

Normal and optimized Python runs are required to have identical receipts.
The code raises explicit exceptions rather than using disabled `assert`s.
These finite controls support the proofs; they do not expand the all-order
statements' declared hypotheses.

## Historical evidence repair

While tracing the old Q≈200 Collatz ridge cutoff, the actual producer
`04-computation/experiments/collatz_h1_20260929_mu1_track.py` was inspected.
It computes with NumPy complex128, floating phases and a finite valuation
cutoff. The corresponding theorem/result used FINITE-EXACT for agreement
to4e-14 and for observed Fourier rates. That status was unsupported.

THM-4519's canonical status, the result's maintained warning and the index
now separate VERIFIED numerical Fourier data from the FINITE-EXACT integer
seed inventory. Fitted limits and an all-seed cutoff near200 remain
EMPIRICAL / OPEN. The exact renewal identity and conditional implication
under H1 are unchanged. The repair does not claim that any numerical value
is false, nor that the full old paper has received an independent audit.
The failure mechanism and source paths are recorded in MISTAKES.md.

## Literature and incoming scope

All new proofs here are given in full; no external-priority claim is made.
The200sqrt(pN) construction and the order200 enumeration are inherited repo
results with file links. The discrepancy theorem behind the former, and
the external OEIS comparison behind the latter, were not re-proved or
re-run. The incoming final Session110 report explicitly distinguishes
portable exact replays from a failed consolidated Python timeout; this
packet does not turn that timeout into a successful verification.

Validation completed: normal and optimized receipts agree byte for byte; all26
local Markdown links resolve to present or tracked incoming sources. The
repository documentation checker and whitespace checks passed before the
checkpoint.
