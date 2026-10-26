# Proof, scope, implementation and integration audit

2026-10-26. **VERIFIED independent proof reviews and exact finite replays.**
All-parameter claims rest on the written proofs. Finite gates test their
implementations, boundary examples and imported data; they are not substitutes
for those proofs. None of the checkers imports a producer's proof as code.

## Proof responsibilities and independent paths

| Lane | Producer | Independent review and computation | Outcome |
|---|---:|---|---|
| Aliquot fibres / unique6 | 19,801 gates | [Chain agent's audit](chains/AUDIT_ALIQUOT.md), Euclid–Euler parity/mod3 review and separately implemented additive divisor sieve,16,005 gates | PASS |
| Powers / doubling / abc | 268,651 gates | [Aliquot agent's audit](aliquot/AUDIT_CHAINS.md), full address and density proofs, direct nondelayed-even census at4^m for m1–7 | PASS |
| Divisor complexes / tournament / marked F8 | 8,494 gates | [Independent audit](tournament_audit/AUDIT_TOURNAMENTS.md), literal divisors rather than exponent boxes, separate feedback subset DP and full permutations, independent F8 multiplication,2,042 gates | PASS |
| Integrated microcosm / shapes / scope repairs | 32 displayed bases, complete targets2–32 | [Integration audit](root_audit/AUDIT_INTEGRATION.md), independent first-base and valuation-address checks, literal seven-vertex feedback witness, correction preimage reconstruction | PASS |

The three producer counts sum to296,946. The16,005 and2,042 independent
counts are separate. Replaying the same producer normally and under `-O`
does not double its mathematical evidence. The integration audit records its
own351 finite gates,100 cyclic-product gates, and the exact reviewed hashes.

The root read the three proof notes and independently checked the inverse
bound, unique6 congruence, operation-address domains, delayed-square census,
perfect-number diagram, primitive abc quotient, multicomplex flag boundary
and marked-field constructions. The root also enumerated all5040 orders of
the H127 witness and obtained feedback minimum5. The closed verifier replays
every supplied executable in normal and optimized Python against frozen
bytes, so test removal by optimization or dependence on another checkout
cannot silently change the result.

## Boundary gates retained

- n=1 is separate in the divisor-complex theorem: the void proper-divisor
  family lacks the empty face of a graph independence complex. The theorem
  domain is n>=2. Prime coordinates absent from the complex are ghosts,
  removed before stating a graph realization.
- The simplex boundary comes from squarefree monomial supports as faces.
  The proper-divisor order complex including1 is instead a cone. Removing
  both endpoints gives its boundary's barycentric subdivision. The cube's
  geometric boundary is a different construction.
- A Frobenius orbit can have full field degree and linearly dependent
  conjugates. Full degree suffices for an individual mark to kill the
  Frobenius stabilizer; a normal basis is additionally needed for the
  binary-word coordinate dictionary. The rank2/rank3 examples test both.
- All-order H omissions, through-eleven tc omissions and open completeness
  are kept separate. Literal cycles on equal vertex supports remain distinct.
- The supplied aliquot examples are inverse fibres.9 maps to4, not3;
  F(1) is infinite; all other listed fibres have a proved complete source
  bound. Zero is a terminal output; s(0) is not silently supplied.
- Unique6 concerns prime-square predecessors. Delayed powers are common,
  and neither a delayed row nor perfection implies that predecessor.
- The abc bound retains an OPEN conjectural antecedent, fixed epsilon and
  unspecified constant. The primitive quotient collapses every perfect
  number, and its discarded gcd is explicitly retained for reconstruction.

The integration review corrected two local wording defects before release
(the earlier base is3, not its square9; graph criterion excludes n1) and
one verifier status mismatch without changing frozen producers. Its review
also separated the full-degree marking requirement from the stronger
normal-basis requirement. These were draft repairs, not counterexamples to
the final theorems. The two inherited-source repairs have their own exact
[lineage](CORRECTION_PREIMAGES.json) and maintained ledger entries.

## Reproduction and provenance

Run `python3 -B 05-knowledge/results/aliquot_power_124_20261026/verify_packet.py`
from the repository. Python's standard library suffices. The verifier fails
on missing, altered or extra packet files, altered pinned source documents,
correction text outside the recorded before/after transition, or changed
replay bytes. It does not rewrite receipts. All JSON and source artifacts
use LF under the packet's `.gitattributes`.

The source snapshot began at61d7cca5d8. The incoming563022f893 marked-port
law was read and its elementary orbit proof checked; the concurrent
fa62e88c5f S60 reservation was recognized as RESERVED and supplies no theorem
dependency. This packet makes no new theorem-ID reservation. Primary-source
locations and literature scope are preserved in the child notes and manifests.
Sparse checkout was used to bypass unrelated unavailable promisor blobs;
every source needed for these proofs and repository document validation is
present and pinned where it supplies mathematics.

The late [S69 integration appendix](INCOMING_S69.md) was added after the
independent integration audit froze. Its elementary proofs and arithmetic
examples were separately checked by the integrating root; it is not covered
by that earlier file hash. The root manifest records this addition explicitly.
