# Audit of precision resets, marked units and witness transport

**INDEPENDENTLY AUDITED as itemized / FINITE-EXACT reproduction /
Collatz OPEN, 2026-09-26.** The symbolic results are elementary proved
claims in the linked notes, not formalized proof-assistant theorems.

## 1. Proof and scope review

The root researcher independently reviewed the
[product note](reset_20260926_inequality.md): both noncollision branches
have the same product law; integer balance gives the stated sharp3/4
reset constant; the Mersenne first-collision word and endpoint are exact;
the unbounded imbalance example is typed as an arbitrary admissible
representation rather than a canonical-reachability theorem; and the
logarithmic-role obstruction retains that same domain restriction.

A separate researcher independently reviewed all of the
[clock note](reset_20260926_clock.md), including the different P5 clock,
the exceptional small source y=5, both sharp constants and families,
the actual177->133->25->19->29 hostile, telescoping through swaps, and
the ternary run-length induction. No repair was needed.

The [colour note](reset_20260926_colours.md) received an independent
proof review and additional exhaustive controls: all binary-mask fibres
through143;523,776 triples enumerated directly by(q,r), independently of
the producer's even-summand bijection; and the complete F10=55 expansion
language with8192 words. The two/three-fibre classification, exactly two
charge-preserving sheets, legal35-blue row, forced-red35 in every36-row,
and restricted floor(y/3) count passed.

One wording repair was made before publication. A pigeonhole obstruction
for preserving EVERY decomposition was initially stated too broadly as
a requirement on any successful inequality. The repaired statement
allows a canonical selector, an arbitrary potential of y alone, or a
quotient preserving descent witnesses without preserving every split.
No numerical result or proved lemma changed.

The root researcher reviewed the [logic note](reset_20260926_logic.md)
and directly checked the primary HF1/HF2/HF3 source. The pure-membership
reformulation derives empty existence from the false-predicate induction
instance. The reverse reflection proof separately derives Pairing, Union,
Power Set, Replacement and Infinity; its existential witness and relation
are reflected together. Choice remains explicit. Finite subset-closed
containers are not confused with power-set-object closure. Classical
HF/PA equivalence and ZF reflection remain CITED, not experimentally proved.

The root researcher independently derived and reviewed the
[swap-lift inequalities](reset_20260926_swaplift.md), including the final
collision R=A and clipped final step P<R<A. The proof compares the actual
lifted endpoint with the ORIGINAL pre-swap source. It never assumes
convergence for arbitrary3q. The frozen171-row bank and the1431 finite
small-source exceptions yield the65 full residue classes. Their density
does not imply that an uncovered orbit must enter them.

A separate researcher independently checked the swap-lift proofs and
ran1,824 first/later core-crossing prefixes for odd q<=1023, with34,060
endpoint checks including101-digit coefficients. The boundary P<R is
essential: q=1,J=2,A=5,P=R=1,v=1 gives actual U^3(3)=1, while incorrectly
extending the clipped formula would give25. The theorem retains the
strict boundary. The final family containing27 was additionally checked
by direct iteration in the root audit after its addition.

## 2. Reproduction and independent arithmetic paths

Run

    python -X utf8 -B 04-computation/experiments/reset_20260926_audit.py

[The audit script](../../04-computation/experiments/reset_20260926_audit.py)
runs each of the five producers normally and with Python optimization,
and compares both outputs with the frozen file. Checks use explicit
exceptions and survive `-O`. Newline and optional UTF-8 BOM normalization
are explicit; no mathematical output is discarded.

The audit then independently uses direct odd Collatz iteration, without
the controller implementation, to verify64 Mersenne episodes and the
P5 hostile. It parses every row of the frozen swap table and reconstructs
its core chronology and affine carry. For each row it replays four lifted
sources, including one with more than512 bits. A separately constructed
low-bit trie recovers the exact union density, independently of the
producer's nested-cylinder removal algorithm.

The [audit output](reset_20260926_audit.out) records agreement of all five
normal/optimized/frozen outputs, all171 independently reconstructed core
certificates,684 direct lifted-source checks and the independent trie
measure. Finite tests audit the implementations and explicit certificates;
the all-height conclusions depend on the displayed algebraic proofs.

## 3. Integration and limits

Incoming correction8287bf6b6 was inspected before integration. Its repaired
Gilbreath wall theorem motivates the invariant-domain question in the
[board](reset_20260926_board.md), without supplying a Collatz implication.
No scarce theorem IDs were reserved and no old theorem was promoted from
RESERVED. Files are named by this session prefix. A separate SHA-256
manifest pins the final artifacts with LF-normalized text.

Still OPEN: a global budget for valuation-one growth and unbalanced resets;
coverage of arbitrary coefficients and insufficient precision; a rule
continuing the supplied coloured diagonal sequence; and Collatz itself.
