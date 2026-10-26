# Session 65: guarded machines and faithful tournament interfaces

**Second independently audited checkpoint — 2026-10-26.** The exact statements
below are separated from cited external results and bounded computation.
No global arithmetic conjecture, new asymptotic machine classification, or
unrestricted tournament-connectivity result is claimed.

## Audited results

### Native machines and a tournament tape observer

The [machine report](session65-cryptids/REPORT.md) identifies the spaghetti
sandwich, bouncer-and-bell, and irregular-bell tables. These three four-state
binary tables are total: their nonhalting is immediate from totality. Their
interesting long-term morphology is a separate question, not a resolved
Cryptid halting problem.

**PROVED:** all-parameter guarded macros on the stated sections preserve the
full finite-support tape, head, control and exact native clock. Encoding a bit
as TT3 or C3 and a word as their ordered join gives

    d_iso(E(u), E(v)) = Hamming(u,v).

The canonical strong-component order decodes the word even when adjacent zero
blocks merge. With window origin, head, control and program retained, native
steps become controlled marked local arc rewrites. This is neither autonomous
tournament dynamics nor the repository's K_P operation. Scalar H and triangle
counts lose word order; an unbounded native clock can accompany one changed
arc. See the [independent audit](session65-audit-cryptids/REPORT.md).

### Signed pair orbits, Frobenius defects, and the two fifteens

**PROVED:** for any finite vertex permutation sigma, let its unordered-pair
orbits have lengths L_O and endpoint-return signs epsilon_O. The complete
orientation-defect enumerator is

    sum_T z^D_sigma(T) = product_O ((1+z)^L_O + epsilon_O*(1-z)^L_O).

The minimum defect is exactly the number of even vertex cycles. An invariant
tournament exists iff sigma has odd order; permutation sign alone retains
only parity. The [general theorem](session65-fifteens/GENERAL_PERMUTATION.md)
and [independent proof](session65-audit-fifteens/AUDIT.md) give all quantifiers,
minimizer counts, involution specialization and hostile examples.

Applying the involution formula to supersingular Frobenius yields minimum
defect (s_p-r_p)/2, which equals g(X_0^+(p)) by the **CITED** arithmetic theorem.
The explicit p=37 calculation independently gives the nonzero-defect hostile.
This is not a canonical arithmetic tournament, a new proof of Ogg's complete
classification, or an identification with Sloane's exceptional exponent 15.

The [fifteens report](session65-fifteens/REPORT.md) distinguishes zero-free
Sloane powers from two-free Erdos powers. Exact finite-depth ternary lifting
retains an exponent-parity coordinate: the coherent all-ones branch has
inverse-limit coordinates (0,-1), so no ordinary integer realizes it.
Multiplication by two reverses the 3-adic tournament orientation, but digit
avoidance is not invariant under the relevant kernel. Neither conjecture is
resolved by these maps or the bounded censuses.

### Binary 3-adic continuation and the supplied sequence

The [weak-Collatz/MB report](session65-weakcollatz/REPORT.md) identifies the
supplied sequence and proves the binary 3-adic survivor-set description,
density-one integer halting, and a sparse exception bound. Its alternating-lex
orientation is transitive and matches the stated tournament bridge. The
[independent audit](session65-audit-weakcollatz/AUDIT.md) separates those proofs
from finite native-machine and external source observations. The
[portable amendment](session65-audit-weakcollatz/PORTABILITY.md) reruns authored
math/native checks without redistributing unlicensed external source; historical
external-source gates are explicitly not rerun by that entry point.

### Goldbach machines, multiply-perfect constraints, and typed reductions

The [Goldbach and concept-lineage report](session65-twentyfive/REPORT.md) tracks
the 2016 27-state, later 25-state, and cited 21-state results without inferring
Goldbach from bounded runs. Licensed Lean sources are inspected, not locally
built. The [independent audit](session65-audit-goldbach/REPORT.md) checks the
exact source closure and transition-table transport. The
[multiply-perfect report](session65-root/multiperfect/REPORT.md) proves scoped
parity-budget consequences, not odd-perfect-number nonexistence.

The [typed reduction atlas](session65-reductions/REPORT.md) specifies carriers,
predicates, preserved structure and losses. A shared universal target is not
an equivalence of fixed problems; finite verification is not uniform halting
decidability. Its [independent audit](session65-audit-reductions/REPORT.md) and
[offline amendment](session65-audit-reductions/PORTABILITY_AMENDMENT.md) check
all eleven pinned inputs, identical mathematical consequences, and missing/
wrong-byte/provenance hostile controls without requiring Git at replay time.

## Evidence and reproduction

The [concept board](CONCEPT_BOARD.md) records preserved predicates, losses and
hostiles. Both lanes have separate independent audits and deterministic
normal/optimized receipts. The cumulative second checkpoint runs fourteen authored mathematical/provenance
checkers and hostile-test drivers in both Python modes, including every shipped
runtime dependency. A C compiler (`gcc`) is required for the independently
authored native MB replay; no external machine simulator is bundled. From this packet directory:

```sh
python3 verify_packet.py
python3 replay_frozen.py . REPLAY_TASKS.json /tmp/session65-replay.json --scratch /tmp/session65-replay
python3 -O verify_packet.py
```

The replay copies the complete packet before executing, clears ambient Python
import overrides, checks saved stdout, and rejects input mutation. The output
receipt and scratch directory must stay outside this packet. These commands
do not require a source checkout, network, Coq, Lean or third-party Python
packages. Python runs in isolated mode (`-I -B`). Historical source-audit entry
points are retained as evidence, not substituted for portable replay tasks. Frozen notes retain their original wording (including a historical
“audit pending” label); the linked final PASS audits resolve that audit state
without altering the sealed proof bytes.

[Provenance and licensing](PROVENANCE.md) distinguishes shipped source from
external observations. [Original-context links](SOURCE_LINK_INDEX.md) preserve
historical links in copied own-repository inputs. No newly built Coq/Lean
formalization is claimed. No scarce theorem identifiers were reserved.
