# Session 65: guarded machines and faithful tournament interfaces

**First independently audited checkpoint — 2026-10-26.** The exact statements
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

## Evidence and reproduction

The [concept board](CONCEPT_BOARD.md) records preserved predicates, losses and
hostiles. Both lanes have separate independent audits and deterministic
normal/optimized receipts. The full first checkpoint runs five authored
mathematical/provenance checkers in both Python modes, including every shipped
runtime dependency. From this packet directory:

```sh
python3 verify_packet.py
python3 replay_frozen.py . REPLAY_TASKS.json /tmp/session65-replay.json --scratch /tmp/session65-replay
python3 -O verify_packet.py
```

The replay copies the complete packet before executing, clears ambient Python
import overrides, checks saved stdout, and rejects input mutation. The output
receipt and scratch directory must stay outside this packet. These commands
do not require a source checkout, network, Coq, Lean or third-party Python
packages. Frozen notes retain their original wording (including a historical
“audit pending” label); the linked final PASS audits resolve that audit state
without altering the sealed proof bytes.

[Provenance and licensing](PROVENANCE.md) distinguishes shipped source from
external observations. [Original-context links](SOURCE_LINK_INDEX.md) preserve
historical links in copied own-repository inputs. No newly built Coq/Lean
formalization is claimed. No scarce theorem identifiers were reserved.
