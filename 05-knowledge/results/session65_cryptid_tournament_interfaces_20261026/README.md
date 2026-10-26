# Session 65: guarded machines and faithful tournament interfaces

**Final independently audited research checkpoint — 2026-10-26.** The exact statements
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

### A common framework, with necessary obstruction tests

The [synthesis](session65-root/SYNTHESIS.md) and
[connection framework](session65-root/CONNECTIONS.md) organize the results by
carrier, dynamics, observation, preserved predicate, and explicit obstruction.
They do not identify the open problems with one another. The
[fifteen-kernel corollary](session65-root/fifteen-kernel/REPORT.md) realizes the
signed-pair formula on the actual fifteen duads and fifteen synthemes:
transposition defects have minima 4 and 6, respectively. The
[independent integration audit](session65-audit-synthesis/REPORT.md) verifies
all actual pair orbits and literal minimizers, not a census of all 2^105
orientations. Its [portable amendment](session65-audit-synthesis/PORTABILITY_AMENDMENT.md)
retains all mathematics, checks 174 shipped dependencies, and explicitly labels
the one omitted historical source scan as not locally verified.

The additional [literal-write obstruction](session65-k-write/REPORT.md) proves
that a finite chain of full-order, fixed-label K operations cannot implement
one literal tape-bit write when n=3m and m is 2 or 3 modulo4. This parity gate
does not prohibit unlabelled or relabel-enabled simulations, local-factor K,
or padded/gadget encodings; bypassing it is not a sufficiency theorem.

The separately checked [palindrome guard](session65-palindrome-guard/REPORT.md)
sharpens the native-color gate: E(w) is pure blue exactly for palindromic w;
otherwise every Hamiltonian word is black. Framing as E(0w1) removes this
particular black-move obstruction, not the remaining K-simulation problem.
The native-K discussion is for positive order; the empty tape has only its
explicit vacuous color convention. This appendix leaves the sealed synthesis
unchanged and has its own frozen normal/optimized replay.

### Preserved incoming work and scope corrections

The incoming analogy packet is historical and is not covered by Session65's
mathematical audits. Read the [incoming-review overlay](session65-incoming-review/REPORT.md)
before using its classification claims: it corrects the common-target/fixed-Pi2
scope, separates JC2 from DC2, supplies the audited LRC finite-check scope, and
repairs Goldbach machine lineage and verification labels. Original incoming
bytes are retained. Later concurrent and pre-audit packets are preserved, not
silently promoted into this session's conclusions.

In particular, the incoming standalone Bell correction concerns the separate
three-state table `1RB0LC_1RC1RA_1LA0RB`; our named bouncer-plus-bell example is
the four-state table `1RB0LA_0RC1LA_1RD0RA_0LB1RB`. Neither correction nor macro
formula is transferred between these different machines.

## Evidence and reproduction

The [concept board](CONCEPT_BOARD.md) records preserved predicates, losses and
hostiles. Both lanes have separate independent audits and deterministic
normal/optimized receipts. The cumulative final checkpoint runs nineteen authored mathematical/provenance
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
