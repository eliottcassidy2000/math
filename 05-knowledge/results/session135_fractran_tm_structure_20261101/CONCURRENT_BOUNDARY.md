# Concurrent conjecture-machine work: a pinned source boundary

**Read-only source interpretation, not an independent audit or a current-status
override.** The concurrent commits below were preserved as ancestors of this
session's checkpoints. Consult the current canonical source if later review
changes its status. None is a proved dependency of our basin or TM results.

## Later source audit checked before final publication

[Audit commit ea0ca603](https://github.com/eliottcassidy2000/math/commit/ea0ca6035588eef8ac41af54d87b8ca077ce14b0)
records **SOUND WITH CORRECTIONS**, followed by close-out30b6a26. We read
[the audit report](https://github.com/eliottcassidy2000/math/blob/ea0ca6035588eef8ac41af54d87b8ca077ce14b0/04-computation/experiments/circulant_fractran_goldbach_20261101/audit/AUDIT.md);
we did not independently replay that other session's audit scripts.

The corrected circulant count is131068, not131071. The all-order arc argument
covers a sufficient subfamily of inherited THM-4730 paths; the broader
three-interval assertion remains only finitely checked. A blanket affine-
isomorphism shortcut is replaced by the required prime-degree tournament
Burnside/normalizer argument, and family provenance is credited. Our elementary
word/path correspondence uses translations only and does not depend on that
shortcut.

The same report accepts the three program equivalences given the cited kernel
contract, while retaining G4/L3 as failed as originally stated. The25-state
Goldbach machine and PRIMEGAME share a trial-division idea; their divisor
orders differ and no predicate-preserving map was established between those
two particular programs. This is distinct from our emitted, audited machine
translations. The historical pre-audit snapshot below is retained for lineage,
not asserted to be the current review status.

## Historical pre-audit snapshot

At [commit1224a73](https://github.com/eliottcassidy2000/math/commit/1224a73dfc37084ea96a775fc8cfbba30e69a879),
the exact-slug [THM-4856 source](https://github.com/eliottcassidy2000/math/blob/1224a73dfc37084ea96a775fc8cfbba30e69a879/01-canon/theorems/THM-4856-goldbachgame-fermatgame-legendregame-halts-iff-counterexample-fractran-programs-built-on-primegames-divisor-loop.md)
labels its package **PRE-AUDIT** and reports three conditional kernel-based
program equivalences: GOLDBACHGAME37 fractions, FERMATGAME36 and
LEGENDREGAME33. Their proposed halt conditions are a Goldbach counterexample,
a sixth Fermat prime, and a prime-free square interval, respectively.
The author reports finite positive controls and literal/dispatch comparisons.
We have read those source claims, not replayed their reported large-step-count runs
or independently promoted the complete program proofs.

The source explicitly retains two failed controls as stated:

- G4 halted during target12, but its check reported11 after the tested unit
  had been consumed. This is a checkpoint/observation failure, not a Goldbach
  counterexample.
- L3 truncated a register that also controlled the next interval jump. The
  source reports a corrected separate-budget control halting at n5. Do not
  label the original L3 run PASS retroactively.

The [circulant source at the same commit](https://github.com/eliottcassidy2000/math/blob/1224a73dfc37084ea96a775fc8cfbba30e69a879/01-canon/theorems/THM-4855-circulant-tournaments-as-s-sequencings-segment-sum-criterion-and-explicit-non-palindromic-hamiltonian-paths-for-every-odd-circulant-of-order-at-least-seven.md)
uses actual connection-set words and distinct modular partial sums, with
finite witness checks and an asserted all-order construction. We recover only
the elementary word/path correspondence in our synthesis; we do not claim
independent verification of its complete all-order witness families.

## The reusable mechanism and its limit

A total exact primality kernel can be wrapped in an outer candidate search.
If the wrapper enumerates every required candidate, resets correctly, and
halts precisely on its failure predicate, the corresponding halts-iff
statement follows. Those are separate proof obligations for each wrapper.
Changing PRIMEGAME's exit1/17 to a new control and replacing its outer reset
changes the program. The original three-basin theorem cannot be transferred
merely because the same eleven inner fractions appear.

Our costed word certificate can validate a particular repeated trace, and our
fixed-horizon cap can certify finite prefixes. Neither supplies an unbounded
success theorem for a conjecture search. Likewise37 fractions versus a
25-state Turing description is not a model-independent complexity comparison:
alphabets, tapes, input encoding, guards, outputs and transition costs differ.

For a concrete encoding boundary,43 is an inert spectator in the original
PRIMEGAME: it reaches2365<->2795 (43 times the55<->65 cycle). In the other
list,43 is GOLDBACHGAME's initial control. A common integer or shared inner
fractions do not identify the programs' input meanings or global dynamics.
