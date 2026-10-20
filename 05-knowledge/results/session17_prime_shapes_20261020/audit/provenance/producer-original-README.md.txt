# Session 17: two exact arithmetic–graph bridges with their marking boundaries

The two proof notes are complementary:

1. [Prime balance and tetrahedral complement-pairs](PRIME_BALANCE_TETRAHEDRAL_BRIDGE.md)
   recovers F, S, U; identifies the same 4-singleton/3-doubleton fibre law
   on divisors and the three-edge tetrahedral graph family; constructs an
   oriented bijection; and proves that full unmarked S4-equivariance is
   impossible.
2. [F5 path pairs and the icosahedral return graph](F5_PATH_PAIR_ICOSAHEDRON.md)
   makes 625→576 an actual zero-row/zero-column deletion, then transports
   twelve signed vector classes to the previously proved regular-five
   two-complement return component. It retains the determinant relation
   and gives the precise SL2(5)/±I action, while exposing dependence on the
   chosen five-state arithmetic marking.

These are scoped structural connections, not an all-order trijection or a
Collatz proof. The old divisor classification, source path decoder,
three-edge graph classification, and regular-return theorem are explicitly
credited. No new canon theorem ID is reserved.

## Inheritance and active board

- Anchor: the user's three-edge four-vertex family and F=S+U.
- Niche: why the same ten-point carriers have inequivalent S4 actions.
- Wildcard, promoted after exact verification: the F5 determinant model of
  the twelve-state regular return component.
- Concepts retained: exponent boxes; actual fibres; complement versus
  opposition; orientation gauges; pointed path products; determinant and
  action; source legality.
- Closest mechanisms: the prior prime-to-extra map, ten-letter simplex
  decoder, balanced path product, and explicit two-pentagon return graph.
- Hostiles: the transposition fixed-point mismatch; divisor complement's
  different cycle type; x↦x^3 changing determinant adjacency; rejected
  divisor profiles re-entering balance under Collatz halving.
- Corrected near miss: equal counts/orbit sizes do not give an equivariant
  map, and 25^2−25 is 600, not 576.
- Least-used sidecars recovered: the arithmetic prime-to-extra fibre map
  and the orientation-parity sheet of labelled regular five-tournaments.

The arithmetic–tetrahedral map preserves a quotient and fibre sizes, while
the field–return map preserves a specified adjacency. Those are different
contracts. Neither is promoted to a common universal theorem merely because
both use small symmetry groups.

## Checks

`verify.py` validates the tetrahedral and divisor note. `verify_field.py`
independently validates the field model against the prior theorem's explicit
bit codes. Both use active exception checks and exact integer arithmetic.
Normal and optimized runs agree byte-for-byte; outputs and hashes are
retained. Finite checks corroborate the printed proofs and do not replace
universal quantifiers.

A bounded startup pass read the current canon routes and guardrails at
f3e658d7cc10. The message processor was unavailable because no `.machine-id`
was configured; no identity was invented and no repository files were
modified by this staged research. All cited historical sources were read
as provenance, not as authorities overruling current corrected statements.

## Open next question

Can the balanced-path marking be characterized by a native operation on
actual Hamiltonian paths, strong enough to select determinant adjacency
without imposing an arbitrary F5 labelling? The x^3 hostile shows that the
pointed five-set plus its involution alone is insufficient. A positive
answer must identify an additional source relation and prove its transport,
not merely recognize another icosahedral count.
