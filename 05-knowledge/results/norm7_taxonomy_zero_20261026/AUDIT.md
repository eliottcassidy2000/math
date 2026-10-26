# Proof audit, exact computation, and publication scope

2026-10-26. Audits distinguish cited structural inputs, deductions, exact
finite universes and open arithmetic identifications. No Lean build is claimed.

## Independent reviews

- **Quadratic atlas: PASS.** `planar_collatz` checked the classical quaternion
  splitting iff, the inherited residue-field embedding, all four minimum-degree
  cases and relative versus absolute torification. The [review](certificates/AUDIT_QUADRATIC_FIBRE_ATLAS.md)
  records an independent117-gate full-triple congruence path in addition to the
  author's16,773-gate replay.
- **Extension transport: PASS.** `arithmetic_shapes` checked E1–E3 and all signs,
  then used different defining polynomials for F9,F27,F49 and literal cubic
  products. Its [review](arithmetic/AUDIT_EXTENSIONS.md) records363 independent
  controls and359 producer gates. In particular it checked both A/B involutions;
  the agreement is not a comparison of two imports of the same evaluator.
- **Taxonomy: PASS by root.** Reviewed the intrinsic block recovery argument,
  strong-module recursion in both directions, carousel primeness, successor
  reconstruction, full coloured stabilizer and its splitting, Möbius inversion,
  and complemented reflection. Reopened BCIZ Theorem24, Proposition51 and
  Corollary57. These justify the full labelled realization fibre; the elementary
  81/27 hostile does not depend on their hard direction. Reviewed the adjacency
  backtracker and the full-permutation oracle, then replayed the producer.
- **Conductor support: PASS by planar_collatz and root.** The
  [independent review](certificates/AUDIT_CONDUCTOR_SUPPORT.md) checks the
  genus-zero integral compactification, terminating odd-pair blowups, double
  cover charts and full boundary, and finite locally constant cohomology over
  the integral base. Both reviewers read Knudsen's actual theorem page and
  the Stacks statements. The argument does not use the false assertion that
  arbitrary generic Galois extensions of unramified representations stay
  unramified. The support assertion is not transferred to the hourglass's
  characteristic-three degeneration.
- **Finite certificates:** the [independent review](taxonomy/AUDIT_CERTIFICATES.md)
  checks the reductions, determinant, actual rank-eight Galois hostile and exact
  GL₂ comparison hypotheses. Three small issues were repaired: nontrivial
  codomain for the late perturbation, positive Collatz inputs via n+1, and an
  explicit q*a84 column check. Final receipts have1,494 gates.

The incoming corrected THM-4803 source is retained unchanged, along with its
stronger independent proof audit. No duplicate canon promotion is made here.

## Computation universes and controls

| Checker | Universe and independent or hostile control |
|---|---|
| `arithmetic/check_quadratic_atlas.py` | Rational squareclasses with numerator−30..30, denominator1..15; squarefree δ−40..40 and D−100..100; explicit local primes and exact quadratic witnesses. Dyadic and real-place omissions are rejected. |
| `certificates/audit_quadratic_atlas.py` | All primitive triples modulo8,9,25,49 for signed squarefree δ−20..20; independent local-embedding and minimum-degree controls. |
| `motives/check_extensions.py` | Six specified graphs; p=3,5,7,11,13,17,19 and F_(p²), plus F27. Literal products at q≤27, elliptic recurrence on base-field parameters, singular fibres, both signs, tagged undefined versus certified zero. |
| `arithmetic/audit_extension_independent.py` | Five independently modelled fields of sizes5,7,9,27,49, literal products, A/B transforms and all-field proof review. |
| `taxonomy/check_taxonomy.py` | All1,099 labelled tournaments through order5 against n! canonization; all pairs of inflated words at lengths3/5; all words at odd lengths through15; maximal-module enumeration, Hamilton-path DP, Paley7, selfconverse-block and noncanonical-cut hostiles. |
| `certificates/check_certificates.py` | Thirteen explicit primes≤103, exact and modular determinants, four independently recomputed raw sums, eta/moment controls, fourteen actual noncubic primes with43 withheld, extension/lattice/isogeny losses, rank-sensitive Newton identities. |

All runtime gates are explicit exceptions and remain active under `python -O`.
Every frozen normal/optimized pair is byte-identical. The packet verifier
checks a closed file manifest before replay, writes temporary outputs, and
never overwrites a failed expectation. The optional PARI replay independently
reconstructs the modular coefficient rows. Hashes establish byte identity,
not mathematical validity; the reviews above supply the separate proof check.

No finite prime bank identifies a pure TT7 motive here. No finite tournament
enumeration proves all-order black-core connectivity. The inherited K3 and
Fano geometry is recovered through explicit cited repository maps, not by
equating numerical occurrences of six and seven.
