# Proof and exact-computation audit

**Same-agent audit with alternate arithmetic routes.** No independent referee,
formal Lean build, or major-open-problem resolution is claimed. Local theorem
labels are qualified by their file paths. Results remain in this research
packet and do not consume reserved canon IDs.

Reproduction from the repository root:

```text
python -B 05-knowledge/results/stabilizer_necklace_kernels_20261028/audit.py
python -B agents/check_docs.py
```

The audit reruns all three mathematical scripts normally and with `-O`,
compares their JSON byte-for-byte, verifies links, and checks the actual
81/27 pair. [audit_receipt.json](audit_receipt.json) records checks and hashes.
The source-recovery census is provenance, not a mathematical test.

| Check | Declared universe | Independent evidence / limit |
|---|---|---|
| Weighted stabilizers | Full recurrence through k=4; literal binary through k=2; literal ternary palette at k=2 | Actual groups and canonical orbits; orbit-stabilizer mass identity |
| Conditional futures | Two binary palettes, k=0..2, A_k and C_k | Each actual first stabilizer evaluated by cycle index |
| Sign recurrence and extrema | k=0..10, all digits through59048 | Recursive blocks versus digit formula; extrema checked on every prefix |
| Integer determinant | Literal matrices through N=81 | Fraction-free Bareiss versus Lucas product recurrence |
| Smith normal form | Literal matrices through N=27 | Euclidean row/column operations versus the CRT-derived invariant factors |
| Finite-field rank | N=1,3,9,27,81,243; p=2,3,19 | Polynomial gcd including characteristic3's repeated roots |
| Frobenius equation | Formal coefficients through N-1 for k<=6 | Direct convolution versus Mahler sections |
| Fixed-clock hostile | All81 A_2 translations; all3 C_2 translations and futures | Same colouring orbit, unequal determinant; complete three-state future separation |
|81/27 | Two literal12-vertex tournaments and every4096 subset | Equal cyclic triples/Hamilton counts; different full stabilizer distributions |

The generic Euclidean Smith run at81 was stopped after poor coefficient-growth
performance; it is not recorded as completed. The successful literal Smith
universe ends at27. The all-depth theorem instead rests on the comaximal-ideal
proof in COKERNELS C1 and the established odd golden-ring normal form.

Proof audit: all recursion bases include k=0; the minimum-set formula states
k>=1 with the trivial case separately; palette multiplicity is not cubed
by cost substitution; regular means setwise stabilizer trivial; finite-field
root arguments handle p=2,3,5 separately; Smith factors count both the
2-primary part and odd part; the cyclic shift's orientation and twist are
explicit. The kernel statistics belong to a specified colouring and cyclic
coordinate, not the unmarked tournament class. Their failure to descend is
demonstrated and repaired by retaining the clock.
