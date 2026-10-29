# Sha-inspired recurrences with marked ports, and an F1 boundary

2026-10-29. **PROVED scoped / CITED classical inputs / FINITE-EXACT /
OPEN arithmetic realization and general Sha finiteness.** Same-agent proof
and alternate-computation audit; no external priority claim.

This continuation incorporates the concurrent [Session96 Sha packet](../session96_sha_pairing_growth_20261029/README.md)
before extending it. The strongest new result is an exact repeated-join
law retaining the port, followed by a sharp two-port escape from it.

* [JOIN_RECURRENCES](JOIN_RECURRENCES.md), J1--J2: for base discriminant
  exponent d and marked-port order h, **every join of r>=2 copies has
  exponent lcm(d,h^2)**. It stabilizes after the first doubling, while
  the group order continues to grow. For base C_p^2 and nonzero port,
  the whole paired group has recurrence G_(r+2)=G_r + H_(p^2).
* J3: one cross channel cannot introduce new discriminant primes. Two
  channels can: explicit eight-vertex tournaments attach two C_3^2 blocks
  to produce **C_5^2 and C_11^2**. Their joint port contractions explain
  both values, and rank two is the sharp threshold in this matrix category.
* [F1_QUADRATIC_BOUNDARY](F1_QUADRATIC_BOUNDARY.md): the quadratic geometry
  used by the classical Selmer model and a symplectic geometry have the
  same thin value 2^n. Their first derivatives differ by n*2^(n-1).
  In dimension four over F_5 their maximal-isotropic counts are 12 versus
  156; over F_11, 24 versus 1464. The ambient form is indispensable.
* [RECOVERY_AND_ARITHMETIC](RECOVERY_AND_ARITHMETIC.md): recovered September
  doubling/squaring forests, October 4 tournament modes, October 13 Monsky
  receipts, and October 14 golden clocks. Reproducing classical descent
  on E_n:y^2=x^3-n^2x proves the **2-primary Sha vanishes for n=5 and 11**,
  with ranks 1 and 0 respectively, using cited inputs. This is a real
  arithmetic calibration, not a new general finiteness result.

The finite-pairing model of Sha via alternating matrix cokernels is
classical (BKLPR); the new synthesis is the marked composition law and its
obstructions, not the invention of that model. Tournament matrices have
odd discriminant and cannot alone encode the twice-square 2-primary defect.

## Reproduction and limits

From the repository root:

```text
python -B 05-knowledge/results/sha_marked_join_f1_20261029/check.py
python -O -B 05-knowledge/results/sha_marked_join_f1_20261029/check.py
python -B 05-knowledge/results/sha_marked_join_f1_20261029/audit.py
```

The producer uses exact inverses, Bareiss determinants and finite-field
ranks. The audit separately computes Smith valuations modulo p^3 and
symbolic polynomial jets. Literal repeated joins are tested through eight
copies (up to 64 vertices), for p=3,5,11; these are selected sources, not
an exhaustive tournament census. All two-dimensional subspaces of F_q^4
and all 2-by-2 graph matrices are exhausted for q=2,3,5,11. General claims
have proofs rather than extrapolated tables. See [AUDIT](AUDIT.md).
