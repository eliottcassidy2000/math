# HYP-2434 - Dimension 72 splits lattice support from binary code support

**Status:** OPEN synthesis; scalar data confirmed, external lattice existence cited.
**Source:** codex-2026-06-11-P4.
**Companions:** HYP-2415, HYP-2425, HYP-2430, HYP-2433, OPEN-Q-066,
T784.

## Statement

Dimension/length 72 is the cleanest support-realization split in this thread:

- an extremal even-unimodular 72-dimensional lattice of minimum 8 exists
  (Nebe);
- an extremal binary doubly-even self-dual `[72,36,16]` code is still open;
- both scalar modular shadows pass their triangular cancellation gates.

Thus the obstruction is not "modular forms at 72." It is which support category
can realize the scalar shadow.

## Evidence

The theta/code atlas computes:

```text
lattice dim 72: kill q^1,q^2,q^3 -> q^4 shell 6218175600
code len 72:    kill weights 4,8,12 -> weight 16 shell 249849
```

Nebe's construction gives the lattice support side; the coding literature still
treats the `[72,36,16]` code as open and support-constrained.

## Prediction

Any proof or disproof of the binary code should be expressible as a retained
support obstruction that the theta scalar shadow forgets: a design layer,
neighborhood constraint, automorphism restriction, matroid/Tutte leakage, or
polarization/frame obstruction.

Scalar modular positivity should not be used as a proxy for code existence.
