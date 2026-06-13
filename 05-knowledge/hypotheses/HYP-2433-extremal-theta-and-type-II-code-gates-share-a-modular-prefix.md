# HYP-2433 - Extremal theta and Type II code gates share a modular prefix

**Status:** COMPUTATION CONFIRMED through stored dimensions/lengths; proof lens OPEN.
**Source:** codex-2026-06-11-P4.
**Companions:** HYP-2425, HYP-2428, HYP-2429, HYP-2434, T784.
**Artifacts:** `04-computation/theta_code_lattice_gate_codex.py`,
`05-knowledge/results/theta_code_lattice_gate_codex.out`.

## Statement

For dimensions/lengths `24m`, extremal even-unimodular lattice theta series and
extremal binary Type II code weight enumerators have the same scalar
cancellation shape:

```text
kill the first m forbidden layers by a modular invariant basis,
then expose a large first shell.
```

For lattices the basis is `E4^(3m-3j) Delta^j`. For codes the basis is the
Gleason ring in `A=x^8+14x^4y^4+y^8` and
`B=x^4y^4(x^4-y^4)^4`.

## Evidence

The stored atlas computes:

```text
dim 24: theta first shell q^2 = 196560
dim 48: theta first shell q^3 = 52416000
dim 72: theta first shell q^4 = 6218175600

len 24: Type II A_8  = 759
len 48: Type II A_12 = 17296
len 72: Type II A_16 = 249849
```

The length/dimension-72 rows are both third-level cancellation gates: three
initial layers are killed, and the first nonzero layer is forced by the scalar
modular form.

## Proof Route

Treat both scalar shadows as the same filtered modular operation:

```text
choose invariant basis -> solve triangular vanishing system -> test support realization.
```

The triangular solve is easy and exact. The theorem-facing content is the last
arrow: realizing the scalar shadow as a lattice or a code support object.
