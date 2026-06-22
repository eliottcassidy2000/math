---
id: THM-567
title: QR/NQR even residue pairing is forced exactly when -1 is a non-residue
status: PROVED
date: 2026-06-22
session: codex-2026-06-22-S94
related:
  - HYP-2854
  - HYP-2857
  - HYP-2861
  - HYP-+2866
  - HYP-2867
  - OPEN-Q-108
---

# THM-567: QR/NQR Even Residue Pairing

## Statement

Let `q` be an odd prime.  Let `Q` be the nonzero quadratic residues in
`F_q`, and let `N` be the non-residues.  For a function

```text
F : F_q^* -> A
```

with values in any additive abelian group, assume `F(r)=F(-r)` for every
`r in F_q^*`.

If `q == 3 mod 4`, then

```text
sum_{r in Q} F(r) = sum_{r in N} F(r).
```

Equivalently, for every even residue mass, the QR shell and NQR shell have
equal total mass.

Conversely, if `q == 1 mod 4`, evenness alone does not force equality.  In
that case there are even functions for which the two shell sums differ.

## Proof

The Legendre symbol gives

```text
(-1 | q) = (-1)^((q-1)/2).
```

Thus `-1` is a non-residue exactly when `q == 3 mod 4`.

If `q == 3 mod 4`, multiplication by `-1` sends every quadratic residue to a
non-residue and every non-residue to a residue.  Since `q` is odd, this map has
no fixed point on `F_q^*`; it is a bijection

```text
Q -> N,    r |-> -r.
```

Using `F(r)=F(-r)`,

```text
sum_{r in Q} F(r) = sum_{r in Q} F(-r) = sum_{s in N} F(s).
```

This proves the equality.

If `q == 1 mod 4`, then `-1` is a quadratic residue, so multiplication by
`-1` preserves both `Q` and `N`.  Evenness therefore does not move mass between
the two shells.  For example, the function `F=1_Q` is even because
`r in Q iff -r in Q`; its QR sum is `|Q|` and its NQR sum is `0`.  Equality is
not forced.

## LRC14 Application

For `q=7`, `q == 3 mod 4`, so every even residue mass on `F_7^*` is split
exactly half between

```text
QR(7)  = {1,2,4},
NQR(7) = {3,5,6}.
```

In Fourier language, real circle indicators have conjugate-symmetric
coefficients, so absolute or energy masses grouped by `n mod 7` are naturally
even under `n -> -n`.  THM-567 therefore explains why the q-uniform
correlation constants in HYP-2854/HYP-2857 balance through the quadratic
character at `q=7`.

The theorem does not by itself bound the signed phase of
`chat(n) * conj(ghat(n))`.  It supplies the exact shell-balance carrier; the
remaining LRC14 work is to control the signed Paley-cut projection of the low
resonance channels.
