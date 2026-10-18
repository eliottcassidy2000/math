---
id: HYP-9291
title: "Cooper's s_7 and Apery's numbers agree at every Frobenius point modulo p^2: for every odd prime p and every r >= 1, with n = (p^r - 1)/2, s_7(n) == (-1)^n A(n) (mod p^2). The case r = 1 is THM-4669; the cases r >= 2 are conjectural"
status: >
  OPEN (CONJECTURE) for r >= 2. PROVED for r = 1 (THM-4669).

  FINITE-EXACT evidence for r >= 2 (13 cases, all hold; none holds modulo p^3):
  - r = 2: p = 5, 7, 11, 13, 17, 19, 23, 29, 31 (n up to 480);
  - r = 3: p = 5, 7, 11, 13 (n up to 1098).
session: opus-2026-10-18-S43
source: 05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
scripts:
  - 04-computation/experiments/s7_apery_halfperiod_extended_20261018.py (+ .out, .json)
related: [THM-4669, THM-4528, Ahlgren-Ono 2000, Cooper 2012]
---

# HYP-9291 — s_7 and Apéry at the Frobenius points

## Statement

`A(n) = Σ_k C(n,k)^2 C(n+k,k)^2` and `s_7(n) = Σ_k C(n,k)^2 C(n+k,k) C(2k,n)`.

**Conjecture.** For every odd prime `p` and every `r >= 1`, with `n = (p^r - 1)/2`,
```
s_7(n) ≡ (-1)^n A(n)   (mod p^2).
```

## Why the proof of THM-4669 does not extend directly

THM-4669 uses three steps:
- the exact identity `s_7(n) = Σ_k [C(n,k) C(n+k,k)] C(2k,k) C(k,n-k)`;
- the exact Clausen identity `Σ_k (-1)^(n-k) 4^k a_k^3 C(k,n-k) = Σ_k a_k^2 a_(n-k)^2`, valid for every `n`;
- two half-period congruences, (C1) and (C2), which need every `j <= n` to be a `p`-adic unit.

For `n = (p^r - 1)/2` with `r >= 2`, multiples of `p` occur among the `j`, so (C1) and (C2) fail termwise. The exact
identities survive.

A natural route is the Lucas/Dwork factorisation of `C(n,k)` along `p`-adic digits. Writing
`n = (p-1)/2 + p (p^(r-1) - 1)/2`, it reduces the conjecture to THM-4669 at each digit, plus a carry analysis modulo
`p^2`. This is the usual mechanism behind "Frobenius-point" supercongruences of Dwork type.

## Tests

- **A proof.** Carry out the Lucas/Dwork factorisation at the Frobenius points.
- **Larger `r`.** Check `r = 3` for `p >= 17`, and `r = 4` for `p = 5, 7`.
