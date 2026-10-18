---
id: HYP-9291
title: "Cooper's s_7 and Apery's numbers agree at every Frobenius point modulo p^2: for every odd prime p and every r >= 1, with n = (p^r - 1)/2, s_7(n) == (-1)^n A(n) (mod p^2). RESOLVED: PROVED in S43 as part of THM-4669"
status: >
  PROVED (resolved the same session): THM-4669 now covers every r >= 1. The steps (C1) and (C2) of its proof hold at
  n = (p^r - 1)/2 by valuation bounds; s7_apery_frobenius_proof_check_20261018.py checks them p-adically.
  History: filed as OPEN for r >= 2 after the r = 1 proof; the obstacle named below dissolved.

  FINITE-EXACT evidence for r >= 2 (13 cases, all hold). None holds modulo p^3 except at the non-ordinary prime
  p = 11, where v_11(A) = v_11(s_7) = r. The congruence then holds mod 11^3 at r = 2, mod 11^4 at r = 3 and mod 11^5
  at r = 4. The
  earlier wording 'none holds modulo p^3' was false there (S43 audit m9).
  - r = 2: p = 5, 7, 11, 13, 17, 19, 23, 29, 31 (n up to 480);
  - r = 3: p = 5, 7, 11, 13 (n up to 1098).
session: opus-2026-10-18-S43
source: 05-knowledge/results/factorization_vertex_supersymmetric_bridge_20261018.md
scripts:
  - 04-computation/experiments/s7_apery_halfperiod_extended_20261018.py (+ .out, .json)
related: [THM-4669, THM-4652, HYP-9274, Ahlgren-Ono 2000, Cooper 2012]
---

# HYP-9291 — s_7 and Apéry at the Frobenius points

## Statement

`A(n) = Σ_k C(n,k)^2 C(n+k,k)^2` and `s_7(n) = Σ_k C(n,k)^2 C(n+k,k) C(2k,n)`.

**Statement (now PROVED as THM-4669).** For every odd prime `p` and every `r >= 1`, with `n = (p^r - 1)/2`,
```
s_7(n) ≡ (-1)^n A(n)   (mod p^2).
```

## Why the proof of THM-4669 seemed not to extend (resolved)

THM-4669 uses three steps:
- the exact identity `s_7(n) = Σ_k [C(n,k) C(n+k,k)] C(2k,k) C(k,n-k)`;
- the exact Clausen identity `Σ_k (-1)^(n-k) 4^k a_k^3 C(k,n-k) = Σ_k a_k^2 a_(n-k)^2`, valid for every `n`;
- two half-period congruences, (C1) and (C2), which need every `j <= n` to be a `p`-adic unit.

For `n = (p^r - 1)/2` with `r >= 2`, multiples of `p` occur among the `j`, so (C1) and (C2) fail termwise. The exact
identities survive.

**Resolution.** No digit factorisation is needed.
- In (C1), compare the ratio `C(n,k) C(n+k,k)/((-1)^k a_k^2) = Π_j (1 - p^(2r)/(2j-1)^2)`. Each factor is
  `≡ 1 (mod p^2)`, because `v_p(2j-1) <= r - 1`.
- In (C2), each `p^r/(2(n-i))` lies in `p Z_p`.

The alternative route, kept for the record, was the Lucas/Dwork factorisation of `C(n,k)` along `p`-adic digits. Writing
`n = (p-1)/2 + p (p^(r-1) - 1)/2`, it reduces the conjecture to THM-4669 at each digit, plus a carry analysis modulo
`p^2`. This is the usual mechanism behind "Frobenius-point" supercongruences of Dwork type.

## Tests

- **Proof: done.** The valuation-bound argument, in THM-4669.
- **Still open: higher powers of `p`.** At `r = 1` no congruence holds modulo `p^3` for any odd `p < 4000`. At
  `r >= 2` the only exception found is `p = 11`, which is non-ordinary. Ask whether a corrected
  form, with an explicit `p^2` term, holds modulo `p^3`.
