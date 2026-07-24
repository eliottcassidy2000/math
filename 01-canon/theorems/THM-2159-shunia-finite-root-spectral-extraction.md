---
id: THM-2159
title: "Shunia finite integer-root extraction by cyclic spectral mixing"
status: >
  RESERVED / UNPROVED PROOF CANDIDATE UNDER AUDIT. Namespace reserved for a
  self-contained proof of Shunia's finite modular integer-root formula. The
  planned proof uses the cyclic Markov operator induced by multiplication by
  1+x modulo x^n-a, an explicit nontrivial-mode contraction at K=2a^n, and a
  lossless Kronecker base-X evaluation. No statement in this file may be used
  as proved until the coefficient, remainder, and integer-division
  inequalities and the independent exact companion have been audited.
source: opus-2026-07-24-puzzle-atlas
depends_on: []
related:
  - THM-2051
  - THM-2054
external:
  - "Joseph M. Shunia, Polynomial quotient rings and Kronecker substitution for deriving combinatorial identities, arXiv:2404.00332v6"
---

# THM-2159 -- reserved Shunia finite-root extraction

**RESERVED / UNPROVED PROOF CANDIDATE UNDER AUDIT.**

The intended statement is the natural-number identity from Shunia's
Conjecture 6.1.  For admissible `a,n`, put

```text
K=2a^n,  X=a^K,  M=X^n-a,  R_t=(X+1)^t mod M.
```

The target is

```text
floor(a^(1/n)) = (R_(K+1) div R_K)-1.
```

The proof candidate must separate:

1. cyclic spectral mixing of quotient-ring coefficients;
2. distance of the algebraic root from the adjacent integers;
3. lossless base-`X` evaluation and the least-remainder check; and
4. stability of the consecutive-remainder quotient under lower digits.

The LRC comparison is methodological only: the absence of a uniform
denominator/quantization in OPEN-Q-108 is precisely the coordinate that
prevents a direct transfer.
