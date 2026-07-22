---
id: THM-2062
title: "Two-anchor hereditary-primitivity determinantal sieve and CRT wheel"
status: >
  RESERVED WITH A PROOF CANDIDATE UNDER ADVERSARIAL AUDIT. In a saturated
  rank-two coefficient template, every rank-two deletion has a finite
  determinantal divisor. Hereditary primitivity should exclude at most two
  projective parameter directions per prime, converting each THM-2058
  longitudinal coprime interval into an explicit squarefree CRT wheel. A
  rank-one deletion should instead force one primitive linear form to equal
  plus or minus one. Exact rank-zero, p-dividing-N, and saturation edge cases
  remain under independent review.
source: codex-2026-07-21-LRC-two-anchor-CRT-wheel
depends_on:
  - THM-2058
  - THM-2060
related:
  - THM-2052
  - THM-2053
  - HYP-8846
---

# THM-2062 -- Two-anchor hereditary-primitivity CRT wheel

Let thirteen coefficient rows `r_i=(a_i,m_i)` generate `Z^2`, and let the
primitive parameter be `d=(N,M)`, so `v_i=r_i.d`. For a deletion `ell` whose
remaining rows have rational rank two, reserve

```text
I_ell=gcd_(i<j; i,j!=ell)|det(r_i,r_j)|.
```

The candidate theorem says that a prime can divide
`gcd_(i!=ell) v_i` only when it divides `I_ell`, in which case the bad
parameters form the unique projective annihilator of the deletion rows.
Across all deletions, at most two such projective directions should occur at
one prime. Along a fixed-`N` THM-2058 interval this would become at most two
forbidden residue classes for `M mod p`, with a separate all-or-nothing rule
when `p|N`.

If a deletion has rational rank one, write its rows as integer multiples of
one primitive row `c`; hereditary primitivity should be possible only when
the coefficient multipliers have gcd one and `|c.d|=1`. No claim is made
until all local and rank-degenerate cases have passed audit.
