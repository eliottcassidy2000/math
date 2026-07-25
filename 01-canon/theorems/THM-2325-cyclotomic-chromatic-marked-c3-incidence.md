---
id: THM-2325
title: "Cyclotomic-chromatic marked c3 incidence"
status: >
  SUPERSEDED / NO INDEPENDENT THEOREM. The reserved divisor-splitting
  candidate was overtaken before promotion: the nonnegative-sector
  strengthening of THM-2323 removes the cyclotomic-degree hypothesis and
  absorbs the whole reduced deepest cofactor into one modulus. Consequently
  no cofactor-25 table or separate chromatic-height theorem is needed.
  This file is not a proved dependency and claims no LRC(14) row exclusion.
source: codex-2026-07-25-cyclotomic-chromatic-incidence
depends_on: []
related:
  - THM-2323-primitive-fixed-colour-cross-correlation-and-same-gauge-word-alignment
---

# THM-2325 -- superseded reservation

**SUPERSEDED / NO INDEPENDENT THEOREM.**

The proposed proof split a reduced deepest cofactor as `r=sR`, paid a
totient-width condition with `s`, and used the chromatic number of
`K_7 x K_13` or `K_6 x K_13` to discharge `R`. Before this candidate was
committed, the THM-2323 proof was strengthened:

```text
fixed-colour correlation coefficients are nonnegative;
Galois transport sends a primitive zero to the standard root;
all surviving standard-root phases lie in |arg|<2*pi/7<pi/2.
```

Their real parts therefore have a strict common sign. The primitive
fixed-colour theorem holds for every modulus divisible by seven, with no
totient-degree condition. Taking the modulus to contain the entire reduced
deepest cofactor then gives the stronger incidence conclusion directly.

No result should cite THM-2325. Use the audited strengthened THM-2323.
