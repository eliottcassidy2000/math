---
id: THM-795
title: Hamming-one full-residue rigidity at every lift height
status: CLAIM STUB — scale-free proof and seven finite atoms independently audited; full proof and Fraction replay in progress
source: codex-2026-07-14-S10
depends_on:
  - THM-770   # bounded full-residue classification and shallow chart
related: [THM-724, THM-769, THM-770, HYP-6775, HYP-6800, HYP-6820]
---

# THM-795 — Hamming-one full-residue rigidity at every lift height

> **Namespace claim stub.**  The analytic reduction and every finite rational
> atom have passed an independent exact audit.  The full proof, replay script,
> and shallow-frontier integration are being written now.

For `r in {1,...,12}` and every integer `k>=1`, put

```text
W_k(r)=({1,...,12}\{r}) union {r+13k}.
```

Then

```text
M(W_k(r)) > 1/13.
```

More generally, let `13` not divide `c` and replace `cr` in
`c*{1,...,12}` by a positive integer `w` congruent to `cr mod 13`.  If the
resulting packet is tight at `1/13`, then `w=cr`.  Thus no proper Hamming-one
perturbation of any shallow arithmetic-progression dilation is tight.

The scale-free mechanism uses the `c` lifted missing-owner splice germs.
Tightness forces all their `w`-phases into one closed arc of length `2/13`, so
the deck order `c/gcd(c,w)` is one and `c|w`.  Dilation descent leaves the
one-coordinate family above.  Eleven-speed core-safe intervals eliminate all
but seven low replacements, which have explicit strict rational witnesses.

