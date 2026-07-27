---
id: THM-2594
title: "The realized theta-slaved contraction: S != 0 at the r = 5 window on a common ancestry base"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT (one common-base table, one
  integral per entry; exact Phi_91 reduction with 5/5 split-prime
  cross-checks; python/-O identical modulo elapsed fields; hostile
  controls beta = 0 and flat-slaving vanish EXACTLY); independent
  hostile audit REQUESTED.  MISTAKE-286 TYPING: this is a positive
  instance of THM-2471's deep-root-sidecar CANDIDATE realization on
  the canonical typed row -- NOT a generic THM-2512 bridge theorem,
  NOT a physical current; no scalar row is excluded; ledger 165;
  LRC(14) OPEN.  MISTAKE-281/293 discipline: all factors lie on one
  finite Boolean ancestry fibre and one base integral before every
  marginalization and DFT, but at linked nodes w_u, X_(u,a), and
  Y_(q,e'), not at one circle point.  The k=6 table is verified
  directly; no claim that k=6 exceeds THM-2449's eventual k_0.
source: opus-2026-07-28 (executing the Stage-2 step opened by
  THM-2575/THM-2577's r = 5 window)
depends_on:
  - THM-2471 (stalk, colours, deep-root sidecar (44)-(47))
  - THM-2449 ((14)-(16) response factors; (15) source-safe deletion)
  - THM-2512 ((12)-(15) contraction and factorization; (17)/(19) ignition)
  - THM-2575 / THM-2577 (the r = 5 window at sigma = {b})
related:
  - THM-2560 (the r = 3 base-only diagnostic this supersedes at {b})
  - MISTAKE-281 / MISTAKE-286 / MISTAKE-293 (typing and discipline, obeyed)
script: 04-computation/lrc14_stage2_theta_contraction_opus_20260728.py
output: 05-knowledge/results/lrc14_stage2_theta_contraction_opus_20260728.out
script_sha256: 377b520761ca64dd5964c60584d177dff70bd2367e12a61af5169f1d03218488
output_sha256: 973e1479867900ad814f77106c4a87f12da6d55e971c5c412860fa1afecb7abf
hash_basis: working-tree bytes (LF)
---

# THM-2594 -- the theta-slaved contraction survives

Row `THM-2309 (25)`, owner `j = 1`, word `sigma = {b}`, packet clock
`K = 2`, exact auxiliary response exponent `k = 6` (class
`k = 6 mod 16380`, `eps = +1`); collision at `d = 13^5` (`r = 5`, the
unique-invariant window).  The exponent is evaluated directly below; it is
not asserted to lie beyond THM-2449's unspecified eventual threshold `k_0`.

## Statement

Build the ONE common-base table
`N(u,q,ell,theta) = int A(y,u) F(y,q) cell_ell(y) DEEP_theta(y) QW_ell(y) dy`
(`U = P_{13^5} f`, `V = P_{13^5} e`, `f = 1_{Q_{1,{b}}} P^2 1_{E_1}`),
where the THM-2449-shaped factors are evaluated on their typed linked stalk
nodes (owner clock at `w_u`; deep probe and word at `X_(u,a)`).  They descend
to the common base because the former is root-blind and the latter two are
sheet-independent, realizing precisely `theta = t - 2u` per THM-2471 (45).
Doubly centre and form the THM-2512 (12)-(14) toothpick contraction with the
output index slaved through `theta`. Then at the fixed primitive quadruple
`(tau, a0, alpha, beta) = (1,1,1,1)`:

```text
S = Psi_{1,1}(1,1) != 0   exactly in Q(zeta_91)
```

(72-coefficient integer numerator over
`DEN_S = 2^4 3^3 5 7^3 11 13^19 53`; split primes 547/911/1093/2003/
2549 give 295, 213, 380, 1405, 883, all nonzero). By the exact
factorization `Psi = K_{alpha tau, beta} dtilde(alpha, -beta a0)`
(verified at 4 quadruples) and rational-table Galois transitivity,
**all 5,184 primitive cut coefficients of the slaved bundle fire.**

## Exact Boolean-fibre meaning

Use THM-2471's stalk coordinates

```text
w_u=(y+u)/13,
X_(u,a)=(w_u+a)/d,
Z_(u,a,b)=(X_(u,a)+b)/13^2,
Y_(q,e')=(w_q+e')/d.
```

Expanding `A(y,u)F(y,q)` as in THM-2471 (31) gives a normalized finite sum
of the Boolean products

```text
1_Q(X_(u,a)) 1_E(Z_(u,a,b)) 1_E(Y_(q,e')).
```

At exponent `k=6`, direct identities give

```text
d_(1,ell)(w_u)=cell_ell(y),
Delta_(2u+theta)(X_(u,a))=DEEP_theta(y),
Q^(+1)_ell(13^6 X_(u,a))=QW_ell(y).
```

All three right sides are independent of the marginalized sheets.  Hence
`N(u,q,ell,theta)` is exactly the normalized sum of the corresponding
Boolean products with these three factors inserted *before* integration and
marginalization.  This proves one common ancestry fibre.  It does not put the
three factors at one circle point and does not identify `N` with THM-2449's
one-variable table `H^R`.

## The signal is the slaving

Both hostile controls vanish IDENTICALLY: `beta = 0` (row-zero) and
the u-independent (flat) slaving. So the nonzero value is carried by
the `theta = t - 2u` coupling itself -- the first exact instance where
the collision root, deep phase, owner clock, word, and cut probe pair
nonzero INSIDE one integral on a live-branch row.

## Structural facts

- The slaving transports deep support with the root: per root `u`,
  only deep labels `t in {2u, 2u+1, 2u+2}` can fire
  (`DEEP_theta` empty for `theta in {3..12}`).
- `ell = 0` row vanishes (word guard-safe factor disjoint from
  `cell_0`); reflection symmetry `(ell,theta) <-> (7-ell, 2-theta)`.
- `34/169` ancestry fibres `(u,s)` carry nonzero interaction (none at
  `s = 0`).
- Gates: `J(0) = I_5 = 48602521488933856/337437093630814766589`
  reproduced from scratch; all colour/zero-sum/Hermitian laws;
  independent route-2 aggregate; toy-model pipeline validation;
  160-point slaving spot-check.

## What follows and what does not

The ancestry/Boolean realization program now has, on one row: the
owner-clock host array (THM-2575 B), the all-clock r = 5 window
(THM-2577), and a REALIZED nonzero theta-slaved contraction (this
theorem). Remaining, per MISTAKE-286/293's honest typing: relating this
candidate realization to THM-2512's generic bridge, to a target-active
`B(q)`-type object, to THM-2449's one-point lawful response table, and to
anything uniform over the 165 rows. NOT implied: a physical current, a
transplant theorem, any row exclusion, or LRC(14).
