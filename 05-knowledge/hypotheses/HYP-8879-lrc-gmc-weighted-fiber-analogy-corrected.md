---
id: HYP-8879
title: "LRC and GMC as differently typed weighted fibers"
status: >
  PARTIAL / CORRECTED BY MISTAKE-235. The Fejer-regularized LRC relation-
  lattice expansion is valid and useful as a strict-bulk diagnostic. It is not
  a fixed GMC moment, and no AP-core or LRC(14) reduction follows.
source: death-star-2026-07-21-S102; codex audit 2026-07-21
depends_on:
  - THM-2047
related:
  - THM-2022
  - THM-2051
  - THM-2058
  - MISTAKE-235
script: 04-computation/lrc_gmc2_resonance_unification_deathstar_S102.py
output: 05-knowledge/results/lrc_gmc2_resonance_unification_deathstar_S102.out
---

# HYP-8879 — corrected weighted-fiber analogy

## Exact survivor

For the strict indicator `g_delta(x)=1[||x||>delta]`, THM-2047 justifies the
Fejer-regularized identity

```text
int_0^1 product_i g_delta(v_i t) dt
  = lim_(H->infinity) sum_(k in ker_Z(v)) w_H(k) product_i hat(g_delta)(k_i).
```

This is an infinite signed lattice fiber with sinc coefficients, a fixed
threshold, and a regularization limit. Its zero mode is `(1-2delta)^n`.
Positive value certifies `M(v)>delta`; value zero conflates the tight boundary
`M(v)=delta` with the subthreshold case.

A fixed GMC moment is instead a finite affine-semigroup fiber

```text
{r in N^s: |r|=m, q dot r=0}
```

with multinomial, factorial, radial, and coefficient weights. GMC quantifies
over masses `m`. The two objects share the abstract schema

```text
Z(M,A,b,w)=sum_(x in M, A x=b) w(x),
```

but differ in monoid, grading, regularization, coefficient ring, and target
predicate. No weight- and predicate-preserving map is known.

## What failed

S102 used a raw box cutoff `|k_i|<=9` on four supports of sizes four and five,
with no tail bound; some printed values visibly differ from direct integration.
Those examples prove no thirteen-speed estimate. The repeated number `6/7`
also does not identify THM-878's primitive-class pair-overlap mean with the
Fourier coefficient `hat(g)(0)`.

The claimed AP reduction is false as stated. The Goddyn--Wong row
`{1,...,11,13,24}` is tight and non-AP, so its strict measure is zero. No
metric or radius was defined for “AP-neighborhood.” Even a unique primitive
kernel direction does not prevent cancellation: for `{1,2}` at `delta=1/3`,
the direction is `+/-(-2,1)`, yet all its multiples cancel the `1/9` zero mode.

## Correct use and next test

Use relation height, support, Fourier decay, and initial channels to schedule
experiments, but prove a uniform tail bound in the actual thirteen-speed class
before drawing an LRC conclusion. At equality, switch to THM-2058's exact
denominator packet and owner data or THM-2047's phase-height/Euler carrier.
THM-2051—not S102—is the proved sparse-relation localization.

Tournament-zeta language remains an analogy: no vertices, orientation gauge,
or weight-preserving map from balanced channels to LRC relation multiples has
been supplied.
