---
id: THM-1147
title: Exact endpoint determinant law for two combs, with the accessible-index truncation guardrail
status: CORRECTED.  The general edge identity and its n=m+1 linear specialization are PROVED.  The claim that a complete linear descent supplies the actual four-comb max/mean nonuniformity is REFUTED: a core interval can expose only a short near-constant slice, and some components use arbitrary endpoint indices or core boundaries.  The exact THM-1148 row has actual max/mean 638/573<4/3 while easily satisfying the sharp r=5 gap target.  Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.70; owner: work the arithmetic alignment bias)
depends_on:
  - THM-1141    # which identified nonuniformity as the lever and asked for its arithmetic cause
  - THM-1140    # the four-comb reconnaissance
related: [THM-1097, THM-1094, THM-1148, MISTAKE-171]
script: 04-computation/alignment_arithmetic_kps_S128c70.py (+ .out)
---

# THM-1147 — the exact two-comb gap law

This note isolates one exact coordinate in the labelled endpoint word.  It is
useful pair-edge telemetry, not a statistic of the full four-comb survivor
word.  THM-1148 and MISTAKE-171 record the failure of the former transfer to
a `4/3` actual-mean inequality.

## (I) The conditional endpoint law

A surviving component is bounded either by tooth edges or by a core boundary.
If its left endpoint is the right edge of tooth `m` of comb `a` and its right
endpoint is the left edge of tooth `n` of comb `b`, the exact general identity is

> **L(a,b;m,n) = n/b − m/a − 1/(14a) − 1/(14b)
> = (an−bm)/(ab) − (a+b)/(14ab).**

The originally displayed law is the specialization `n=m+1`, `a<b`, `j=m`,
and `d=b-a`.  Not every component has that index pattern.  The standing worst
case — core [1,3,5,6,7,8,11,12], killers 371/374/377/379 — does:

> [1373/5278, 1385/5306], and 5278 = 14·377, 5306 = 14·379,

i.e. it runs from the **right edge of tooth j = 98 of comb 377** to the **left edge of tooth
j = 99 of comb 379**.  In this specialization:

> **L(a,b,j) = (j+1)/b − j/a − 1/(14a) − 1/(14b) = (a − j·d)/(a·b) − 1/(14a) − 1/(14b).**

Check: (377 − 98·2)/(377·379) − 1/5278 − 1/5306 = **127/142883**, and the measured longest
component is **127/142883** exactly.

## (II) Linear in `j` on one specialized branch

The gap descends linearly from ≈ 1/a at j = 0 to zero at j ≈ a/d:

| (a,b) = (371,379), d=8 | j=0 | 5 | 10 | 20 | 40 |
|---|---|---|---|---|---|
| usable gap × b | 0.856 | 0.748 | 0.640 | 0.424 | ≈0 |

| (a,b) = (371,372), d=1 | j=0 | 50 | 150 | 300 | 370 |
|---|---|---|---|---|---|
| usable gap × b | 0.857 | 0.722 | 0.453 | 0.048 | ≈0 |

Small `d` stretches the complete algebraic branch over many more indices.  This
does not say which indices a fixed core interval exposes, or which candidate
gaps survive the other combs.

## (III) Accessible-index truncation refutes the max/mean inference

The mean of a complete linear descent is about half its maximum.  A core window
can expose only a short, almost constant slice.  THM-1148's exact legal row is

```text
P={1,...,8},  I=[1/14,13/112],  K=(108,109,110,111).
```

Its five final gap lengths are

```text
319/55944, 305/55944, 291/55944, 277/55944, 13/3024.
```

Consequently

```text
mu_actual=191/37296,
L_max/mu_actual=638/573<4/3,
7*111*L_max=319/72>1.
```

The first four gaps obey the specialized law for `(a,b)=(108,111)`:

```text
L_j=(1293-42j)/167832,       j=8,9,10,11.
```

Their own max/mean is `319/298<4/3`.  The full positive branch runs through
`j=30`, but the core sees only this four-index slice.  All four combs remain
essential despite the visible endpoints belonging to 108 and 111.  The
identity is correct; the proposed universal nonuniformity consequence is false.

Three means had been conflated: the complete pair-branch mean, the actual
component mean `mu_actual`, and the uniform-interleaving benchmark
`m0=3/(7 sum k_i)`.  Put `D=L_max/mu_actual` and `B=mu_actual/m0`.  Then

```text
7 k4 L_max>1  iff  D B > (sum k_i)/(3 k4).             (1)
```

The counterexample wins through a large baseline factor `B`, not through
dispersion `D`.  An analytic tail must control their product.

## (IV) The predictor I proposed is refuted

I expected small j·d/k to predict large surviving gaps. It does not. Measured over 160
clustered quadruples:

| | median j*·d_min/k₄ |
|---|---|
| worst half by L·k₄ | **0.345** |
| best half by L·k₄ | **0.758** |

The correlation is **inverted**. The proxy is mis-specified: it uses k₄ and d_min, while the
law involves the actual bounding pair (a,b) and *its* difference, which need not involve k₄
at all.  The endpoint identity and branchwise linearity stand; both the summary
statistic and the former universal max/mean inference do not.
Recording it rather than quietly dropping it, since an inverted correlation is the kind of
thing that reads as noise if you do not say you predicted the opposite.

## Named next

- Extract the actual endpoint-owner word `(a,m;b,n)` plus core-boundary tokens
  for every final component.  The determinant `an-bm`, not a presumed common
  `j`, is the exact discrete coordinate.
- Bound the coupled product `D*B` in (1).  Multiplier/mass/overlap arguments
  naturally control `B`; endpoint dispersion or autocovariograms control `D`.
- Treat obstruction by the other combs and accessible-index truncation as part
  of the object.  Six independent latent full descents are not a faithful model.
