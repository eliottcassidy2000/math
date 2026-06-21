---
id: HYP-2736
title: LRC14 L7 torus-line discrepancy has a 7pq integer-grid defect form
status: OPEN sharp refinement; exact integer-grid reduction plus q<=160 scout
source: codex-2026-06-21
depends_on:
  - HYP-2730
  - HYP-2729
  - THM-562
related:
  - HYP-2731
  - HYP-2733
  - HYP-2734
  - HYP-2726
  - THM-546
  - OPEN-Q-108
---

# HYP-2736: Torus-Line Discrepancy Integer Grid

## Claim

Incoming HYP-2730 already supplies a usable elementary L7 tail proof via
`D_{p,q} <= 14/p`, and HYP-2733 proves the apex-prime zero law
`D_{p,q}=0` iff `7 | p*q`.  HYP-2736 is the sharper arithmetic refinement:
the observed `sup D_{p,q}*q = 12/7` target can be converted into a concrete
integer inequality.

For coprime `1 < p/q <= 43/20`, let `c_ij` be the number of subintervals in
the common breakpoint grid of size `7*p*q` on which

```text
v -> (sector(qv), sector(pv))
```

lands in the 7-by-7 cell `(i,j)`.  Then the torus-line cell discrepancy is

```text
D_{p,q} = sum_ij |49*c_ij - 7*p*q| / (343*p*q).
```

Thus the sharp `D*q` tail constant

```text
D_{p,q} <= 12/(7q)
```

is exactly the integer defect inequality

```text
sum_ij |49*c_ij - 7*p*q| <= 588*p.
```

This is not needed to close the current L7 tail, but it is the cleanest
formalization target for the sharp constant and a bridge to the HYP-2733
mod-7 zero structure.

## Exact Evidence

Script:
`04-computation/lrc14_torus_line_discrepancy_integer_grid_codex_20260621.py`.

Stored output:
`05-knowledge/results/lrc14_torus_line_discrepancy_integer_grid_codex_20260621.out`.

The script verifies the integer-grid formula against the older Fraction
breakpoint computation on the `q<=12` bounded-ratio atlas (`52` ratios), then
scans all `8977` primitive ratios with `q<=160`.

Exact scan summary:

```text
max D*q = 12/7 at p/q=3/2
D <= 12/(7q) violations in scan: 0
largest q with D >= 21/100 is 4
worst q>=5 row: p/q=9/5, D=32/315, D*q=32/63
```

The small-denominator table shows the integer defect inequality is tight only
at `3/2` among the audited leaders:

```text
q=2, p=3: defect_sum=1764 = 588*p
```

For `q>=15`, the worst observed residue class is `(p mod 7, q mod 7)=(4,3)`,
witnessed by `18/17` with `D*q=22/63`.  When `p` or `q` is `0 mod 7`, the
discrepancy is zero in the scan, suggesting that the apex-prime grid splits
the problem into a finite mod-7 block proof plus a monotone quotient argument
inside each nonzero residue class.

## Proof Target

Prove the integer defect inequality

```text
sum_ij |49*c_ij - 7*p*q| <= 588*p
```

for every coprime bounded-ratio pair.  Together with exact finite checks for
`q<=8`, this gives a sharper alternate HYP-2730 non-resonant L7 tail:

```text
q>=9: D_{p,q} <= 12/(7q) <= 4/21 < 21/100 <= margin.
```

The proof should avoid a lossy analytic discrepancy black box if possible,
but it should be treated as a refinement of the already-proved `14/p` route,
not as a prerequisite for the current L7 closure.
The integer-grid form suggests three possible routes:

1. Partition by `(p mod 7, q mod 7)` and prove each residue block by a floor
   sum identity.
2. Interpret `c_ij` as a Sturmian/Christoffel word balance count and use the
   one-dimensional balanced-word discrepancy theorem.
3. Compress the 49 cells into seven row-conditionals; each row is a rotation
   average of a length-`p/(7q)` interval over a `q`-point orbit.

## Assumption Challenge

The useful tournament vertices here are not runners, arcs, or even ratio
channels.  They are residue-pair proof obligations `(p mod 7, q mod 7)`.
This quotient preserves the exact discrepancy constant needed by HYP-2730,
but it destroys the base set `B` and all Delsarte packet data.  That is
acceptable only because HYP-2730 has already separated the universal
torus-line discrepancy bound from the base-dependent bounded function
`0 <= g_B(i,j) <= 1`.
