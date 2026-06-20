---
id: HYP-2691
title: LRC14 sector-state transfer DP as the next proof router
status: OPEN
source: codex-2026-06-20-S58
tangent: T927
depends_on:
  - THM-554
  - THM-551
  - HYP-2675
  - HYP-2683
  - HYP-2684
  - HYP-2690
related:
  - HYP-2676
  - HYP-2677
  - HYP-2680
  - HYP-2681
  - HYP-2682
  - OPEN-Q-108
---

# HYP-2691 - LRC14 Sector-State Transfer DP

## Claim

The useful dynamic-programming address for the LRC14 sector route is the exact
missed-sector transfer kernel, not the row as a scalar `p0` value.  Insert the
speeds in an address order, retain the wall-atom missed-state distribution, and
use THM-554's one-sector deletion recurrence:

```text
p0(P union {e}) - p0(P)
  = sum_s meas{P misses exactly sector s and e lands in sector s}.
```

The next proof target is a transfer inequality:

```text
landing_mass(P,e)
  <= finite_address_templates(P,e) + decorrelation_error(P,e),
```

where low-state / low-growth prefixes route to AP, dyadic, cube-root, and
Ruzsa/Freiman atlases, while high-state prefixes route to Weyl/BV
decorrelation.

## Evidence

The script `04-computation/lrc14_state_transfer_dp_codex_s58.py` computes exact
prefix transfer packets over the common wall refinement and stores the run in
`05-knowledge/results/lrc14_state_transfer_dp_codex_s58.out`.

The exact theorem-level recurrence is asserted internally: after one
insertion, the missed-sector state only stays fixed or deletes one sector, and
the total closure mass equals the exact `p0` increment.

Order choice materially changes the intermediate DP size even though final
`p0` is invariant.  In the named row bank, the global order-strategy averages
were:

```text
greedy-min-support    avg_max_support=45.88  avg_area=253.00
increasing            avg_max_support=46.75  avg_area=265.50
dyadic-tower          avg_max_support=48.12  avg_area=260.88
decreasing            avg_max_support=55.50  avg_area=329.25
residue-layer         avg_max_support=55.14  avg_area=321.43
```

The strategy winner depends on the row's structure:

- AP-like rows prefer the greedy minimum-support order.
- The HYP-2671 dyadic block prefers the dyadic-tower order.
- The three-cluster true-wide row prefers greedy minimum support.
- The AP-triple phase row prefers increasing order.

This is a proof signal rather than a speed trick: the right address order is
part of the finite-resonant classifier.

## Residual-Budget Addendum

The companion script `04-computation/lrc14_transfer_residual_budget_codex_s58.py`
turns each insertion into the local one-far residual

```text
residual(P,e) = [p0(P union {e})-p0(P)] - p1(P)/7.
```

Here `p1(P)/7` is the decorrelated transition value.  The script compares the
exact residual with the THM-546 style component budget

```text
|residual(P,e)| <= (6/49) * V(P) / e,
```

where `V(P)` is the circular number of one-missed-sector runs by sector.
The bound held with exact assertions on every audited prefix and schedule.

The largest observed pressures were still well below `1`:

```text
boundary-leader increasing: max_pressure=27/65 ~= 0.415385
AP9 increasing:             max_pressure=27/65
dyadic-block dyadic-tower:  max_pressure=329/864 ~= 0.380787
direct-risk increasing:     max_pressure=943/2970 ~= 0.317508
three-cluster increasing:   max_pressure=2478259/11606400 ~= 0.213525
AP-triple increasing:       max_pressure=392675/1945944 ~= 0.201792
```

This says the DP-local version of THM-546 has large cap-level slack on the
named banks.  The hard steps are not generic high-support steps; they are
finite structured prefixes, exactly the rows that should be routed through
AP/dyadic/cube-root/Ruzsa address templates before decorrelation.

The finite AP-prefix template was isolated by
`04-computation/lrc14_ap_prefix_transfer_template_codex_s58.py`.  For

```text
d*{0,1,...,m-1} -> add d*m,   1<=m<=13,
```

the exact pressure peak occurs at `m=6`:

```text
p1=11/42, delta=31/210, residual=27/245,
V=13*d, bound=13/49, pressure=27/65.
```

Thus the largest pressure in the named residual bank is explained by a finite
AP6 append template and its integer dilates, not by a true-wide phenomenon.

Concrete exact rows from the scout:

```text
direct-risk-leader (0,4,6,8,10,12,14,15,16):
  p0=321/980, cap_9-p0=11681/70070.
  best schedules keep max_support=55; decreasing reaches 61.

dyadic-block (0,1,2,4,8,12,16,20,24):
  dyadic-tower is the Hamiltonian order-strategy leader.

three-cluster-wide (0,1,2,30,31,32,60,61,62):
  p0=141899/635376, cap_9-p0=24615911/90858768.
  greedy-min-support beats increasing on max_support 55 vs 60.
```

## Proof Route

The live LRC14 proof obligation from HYP-2675 is still

```text
span(E)>14, 8<=|E|<=12  =>  p0(E)<=cap_|E|.
```

The DP route sharpens it to a sequence of exact insertion obligations.  At
each prefix `P`, the only way to increase `p0` is through one-missed-sector
landing packets.  Thus a proof can aim to show that every insertion step is in
one of two regimes:

1. **Finite addressed regime:** state support is small or low-growth.  Then the
   row is classified by an exact finite address, such as AP, dyadic tower,
   AP-triple/cube-root phase, same-sign Ruzsa packet, or residue-private
   template.
2. **Decorrelated regime:** state support is large or high-entropy.  Then the
   one-missed-sector landing mass is bounded by Weyl/BV decorrelation and
   cannot consume the cap margin.

This reframes HYP-2683's successful coarse `state_mass` separator as a
dynamic transfer theorem: prove the separator at each insertion, not only at
the final row.

## Assumption Challenge

The tournament vertices used by this scout are insertion schedules and proof
state addresses, not runners or arcs.  The quotient preserves the LRC predicate
`p0(E)<=cap_k` and the exact insertion increments.  It destroys the continuous
location inside a wall atom except for the missed-sector state and measured
transition mass.

## Status

No LRC14 proof is claimed.  The concrete next sharp target is the transfer
inequality for one-missed-sector landing mass, with the AP6 append template as
the first finite low-state atom and a high-state Weyl/decorrelation bound for
the complement.
