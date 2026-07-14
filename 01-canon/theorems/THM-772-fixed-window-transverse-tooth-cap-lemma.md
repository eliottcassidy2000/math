---
id: THM-772
title: Fixed-window transverse-tooth cap lemma
status: PROVED (elementary circle geometry plus THM-755)
source: codex-2026-07-14-S2
depends_on:
  - THM-755   # capped-envelope criterion
related: [THM-752, THM-761, THM-766, THM-771, HYP-6815, HYP-6830]
---

# THM-772 - Fixed-window transverse-tooth cap lemma

## Statement

Let `B` be a nonempty finite set of distinct positive integer speeds. Suppose
its closed `1/14`-good set contains a proper circle interval `J` of length
`L>0`, and put `S=sum(B)`. For an integer `N` not in `B`, define

```text
P_N = B union {N}.
```

Write `mu_N=|G'_{P_N}|` and let `r_N` be the number of components of that good
set. Then

```text
mu_N >= 6L/7 - 2/(7N),                 (1)
r_N  <= N+S.                            (2)
```

Consequently, for every integer `a>=2` such that `aN` is not in `P_N`,

```text
r_N/(aN mu_N)
 <= (N+S)/(aN(6L/7-2/(7N))).           (3)
```

In particular,

```text
limsup_(N->infinity) r_N/(aN mu_N) <= 7/(6aL).       (4)
```

There is also a completely rational THM-755 threshold. Put `q=333/106<pi`.
If

```text
q*a*(6L/7) > 1                         (5)
```

and

```text
N > (q*2a/7+S)/(q*6aL/7-1),            (6)
```

then the capped-envelope inequality holds for the peel `v=aN`:

```text
pi*(aN)*mu_N > r_N.
```

Hence THM-755 closes `P_N union {aN}` at the `1/14` threshold. Covering and
primitivity are not hypotheses of the analytic lemma; when they hold, this is
a terminal certificate in the LRC14 covering branch.

## Proof

The `N`-runner's strict danger set is the union of `N` disjoint open teeth,
each of length `1/(7N)`, centred at `k/N`. A tooth can meet `J` only if its
centre lies in the circular `1/(14N)`-neighbourhood of `J`. This expanded
interval has length `L+1/(7N)`. A grid of mesh `1/N` meets any proper circle
interval of length `X` in at most `NX+1` points, and therefore fewer than

```text
N(L+1/(7N))+1 = NL+8/7 < NL+2
```

centres can contribute. Bounding each intersection by the full tooth length
gives

```text
mu_N >= L-(NL+2)/(7N)
     = 6L/7-2/(7N),
```

which is (1). No independence or equidistribution is used.

For (2), runner `p` contributes `p` open danger arcs. Their union therefore
has at most `sum(P_N)=N+S` components. Whenever (1) is positive its complement
is nonempty, and the complement of `m` open arcs on the circle has at most `m`
components. Thus `r_N<=N+S`. Division gives (3), and taking the limit gives
(4).

For the exact cap, multiply (1) by `q*aN` and use (2):

```text
q*aN*mu_N-r_N
 >= (q*6aL/7-1)N-(q*2a/7+S).
```

Conditions (5) and (6) make the right-hand side strictly positive. Since
`q<pi`, this implies `pi*aN*mu_N>r_N`, exactly THM-755's strict capped-envelope
criterion. This proves the claim. ∎

## Structural meaning

A new high frequency can create `Theta(N)` endpoint walls inside an old good
component, so raw fragmentation can diverge with no divisor packet. But its
danger duty cycle removes asymptotically only one seventh of the retained base
window, while a proportional peel supplies `Theta(N)` sampling resolution.
The theorem preserves the three quantities that must be compared together:

```text
base-window mass, endpoint-wall rate, peel rate.
```

None of the three alone is a compactness coordinate. Formula (3), rather than
`r_N`, is the theorem-facing quotient.

## Assumption challenge and Tournament Analysis

The useful vertices are not runners. Consider base safe components, new danger
teeth, and the named peel obligation. The exact carrier is the bipartite
component-tooth incidence with interval lengths and the peel rate as sidecars.
It preserves the measure loss, component upper bound, and cap predicate.

A runner-pair tournament destroys how many teeth land in one base component
and cannot reconstruct (1). A diagnostic tournament may orient candidate
carriers by predicate retention, with ties resolved

```text
component-tooth incidence
  -> owner-coloured endpoint word
  -> peel-relative scalar
  -> divisor-support profile
  -> raw runner tournament.
```

This path records a declared gauge, not a proof tournament. The theorem itself
is the metric incidence inequality above.

## Scope

THM-772 proves eventual closure only when a fixed positive safe window and a
large enough proportional peel are both present. It does not give a uniform
lower bound for `|G'_B|` over arbitrary cores, classify families with a
sublinear peel, or prove the global peel-relative splice of HYP-6830. It turns
one broad transverse degeneration into a terminal analytic face and makes the
remaining obstruction more precise: simultaneous collapse of base-window
mass relative to endpoint and peel rates.
