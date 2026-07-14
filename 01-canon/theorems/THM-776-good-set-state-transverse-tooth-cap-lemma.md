---
id: THM-776
title: Good-set-state transverse-tooth cap lemma
status: PROVED (elementary circle geometry plus THM-755)
source: codex-2026-07-14-S2
depends_on:
  - THM-755   # capped-envelope criterion
related: [THM-752, THM-761, THM-766, THM-775, HYP-6815, HYP-6830]
---

# THM-776 - Good-set-state transverse-tooth cap lemma

## Statement

Let `B` be a nonempty finite set of distinct positive integer speeds. Write

```text
mu  = |G'_B|,
r_B = number of components of G'_B.
```

For a positive integer `N` not in `B`, put `P_N=B union {N}`, and write
`mu_N=|G'_{P_N}|` and `r_N=r(P_N)`. Then

```text
mu_N >= 6mu/7 - 2r_B/(7N),              (1)
r_N  <= N+r_B.                           (2)
```

This is the invariant-level form: the state `(|G'_B|,r_B)` controls a
transverse frequency insertion without retaining every endpoint.

There is a useful marked-window form. If `G'_B` contains a connected closed
circle interval `J` of length `L>0`, and `S=sum(B)`, then

```text
mu_N >= 6L/7 - 2/(7N),                  (3)
r_N  <= N+S.                             (4)
```

If `mu>0`, then for any integer `a>=2` with all speeds in
`B union {N,aN}` distinct,

```text
r_N/(aN mu_N)
 <= (N+r_B)/(aN(6mu/7-2r_B/(7N)))       (5)
```

whenever the denominator is positive. In particular,

```text
limsup_(N->infinity) r_N/(aN mu_N) <= 7/(6a mu).   (6)
```

The marked-window version replaces `(mu,r_B)` in (5)-(6) by `(L,S)`.

There is a completely rational THM-755 threshold. Put `q=333/106<pi`. If

```text
q*a*(6mu/7) > 1                                                (7)
```

and

```text
N > r_B*(q*2a/7+1)/(q*6a*mu/7-1),                             (8)
```

then the peel `v=aN` satisfies

```text
pi*(aN)*mu_N > r_N.
```

Alternatively, one marked safe interval gives the sufficient conditions

```text
q*a*(6L/7) > 1,
N > (q*2a/7+S)/(q*6aL/7-1).                                  (9)
```

In either case THM-755 closes `P_N union {aN}` at the `1/14` threshold.
Covering and primitivity are not analytic hypotheses. To call the result a
13-speed LRC14 row, additionally require `|B|=11` and check those arithmetic
properties.

## Proof

The `N`-runner's strict danger set is the union of `N` disjoint open teeth,
each of length `1/(7N)`, centred at `k/N`.

First decompose the positive-length part of `G'_B` into closed circle intervals
`J_i` of lengths `L_i`; there are at most `r_B` of them and
`sum_i L_i=mu`. A tooth can meet `J_i` only if its centre lies in the circular
`1/(14N)`-neighbourhood of `J_i`. That expanded interval has length
`L_i+1/(7N)`. A mesh-`1/N` grid meets it in fewer than `NL_i+2` points. Bounding
each intersection by a full tooth gives

```text
mu_N >= sum_i (L_i-(NL_i+2)/(7N))
     >= 6mu/7-2r_B/(7N),
```

which proves (1). Applying the same argument only to a named interval `J`
proves (3).

Remove the `N` open teeth from `G'_B` one at a time. Removing one connected
open arc from a finite union of closed arcs and points can increase its number
of components by at most one: it either trims, deletes, or splits one retained
component. Hence `r_N<=r_B+N`, proving (2). Independently, the union of all
danger arcs for `P_N` has at most `sum(P_N)=N+S` components, so its complement
has at most that many components; this proves (4). These bounds count singleton
boundary components as well as nondegenerate interval components.

Division gives (5), and the limit gives (6). For the full-state exact cap,
multiply (1) by `q*aN` and use (2):

```text
q*aN*mu_N-r_N
 >= (q*6a*mu/7-1)N-r_B*(q*2a/7+1).
```

Conditions (7)-(8) make this strictly positive. The same calculation with
(3)-(4) gives (9). Since `q<pi`, either inequality implies
`pi*aN*mu_N>r_N`, exactly THM-755's strict capped-envelope criterion. ∎

## Iterated transport and the order gauge

The state bound composes. Let `B_0=B`, and successively adjoin distinct positive
frequencies `N_1,...,N_k` not already present. Put

```text
B_j  = B_{j-1} union {N_j},
mu_j = |G'_{B_j}|,
r_j  = r(B_j).
```

Repeated application of (1)-(2) gives

```text
r_j <= r_0 + sum_(h<=j) N_h,                                  (10)

mu_j >= (6/7)^j mu_0
        - (2/7) sum_(i=1)^j (6/7)^(j-i)
            * (r_0+sum_(h<i) N_h)/N_i.                        (11)
```

Indeed, (10) is immediate by induction. Substituting its `j-1` case into the
one-step mass bound and iterating proves (11).

The final set `B_k` is independent of the insertion order, but the certified
lower bound (11) need not be. Therefore every permutation gives a valid lower
bound for the same `mu_k`, and their maximum is valid as well.

This maximum has a canonical order. On an enclosure state `(m,r)`, write

```text
T_N(m,r)=(6m/7-2r/(7N),r+N).
```

For `x<y`, direct subtraction of the mass coordinates gives

```text
[T_y(T_x(m,r))]_mass - [T_x(T_y(m,r))]_mass
  = (2/7)*(y-x)/(xy)*(x+y+r/7) > 0.                 (12)
```

The component coordinate is the same after either two-step order, so every
later transition preserves this strict improvement, multiplied by a positive
power of `6/7`. An adjacent-exchange argument therefore shows that (11) is
uniquely maximized by inserting the distinct frequencies in increasing order.
In a fixed-core four-far chart, a one-peel capped-envelope search starts with 24
conceivable gauges,

```text
4 choices of peel * 3! insertion orders = 24
```

but (12) reduces them to four canonical state certificates, one per peel. This
is a finite proof calculus, not a finite classification of the chart: all four
lower bounds may be nonpositive, and the compression deliberately forgets
endpoint owners and correlations.

## THM-775 specialization

For

```text
B={1,...,9,15,110},  L=1/1540,  S=170,  a=1092,
```

the marked-window threshold (9) has exact crossing

```text
11734415/9278 < 1265.
```

Thus THM-776 supplies the elementary infinite tail in THM-775. THM-775's
divisor-packet statement, fragmentation lower bound, and 176-prime finite base
remain separate exact inputs.

## Structural meaning

A high frequency can create `Theta(N)` endpoint walls inside an old good set,
so raw fragmentation can diverge with no divisor packet. But its danger duty
cycle removes asymptotically only one seventh of each retained interval, while
a proportional peel supplies `Theta(N)` sampling resolution. The minimal state
used by this proof operation is

```text
(safe mass, safe component count, new wall rate, named peel rate).
```

No single coordinate in that tuple is a compactness invariant. The quotient
`r_N/(aN mu_N)`, rather than `r_N`, is theorem-facing.

## Assumption challenge and Tournament Analysis

The useful vertices are base safe components, new danger teeth, and the named
peel obligation, not runners. The exact carrier is the bipartite
component-tooth incidence with lengths and peel rate as sidecars. For the
upper bounds above, it compresses further to `(|G'_B|,r_B,N,a)` without loss.
For several insertions, the remaining-frequency multiset is also retained;
insertion order is a proof gauge rather than object data, and (12) fixes its
canonical increasing representative.

A runner-pair tournament destroys how many teeth land in one component and
cannot reconstruct (1). A diagnostic tournament may orient candidate carriers
by predicate retention, with ties resolved

```text
component-tooth incidence
  -> (safe mass,component count,wall rate,peel rate)
  -> owner-coloured endpoint word
  -> peel-relative scalar
  -> divisor-support profile
  -> raw runner tournament.
```

This Hamiltonian path is a declared gauge, not the proof. The proof is the
metric incidence inequality.

## Scope

THM-776 proves eventual closure only when positive safe mass persists and the
peel rate is large enough relative to it. It does not give a uniform lower
bound for `|G'_B|` over arbitrary cores, classify sublinear peels, or prove the
global HYP-6830 splice. It makes the remaining obstruction more precise:
simultaneous collapse of safe mass relative to component and peel rates.
