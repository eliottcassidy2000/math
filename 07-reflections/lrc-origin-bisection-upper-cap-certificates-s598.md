---
source: codex-2026-06-03-S598
status: proved aggregate certificate + computation; refines THM-398 Lemma H and Lemma F
tags: [LRC, Cprime, n14, cap-face, endpoint-owner, origin-bisection, THM-398, HYP-2140, tournament-analysis]
---

# Origin-bisection laws become upper cap certificates

The user's phrase landed cleanly inside the Cprime cover problem.

For `S=S' union {v=nw}`, the old dual cap face used the whole danger capacity in
one primary `n`-clock cell:

```text
total cap capacity = 2/n^2.
```

But each danger cap has an origin: its center `j/(nw)`.  That center bisects the
cap.  If a small-owner congruence pins a left endpoint of a safe component to
the center, then any cover must put the component in the upper half of that cap.
There are `w` upper half-caps per primary cell, each of length `1/(n^2 w)`, so
the upper capacity is only

```text
upper capacity = 1/n^2.
```

Right endpoints give the lower-cap dual.

## Certificate

Let `C=(a,b)` be a component of `G(S')`.

If the left endpoint has small owner `u<n`,

```text
a = (k n + 1)/(n u),
```

and the `v=nw` cover congruence pins it to a center,

```text
w(k n + 1) == 0 mod u,
```

then coverability forces `C` into the upper half-cap starting at `a`.  Therefore,
if the total length of all such upper-forced components in one primary cell
exceeds `1/n^2`, the cover is impossible and `S` is loose.

This is Lemma F aggregated by side.  Lemma F says one pinned component must fit
one half-radius.  S598 says many pinned components in the same cell must fit the
whole one-sided half-cap budget.

## Computation

`04-computation/lrc_origin_bisection_upper_caps_s598.py` audits deterministic
multiple-of-`n` samples.  After the total-cell cap test fails, the bisection cap
routes:

```text
n=6:  289
n=8:  412
n=10: 143
n=12: 84
n=14: 23
```

The named calibration row is the point:

```text
apex_AP_replace_13_n14
total cell surplus = -5/1176
upper cap surplus  =  1/1176
```

So it is not a total-mass overload, but it is an upper-cap overload once the
origin-bisection law is used.

## Rebase integration: HYP-2142

After this note was committed locally, origin added Opus S595/HYP-2142: the
large-owner cover residual is a 2x2 determinant / bounded CRT problem.  This
fits the layer stack exactly.  Total-cell caps and origin-bisected one-sided
caps are aggregate capacity exits; rows under all those caps should then pass
to the HYP-2142 determinant language, where the remaining question is whether
one bounded multiplier `w` can satisfy all owner/slack blocks at once.

## Tournament Analysis

Vertices are proof certificates: total cell cap, upper cap, lower cap, Lemma F,
Lemma E, Lemma C, Bprime, and residual.  The pair observable is
`(certificate_rank, sampled_route_count, name)`, with the switch favoring
stronger aggregate certificates before local component criteria.  The sample
fingerprint is transitive: score histogram `{0:1,...,7:1}`, no directed
3-cycles, singleton SCCs, and one Hamiltonian path.

Assumption challenged: cap vertices are not just primary cells.  A proof may
need side-labelled half-cells whose labels come from endpoint-owner congruence
laws.  The quotient preserves "forced into upper/lower capacity"; it destroys
exact phase order inside a cell, so mixed fibers must be lifted by endpoint
owner, side, and residue labels.
