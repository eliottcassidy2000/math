# LRC14 Smoothing Switchboard

The creative move in this pass was to stop treating smoothing as a global
choice.  Fejer, Ramanujan, large-sieve, boundary-moment, and Kaczynski kernels
are not rivals.  They are admissible only in the packet families where their
forgotten labels are either irrelevant or explicitly routed elsewhere.

The switchboard script makes that concrete on a small but representative bank:
AP/GW, K33, petals, covering rows, q-witness rows, and the selected interval
Fejer certificates.

## Main Lesson

The right lemma shape is:

```text
packet route first, smoothing kernel second.
```

This is the same lesson as HYP-2982/HYP-2983, but more operational.  A
large-sieve bound without endpoint owners is too lossy.  A Ramanujan projector
without safe-open/boundary status is too lossy.  A Fejer negative without an
interval and packet key is too lossy.  A Kaczynski boundary phrase without a
named approach class is too lossy.

The switchboard says what each route is allowed to forget.

## Readout

The finite audit classifies `16` packets:

```text
AP/GW boundary atoms                  2
Fejer interval certificates           5
Ramanujan pre-split handoffs          6
covering/lift boundary-moment rows    2
K33 state-lift debt                   1
```

The important negative result is that no audited row needs a raw scalar
smoothing route.  Every useful kernel is attached to a labelled side channel:
endpoint zero-credit pairs, rational Fejer center/degree/sign, first strict
Ramanujan `q`, lift chart, or K33 state debt.

## Why This Helps

OPEN-Q-108 has accumulated many almost-proof routes.  The failure mode is
usually the same: a route proves something after forgetting a coordinate that
later turns out to be load-bearing.

The switchboard is a way to phrase the remaining lemma without that mistake.
It should become a familywise proof obligation:

```text
for each packet family, prove that its selected kernel either certifies a
strict safe interval or routes the packet to one named residual bucket.
```

The next useful stress test is full-bank scaling: run this switchboard over the
HYP-2963 packet bank, attach HYP-2981 Fejer records where available, attach
HYP-2979 first strict q data, and count the packet families not classified by
the current routes.
