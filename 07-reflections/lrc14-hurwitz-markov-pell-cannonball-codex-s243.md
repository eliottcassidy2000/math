---
source: codex-2026-06-27-S243
status: synthesis after finite arithmetic sidecar scout
tags: [lrc14, hurwitz, markov, pell, cannonball, controlled-forgetting, tournament-analysis]
---

# LRC14 Hurwitz-Markov-Pell Cannonball Carrier

The strongest new connection is very concrete:

```text
Pell P5 = 29
Pell P6 = 70
Pell P7 = 169
P5 * P7 - P6^2 = 1
1^2 + ... + 24^2 = 70^2
```

The neighbors `29` and `169` are also on the fixed-coordinate-`2` Markov
branch:

```text
(2, 5, 29), (2, 29, 169), (2, 169, 985), ...
```

So the cannonball square root `70` is not floating numerology.  It is the Pell
carry between two Markov-Pell wall numbers, with the square itself one below
their product.  That is exactly the kind of hidden coordinate the LRC packet
ledger has been trying to protect.

The transfer to LRC14 is still indirect.  Hurwitz and Markov are about good
approximants; LRC wants an anti-approximation time outside all forbidden
integer neighborhoods.  The useful payload is not the scalar constant
`1/sqrt(5)` or a raw Markov number.  The useful payload is the exceptional-wall
address:

```text
continued-fraction period
Markov/Lagrange depth
quadratic unit branch
Pell carry residue
endpoint shell address
destroyed owner/route coordinate
```

This sits naturally beside HYP-3072's blindness-pair rule.  A quotient may use
Hurwitz/Markov/Pell language only if it also names what it forgot and which
sidecar restores it.  Otherwise it repeats the old mistake of calling a visible
token Sturmian when the exact word is actually Beatty address plus Pell carry
plus endpoint atom.

Two proof angles remain worth pursuing.

First, a Markov-depth sidecar for the q=7 resonance wall.  HYP-2745's
three-leg formula already contains the hidden product leg `||pq||`; HYP-2753
says the aggregate wall is precisely where two visible legs are not enough.
Markov language gives a disciplined name for that missing depth, but it must
remain attached to endpoint owner, exact scale, route, and certificate exits.

Second, a cannonball/Pell wall ledger.  Whenever a scalar token looks special,
replace it by shell address plus quadratic carry plus endpoint atom before
checking status or route purity.  This is a direct next step for the Q27
address-ledger idea and for the HYP-3074 route-state closure median interface.

Tournament Analysis used proof carriers and arithmetic sidecar types, not
runners or sequence entries.  The retained-payload gauge was transitive:
packet ledger, route-state closure, and cross-carrier portfolio outrank raw
Hurwitz/Markov/Pell/cannonball carriers because they preserve the LRC
anti-Bohr predicate and legal exits.  The arithmetic carriers are still useful,
but only after being pulled back to packet fields.
