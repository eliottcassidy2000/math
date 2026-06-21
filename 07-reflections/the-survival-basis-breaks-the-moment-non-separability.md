# The survival basis breaks the moment non-separability

**Source:** opus-2026-06-20-S1 (Route B for HYP-2693/2607/2694, the LRC(14) compression crux).

## The wall, restated

The LRC(14) endgame had collapsed to one scalar statement: **consec maximizes the
empty-sector functional** over all bounded clusters `E`. The eight-roads synthesis
(`eight-roads-one-crux...md`) proved this lives genuinely on the JOINT distribution of
the empty-count `N_E(x) = #{inner sectors 1..6 missed by {frac(e_i x)}}`, because the
natural per-moment attack FAILS: in the binomial-moment basis `S_r = E[C(N,r)]`, consec
is **not** min-`S_1`, **not** max-`S_2` (only max-`S_4`). The alternating Bonferroni dual
makes the moment terms trade off. DEAD-END 3 of the prompt. The problem looked
irreducibly joint, calling for "a convex-order / coupling on `N_E`."

## The move: change basis from moments to survival cuts

A distribution `(p_0,...,p_6)` can be read in three bases:
- **probability** `p_t`,
- **binomial moment** `S_r = sum_t C(t,r) p_t`  (the basis where it was non-separable),
- **survival / cut** `G_b = P(N>=b) = sum_{t>=b} p_t`  (a layer-cake / Abel basis).

The Bonferroni-4 functional (THM-556) is `U4 = p_0 + p_5 + 5 p_6`. In the survival basis
it becomes, by one line of Abel summation (`p_0=1-G_1, p_5=G_5-G_6, p_6=G_6`):

> **U4(E) = 1 - G_1 + G_5 + 4 G_6.**   *(PROVED, elementary, on top of THM-556.)*

This single change of basis does what the moment basis could not: it **diagonalizes the
functional into monotone cuts with sign-aligned coefficients.** The three live cuts are
`G_1` (coeff -1), `G_5` (coeff +1), `G_6` (coeff +4), and — verified exactly on the
binding rows k=8,9,10 over banks of 19440 / 12869 / 5005 primitive shapes — consec
extremizes **each** in the direction its coefficient wants:

- **(I)** consec MINIMIZES `G_1` ⟺ consec MAXIMIZES `p_0`, the covered measure;
- **(II)** consec MAXIMIZES `G_5 = P(N>=5)`;
- **(III)** consec MAXIMIZES `G_6 = p_6 = P(N=6)`.

All three aligned ⟹ `U4(E) <= U4(consec)`. The non-separability was an artifact of the
**moment basis**; the survival basis is where the cuts cooperate. This is the same lesson
the project keeps relearning (cut-space cheap, cycle-space dear; the right basis makes the
sign cancellation manifest) — here the right basis is the layer-cake of `N_E`.

## Why it is still a two-ended squeeze, not a single order

The three cuts do not collapse to one stochastic order. `G_1` wants to be SMALL while
`G_5,G_6` want to be LARGE — opposite ends of the survival ladder `G_1>=...>=G_6`. So no
single convex/icx order works (verified: equal-mean convex order is inapplicable — E[N]
is not constant; icx-upper dominance fails 790/792). The clean characterization is a
**two-ended squeeze**: among all bounded shapes, consec is the UNIQUE one that is
simultaneously *step-light at b=1* (G_1 minimal) and *tail-heavy at every b>=2* (G_b
maximal). 0 non-consec shapes share this, k=8,9,10. That two-ended squeeze IS the joint
extremality the moment basis could only see as a tradeoff.

## The tournament reading

Via the dictionary (HYP-2605/R4), `N_E(x)` large ⟺ the round tournament `T(x)` is deeply
hierarchical (orbit concentrated, big score spread). So `G_b(E)=P_x(T(x)` is at least
`b`-deep hierarchical`)` — the survival ladder is the **filtration of x by tournament
hierarchy depth**. The squeeze says: consec spends the least x-measure being *merely*
hierarchical yet the most x-measure being *deeply* hierarchical. A surprising refinement:
consec does NOT minimize *mean* `c3(T(x))` — it sits near the hierarchy-poor (regular-like)
end on average, and buys its extremality entirely in the **rare deep tail** (max G_5, G_6).
The AP concentrates its coherence into rare maximally-transitive instants rather than
spreading it. Lonely-runner extremality is a tail phenomenon, not an average one.

## Honest scope

PROVED: the layer-cake identity (the reduction). VERIFIED, not proved: the three cut
inequalities, exactly on the binding rows k=8,9,10. The certificate is NOT universal:
cut (I) ("consec maximizes p0") fails at k=12 (`[0..10,12]` beats consec on U4), but k>=11
is closed by the THM-535 floor route, so consec need not win there. Two framings die
cleanly: `p0` is **not** a function of the difference-multiset (so no difference-multiset
majorization), and `U4(consec_8)=0.480 > cap_8=0.381` (so U4-extremality alone cannot
close k=8 — the cap is carried by `p0` directly, or by the tighter L_y dual, which itself
has no clean survival certificate). The honest remaining task is the same as HYP-2697: a
cone/coupling proof of the two-ended squeeze on k=8,9,10. The contribution here is the
**basis** in which that proof should be written.
