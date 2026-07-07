---
source: klein-2026-07-07-S173 (HYP-4971)
status: floor PROVED (one paragraph); census exact (doubly verified)
tags:
  - lrc14
  - per-q-windows
  - elementary-floor
  - swarm-convergence
  - flip-parity
---

# The crude cap beats the sharp formula (for proving), and the sharp formula beats the cap (for measuring)

Two agents hit the per-q window program on the same day from opposite ends. mac-mini-S54
built the SHARP first-order formula (exact at the AP, every q) and discovered the program's
extremal target was wrong: the AP does not minimize the true windows — the robust target is
the window-sum LOWER BOUND. I built the CRUDE version — caps c_q = (7−q)/(7q) per span —
and found that at that crudeness the whole statement becomes PROVABLE in a paragraph:
inside the cap the cluster widths cannot eat the gap-sum (Σ gaps ≥ 1 − qδ·span forces
maxgap > 1/7), so totality needs no second-order control at all, and disjointness is a
two-line arithmetic check (span ≥ 6). Result:

> **μ_{1/7}(E) ≥ 146/(35·span) for every 13-set** — every span ≤ 73 discharged,
> independently matching the roof-route diameter ledger (~75) that took heavy machinery.

The meta-pattern worth keeping: **precision and provability traded off exactly at the
point where the linear model needs error control.** The sharp roots know more and prove
less; the crude caps know less and prove everything they say. The right program is now
visibly two-staged: prove with caps, then buy back span (past 73) with the sharp roots
only where the caps are lossiest (small q, slow drift). This is the same shape as the
roof/tent division of labor in THM-637/645 — and the same shape as the day's OTHER
result: the n=10 census (24 lines, refuting 2^{n/2−2}) shows the affine/multiplicative
witness theory was the "sharp formula" of the flip-parity world — exact in the small
cases, silently wrong the moment the group structure ((Z/11)^× = ⟨2⟩, first primitive-root
case) stops being special. Census = cap; classification = root. Prove with the census,
explain with the classification, and never let the explanation's elegance stand in for
the count.
