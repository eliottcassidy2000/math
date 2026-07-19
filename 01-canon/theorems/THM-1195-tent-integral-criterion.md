---
id: THM-1195
title: WITHDRAWN TENT-INTEGRAL CRITERION — the exact integral computation survives, but the claimed cell upper bound has its inequality reversed; zero-cell lower envelopes are concave and lie above their maximum tent
status: WITHDRAWN / PROOF REFUTED by THM-1201 and MISTAKE-174. The implication integral g >= 1/28 => LRC(14) is not established. Numerical integral values remain exact, but “CERTIFIES” labels and the deterministic/probabilistic threshold interpretation are invalid
source: opus-2026-07-17-S388 (owner: work a new creative angle on the LRC 14 open math)
depends_on: [THM-1185 (the measure-vs-pointwise triage that motivated a hybrid), THM-1170 (the critical-point structure supplying the arrangement), THM-1120 (the tight families that bound what any criterion can do)]
scripts: 04-computation/tent_integral_opus_S388.py -> 05-knowledge/results/tent_integral_opus_S388.out; correction audit 04-computation/lrc14_tent_cell_direction_audit_codex_S78.py -> 05-knowledge/results/lrc14_tent_cell_direction_audit_codex_S78.out
---

# THM-1195 — withdrawn criterion

> **Correction (codex-S78):** the central geometric inequality below is false.
> Between consecutive combined zeros, `g` is concave and lies **above** the
> two chords through a maximum.  Thus cell area is at least `max*L/2`, not at
> most that quantity.  THM-1201 gives the proof and an exact counterexample
> inside the 13-speed tight family `{1,...,11,13,24}`.  The remainder of this
> file is retained as a record of the withdrawn argument.

## The criterion

g(t) = min_v ‖vt‖ vanishes at every k/v and rises to a local maximum in
each gap between consecutive zeros. On a cell of length L the graph lies
under the tent of height max, so its area is at most max·L/2. Hence if
*every* local max were < 1/14,

> ∫₀¹ g < (1/14) · (Σ L_i)/2 = (1/14)(1/2) = **1/28**.

Contrapositive:

> **∫₀¹ min_v‖vt‖ dt ≥ 1/28  ⟹  some local max ≥ 1/14  ⟹  LRC(14) holds.**

A single integral, computed exactly over the arrangement of zeros, peaks
and pairwise crossings.

## What it delivers

| family | ∫g | ratio to 1/28 | verdict |
|---|---|---|---|
| {1,…,13} | 0.03051 | 0.854 | fails (must) |
| {1,…,11,13,24} | 0.03047 | 0.853 | fails (must) |
| 12 random | — | min 0.921, **median 1.012**, max 1.102 | fires 9/12 |

It is **sufficient only**. It fails on the tight families — necessarily,
since their maxima equal 1/14 exactly, so the tent bound is saturated — and
on roughly a quarter of random families. It cannot prove LRC(14).

## The real content: why 1/28 appears twice

The tent bound's threshold is λ·(1/2). The independence heuristic — treat
the n values ‖vt‖ as independent uniforms on [0,1/2] — predicts

> E[min of n uniforms on [0,1/2]] = (1/2)/(n+1).

These two thresholds coincide precisely when

> **λ·(1/2) = (1/2)/(n+1),  i.e.  λ = 1/(n+1)**

which at n = 13 speeds is λ = 1/14 — **exactly the conjectured gap**. So the
observed knife-edge (median ratio 1.012, essentially 1) is not numerical
accident. It is LRC(n)'s own tightness at λ = 1/n, showing up in two
computations that know nothing about each other: a deterministic geometric
bound and a probabilistic expectation.

That explains something the programme has kept running into without naming
it. LRC(14) sits at a threshold rather than having slack **because the
conjectured λ is exactly the value at which the geometric and probabilistic
estimates cross**. Any method whose accuracy is of the order of the gap
between those two estimates will land on the boundary — which is what
S₁ = 13/7 > 1, the k = 7 arity ceiling (THM-1155), and now this criterion
all do.

## Status

Withdrawn.  The arrangement script computes the integral exactly, but the
integral is not a proved certificate at threshold `1/28`.  See THM-1201 and
MISTAKE-174.
