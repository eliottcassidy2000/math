---
id: THM-4439
title: "All-node two-jet metric precision by terminal clusters"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 fourth research wave
depends_on:
  - THM-4435-four-node-metric-blindness-and-universal-hermite-precision
proof: 05-knowledge/results/hermite-terminal-cluster-precision-overnight-hexagon-sep05.md
script: 04-computation/hermite_metric_loss_probe_overnight_hexagon_sep05.py
independent_script: 04-computation/hermite_terminal_precision_referee_overnight_hexagon_sep05.py
independent_output: 05-knowledge/results/hermite_terminal_precision_referee_overnight_hexagon_sep05.out
---

# THM-4439 -- All-node two-jet metric precision by terminal clusters

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[full proof and independent full-Smith audit](../../05-knowledge/results/hermite-terminal-cluster-precision-overnight-hexagon-sep05.md)
are part of this theorem. No external priority is claimed.

Let X be any finite set of n>=2 distinct integers and p any prime. Put
S_x=sum_(y!=x)v_p(x-y). A terminal cluster C of depth f is a nonsingleton
ball intersection C={y in X:v_p(y-x)>=f} whose distinct points all have
pairwise depth exactly f. Its children in the distance tree are all leaves;
its size m is between2 and p, and S_x is constant on C, denoted S_C.

The complete values/first-Hasse-derivatives observer on polynomials of
degree<2n has largest p-Smith exponent

```text
L_p=max_(terminal C) [2S_C+max(0,f_C-[|C|=p])].
```

For n=1 it is zero. For every N>=1, data modulo p^(N+L_p) determine all
coefficients modulo p^N, and one less digit fails uniformly. Thus the
sharp two-jet loss depends only on the valuation tree at **every** node
count. Partial observations, moving modules and higher jets are excluded.

The exact inverse from [THM-4435](THM-4435-four-node-metric-blindness-and-universal-hermite-precision.md)
gives candidates2S_x and2S_x-v_p(2q_x), q_x=sum_(y!=x)1/(x-y).
A nonterminal leaf can be replaced by a deeper sibling, increasing2S+f
by at least three. At a terminal cluster simultaneous cancellation loses
at most one digit. Therefore a nonterminal leaf cannot control the maximum.

At odd p with m<p children, the normalized root polynomial G has G''
nonzero modulo p, so some reciprocal sum is a unit. With m=p children,
G is U^p-U modulo p; G''/p has degree p-2 with unit leading coefficient.
The outside terms contribute only a constant modulo p after first division.
They cannot cancel that nonconstant polynomial at all p residues: exactly
one digit is lost. At p=2 the two-child reciprocal sum is a unit, and the
ordinary derivative's factor two instead supplies the one-digit loss.
Maximizing over the complete cluster before discarding unit data is essential.

Two nodes at gap d recover3v_p(d)-min(v_p(d),v_p(2)); depth-zero and
singleton cases are retained. THM-4435's isometric four-node families still
have different intermediate factors for every dyadic depth>=3, while both
largest losses are7e+5. The full metric-only partition remains REFUTED.

Primary exact tests cover220,616 complete head rows and60,000 seeded
high-unit isometries at six primes with up to30 nodes. The independent
referee gives7,944 fresh reciprocal/tree comparisons and176 literal integer
Smith comparisons, including full terminal p-clusters with nearby outsiders,
signed/deep nodes, n=1/2, and the metric twins. Full proof and boundary
reviews passed. The unbounded theorem is the cluster argument, not a
finite extrapolation. General higher-Hasse-jet precision remains OPEN.
