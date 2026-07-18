---
id: THM-1085
title: THE MULTIPLICATIVE-MINIMUM BOUND IS PROVABLE BUT NOT UNIFORM, AND SUM-FREENESS DOES NOT CERTIFY — the absolute-value route CONVERGES (my log-divergence heuristic was wrong: at height H a rank-(k−1) lattice has ~H^(k−2) vectors each ~H^(−k), giving Σ H^(−2)), so |δ| ≤ Σ ∏1/(π|nᵢ|) is a legitimate bound, but its constant is NOT uniform — abs-sum·m₇·π^k ranges 2.4→19.0 at k=3 and 5.4→50.7 at k=4, growing with k, while the SIGNED constants stay ≲2 and ≲4; so m₇ alone is insufficient and the count of near-minimal vectors matters. Separately, the Schur/sum-free direction is REFUTED as a certificate: the perfectly sum-free family {1,3,…,25} has BONF5 = −4.96, the worst in a 26-family sample, because sum-freeness controls δ(S) for |S|≥3 while BONF5 is dominated by S₂/S₃
status: convergence and non-uniformity measured at k=3,4 over 20 families; the Schur/BONF5 correlation measured over 26 random 13-families (0/26 BONF5-positive) plus four named families; the refutation of "sum-free ⟹ certifiable" is exhibited by an explicit counterexample
source: opus-2026-07-17-S370 (owner: work the multiplicative minimum bound and the Schur/sum-free direction)
depends_on: [THM-1080 (m₇ and the Schur equivalence this tests), THM-1070 (the same constant-growth failure mode), THM-1065 (level 5 insufficient in general), MISTAKE-154/156 (why dilation-invariance of a predictor matters)]
scripts: 04-computation/multminimum_bound_opus_S370.py, multmin_constant_opus_S370.py, schur_bonf5_opus_S370.py -> 05-knowledge/results/
---

# THM-1085 — how far each direction actually goes

## 1. The multiplicative-minimum bound: provable, not uniform

I expected the absolute-value route to diverge logarithmically. **It does
not** — increments at doubling N fall geometrically (0.0131 → 0.0087 →
0.0054 for (2,3,5)), and the corrected heuristic agrees: at height H a
rank-(k−1) lattice has ~H^(k−2) vectors each contributing ~H^(−k), giving
Σ_H H^(−2), convergent for every k. So

> |δ(S)| ≤ Σ over full-support 7-free n ∈ Λ of ∏ 1/(π|nᵢ|)

is a legitimate, provable bound. The problem is its **constant**:

| | k=3 | k=4 |
|---|---|---|
| abs-sum · m₇ · π^k | **2.4 → 19.0** | **5.4 → 50.7** |
| \|δ\| · m₇ · π^k (signed) | 0.27 → 1.94 | 0.42 → 3.83 |

The absolute constant is not uniform and grows with k; the signed one
stays small. The gap is largest exactly where many near-minimal vectors
exist — (6,10,15) has m₇ = 60 yet an absolute constant of 19.0, because
its large pairwise gcds create a thicket of comparable vectors.

**Conclusion: m₇ alone does not suffice.** The bound needs a second
invariant counting near-minimal vectors, or it needs the cancellation
that the signed sum enjoys and absolute values discard. Note this is the
**same failure mode as THM-1070** — a valid bound whose looseness grows
with k — arriving now by a different route. That is twice that discarding
sign has cost an order of magnitude per level.

## 2. The Schur direction: a real invariant, but NOT a certificate

The Schur count (number of triples a+b=c inside V) has the property
min-speed lacked: it is **dilation-invariant**, so unlike the threshold
of MISTAKE-154/THM-1055 it is an admissible predictor of BONF5. It also
identifies the extremal family sharply — {1,…,13} has Schur count **36**
against 0–6 for everything else, and it is the unique tight family.

But it does not certify. Over 26 random 13-families:

| class | median BONF5 | BONF5 > 0 |
|---|---|---|
| Schur ≤ 1 | −0.753 | 0/4 |
| Schur ≥ 3 | −0.864 | 0/18 |

A weak trend in the predicted direction, and **0/26 positive**. The
decisive counterexample is explicit:

> **{1, 3, 5, …, 25} is perfectly sum-free** (a sum of two odds is even,
> so never in the set) **and has BONF5 = −4.96** — the worst value in the
> entire sample — while being genuinely lonely (uncovered = 0.1159).

**Why:** sum-freeness controls δ(S) for |S| ≥ 3, which is exactly what
THM-1080 established. BONF5 is dominated by S₂ and S₃. Controlling the
3-body terms while leaving the 2-body terms untouched does not move the
level-5 certificate. So "sum-free ⟹ certifiable" is **refuted**, and the
fleet should not pursue it.

Note also that AP-ness and sum-freeness are not opposites: the AP
{1, 9, 17, …, 97} (d = 8) has Schur count **0**. Only the a = d APs — the
dilates of {1,…,13} — are additively rich. THM-1080's remark that "APs
are the additively richest 13-sets" is therefore too broad and should be
read as applying to that diagonal only.

## What survives

THM-1080's structural claim stands: additive richness drives the k-body
terms, and m₇ is the leading invariant. What this file removes is the
hope that either fact yields a certificate on its own. The remaining
route is unchanged and now doubly indicated: **a bound that keeps the
signs**, since both attempts that discarded them (containment/
fragmentation in THM-1070, absolute values here) lost an order of
magnitude per level in exactly the same way.
