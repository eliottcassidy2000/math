---
id: THM-1085
title: THE MULTIPLICATIVE-MINIMUM SERIES CONVERGES ABSOLUTELY, BUT ITS m7-NORMALIZED SIZE VARIES SHARPLY, AND SUM-FREENESS DOES NOT CERTIFY — THM-1092 proves |δ| ≤ Σ ∏1/(π|nᵢ|) by a coefficient-uniform dyadic majorant; experiments show abs-sum·m₇·π^k ranging 2.4→19.0 at k=3 and 5.4→50.7 at k=4 while the signed values stay much smaller, so m₇ alone is not a sharp predictor and near-minimal-vector counts or cancellation matter. Separately, {1,3,…,25} is sum-free but has BONF5 = −4.96, refuting sum-free as a certificate
status: absolute convergence is proved in THM-1092 for every support and every nonzero coefficient vector; the m7-normalized variation is measured at k=3,4 over 20 families, not proved unbounded; the Schur/BONF5 correlation is measured over 26 random 13-families plus four named families; the implication "sum-free ⟹ certifiable" is refuted by an explicit counterexample
source: opus-2026-07-17-S370 (owner: work the multiplicative minimum bound and the Schur/sum-free direction)
depends_on: [THM-1092 (absolute convergence), THM-1080 (m₇ and the Schur equivalence this tests), THM-1070 (the same constant-growth failure mode), THM-1065 (level 5 insufficient in general), MISTAKE-154/156 (why dilation-invariance of a predictor matters)]
scripts: 04-computation/multminimum_bound_opus_S370.py, multmin_constant_opus_S370.py, schur_bonf5_opus_S370.py -> 05-knowledge/results/
---

# THM-1085 — how far each direction actually goes

## 1. The multiplicative-minimum bound: convergent, but not sharp from `m7` alone

I expected the absolute-value route to diverge logarithmically. **It does
not.**  The original doubling experiment (0.0131 → 0.0087 → 0.0054 for
(2,3,5)) suggested convergence, and THM-1092 now proves it: on a dyadic
shell, choose a largest coordinate, solve for it from the resonance equation,
and sum the other reciprocal coordinates harmonically.  The shell is at most
`O_k(2^(-r)(r+1)^(k-1))`, a summable bound uniform in the nonzero speed
coefficients.  Therefore

> |δ(S)| ≤ Σ over full-support 7-free n ∈ Λ of ∏ 1/(π|nᵢ|)

is a legitimate, provable bound.  The issue is not convergence but how much
coefficient-sensitive information is lost when signs are removed:

| | k=3 | k=4 |
|---|---|---|
| abs-sum · m₇ · π^k | **2.4 → 19.0** | **5.4 → 50.7** |
| \|δ\| · m₇ · π^k (signed) | 0.27 → 1.94 | 0.42 → 3.83 |

Across this sample the `m7`-normalized absolute value varies sharply and is
larger at `k=4`; this is evidence, not a proof of unbounded growth in `k`.
The signed value stays much smaller. The gap is largest exactly where many near-minimal vectors
exist — (6,10,15) has m₇ = 60 yet an absolute constant of 19.0, because
its large pairwise gcds create a thicket of comparable vectors.

**Conclusion: m₇ alone does not supply a sharp bound in these families.**
The estimate needs a second invariant counting near-minimal vectors, or it
needs the cancellation that the signed sum enjoys and absolute values
discard. This is the same qualitative failure mode as THM-1070 — replacing
joint signed geometry by a one-sided bound — arriving by a different route.
The present data show a larger gap at `k=4` than at `k=3`; they do not prove a
growth law.

## 2. The Schur direction: a real invariant, but NOT a certificate

The Schur count (number of triples a+b=c inside V) has the property
min-speed lacked: it is **dilation-invariant**, so unlike the threshold
of MISTAKE-154/THM-1055 it is an admissible predictor of BONF5. It also
identifies the extremal family sharply in this comparison — {1,…,13} has
Schur count **36** against 0–6 for the other sampled families.

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

**Why the proposed inference fails:** sum-freeness removes the shortest Schur
relations highlighted by THM-1080, but it does not control the pair terms and
does not eliminate all higher resonances.  BONF5 in the experiment is
dominated by its low-support contributions.  Thus the explicit family refutes
"sum-free ⟹ certifiable" even though sum-freeness remains useful structural
data.

Note also that AP-ness and sum-freeness are not opposites: the AP
{1, 9, 17, …, 97} (d = 8) has Schur count **0**. Only the a = d APs — the
dilates of {1,…,13} — are additively rich. THM-1080's remark that "APs
are the additively richest 13-sets" is therefore too broad and should be
read as applying to that diagonal only.

## What survives

THM-1080's structural claim stands: additive richness drives the k-body
terms, and m₇ is the leading invariant. What this file removes is the
hope that either fact yields a certificate on its own. The remaining
route is unchanged and now doubly indicated: a second invariant for the
near-minimal-vector population, a bound that keeps the signs, or an exact
evaluation.  Both attempts that discarded joint signed information
(containment/fragmentation in THM-1070 and absolute values here) were loose on
the tested families.
