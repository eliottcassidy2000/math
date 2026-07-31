---
id: THM-3004
title: "Circuit sign-change cluster law, and refutation of the two-sign classifier"
status: REFUTED (the classifier) + VERIFIED-EXACT (the cluster law) / AWAITING INDEPENDENT HOSTILE AUDIT
source: klein-S428
depends_on:
  - THM-3000-fixed-edge-cumulant-curvature-universality-and-bounded-jet-transfer
  - THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law
  - THM-3003-antipodal-circuit-rigidity-and-the-multipole-spread-criterion
related:
  - THM-2991-pf-infinity-arbitrarily-delayed-newton-ratio-return
script: 04-computation/gmc_circuit_sign_change_cluster_law_and_classifier_refutation_thm3004.py
output: 05-knowledge/results/gmc_circuit_sign_change_cluster_law_and_classifier_refutation_thm3004.out
script_sha256: 4e7e6a6c6c273ee6972c2035ecff489c939ecc6c3b2fba69c111733007b6fc8f
output_sha256: 1309387ee4ede41258fdd9235cd07ad17088d702a377261d16e2d06e3a20d8a8
hash_basis: LF-normalized bytes
---

# THM-3004 -- circuit sign-change cluster law, and refutation of the two-sign classifier

**REFUTED (the classifier) + VERIFIED-EXACT (the cluster law).**

This file **retracts** the hypothesis of
[THM-3001](THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law.md)
section 6 and replaces it with the correct law.  It also records the census
defect that let the false hypothesis stand, as MISTAKE-337.

## 1. The refutation, with a degree-five witness

Notation is THM-3000 section 1; `c_k=log(R_k/R_(k-1))`, `k=2..d-1`.

**Refuted claim** (THM-3001 section 6): for real-rooted positive `N` with
bounded jets, the pair `(sign C(mu), sign C(mu*))` determines the global shape of
`R`; in particular `C(mu)>0` and `C(mu*)>0` gives an interior maximum, i.e.
**exactly one** sign change of `c`.

**Witness** (exhaustive search, smallest degree in the search pool):

    N(n)=(n+1)^2 (n+3)^2 (n+8),   d=5,

with exact ratios

    R_1=256/215,  R_2=1849/1600,  R_3=10000/8643,  R_4=4489/4000,

numerically `1.190698, 1.155625, 1.157006, 1.122250`: **down, up, down**.  The
circuit sign pattern is `-+-`, so `c` has **two** sign changes.  Meanwhile

    C(mu)=+0.14135742...,   C(mu*)=+0.41294497...,

both positive, so the classifier predicts one sign change.  It is wrong.

All roots are real and positive, so `N` is PF-infinity, Hurwitz stable and
strictly ULC (`R_k>1` for every `k`, checked exactly).  The witness therefore
lives in the interior of every hypothesis class the lane cares about; it is not
a boundary artefact.

**Corollary.**  Unimodality of the Newton-ratio sequence is FALSE for
real-rooted positive polynomials, from degree five on.

## 2. The mechanism, and why the witness was predictable

This was found by the multipole picture of THM-3003 section 3, *before* the
search, not by luck.  For well-separated clusters `e_k` is dominated by the
product of the `k` largest roots, so

    log e_k = sum_(i<=k) log r_(i) + O(1),
    c_k = -Delta^2(sorted log-root step function)_k + binomial term.   (1)

A step function with `m-1` kinks has second difference supported on `2(m-1)`
positions with alternating signs.  Three clusters therefore already force
several sign changes -- and unequal cluster **sizes** are what move the kinks
into the interior where they are visible.

## 3. The cluster law (VERIFIED-EXACT, attained)

**Law.**  For `m` well-separated root clusters the maximum number of sign
changes of the circuit is exactly

    0            for m=1,
    2m-3         for m>=2.                                             (2)

Observed maxima over randomised size profiles with separation `1000`:

| `m` | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| max sign changes | 0 | 1 | 3 | 5 | 7 | 9 |
| `2m-3` | -- | 1 | 3 | 5 | 7 | 9 |

The `2(m-1)` spikes of (1) would give `2m-3` gaps between them, and the boundary
trims exactly one; both the upper bound and its attainment are confirmed on
every `m<=6`.

**This is the honest replacement for the classifier.**  The two end curvatures
of THM-3001 control the two *ends* of the circuit and nothing else; the number
of interior reversals is a **cluster count**, i.e. a property of the root
measure's support structure, invisible to any fixed finite set of low moments.

## 3a. Localization: WHERE the reversals sit, and the free-fermion reading

`e_k` is exactly the canonical `k`-fermion partition function of a system whose
grand partition function is `prod_i(1+r_i t)`: single-particle Boltzmann weights
`r_i`, fugacity `t`.  So `log h_k` is a **free energy** in the binomial gauge and
the circuit is its **third derivative in particle number**.  Well-separated root
clusters are **bands**; filling a band is a first-order transition in `k`.

**Localization law (VERIFIED-EXACT).**  Sort the clusters by root size,
descending, with sizes `s_1,...,s_m`.  Then every partial sum
`s_1+...+s_j` (`j=1..m-1`) -- the band-filling numbers -- is a sign-change site
of the circuit, and exactly one further site lies strictly between each pair of
consecutive band-filling numbers.  Hence

    sign changes = (m-1) band boundaries + (m-2) inter-band sites = 2m-3,   (3)

which is a complete structural account of (2), not just a count.  Confirmed on
`(4,5,3)`, `(3,3,3)`, `(6,2,4)`, `(2,7,3)`, `(5,4,2,3)`, `(3,3,3,3)`,
`(2,3,4,2,3)`; e.g. `(2,3,4,2,3)` has band-filling numbers `[2,5,9,11]` and
reversal sites `[2,3,5,7,9,10,11]`, exactly `2*5-3=7` of them.

**Exact zeros are forced by rigidity, not extra reversals.**  Equal cluster
sizes under geometric separation make the log-root measure symmetric, so
THM-3003 section 1 gives `c_k=-c_(d+1-k)`; when `d` is odd, `k=(d+1)/2` is a
fixed point of that involution and the circuit vanishes there **exactly**.
Checked: `(3,3,3)` has `d=9` and an exact zero at `k=5`, precisely as predicted,
while `(2,2,2)`, `(4,4,4)` and `(2,2,2,2)` have even `d` and no zero.  Sign-change
counts must therefore be taken on the zero-filtered sign word; counting a tie as
two reversals inflates `(3,3,3)` from `3` to `4` and would falsify (2)
spuriously.

## 4. What survives: the classifier is true for at most two clusters

Exhaustive over `936` two-cluster configurations -- `d=4..16`, ratios
`{1/3,2,3,5,7,10,100,10^4}`, every multiplicity split -- the circuit has **at
most one** sign change, with zero violations.  Consistent with (2) at `m=2`.

So the correct scoping is:

    m<=2 : the two-sign classifier holds (verified exhaustively);
    m>=3 : it fails, and the failure rate is not small.

## 5. Census defect (MISTAKE-337)

Among all three-cluster configurations with `d=6..12` over the root pool
`{1,2,3,5,10}`, `51/2100 = 2.43%` have two or more sign changes -- not a rare
pathology.  Restricted to **equal** cluster sizes (`d//3` each), the failure rate
is `0/30 = 0%`.

THM-3001 section 6's `42/42` census varied the root ratios, the number of
clusters, and the degree, but held the cluster **sizes** equal in every
three-cluster row.  That is precisely the axis the failure lives on.

> **Rule.**  A census is evidence only about the axes it varies.  Sample size
> does not substitute for an un-varied coordinate.  Before quoting an `n/n`
> census, list the coordinates of the configuration space and mark which ones
> were actually moved.

## 6. Consequences for the lane

- THM-3000, THM-3001 sections 1--5 and THM-3003 are **untouched**: they concern
  the two ends and the reversal involution, all of which remain exact.  Only
  THM-3001 section 6 is retracted.
- THM-3001's proved necessary condition `C(mu)>=0>=C(mu*)` for asymptotic global
  no-return still stands, and is still the cheap two-scalar screen.  What is now
  known is that it is **not** sufficient, and no bounded set of moments can be.
- This strengthens THM-3001 section 2 in spirit: not only can no reversal-closed
  class prove no-return, but the obstruction is not even a low-moment condition.
  Any successful family-specific proof must control the **cluster structure** of
  the core's root measure, not just its first few multipole moments.
- For the first-gap family this is a concrete new question: how many separated
  scales does the wall-stripped core's root measure have?  The wall itself has
  roots at `-1,...,-M` with a five-band multiplicity profile (THM-2997 (9)),
  which is a *continuum* of scales rather than a few clusters, so (2) is not
  immediately alarming -- but it must be checked rather than assumed.

## 7. Reproduction

    python3 04-computation/gmc_circuit_sign_change_cluster_law_and_classifier_refutation_thm3004.py

Five parts, all reporting `True`: the exhaustive minimal-witness search with
exact rationals and a Newton `R_k>1` control, the cluster law through `m=6`, the
exhaustive two-cluster scope, the census audit, and the band-localization law
with its forced-zero control.
