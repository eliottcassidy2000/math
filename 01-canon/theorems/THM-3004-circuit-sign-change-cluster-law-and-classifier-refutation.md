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

**Scope clarification, 2026-09-06 continuing8.** The general cluster wording
below records finite evidence. A new
[independently audited theorem](../../05-knowledge/results/continuing8_20260906_newton_clusters.md)
proves exactly `2K-3` changes for every fixed profile with multiplicities at
least two, under explicit narrow-block and separation bounds. It does not
promote the arbitrary-separation candidate. The complementary
[circuit surjectivity theorem](../../05-knowledge/results/continuing8_20260906_newton_universality.md)
realizes every positive circuit-ratio vector with distinct negative roots
and fixed root-parameter sum. Thus real-rootedness alone does not restrict
the circuit word; the cluster geometry or further moments must be retained.

This file **retracts** the hypothesis of
[THM-3001](THM-3001-newton-circuit-reversal-involution-and-two-end-curvature-law.md)
section 6 and replaces it with the scoped candidate below. It also records the census
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

All root parameters are real and positive (the polynomial roots are negative), so `N` is PF-infinity, Hurwitz stable and
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

**Finite-supported candidate.** For `m` well-separated root clusters the maximum number of sign
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


## 3b. UPGRADE: the bound is `2K-3` for `K` DISTINCT roots, and attainment is PROVED

An independent adversarial pass (klein-S428 second wave) strengthens (2) in two
ways, both re-derived here:

**(i) Upper bound for arbitrary distinct roots.**  The bound is not about
separation: over `375,698` exact tests the circuit has at most `2K-3` sign
changes where `K>=2` is the number of **distinct** roots, with zero violations.
Separation is what makes it *attained*; THM-3005 shows the other extreme, a
near-continuum profile with `K=15` distinct scales and **zero** sign changes.

**(ii) Attainment is PROVED, not merely observed.**  For `K>=2`, `m>=2`, `T` large,

    N_(K,m,T)(n)=prod_(j=0)^(K-1)(n+T^j)^m,   d=Km,

has at least `2K-3` sign changes, and exactly `d-3` -- the **maximum** a length
`d-2` vector admits -- when `m=2`.  Tropical/Newton-polygon proof: with
`k=qm+s`, `e_k` is dominated by the unique top-down filling, so
`e_k=C(m,s)T^(D(k))(1+O(1/T))` with `D` concave piecewise linear of slope
`K-1-q`; hence `Delta^2 D(k-1)=-1` exactly when `m | k` and `1<=k/m<=K-1`, and

    log R_k = 1_(m|k, 1<=k/m<=K-1) log T + B_k + O(1/T),

with `B_k` bounded and `T`-free.  So `c_(qm)=+log T+O(1)` and
`c_(qm+1)=-log T+O(1)` for `q=1..K-1`, giving an alternating subsequence of
length `2(K-1)`.  Deleting entries cannot increase sign changes, so the
subsequence count is a lower bound.  For `m=2` the indices `{qm,qm+1}` exhaust
`[2,d-1]` and the whole circuit alternates.

Verified independently: `(n+1)^2(n+3)^2(n+9)^2` gives
`R=(65/57, 4693/4005, 71289/61009, 4693/4005, 65/57)`, word `+-+-`, `3=d-3`
changes -- already at `T=3`, not merely for large `T`.  Note this is a geometric
(log-symmetric) profile, so THM-3003 section 1 correctly predicts the palindrome:
the palindrome is **W-shaped**, not single-humped.

**(iii) The necessary condition is demonstrably not sufficient.**  The hostile
named in THM-3001 section 6 exists: roots `1^u 3^(8u) 100^(51u)` (`d=60u`) have

    C(mu)=9523726464503/27595322265625=+0.34512104...,
    C(mu*)=-4425627648481/2464928260081=-1.79543872...,

so THM-3001's limiting screen holds with a fixed strict margin, the root measure has fixed
compact support (bounded jets at both ends), and yet the circuit is
`+...+-----+` with a strict interior dip that persists and widens proportionally
to `d`.  Verified exactly at `u=1`.

**(iv) Minimality is only over the search pool.**  Section 1's witness is minimal
over an **integer** root pool.  Degree five is genuinely minimal (a length-`d-2`
circuit needs `d>=5` to admit two sign changes), but small *spread* is not:
`(n+1)(n+3/2)^2(n+8/3)(n+11/4)` is real-rooted with word `+-+` and max/min spread
only `11/4`, well below the integer witnesses' spread `8`.  Do not read
"exhaustive" over an integer grid as exhaustive over real roots.

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

**Second, independent defect (amended after an adversarial pass).**  The census's
`shape_of` decides `INTERIOR-MAX`/`INTERIOR-MIN` from `R[2]>R[1]` and
`R[d-1]<R[d-2]` alone -- the two **end** circuits -- so it is blind to interior
oscillation whenever the ends happen to match.  Exact split over the `51` failing
configurations: `46` would have been caught (`shape_of` returns `MIXED`), `5`
slip through.  So the census was **not** vacuous, and an audit note claiming its
agreement was "logically compatible with arbitrary interior oscillation" is too
strong; the blind spot is precisely the W-shaped palindromes.

> **Rule.**  A census is evidence only about the axes it varies.  Sample size
> does not substitute for an un-varied coordinate.  Before quoting an `n/n`
> census, list the coordinates of the configuration space and mark which ones
> were actually moved.

## 6. Consequences for the lane

- THM-3000, THM-3001 sections 1--5 and THM-3003 are **untouched**: they concern
  the two ends and the reversal involution, all of which remain exact.  Only
  THM-3001 section 6 is retracted.
- THM-3001's repaired quantitative necessary condition
  `C(mu_d)>=-O(1/d), C(mu_d*)<=O(1/d)` for asymptotic global no-return still
  stands and is still the cheap two-scalar screen.  What is now
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
