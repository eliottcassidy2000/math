---
id: THM-1061
title: THE r=3 CLUSTERED CASE IS UNIFORMLY CLOSED BY THE EXACT TWO-COMB COMPONENT THEOREM; THE ORIGINAL SAMPLED SCALING ARGUMENT REMAINS WITHDRAWN, AND THE SOFT-WEYL LADDER IS REFUTED
status: PROVED uniformly for r=3 by THM-1094. Its exact 9,246,070-pair bank plus analytic tail proves the strong component theorem L>1/(3k2); independently, its sharp periodic-discrepancy corollary proves the weaker sufficient L>1/(7k2) from only the 66-core atlas and elementary inequalities, so uniform r=3 closure does not depend on the large bank. The older bounded interval computation and 3,408,751-row finite horn remain exact independent checks. MISTAKE-163 remains in force against the original sampled all-scale inference; THM-1094 supplies a new proof rather than rehabilitating that inference. The Fourier coefficient formula and elementary sign law are exact; the relation-support truncation is computationally refuted
source: kind-pasteur-2026-07-18-S128 (cont.60; owner: run the r=3 finite horn, work the remaining LRC(14) mathematical pieces, think soft Weyl bounds)
depends_on:
  - THM-1051         # the r=2 dichotomy this extends, and whose scan-window gap (III) repairs
  - THM-1094         # exact two-comb component theorem; the uniform all-scale bridge
related:
  - THM-930, THM-935 # the Bonferroni / relation-mass ladders that (IV) explains Fourier-side
  - THM-1026         # opus: the 13/7 overshoot; (V) is the same failure in Fourier clothing
  - MISTAKE-163      # the old sampled proof was invalid; the new proof does not erase that audit
script: 04-computation/r3_measure_horn_kps_S128c60.py, r3_pair_removal_kps_S128c60.py, r3_finite_horn_kps_S128c60.py, horn_scaling_check_kps_S128c60.py, soft_weyl_kps_S128c60.py, weyl_support_ladder_kps_S128c60.py, lrc14_r3_two_comb_component_exact_codex_S73.cpp, lrc14_r3_two_comb_extremal_replay_codex_S73.py (+ .out)
---

# THM-1061 — r=3 uniformly closed, and the soft-Weyl negative

> **Audit correction and repair (codex-S73; MISTAKE-163; THM-1094).**  The
> finite horn and bounded exact interval scan below were valid.  The former
> all-scale inference from sampled ratios was not and remains withdrawn.
> THM-1094 now closes the gap by a different proof: an exact 9,246,070-pair
> endpoint bank plus an elementary mass/component tail covering every omitted
> scale.

## (I)–(II) The r=3 clustered case

Same two-horn scheme as THM-1051, one killer further. The core now has 10 speeds, so S(P)
is roomier: over the 66 ten-speed cores the largest component of S(P) ranges 0.0089–0.0536
(median 0.0211) and the total measure 0.084–0.207.

**Measure horn.** Remove the two smaller killers from S(P) exactly; the last killer needs
only k₃ > 1/(3L), where L is the largest surviving component. Note the threshold does *not*
depend on r — THM-1051's crude 2r/(L(7−r)) formula was weaker than necessary. Worst case
over all 66 cores and all killer pairs below 900:

> **L = 0.00077446** at core [1,2,4,5,6,7,8,9,10,11], killers (864, 897)
> **threshold 1/(3L) = 430.4**

**Finite horn.** Covering forces a multiple of 13 *and* a multiple of 14 among the killers
(a core inside {1,…,12} supplies neither), which prunes the triple enumeration hard.
Exhaustive check by bitmask intersection against the small-modulus criterion at q ≤ 40:

> **3,408,751 covering families with all three killers < 431 — 3,408,751 certified, 0 uncertified.**

430.4 < 431, so these two original horns overlap throughout the exact interval
scan's `k2<900` window.  They are now retained as independent bounded checks;
uniform closure comes from THM-1094.

## (III) The exact all-scale bridge — THM-1094

For every ten-speed core `P subset {1,...,12}` and every
`13 max(P)<k1<k2`, THM-1094 proves

```text
L(P;k1,k2)>1/(3k2),
```

where `L` is the longest component after removing both danger combs from the
core-safe set.  Its proof has two exhaustive pieces in the logical, not
heuristic, sense:

1. Every core has a component of length `ell>=1/112`.  Tooth incidence and
   occupied-length bounds prove the desired inequality whenever

   ```text
   56ell k2-49ell k1-24(k2/k1)-185>0.
   ```

   This includes every `ell k1>=209/7`, hence every `k1>=3344`, and gives an
   explicit `k2` tail below that height.
2. The exact rational endpoint referee exhausts the finite complement:
   66 cores, 9,246,070 guarded pairs, zero failures.  Its hardest finite-bank
   row has `3k2L=158/119` at core
   `{1,2,3,5,6,7,8,9,10,11}` and `(k1,k2)=(153,159)`.

The strong theorem is more than closure needs.  THM-1094 also proves directly,
using the sharp periodic bound

```text
|I intersect D_k| <= |I|/7+6/(49k),
```

that `L(P;k1,k2)>1/(7k2)`.  For 65 cores this follows from the exact atlas
inequality `ell(P)(13max(P)+1)>=1727/1008>5/3`.  On the unique smaller core,
the same inequality handles `k1>=187`; the residual
`157<=k1<k2<=187` has at most five components and closes by the elementary
positive polynomial

```text
5k1k2-96k2-656k1 >= 5870.
```

For a third killer `k3>k2`, the sharp periodic bound is then strictly less
than `L`.  This proves the uniform three-killer clustered case from the
66-core atlas alone, without a finite/measure split or covering hypothesis.
The 9,246,070-row bank remains the proof of the stronger `1/(3k2)` structural
statement and an independent stress test of the endpoint carrier.

The historical scaling samples remain useful diagnostics, but not proof.
In particular, the exact hardest finite-bank ratio is `R=119/158=0.753...`,
well above the old sampled `0.4341`; this quantitatively confirms why
MISTAKE-163 had to be recorded.  THM-1051 remains uniform by its own
finite/mixed/explicit-large split.

## (IV) Soft Weyl: two structural facts

Expanding the safe indicator, 1_{‖x‖≥λ} = Σ_m c_m e(mx) with c_0 = 1−2λ and
c_m = −sin(2πmλ)/(πm). At **λ = 1/14** this is c_m = −sin(πm/7)/(πm), so:

1. **The Fourier support avoids every multiple of 7.** c_{7k} = 0 identically — a structural
   zero tied to the 7 in 1/14 = 1/(2·7). Relations supported on multiples of 7 are invisible.
2. **Sign law:** c_m < 0 for m mod 14 ∈ {1,…,6} and c_m > 0 for m mod 14 ∈ {8,…,13}
   (checked for all m < 60). Hence a relation with all |m_j| ≤ 6 contributes with sign
   **(−1)^s**, s its support.

By orthogonality μ = Σ_{m·v = 0} ∏_j c_{m_j}, so the main term is (1−2λ)¹³ = (6/7)¹³ =
0.13480 and everything else is a relation correction. The sign law makes that correction an
**alternating ladder in relation support** — the Fourier-side explanation of the Bonferroni
alternation in THM-930/935, with relations in place of events. Measured corrections:

| family | μ | μ/main |
|---|---|---|
| tight {1,…,13} | **0 exactly** | **0.0000** |
| {1,…,12} ∪ {169} | 0.02396 | 0.1777 |
| {2,…,12} ∪ {169,182} | 0.07462 | 0.5535 |
| {2,…,12} ∪ {1000,2000} | 0.08304 | 0.6160 |
| {1,…,11} ∪ {312,364} | 0.04247 | 0.3150 |

The tight family sits at exactly μ/main = 0 — the correction cancels the main term precisely,
which is a sharp consistency check on the whole expansion.

## (V) The ladder diverges — the soft-Weyl route is refuted

The alternation is real but **useless**, because the terms grow instead of decaying:

| family | w₂ | w₃ | w₄ |
|---|---|---|---|
| tight {1,…,13} | +1.120 | −5.227 | +12.056 |
| {2,…,12} ∪ {169,182} | +0.938 | −3.135 | +5.834 |
| {2,…,12} ∪ {1000,2000} | +1.017 | −3.135 | +5.809 |

Partial sums 1+w₂, 1+w₂+w₃, 1+w₂+w₃+w₄ run 2.12, −3.11, +8.95 against a true value of 0.
Consecutive truncations do "bracket" μ/main, but **vacuously** — the bracket is [−3.11,
2.12] around a quantity that is trivially in [0,1]. I am recording this explicitly because
the bracketing check returns `True` and could easily be mistaken for a result.

The absolute (soft) bound fails by the same margin: the relation-weight sum at support ≤ 3
is 4.07–6.35 against the 1 it would need, dominated by **doubling relations** 2v_i − v_j = 0
at weight 0.04678 apiece — and those live in the *core* ({1,…,12} contains 2·2=4, 2·3=6,
2·4=8, 2·5=10, 2·6=12), so no choice of killers can reduce them.

**Diagnosis:** the number of relations of support s grows faster than ∏|c| decays. This is
opus's 13/7 overshoot (THM-1026) in Fourier clothing — the same wall, reached by a different
road. Soft Weyl does not open it.

## Named next
- THM-1097 has now completed the analogous three-removal all-scale bridge for
  `r=4`, using THM-1094's sharp periodic estimate
  `|I intersect D_k|<=|I|/7+6/(49k)` and the target
  `L>1/(7k_max-removed)`.  The next clustered problem is the four-removal
  `r=5` bridge, where the coarse mass/component slope crosses the wrong way.
- The doubling relations of (V) are the concrete obstruction; a bound that exploits their
  *sign coherence* rather than their absolute size is the only Fourier route left open.
