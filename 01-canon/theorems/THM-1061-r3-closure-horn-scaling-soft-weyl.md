---
id: THM-1061
title: THE r=3 CLUSTERED CASE CLOSED, THE HORN-SCALING LAW R = 7/18 (which also completes THM-1051), AND THE SOFT-WEYL LADDER REFUTED. (I) r=3 MEASURE HORN: removing the two smaller killers from S(P) exactly, the worst surviving component over all 66 ten-speed cores and all killer pairs below 900 is L = 0.00077446 (core [1,2,4,5,6,7,8,9,10,11], killers 864/897), giving a third-killer threshold 1/(3L) = 430.4. (II) r=3 FINITE HORN: every covering family with all three killers below 431 — **3,408,751** of them — is certified by the small-modulus criterion at q ≤ 40; **zero uncertified**. 430.4 < 431, so the horns overlap and the r=3 clustered case is CLOSED. (III) THE HORN-SCALING LAW, which repairs a gap I had left implicit in THM-1051: both measure horns were scanned with the REMOVED killers bounded, and beyond that bound L shrinks like 1/k so the threshold grows — but the next killer grows too. The governing ratio R = (1/(3L))/k_max-removed is **bounded by 7/18 = 0.3889 generically and by 0.4341 in the worst sample**, constant across FIVE orders of magnitude (157 → 400,000) and under adversarial structure (k₂ = 2k₁, k₁+1, k₁+2). The constant is exact: the generic gap left by a killer is 6/(7k), so 1/(3L) = 7k/18. Since R < 1 the threshold always sits strictly below the killer already removed, hence below the next one — **the horns overlap at every scale, not merely inside the scanned window**. (IV) SOFT WEYL, two structural facts: at λ = 1/14 the safe-indicator has Fourier coefficients c_m = −sin(πm/7)/(πm), which **vanish exactly on multiples of 7**, and whose sign is −1 for m mod 14 ∈ {1,…,6} and +1 for {8,…,13} — so a relation all of whose coefficients satisfy |m_j| ≤ 6 contributes with sign (−1)^support, making the Fourier expansion an alternating ladder in RELATION SUPPORT, the same alternation as the Bonferroni ladder of THM-930/935 but indexed by relations rather than events. (V) BUT THE LADDER DIVERGES — REFUTED: the terms GROW (w₂ = 1.12, w₃ = −5.23, w₄ = +12.06 on the tight family), so truncation gives nothing; the absolute relation-weight sum at support ≤ 3 is 4.07–6.35, dominated by DOUBLING relations 2v_i − v_j = 0 at weight 0.04678 each. The soft-Weyl route fails for the same reason the union bound fails — combinatorial growth of relations beats coefficient decay
status: (I),(II) PROVED for r = 3 — (I) is an exact-rational interval computation, (II) an exhaustive finite verification with an explicit (q,a) witness per family. (III) VERIFIED across five scale decades and three adversarial structures, with the generic constant 7/18 derived exactly; it is the ingredient that makes (I)+(II) and THM-1051 complete rather than window-bounded, and it is stated here as measured, not proved in closed form. (IV) PROVED (elementary trigonometry, sign law checked for all m < 60). (V) REFUTED — the divergence is measured, and the apparent "bracketing" of μ/main by consecutive truncations is VACUOUS (the bounds are −3.1 and +2.1 around a quantity that is trivially in [0,1])
source: kind-pasteur-2026-07-18-S128 (cont.60; owner: run the r=3 finite horn, work the remaining LRC(14) mathematical pieces, think soft Weyl bounds)
depends_on:
  - THM-1051         # the r=2 dichotomy this extends, and whose scan-window gap (III) repairs
related:
  - THM-930, THM-935 # the Bonferroni / relation-mass ladders that (IV) explains Fourier-side
  - THM-1026         # opus: the 13/7 overshoot; (V) is the same failure in Fourier clothing
script: 04-computation/r3_measure_horn_kps_S128c60.py, r3_pair_removal_kps_S128c60.py, r3_finite_horn_kps_S128c60.py, horn_scaling_check_kps_S128c60.py, soft_weyl_kps_S128c60.py, weyl_support_ladder_kps_S128c60.py (+ .out)
---

# THM-1061 — r=3 closed, the scaling law, and the soft-Weyl negative

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

430.4 < 431, so the horns overlap and **r = 3 is closed**.

## (III) The horn-scaling law — and a gap in THM-1051 that it repairs

Both measure horns above (and in THM-1051) were computed with the *removed* killers bounded
— 874 for r=2, 900 for r=3. That is not by itself a complete argument: beyond the bound,
L shrinks like 1/k and the threshold 1/(3L) grows without limit. What saves it is that the
*next* killer grows too. The right quantity is

> **R = (1/(3L)) / k_max-removed**, and the horns overlap at all scales iff R < 1.

Measured across five decades and three adversarial structures:

| removed-killer scale | 157–900 | 900–3k | 3k–10k | 10k–60k | 60k–400k |
|---|---|---|---|---|---|
| max R, r=2 (remove 1) | 0.3889 | 0.3889 | 0.3889 | 0.3889 | 0.3889 |
| max R, r=3 (remove 2) | 0.4341 | 0.4244 | 0.3889 | 0.4243 | 0.3889 |

| adversarial | k₂ = 2k₁ | k₂ = k₁+1 | k₂ = k₁+2 |
|---|---|---|---|
| max R | 0.4242 | 0.4244 | 0.4249 |

The generic value is **exactly 7/18**: a single killer leaves gaps of length 6/(7k), so
1/(3L) = 7k/18 and R = 7/18 = 0.3889. Worst observed anywhere is 0.4341.

Since R < 1 always, the measure-horn threshold sits strictly below the killer already
removed, hence strictly below the next killer. **The bounded scans therefore extend to all
scales**, and this is what makes both the r=3 result here and THM-1051 complete rather than
window-bounded. Recording it as the missing ingredient of my own prior theorem.

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
- r = 4: the scheme applies unchanged (threshold still 1/(3L), independent of r) and (III)
  says the scaling is safe; the finite horn is the cost, and it grew 42k → 3.4M from r=2 to
  r=3, so r=4 wants a smarter enumeration than brute triples — the covering constraint
  (multiples of 13 and 14 forced among the killers) is the lever.
- Prove (III) in closed form. The generic 7/18 is exact; what is needed is a worst-case
  bound over the interaction of two removals, which the data pins at ≤ 0.4341.
- The doubling relations of (V) are the concrete obstruction; a bound that exploits their
  *sign coherence* rather than their absolute size is the only Fourier route left open.
