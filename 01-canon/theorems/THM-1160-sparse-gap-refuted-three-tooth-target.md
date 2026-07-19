---
id: THM-1160
title: The conditional sparse-gap reduction is sound; the sampled sparse-gap heuristic and the proposed three-tooth factor both fail
status: CORRECTED — the m=2 implication is proved; the 800-row observation is finite only; THM-1161 exactly refutes the universal m=3 spacing factor. Uniform r=5 remains OPEN
source: kind-pasteur-2026-07-18-S128 (cont.71; owner: prove the six-linear-functions statement for the four-comb tail)
depends_on:
  - THM-1147    # the exact conditional two-comb gap law
  - THM-1141    # the nonuniformity lever
  - THM-1140    # the four-comb frontier (as corrected by THM-1137/MISTAKE-169)
related: [THM-1137, THM-1097, THM-1148, THM-1161, MISTAKE-169, MISTAKE-171]
script: 04-computation/sparse_gap_lemma_kps_S128c71.py (+ .out)
---

# THM-1160 — the sparse-gap experiment and its corrected local conclusion

## Correction from THM-1161

The conditional `m=2` reduction in section I is sound.  Sections II--IV
originally made two further leaps which are not theorems:

1. absence of a 2-sparse gap in 800 sampled rows does not prove universal
   absence;
2. the proposed universal spacing bonus for `m=3` is false.

THM-1161 supplies an exact infinite legal family with one tooth from each
foreign comb and

```text
(longest piece)/(THM-1160 equal-split bound)
  = 1 + 20/(3K+6) -> 1.
```

Already at `K=132` the ratio is `211/201<1.295<4/3`.  Hence the sharp
universal local factor is one.  Also, `4/3` and `1.295` were conflated in
the original prose: the former is a strictly stronger factor, not an
equivalent restatement of the latter.  The historical computation below
remains useful telemetry, but it does not reduce the global four-comb
theorem to a valid one-gap spacing lemma.

## (I) The reduction is sound

Inside a core-safe component, the smallest killer k₁ cuts gaps of length **exactly**
6/(7k₁). If one such gap contains at most m teeth of k₂,k₃,k₄ — each of width ≤ 1/(7k₂) —
then it splits into ≤ m+1 pieces of total length ≥ 6/(7k₁) − m/(7k₂), so

> **L ≥ ( 6/(7k₁) − m/(7k₂) ) / (m+1).**

At **m = 2** the four-comb requirement L > 1/(7k₄) becomes

> 6/k₁ − 2/k₂ > 3/k₄,

which for clustered killers reads 4/k > 3/k — true with 33% room. Verified on every row:

| killers | 6/k₁ − 2/k₂ | 3/k₄ | holds |
|---|---|---|---|
| 157,158,159,160 | 0.02556 | 0.01875 | ✓ |
| 371,374,377,379 | 0.01082 | 0.00792 | ✓ |
| 550,553,554,558 | 0.00729 | 0.00538 | ✓ |
| 157,314,628,1256 | 0.03185 | 0.00239 | ✓ |

So **a 2-sparse gap would prove the four-comb theorem.**

## (II) No 2-sparse gap occurs in the frozen 800-row sample

| regime | trials | min foreign teeth over all k₁-gaps | trials with a 2-sparse gap |
|---|---|---|---|
| consecutive | 160 | **3** | **0** |
| step ≤ 3 | 160 | **3** | **0** |
| step ≤ 8 | 160 | **3** | **0** |
| step ≤ 30 | 160 | **3** | **0** |
| spread ×1.3 | 160 | **4** | **0** |

The first four regimes have sampled minimum three, the spread regime has
sampled minimum four, and none of the 800 rows has minimum two.  This
refutes the proposed averaging heuristic as a proof route; it does not prove
that every legal configuration has at least three. My counting heuristic was that each k₁-period holds one
tooth of each foreign comb, of which ~1/7 land inside the k₁-tooth and are harmless, giving
an average 18/7 ≈ 2.571 < 3 and hence a 2-sparse gap by pigeonhole. **That inference is wrong**: a
tooth straddling the k₁-tooth boundary *still cuts the gap*, so it is not harmless. All
three foreign teeth per period can count in the sampled clustered rows, and
this pigeonhole has nothing to bite on without a separate global phase argument.

## (III) At m = 3 the crude bound falls 23% short

On the standing worst case — core [1,3,5,6,7,8,11,12], killers 371/374/377/379 — the
component is [71/154, 41/84], it contains ten k₁-gaps, and every one holds exactly three
foreign teeth. The equal-split bound is

> (6/(7·371) − 3/(7·374))/4 = 1131/3885112 = 0.00029111,  giving 7·k₄·L ≥ **0.7723**.

Short of 1 by 23%.

## (IV) The standing row suggests tooth-position leverage, but not a theorem

The reported longest piece on this row is 0.0008888 — **3.05× the equal-split bound**. The crude bound
assumes three teeth split the gap into four *equal* pieces. THM-1147's
exact law says why (the gap from tooth j of a to tooth j+1 of b is (a − jd)/(ab) minus
radii, linear in j, so consecutive spacings differ systematically).

The original note proposed the following statement about **three numbers in one interval**:

> **Three foreign teeth inside a k₁-gap cannot split it into four pieces all shorter than
> 1.295 × the equal-split value.**

The `3.05` ratio is a measurement on one row, not uniform room.  THM-1154
refutes the displayed `1.295` target by an exact legal family whose factor
tends to one.  A successful four-comb proof must therefore choose among
multiple gaps/components or use a global phase-coupling invariant; it cannot
certify an arbitrary three-tooth gap by a fixed multiplicative bonus.

## Honest status

The sampled sparse-gap heuristic is not a theorem, the universal three-tooth
factor is exactly refuted, and the four-comb theorem is not proved. **Uniform
r=5 remains open.**  The sound remnant is the conditional `m=2` implication
and the exact wall-position language supplied by THM-1147.

## Named next
- Seek a global selection lemma showing that *some* gap has enough absolute
  survivor length.  THM-1154 proves that local relative dispersion alone
  cannot supply it.
- Retain the eight wall events of a gap as the faithful local carrier, then
  compare their phase stalks across consecutive `k1`-gaps.  Runner-only and
  single-gap piece tournaments forget exactly this global chronology.
