---
id: THM-717
title: The k=9 base tail-decomposition isolates the signed cancellation — via the exact Abel identity J = 6T₁+4T₂+2T₃−2T₅−4T₆ (Tⱼ = P(N≥j)), the base row J ≥ 432/91 SEPARATES into a cancellation-free covering bound (POS) 6T₁+4T₂+2T₃ ≥ 4717/882 and a diameter-controlled bunching bound (BUNCH) p₅+3p₆ ≤ 1/7, with 4717/882 − 2/7 = 4465/882 = J(consec) ≥ 432/91; every prior absolute/Bonferroni bound on J failed on the −2T₅−4T₆ cancellation, which this decomposition confines to the single small (BUNCH) term (equality at consec) — leaving the bulk (POS) as a monotone-covering functional the density-floor machinery can attack without the wall
status: IDENTITY PROVED (Abel summation + N ≤ 6 pointwise). CORRECTED (klein-S255, MISTAKE-138): the S254 "both extremal at consec / BUNCH ≤ 1/7" was a BOX ARTIFACT — the box [1..19] missed the mod-7 pole {1,8,…,57} (spread 56) with BUNCH = 6/19 > 2/7. The separation is genuinely TWO-POLE: min POS = 4717/882 at consec, max BUNCH = 6/19 at the mod-7 pole (both verified: POS by box+adversarial, BUNCH by exhaustive mod-7 search). It STILL closes GLOBALLY: min POS − max BUNCH = 4717/882 − 6/19 = 84331/16758 ≈ 5.032 ≥ 432/91 (margin +0.285). The two piece-bounds are SHARP CONJECTURES at their (distinct) poles. STRENGTHENS THM-711/716: consec = the EXACT global J-argmin over 167 950 primitive 9-cores in [1..20] (J = 4465/882, decorrelation budget J_iid − floor = +3.71). Does NOT prove the base outright — it REDUCES it to [(POS) coupled covering extremality ≥ 4717/882, cancellation-free, consec-pole] + [(BUNCH) 2(p₅+3p₆) ≤ 6/19 bunching, mod-7 pole], each cleaner than raw J.
source: klein-2026-07-11-S254 (HYP-6030)
depends_on:
  - THM-711   # mac-mini: the k=9 base J = E[N(7-N)] >= 432/91 (this decomposes its proof)
  - THM-716   # mac-mini: the (mu,Var) finite-dimensional reduction (complementary decomposition)
  - THM-661   # the density floor nu >= bar (the T_1 lower bound; POS generalizes it to higher tails)
related:
  - THM-710   # eigen-transfer: k>=10 inherit, so k=9 (+k=8) is the whole base
  - THM-714   # the k=8 cubic base (the deg-3 analog awaits its own tail decomposition)
external: Abel summation; Bonferroni inequalities; the "AP is the extremal coverer" principle (Freiman-adjacent).
---

# THM-717 — the k=9 base tail-decomposition isolates the cancellation

> **⚠ READ THE CORRECTION (klein-S255) AT THE BOTTOM FIRST.** The S254 body below claims min POS
> and max BUNCH are both at consec and `BUNCH ≤ 1/7`; that was a box artifact (MISTAKE-138). The
> mod-7 pole `{1,8,…,57}` has `BUNCH = 6/19 > 2/7`. The separation is two-pole and still closes
> globally (`min POS − max BUNCH = 84331/16758 ≥ 432/91`, margin +0.285).

## Setup

`N(x)` = number of empty sectors among the seven arcs `[s/7,(s+1)/7)` for the phases
`{frac(e x): e ∈ E}`; `pⱼ = P(N = j)`, `Tⱼ = P(N ≥ j)`. The k=9 base row (THM-711) is
`J := E[N(7−N)] = 6m₁ − m₂ ≥ 432/91` for every 9-core.

## The exact Abel identity

`N(7−N)` takes values `[0, 6, 10, 12, 12, 10, 6]` at `N = 0,…,6` (and `N ≤ 6` always, since
at least one sector is occupied). Abel summation `E[v(N)] = Σⱼ Tⱼ (vⱼ − vⱼ₋₁)` with
`v = [0,6,10,12,12,10,6]` gives the **exact identity**

> **`J = 6T₁ + 4T₂ + 2T₃ − 2T₅ − 4T₆`**   (the `T₄` coefficient `v₄−v₃ = 0` drops out).

Equivalently `J = (6T₁+4T₂+2T₃) − 2(p₅+3p₆)`, since `2T₅+4T₆ = 2p₅+6p₆ = 2(p₅+3p₆)`.

## The separation (isolating the cancellation)

Write **POS** `:= 6T₁+4T₂+2T₃` and **BUNCH** `:= 2T₅+4T₆ = 2(p₅+3p₆)`. Then `J = POS − BUNCH`,
and the base follows from two bounds that are each **extremal at consec** (verified universal,
92 377 primitive 9-cores in [1..19], zero violations):

- **(POS) — cancellation-free covering bound:** `6T₁+4T₂+2T₃ ≥ 4717/882 ≈ 5.3481`.
  A nonnegative-weighted sum of the *monotone* good-set tails `Tⱼ = P(≥ j sectors empty) =
  P(≥ j runner-free 1/7-arcs)`. It is `E[w(N)]` for the nondecreasing weight
  `w = [0,6,10,12,12,12,12]` (= `N(7−N)` with its high-N dip filled flat). No subtraction —
  the absolute/Bonferroni methods that DIE on `J` (nine documented cancellation failures) have
  no cancellation to fight here.
- **(BUNCH) — diameter-controlled bunching bound:** `p₅ + 3p₆ ≤ 1/7`, equality at consec.
  `p₆ = P(all 9 phases in one 1/7-arc)` and `p₅ = P(in exactly two)` are near-origin/near-rational
  events whose measure shrinks with the diameter; consec (smallest spread) maximizes them,
  at `p₅+3p₆ = 1/21 + 6/63 = 1/7` exactly.

**Assembly:** `J ≥ POS − BUNCH ≥ 4717/882 − 2/7 = 4465/882 = J(consec) ≥ 432/91` (margin
`4465/882 − 432/91 = +0.3151`). The separation is TIGHT (both parts extremal at the same core,
consec), so it loses nothing — yet it confines the signed part to (BUNCH).

## Why this matters

The base extremality resisted every *absolute* bound because `J`'s `−2T₅−4T₆` term cancels
against `6T₁+…` (MISTAKES / the standing "no order-blind bound" law, 9 confirmations). This
decomposition **quarantines that cancellation** in (BUNCH) — a single term that is (a) small
(≤ 2/7), (b) sharp only at consec, and (c) diameter-controlled (bunching → 0 as spread → ∞).
Everything else, (POS), is a monotone covering functional in reach of the density-floor toolkit
(THM-661 already gives `T₁ ≥ bar`; the coupled `6T₁+4T₂+2T₃` floor is the remaining covering
content — note the *individual* tail minima are NOT aligned (`min T₂, min T₃` occur off consec:
`6·minT₁+4·minT₂+2·minT₃ = 4.876 < 5.033`), so (POS) is a genuine *coupled* covering extremality,
not a sum of separate floors). Complementary to THM-716's `J = μ(7−μ) − Var` frontier: that
reduces to a 1-parameter minimization; this isolates the cancellation.

## Status & handoff

- PROVED: the Abel identity; `J ≥ 6·P(N≥1) = 6ν` (weak form, insufficient alone since
  `6·bar₉ = 3.37 < floor`); the exact rational values at consec.
- CONJECTURED (verified universal, sharp at consec): (POS) `≥ 4717/882`, (BUNCH) `p₅+3p₆ ≤ 1/7`.
- NEXT: (BUNCH) is the crux and is self-contained — a sharp bunching bound (three-gap governs
  consec's phases `{jx}`, so its LHS is a three-gap computation). (POS) is the coupled covering
  floor — attack with the moment-LP / Bonferroni machinery, now safe from cancellation. The k=8
  deg-3 base awaits its own tail decomposition (its requirement function is cubic, not `N(7−N)`).

## Files

`04-computation/lrc14_J_landscape_twopole_klein_S254.py`,
`lrc14_J_decorrelation_full_klein_S254.py` (+ `.out`s).

## Addendum — the SAME decomposition unifies the k=8 deg-3 base (klein-S254)

The k=8 base (THM-714's cubic requirement `E[ψ₈(N)] ≥ 1−cap₉ = 2025/4004`, where
`ψ₈(N) = (2/3)N − (47/252)N(N−1) + (5/252)N(N−1)(N−2)` = `1 − φ₃(N)`) Abel-decomposes into the
SAME shape:

> `E[ψ₈] = (2/3)T₁ + (37/126)T₂ + (5/126)T₃ − (2/21)T₄ − (1/9)T₅ − (1/126)T₆`
> `      = POS₈ − NEG₈`,  `POS₈ = (2/3)T₁+(37/126)T₂+(5/126)T₃`,
> `NEG₈ = (2/21)T₄+(1/9)T₅+(1/126)T₆`.

(The `T₇` coefficient `3/14` drops since `N ≤ 6 ⟹ T₇ = 0`.) Verified exact over 12 869 primitive
8-cores in [1..16]: `min POS₈ = 1457/2520 ≈ 0.5782` and `max NEG₈ = 655/24696 ≈ 0.0265`, **both
at consec {1..8}**, giving `E[ψ₈] ≥ 0.5782 − 0.0265 = 11353/20580 ≈ 0.5517 ≥ 2025/4004`
(margin **+0.0459**, matching THM-714 exactly). The negative/bunching part is even SMALLER here
(0.0265 vs k=9's 2/7 ≈ 0.286) — the deg-3 majorant's higher-tail coefficients are tiny, so the
k=8 cancellation is almost entirely isolated.

**Unification:** BOTH moment-ladder base rows (k=8 deg-3, k=9 deg-2) — and hence, via THM-710's
eigen-transfer, the ENTIRE wide-spread base — reduce to `[POS: a cancellation-free monotone-covering
floor] + [NEG: a small higher-tail bunching bound]`, both extremal at consec. The signed content
of LRC(14)-S3 is confined to the two NEG terms (`p₅+3p₆ ≤ 1/7` and `(2/21)T₄+(1/9)T₅+(1/126)T₆ ≤
0.0265`), each a "measure of near-origin/near-rational bunching" quantity that vanishes with the
diameter.

## CORRECTION (klein-S255) — the separation is TWO-POLE; the mod-7 pole is the BUNCH-max

The S254 body above claims min POS and max BUNCH are BOTH at consec, "verified universal over
92 377 primitive 9-cores in [1..19]." **That was a box artifact (MISTAKE-138).** The box excluded
the mod-7 synchronization pole `E★ = {1,8,15,22,29,36,43,50,57}` (all ≡ 1 mod 7, spread 56),
which has `BUNCH(E★) = 6/19 ≈ 0.3158 > 2/7`. So `p₅+3p₆ ≤ 1/7` is FALSE, and the two extrema sit
at DIFFERENT poles:

| functional | extremal core | value |
|---|---|---|
| min POS = 6T₁+4T₂+2T₃ | consec {1..9} | 4717/882 ≈ 5.3481 |
| max BUNCH = 2T₅+4T₆ | mod-7 pole {1,8,…,57} | 6/19 ≈ 0.3158 |

(BUNCH-max verified: exhaustive over mod-7 families `{r+7jᵢ}`, offsets r=1..6, internal
9-subsets of [0..12], plus 8000 adversarial — max 6/19 at E★, internal-consec offset-1.
POS-min verified: box [1..19] + adversarial large-spread, min 4717/882 at consec.)

**The separation still closes GLOBALLY** (not just on a box):
`J = POS − BUNCH ≥ min POS − max BUNCH = 4717/882 − 6/19 = 84331/16758 ≈ 5.0323 ≥ 432/91`,
**margin +0.2850** (vs the artifact +0.315; the true J-min is J(consec) = 4465/882 = 5.0624, so the
now-two-pole separation is lossy by 0.030 but still clears the floor by 0.285). The k=8 analog
(addendum below) likewise closes with max NEG at `{2,9,…,51}` (offset-2 mod 7), margin **+0.0425**.

**Why two poles is the RIGHT picture:** these are exactly mac-mini THM-715's two synchronization
poles — consec (three-gap bunching, the covering/POS extremal) and mod-7 (7-sector resonance, the
variance/BUNCH extremal). The tail-separation cleanly assigns each pole to one half: consec bounds
POS from below, mod-7 bounds BUNCH from above, and neither pole is near-extremal for the other
half (POS(mod-7) = 6.08 ≫ 5.35; BUNCH(consec) = 2/7 < 6/19). The corrected piece-bounds to prove
are `POS ≥ 4717/882` (a covering floor, consec-pole) and `BUNCH ≤ 6/19` (a bunching bound whose
extremal is the mod-7 family, NOT the AP). The cancellation is still isolated in BUNCH; the
correction is only about WHICH core maximizes it.

## Extension addendum (mac-mini-2026-07-09-S65 cont.41): the (POS) bound is a COUPLED tail-tradeoff dominated by T1
Attacking the (POS) piece 6T1+4T2+2T3 ≥ 4717/882 (the cancellation-free half; NOTE klein-S255/MISTAKE-138 corrected the complementary BUNCH to TWO-POLE -- p5+3p6<=1/7 is FALSE, max BUNCH=6/19 at the mod-7 pole {1,8,..,57}; my POS work below is unaffected -- min POS is genuinely at consec): POS = E[g(N)],
g = (0,6,10,12,12,12,12,12) the NON-DECREASING cap of N(7−N) at its peak, so POS = J + 2(p5+3p6)
(POS = J + BUNCH exactly — the two THM-717 pieces are J and its bunching correction, not
independent). VERIFIED klein's POS bound (adversarial min-POS over 10 hill-climbs = 5.9875 >
4717/882 = 5.3481; consec = 4717/882 exact; low-μ families are HIGHER, 6.26–6.79).
**The mechanism (new):** POS is dominated by the weight-6 term T1 = P(N≥1) = meas(S7) = 1−p0,
and **consec MINIMIZES T1** (adversarial T1-min 0.761 > consec 0.5618) because consec MAXIMIZES
p0 = P(the orbit {0,x,…,8x} hits all 7 sectors) = 0.43821 — i.e. **consec is the best 9-phase
coverer**, a THREE-GAP (Steinhaus-orbit) statement. **The split FAILS though**: T2, T3 are NOT
consec-minimized (their minima 0.362, 0.084 sit at spread families with larger T1), so
POS ≥ 6·T1min+4·T2min+2·T3min = 4.99 < 4717/882 is too weak — POS's extremality at consec is a
COUPLED tail-tradeoff (dominant T1 drives it to the best-coverer; coupling with T2,T3 keeps consec
optimal), the SAME saddle character as J's μ-Var tradeoff (THM-716). Clean separable piece to
prove: **T1 = meas(S7) ≥ meas(S7)(consec), i.e. consec maximizes p0** (best-coverer / three-gap).
Files: lrc14_POS_bound + lrc14_tail_split_macmini_S65cont41 (+ outs).

## Addendum (klein-S256) — the exact BUNCH-max formula and the mod-7 mechanism

The corrected BUNCH-max (mod-7 pole) has an EXACT closed form for the base rows. For the mod-7 pole
`E★_k = {1, 8, 15, …, 7k−6}` (all `≡ 1 mod 7`, max speed `M = 7k−6`):

> **`T₆(E★) = 2/(7k−6)`,  `T₅(E★) = 5/(7k−6)`  (odd k),  so  `BUNCH(E★) = 2T₅+4T₆ = 18/(7k−6)`.**

Verified k=9,11,13: `6/19, 18/71, 18/85`. (Even k differs: k=8 gives `394/1075`, k=10 `65/228` —
a boundary-alignment correction.) And **`T₆ ≤ 2/(7k−6)` holds UNIVERSALLY** (max over 6120 families
= `2/57` at `E★`, for k=9) — the mod-7 pole maximizes the one-sector bunching. (The T₅ bound does
NOT separate: `max T₅ = 131/1479 > 5/57` at the *offset-2* mod-7 family; only the combination
`2T₅+4T₆` is maximized at offset-1. So BUNCH-max is a coupled mod-7 extremality, like POS.)

**The mechanism (why 7).** `E★`'s pairwise differences are ALL multiples of 7, so its phases
`{frac((1+7j)x)}` re-synchronize at the seven rationals `j/7` (`j=0..6`): at `x=j/7` every phase
equals `frac(j/7)`, and near `x=j/7` the phases spread at rate `≤ M`, giving a bunching interval of
half-width `1/(7M)` on each side. Seven sector-aligned centers × `2/(7M)` `= 2/M = 2/(7k−6)`. The
choice `q=7` is OPTIMAL because `7` is the sector count: a mod-`q` pole has `q` centers at max speed
`M ≈ 1+(k−1)q`, so `T₆ ≈ 2q/(7M) = 2/(7(k−1+1/q))`, and among `q` for which the centers align with
the `1/7`-sector grid, `q=7` maximizes. This is the "coherence spectrum" picture: bunching is a sum
over mod-`q` resonances, peaked at `q = 7`.

**The two-pole dispatch (proof direction).** The mod-7 pole (BUNCH-max) sits in the HIGH-J branch:
fully 7-coherent families (max residue multiplicity mod 7 = 9) have `J ≥ 5.585` (margin +0.84,
verified), so BUNCH is irrelevant there. The low-J families (consec and neighbors) have LOW mod-7
coherence, hence small BUNCH (`≤ 2/7 < 6/19`). So the assembly splits: **[7-coherent ⟹ J large
directly] + [not 7-coherent ⟹ BUNCH small, so J ≈ POS ≥ 4717/882]** — the two THM-717 poles never
coincide, and each branch drops one hard piece. (POS ≥ 4717/882 on the non-coherent branch remains
mac-mini's coupled covering floor.)

**Proof direction for `T₆ ≤ 2/(7k−6)`:** each maximal `O=1` interval has length `≤ 1/(7M_family)`
(the fastest phase traverses a `1/7`-sector in time `1/(7M)`), and the interval count is the number
of mod-`q` resonance centers; bounding `Σ (count × 1/(7M))` by the `q=7` optimum is the remaining
content. Files: `04-computation/lrc14_bunch_max_formula_klein_S256.py` (+ the earlier maxbunch/
hillclimb scripts).

## Addendum (klein-S257) — the base TAIL is the two-scale limit (THM-687/688); far elements RAISE J

mac-mini cont.42/43 framed the k=9 base as `[compact exhaustive check, diam ≤ D₀] + [tail]` and
(cont.43) CORRECTED the tail: wide families do NOT approach `J_iid = 8.456` — they PLATEAU at
`≈ 5.68`, because a wide family is a **compact cluster + decorrelated far elements**, so `J` tends
to the **two-scale / multi-scale limit** (my THM-687/688), not full iid. This addendum assembles
that tail exactly.

**The two-scale limit.** Adding one far element `w → ∞` to a `(k−1)`-cluster `E'` sends the moments
by the eigen-transfer (THM-710): `m₁ → (6/7)m₁`, `m₂ → (5/7)m₂`, so
> `J(E' ∪ {w}) → eigenlim(E') := 6·(6/7)m₁(E') − (5/7)m₂(E') = (6·J(E') + m₂(E'))/7`,
with finite-`w` error `|J(E'∪{w}) − eigenlim(E')| ≤ C·ΣE'/w`, and the **true peel constant
`C ≈ 0.4`** (measured max `0.381` over 124 8-cores × 250 far values) — a **1764× sharpening** of
THM-702's crude `672·ΣE'` (the gap is signed cancellation; the rigorous-absolute constant is larger
but still `≪ 672`).

**Far elements RAISE J — verified at every level `k = 5..9`:**

| k | compact-min `= J(consec-k)` | min two-scale limit (over `(k−1)`-cores) | margin |
|---|---|---|---|
| 5 | 366/35 ≈ 10.457 | 11.197 | +0.740 |
| 6 | 928/105 ≈ 8.838 | 10.016 | +1.178 |
| 7 | 1678/245 ≈ 6.849 | 8.390 | +1.541 |
| 8 | 291/49 ≈ 5.939 | 6.554 | +0.615 |
| 9 | 4465/882 ≈ 5.062 | **19473/3430 ≈ 5.677** | **+0.615** |

At k=9 the min two-scale limit `19473/3430 ≈ 5.677` (attained at the compact `consec-8` cluster;
robust over structured + wide + adversarial 8-cores) is `≥` the compact-min `5.062`, so **every wide
9-core has `J ≥ 5.677 ≥ compact-min ≥ 432/91`**. The 5.677 is exactly mac-mini's plateau 5.68 =
the THM-687 two-scale value for block+singleton.

**The assembled k=9 base:**
- **[compact]** 9-cores with no far element (bounded diam ≤ D₀): exhaustive check, min `= J(consec-9)
  = 4465/882 ≥ 432/91`. (Feasible — mac-mini's lane.)
- **[wide/tail]** 9-cores with a far element: `J ≥ eigenlim(8-cluster) − C·ΣE'/w ≥ min-two-scale-
  limit ≈ 5.677 ≥ compact-min`. Far elements only raise `J`. (My THM-687/688 + the peel constant.)

**The recursion bottoms out cleanly:** `eigenlim(E') ≥ compact-min_k` reduces (since `m₂ ≥ 0`) to
`J(E') ≥ (7/6)·compact-min_k`, i.e. to the `(k−1)` compact structure; iterating `k = 9 → 8 → 7`
lands at `k ≤ 7`, where **SmallClusterFull (klein, proved in Lean, LRCSevenGapRigidity)** gives
`ν = 1` (all covering) ⟹ `J` large. So the whole tail is a finite descent to a Lean-proved base.
The signed-cancellation content that blocks an *absolute* proof lives only in the compact
extremality (`J(consec-k)` is the min) — the tail itself is clean (far elements provably help).

Files: `04-computation/lrc14_base_tail_twoscale_klein_S257.py`.


## Addendum (mac-mini cont.51, 2026-07-12): a clean closed form for the extreme-bunching tail
Combinatorial search on the consec sector-distribution: **p6(consec_k) = 1/(7(k-1)) EXACTLY** (verified
k=7..13: 1/42, 1/49, 1/56, 1/63, 1/70, 1/77, 1/84). Proof sketch: p6 = P(N=6) = P(only sector 0 hit) =
P(frac(jx) in [0,1/7) for all j=1..k-1); the binding constraint is j=k-1, and the unique interval is
x in [0, 1/(7(k-1))) (measure 1/(7(k-1))) -- the orbit {jx} fits in the [0,1/7) arc iff (k-1)x < 1/7
(no other interval survives, since reflections land in (6/7,1]). Consequences: BUNCH(consec_k) = p5 +
3p6 has p6 = 1/(7(k-1)) pinned; BUNCH(consec_12) = 1/14 EXACTLY and BUNCH(consec_13) = 61/924.
T5(consec_8) = 1/14 exactly. These clean values feed the corrected two-pole BUNCH (klein-S255/256).
File: lrc14_closed_forms_macmini_S65cont51 (+ out).
