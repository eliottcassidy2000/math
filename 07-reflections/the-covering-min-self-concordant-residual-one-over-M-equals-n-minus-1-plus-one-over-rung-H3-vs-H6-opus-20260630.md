# The covering-min's SELF-CONCORDANT RESIDUAL: 1/M = (n−1) + 1/rung — the loneliness reciprocal splits EXACTLY into the barrier parameter ν=n−1 (the danger-arc count) plus an interior-point residual 1/rung; the tight AP is the analytic center (residual 1, rung 1), coverings are pushed off it (residual ≤ 1/2, rung ≥ 2), and the covering-min is the LARGEST-residual reachable covering (smallest rung a(n)); H3 (residual stays Θ(1), beaters persist for all n) vs H4 (residual → 1/n, construction wins at n≥12) is exactly whether the interior-point residual stays bounded or vanishes — I confirm H3 independently at n=7,8 (covering-min = 2/13, 2/15, below the construction), show the known beaters n=7..10 VIOLATE the lowness lemma HYP-3747 (miss core, yet M<n/Φ₆ ⇒ HYP-3747 is n-dependent, holds n=14, fails n≤10), verify the H6 window (margin Θ(1/n²), cyc/hex ratio→2, and the construction-rung margin IS mac-mini's Dedekind sum HYP-3768), and a structured n=12,13,14 band search finds no beater (honestly uninformative — search-resistant, as klein warned)

*opus-2026-06-30. Owner: study the self-concordant barrier residual and finite-periodic-object work; work H3 and
H6 back and forth. The self-concordant lens on the covering-min turns out to be exact algebra (1/M = ν + 1/rung),
and it reframes the H3/H4 open edge as "does the interior-point residual stay Θ(1) or vanish."*

## 1. The self-concordant residual identity (EXACT — the centerpiece)
Every primitive covering set sits on the Stern–Brocot ray `M = [0; n−1, k] = k/(k(n−1)+1)` at rung `k`
(klein HYP-3732). Invert:
> **1/M = (n−1) + 1/k.**
This is a **self-concordant "barrier parameter + residual" split.** The canonical self-concordant barrier is
`−log` (Nesterov–Nemirovskii, ν=1); the covering LOG-BARRIER `F = −∫₀¹ log m_w(t) dt` (my HYP-3765 object,
`m_w` = danger cover) is a superposition of `−log`s over the **n−1 danger arcs**, so its natural barrier
parameter is `ν = n−1`. The identity says the loneliness reciprocal `1/M` equals that parameter plus a residual
`1/k`. Verified exactly (fractions):

| n | covering-min M | 1/M | = (n−1) + residual | rung k |
|---|---|---|---|---|
| 7 | 2/13 | 13/2 | 6 + **1/2** | 2 |
| 8 | 2/15 | 15/2 | 7 + **1/2** | 2 |
| 9 | 4/33 | 33/4 | 8 + **1/4** | 4 |
| 10 | 4/37 | 37/4 | 9 + **1/4** | 4 |
| 11 | 3/31 | 31/3 | 10 + **1/3** | 3 |
| 14 construction | 14/183 | 183/14 | 13 + **1/14** | 14 |

**The AP is the analytic center.** The LRC floor `1/n` is rung 1: `1/M = (n−1) + 1 = n`, residual **1** — the
deepest interior point (the tight AP, unconstrained loneliness optimum). Coverings are pushed OFF it: `M > 1/n`
strictly (THM-523) ⟹ rung ≥ 2 ⟹ residual ≤ **1/2**. The **covering-min is the largest-residual reachable
covering** — the point closest to the AP center that the covering constraint allows — at residual `1/a(n)`,
`a(n)` = the smallest realizable rung.

## 2. H3 confirmed independently (n=7,8) + the known beaters BREAK HYP-3747
Independent exact search (min `M(S) > 1/n` over (n−1)-subsets, rational-t sweep): **C(7) = 2/13, C(8) = 2/15**
— both **below** the construction (7/43, 8/57). So beaters exist; **my own earlier reflection's "C(n)=n/Φ₆
exact for n≥7" was WRONG** (bounded search missed them, exactly the mirage H3 names). klein's H3 vindicated at
n=7,8.

**Sharper:** the four known beaters (n=7,8,9,10 = {1,2,5,6,7,8}, {1,4,5,6,7,11,16}, {1,3,4,5,7,11,18,32},
{1,2,3,5,6,7,8,9,30}) all **miss core elements** {1,…,n−2} (missing {3,4},{2,3},{2,6},{4}) yet have
`M < n/Φ₆`. So they are **counterexamples to the lowness lemma HYP-3747** ("M ≤ n/Φ₆ ⟹ {1,…,n−2}⊂S") at
n=7,8,9,10. HYP-3747 is **n-dependent**: verified at n=14 (klein, single-core-drop), **false at n≤10**
(multi/single-drop spread beaters). This is precisely the "challenges the evidential basis of HYP-3737/3747"
that H3 flags — now made concrete: the lemma's conclusion literally fails below the (alleged) transition.

## 3. H6 window — and the margin IS a Dedekind sum (mac-mini HYP-3768)
The margin above the floor at rung k:
> **M_k − 1/n = (k−1)/(n(k(n−1)+1))**, residual-form `= (1−1/k)·(1/M_k)/n`.
Pinned between (verified n=7..20, ratio → 2, "never closes, never opens"):
- **hexagonal edge** (rung 2, max residual 1/2): `1/(n(2n−1))` — the smallest margin any covering can have (H5's
  LRC floor-margin);
- **cyclotomic edge** (rung n, the construction): `(n−1)/(n·Φ₆)`.
Both `Θ(1/n²)`; `cyc/hex = (n−1)(2n−1)/Φ₆ → 2`. **New tie-in (mac-mini HYP-3768):** the *cyclotomic-edge*
(construction) margin `(n−1)/(nΦ₆) = −12·s(n,Φ₆)/n²` is the **observer's Dedekind sum**, positive exactly
because n has order 6 mod Φ₆ (hexagonal/Eisenstein); an order-4/square would give `s=0`, zero margin. So the
two endpoints of the H6 window are the two figurate arithmetics — **hexagonal** `n(2n−1)` and **cyclotomic**
`Φ₆=n²−n+1` — and the residual `1/rung` interpolates between them. The self-concordant residual (§1) and the
Dedekind margin (HYP-3768) are the reciprocal and additive readings of the same `Θ(1/n²)` gap.

## 4. H3 ↔ H6, back and forth (the synthesis)
- **From H6 to H3:** the window `[hexagonal, cyclotomic]` is a window of *residuals* `1/k ∈ [1/n, 1/2]`. The
  covering-min sits at residual `1/a(n)`. H3 says `a(n)` stays small (2,2,4,4,3,…) ⟹ residual stays **Θ(1)**
  (≈1/2–1/4, near the AP center 1); H4 says `a(n)→n` at n≥12 ⟹ residual **→ 1/n** (the construction, far from
  center). So **H3 vs H4 is exactly: does the interior-point residual stay bounded away from 0, or vanish
  like 1/n?**
- **From H3 to H6:** if H3 holds, the covering-min hugs the **hexagonal** end (small rung), and its margin is a
  *small-rung* Dedekind-type sum `(a−1)/(n(a(n−1)+1)) ~ (1−1/a)/n²`, NOT the full construction Dedekind sum
  `(n−1)/(nΦ₆)`. So H3 predicts the true LRC14 margin is **strictly smaller** than mac-mini's `13/2562`
  (it would be `(a(14)−1)/(14(13a(14)+1))`, e.g. `1/378` if `a(14)=2`) — the construction's Dedekind sum is
  only the *loose* (cyclotomic) end. **The Dedekind-margin and the covering-min disagree by exactly the
  rung gap `a(14)` vs `14`** — the same open edge, now in Dedekind-sum language.

## 5. The band search (honest null) and finite periodic objects
A **structured** search (sample from the rung-k middle band [k, D−k] mod D=k(n−1)+1, the actual beater
geometry, covering-checked) over ~20k covering candidates each at n=12,13,14 found **0 beaters**. But this is
**uninformative**: random-in-band hits a specific search-resistant beater with probability ~1/C(D−2k, n−1)
(≈1/C(26,8) for rung-2 n=14) — my method, like klein's random/hillclimb/drop, cannot reliably surface them
(it recovers the M-values of the *known* beaters when handed them, but not by sampling). So the n≥12 null is
consistent with **both** H3 (search-resistant) and H4 — no evidence either way, exactly as klein argued.
**Finite-periodic-object tie-in:** the danger cover `m_w(t)` is a finite periodic object; its log-barrier
(§1) is the self-concordant barrier, its Fourier/circulant spectrum is the Riesz-factorization view (HYP-3765),
and the covering-min is its extremal packing — the rung `a(n)` is the integer-realizability obstruction that a
*continuous* self-concordant program would not see (the continuous optimum is the loose ceiling `1/(n−1)`,
residual→0; the integer/covering constraint forces residual ≥ `1/a(n) > 0`).

## Status
- **Confirmed (opus, exact):** the self-concordant residual identity `1/M = (n−1) + 1/rung`; the covering-min
  at n=7,8 is `2/13, 2/15` (below construction ⟹ H3 holds there, my old n≥7 claim corrected); the H6 window
  (margin Θ(1/n²), cyc/hex→2).
- **Sharp (opus):** the known beaters n=7,8,9,10 violate the lowness lemma HYP-3747 (miss core, M<n/Φ₆) — it is
  n-dependent (holds n=14, fails n≤10), concretely challenging the HYP-3737/3747 evidential basis (H3).
- **Framing (opus):** H3 vs H4 ⟺ the interior-point residual `1/a(n)` stays Θ(1) vs vanishes like 1/n; the H6
  endpoints are the hexagonal (residual 1/2) and cyclotomic (residual 1/n, = mac-mini's Dedekind sum HYP-3768)
  arithmetics; if H3, the true LRC14 margin is smaller than the construction's Dedekind `13/2562`.
- **Open (honest):** a(n) for n≥12 (the open edge); the structured band search null is uninformative
  (search-resistant). The self-concordant barrier is an exact-identity LENS, not (yet) a computed program;
  the residual is the CF/Farey tail (klein HYP-3732), rigorous; the "ν = n−1 barrier parameter" reading is
  the interpretive layer.

Related: HYP-3764/klein-S53 (the open edge; H3/H5/H6), HYP-3768/mac-mini-S64 (margin = Dedekind sum, the
cyclotomic endpoint), HYP-3737 (the alleged n≥12 transition, H4), HYP-3747 (lowness lemma — shown n-dependent),
HYP-3765 (the log-barrier = self-concordant barrier; Riesz/circulant spectrum), HYP-3732 (Stern–Brocot ray),
HYP-3726 (hexagonal margin), THM-523 (rung ≥ 2), OPEN-Q-108. HYP-3769 (this).
Scripts: 04-computation/lrc_covmin_selfconcordant_residual_opus_20260630.py.
