---
source: opus-2026-07-23-S4 (decoding the owner's external artanh snippet, session 4)
status: >
  The snippet's object is now PINNED EXACTLY and uniquely. This reflection reconciles
  the opus-S3/mac-mini-S168 tension, records the structural anchor (arctanh uniqueness),
  logs the ruled-out homes (figurate masses, named constants, hypergeometric, and BOTH
  candidate LRC papers by direct fetch), and adds a new coefficient fingerprint that ties
  BOTH tight LRC(14) families. Home narrowed to the log-energy-beats-floor family with
  LRC(14) wider-gap (n=13) leading; irrationality-measure not fully killed. Honest scope:
  we identify the FAMILY, not a citation; the integers are source-generated.
tags: [decode, artanh, rapidity, linear-form-in-logs, lrc14, wider-gap, log-energy, thm-252, hook-B, second-moment, honest-uncertainty]
related: [HYP-9023, THM-2000, THM-252, HYP-1992, HYP-1892, HYP-3114, CONSTANTS-INDEX]
---

# The two-log artanh form, pinned and homed (opus-S4)

## 1. The object is now exact and unique

Independent re-verification (`04-computation/snippet_eq27_reconcile_opus_S4.py`,
output saved) confirms mac-mini-S168 to all 51 digits and *upgrades* it from "a fit"
to a **uniqueness proof**:

```
RHS(27) = (2457/6592)·log(8847357/2974400) − log(1285/896),   certified > 1/25.
(2457/6592)·L_B − U_A − 1/25 == G   (residual EXACTLY 0), where
  U_A = upper artanh-sandwich at t_A=389/2181, L_B = lower sandwich at t_B=5872957/11821757.
Integer form:  2457·log_B − 6592·log_A > 6592/25 = 263.68.
```

**Why it is unique (resolves opus-S3's "c unpinnable").** Set the target
`T = 1/25 + G`. For each integer `d` on `log_A`, solve `c=(T + d·U_A)/L_B` and measure
the height of `c`:

| d | −2 | −1 | 0 | **+1** | +2 | +3 |
|---|----|----|----|--------|----|----|
| log₁₀ height(c) | 52.5 | 52.7 | 53.0 | **3.8** | 53.0 | 52.7 |

Only `d=+1` collapses the height from ~10⁵³ to 10³·⁸ (c = 2457/6592). A 49-order-of-
magnitude cliff is not a coincidence: it **pins** `(c, d, r) = (2457/6592, 1, 0)`. My
S3 "unpinnable" was an artifact of letting the rational part `r` float; with `r=0`
enforced by the height cliff, the coefficient is forced. The `31⁵·381347⁵` I saw in S3
is (mac-mini is right) an **artifact of the t⁵/5 truncation order**, not the problem.

## 2. The structural anchor: arctanh is not incidental

`03-artifacts/substack/hook-B-arctanh-unique.tex` (repo canon):

> **arctanh is the unique odd formal power series `f` with `e^{f}` a degree-(1,1)
> rational** `(1+at)/(1−at)`.

The snippet's ratio `(1+t)/(1−t) = S_n/S_{n−1}` is *exactly* a degree-(1,1) rational, so
the arctanh is **forced** by the "ratio of two linear quantities" structure — it is not a
stylistic choice. This dissolves the three-way naming dispute: "rapidity" (THM-252),
"Abel–Dini log-edge" (THM-2000 §3.1), and "logit/entropy" (opus) are the SAME unique
object. THM-252 places `2·artanh(rational) = log((m−k)/k)` in the `⊕_p Q·log p` rapidity
lattice; the snippet's two logs live in `{2,5,7,257}` (side A) and `{2,3,5,11,13,2949119}`
(side B). The Collatz-rapidity thread (`collatz-rapidity-defect.md`) shows the repo already
certifies linear forms `K log2 − L log3 − D > 0` with exactly this machinery — the snippet
is one more instance of the repo's "arctanh = universal linearizer of the Cayley formal
group" thesis.

## 3. What it is NOT (ruled out this session)

- **Not figurate masses (THM-2000).** No pairwise difference of the low masses `M(s,d)`
  hits `log_A=0.3606`, `log_B=1.0901`, or the LHS `0.04572`; both logs sit *below* the
  entire mass band `[1.16, 2.0]`. Answers mac-mini's "is (A,B) a ratio of two masses?": **no.**
- **Not a named constant / its convergents.** `A,B` are not close to `√2, ∛3, e, π, φ`.
- **Not hypergeometric.** `896,1285,2974400,8847357` are not binomials, central binomials,
  factorials, or `lcm(1..n)` — so **not** a standard Apéry/Padé (binomial-based) construction.
  This is the main count *against* the irrationality-measure family: those numbers are almost
  always binomial/lcm products; here they are generic (one 7-digit prime `2949119` on side B).
- **Not either candidate LRC paper (checked by direct fetch).** Bedert 2511.16636 ("Riesz
  products & LRC: a wider gap") has eq (27) = `ML(V) ≥ ML(V″) − 1/mᵉ − q/p` (model reduction,
  **no logs**); it is asymptotic `1/(2n)+1/n^{5/3}`, carries no 51-digit certificate. The
  "eleven, twelve, thirteen runners" paper (2604.23906) uses thresholds `1/13`, never `1/25`,
  and has no artanh sandwich. **We identify the method-family, not a source citation.**

## 4. New fingerprint: the weight ties BOTH tight LRC(14) families

klein-S404 noted the numerator `2457 = 3·S₂({1..13}) = 3·819` is the AP core's second
moment (unique to n=13). I add the **denominator** (`snippet_eq27_lrc_families_opus_S4.py`):

```
6592 = 2⁶ · 103,   and   103 = S₁(GW) = Σ{1..11,13,24}   (the OTHER tight family's total speed).
Moreover 103 | x_n(B) = 5872957 = 19·103·3001.
So   2457/6592  =  3·S₂(AP) / (2⁶·S₁(GW)).
```

Both tight LRC(14) extremals (AP `{1..13}` and Goddyn–Wong `{1..11,13,24}`) appear in the
weight — AP through its 2nd moment, GW through its 1st moment. This is a candidate answer to
klein's open "why is `3·Σv²/denominator` the natural balance weight." **Caveat (mac-mini's
meta-warning):** the prefactors `3` and `2⁶` are unexplained, and `2457=3·819` is equally the
figurate fact `3·P₃(13)`; this could be a coincidental factorization of an optimizer output.
It should be *derived* (klein's task) or *reproduced* by a 13-speed Riesz optimization
(mac-mini's task), not asserted.

## 5. Home: two live families, LRC(14) wider-gap leading

Converged cluster taxonomy (mac-mini's discriminators + klein-S404 + this session):

- **(A) log-energy-beats-floor** [loneliness / min-cosine / merit-factor / Cohn–Elkies LP].
  `X=0.04572` is a genuinely-positive floor-beating quantity (no near-cancellation:
  `log_A/log_B=0.331 ≠ c=0.373`), the arctanh odd series is a per-mode log-energy sum, and the
  LRC tent `‖x‖` has an odd-harmonic Fourier series. **Leading sub-home: LRC(14) wider-gap,
  n=13** (klein: 2nd-moment fingerprint + `1/25=1/(2n−1)` + value-band `1/26<0.0457<1/14`;
  reinforced here by `6592=2⁶·S₁(GW)`).
- **(B) irrationality-measure via clean Padé/convergent** — *not* killed (mac-mini: a measure
  is itself a clean-log rate comparison, and the big-prime fingerprint fits), but weakened by
  the non-binomial integer structure (§3).

**Honest uncertainties I will not paper over.** (i) The threshold reading is soft: `X≈0.0457`
also exceeds `1/24=1/(2·12)` and `1/22`, so "1/25" as `1/(2n−1)|₁₃` is a *choice*, not forced —
mac-mini's `X~1/21.9 ⇒ n~11` is an equally literal reading. (ii) We are 4-for-4 mapping onto
repo pets (rapidity → Riesz-LRC → Abel–Dini → tight-family moments); the owner warned it may be
*any* problem. The coefficient fingerprints are the strongest LRC evidence; the arguments `A,B`
do **not** come from these families' power sums, which is unexplained.

## 6. Forward (turning decode into progress on the prize)

klein's template idea is the right one and is the repo-flagship payoff: use the certified
artanh log-energy bound (`certified_logratio_abeldini_opus_S2.py`, already float-free) as a
**Lean-ready certificate style for `inf L>0` / a `>1/(2n−1)` wider-gap bound** on the 13-speed
tight cores `d·{1..13} ∪ {r}`. Decisive discriminator (mac-mini's sub-task): do the amplitudes
`p_A=1285/2181, p_B=8847357/11821757` (or `A,B`) arise as Riesz-integral ratios `∫ΦR/∫R` for a
13-speed core? If yes, the ID is confirmed and we inherit a proof-template; if no after a fair
search, family (A) via LRC weakens and (B) or an off-repo log-energy problem (merit-factor,
Cohn–Elkies) should get a real look.

Artifacts: `snippet_eq27_reconcile_opus_S4.py`, `snippet_eq27_fingerprint_opus_S4.py`,
`snippet_eq27_lrc_families_opus_S4.py` (+ outputs); updates HYP-9023.
