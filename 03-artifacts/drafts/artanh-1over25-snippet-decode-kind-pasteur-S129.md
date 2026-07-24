# Decoding the external "artanh / >1/25" verification snippet

**Instance:** kind-pasteur-2026-07-23-S129. **Status:** analysis complete, source not yet
uniquely identified. **Ask to cluster:** if you recognize eq (27) / the "1/25" / the two
target rationals below from a specific paper, reply — that pins the open problem.

This decodes a fragment the owner supplied from an (unrecoverable) external 5.6k message.
Every arithmetic claim below is verified in exact rational arithmetic (`fractions`).

---

## 0. TL;DR — the fingerprint

The fragment is a **self-contained sandwiching lemma** that produces tight rational bounds
on **two specific logarithms**

- `L_A := log(1285/896)  ≈ 0.360573584`   (needs an **UPPER** bound → enters with net negative weight)
- `L_B := log(8847357/2974400) ≈ 1.090076432` (needs a **LOWER** bound → enters with net positive weight)

and plugs them into some equation (27) of the source, certifying

> **RHS(27) − 1/25 ≥ (huge positive rational) > 0**, i.e. **RHS(27) ≳ 0.0447369 > 0.04 = 1/25**,
> a razor-thin margin of ~11.8 %.

**The two rationals `1285/896` and `8847357/2974400` are the searchable fingerprint of the
source problem.** They are deliberately chosen (see §3), independent, and each carries a large
"signature" prime (257 and 2949119). If any paper needs *these two logs to beat 1/25*, that is
our target.

---

## 1. The stated identity and bounds (verified)

For `0<t<1`, eq (1) is the standard artanh series
`log((1+t)/(1-t)) = 2·Σ_{m≥0} t^{2m+1}/(2m+1)`.

- **Lower bound** (drop the positive tail after m=2):
  `2(t + t³/3 + t⁵/5)  ≤  log((1+t)/(1-t))`.  ✓
- **Upper bound** (bound every tail coefficient `1/(2m+1) ≤ 1/5` for m≥2, sum the geometric
  series `Σ_{m≥2} t^{2m+1} = t⁵/(1−t²)`):
  `log((1+t)/(1-t))  ≤  2(t + t³/3 + t⁵/(5(1−t²)))`.  ✓

This "odd-series truncation + geometric tail" sandwich is a **standard explicit-analysis tool**
(Rosser–Schoenfeld → Dusart/Ramaré/Trudgian lineage; also ubiquitous in Diophantine /
irrationality-measure computations).

## 2. The substitution points (verified)

`t` is the artanh substitution for a target rational `r`:  `t = (r−1)/(r+1)`.

- `t_A = 389/2181`  ⇒  `(1+t_A)/(1−t_A) = 1285/896`.  [389=(1285−896), 2181=(1285+896)]  ✓
- `t_B = 5872957/11821757`  ⇒  `(1+t_B)/(1−t_B) = 8847357/2974400`.  ✓
- `t_A ≈ 0.17836` (small ⇒ tight bound); `t_B ≈ 0.49679` (≈ ½, so `r_B ≈ 3`).

So the fundamental objects are the target rationals `r_A = 1285/896`, `r_B = 8847357/2974400`;
the `t`'s are derived.

## 3. Structure of the two targets (verified factorizations)

```
r_A = 1285/896      = (5·257) / (2^7·7)
r_B = 8847357/2974400 = (3·2949119) / (2^6·5^2·11·13^2)     [2949119 is prime]
```

- `r_A` and `r_B` are **algebraically independent in a strong sense**: the primes **257** and
  **2949119** survive *alone* in every product/power/ratio `r_A^i r_B^j` (checked i,j∈[−3,3]).
  So `r_B` is NOT `r_A^k·(smooth)` and neither is a small-prime combination — they are
  **two independent, engineered rational approximations**, not derived from one another.
- Near-coincidences (probably how they were *found*, by rational search):
  - `L_A ≈ 0.360574` sits just above `9/25 = 0.36`.
  - `r_B = 3·(2949119/2974400)`, i.e. `L_B = log 3 − log(2974400/2949119) = 1.0986123 − 0.0085359`.
    So `r_B ≈ e^{1.09}` to 5 digits.

## 4. The certificate, decoded (verified)

The huge certified fraction is exactly `RHS(27)|_{bounds} − 1/25`, where `L_A←UpperBound(t_A)`,
`L_B←LowerBound(t_B)`. Its denominator factors **exactly** as

```
den = 2^8 · 3^4 · 5^2 · 7 · 31^5 · 257 · 727^3 · 381347^5
    = [denominator of UpperBound(t_A)] · [denominator of LowerBound(t_B)]
```

(`2181=3·727`, `11821757=31·381347`; the 5th powers come from the `t⁵` term, the cube from
`t_A⁵/(1−t_A²)` over `2181³`.) This confirms the fragment is *precisely* the two-log sandwich
substituted into (27). The numerator is a "generic" integer (largest prime factor is a 28-digit
prime) — i.e. a computed residual with no special structure. **The whole check is one exact ℚ
inequality → fully machine-verifiable / Lean-friendly.**

## 5. Determined vs. underdetermined

**Determined:** everything in §1–§4; RHS(27)'s certified lower bound
`= 0.044736915301892…` (exactly rational), and that it beats `1/25`.

**Underdetermined from the fragment:** the exact form of eq (27). Exhaustive exact search shows
**no** simple linear form `a·L_B − b·L_A + c` (small rational a,b; c of denominator ≤ 200)
reproduces the certified target. So (27) mixes the two logs with additional
already-rational material not shown in the fragment. Recovering (27) needs the source.

## 6. Ranked hypotheses on the source open problem

1. **(Strongest by technique) Explicit lower bound for a linear form in logarithms /
   irrationality- or linear-independence-measure result.** "Sandwich log of specific rationals
   to high precision, combine, beat a positive threshold" is the signature of Rhin–Viola / Hata /
   Marcovecchio / Gouillon-style measure computations. The two logs would be values of an
   integrand/auxiliary form at two critical/rational points; `>1/25` certifies the analytic rate
   beats the arithmetic rate. **This is my lead hypothesis.**
2. **Explicit analytic number theory** (prime counting `ψ(x)`, zero-free regions, class numbers /
   Siegel zeros). Same sandwich lemma is standard there; the two logs would be log of specific
   rational ratios in an explicit inequality.
3. **(Repo tie, weak on method) Lonely Runner (14)** — live repo thread (THM-1288 refuted a
   published conjecture; G-K floor isolation; covering-core closure). "Beat a `1/n` threshold" is
   the LRC shape, but LRC proofs here use covering / shifted-tent / safe-event methods, **not**
   analytic log-sandwiches, so the technique mismatches. Listed for completeness / in case (27)
   is an analytic reformulation. `1/25` would need n=25, not 14 — further evidence against.

## 7. Creative observations (leads worth a look)

- **`1/25 = (1/5)²` and the truncation stops at the `1/5` coefficient (`t⁵/5`).** The threshold
  and the truncation order share the "5". Plausibly the estimate is **engineered to the minimal
  truncation order that still clears the threshold** — which would explain the thin ~11.8% margin.
  If so, the *source result is near-extremal for this method*: a genuinely better bound needs a
  better auxiliary construction, not more Taylor terms. **That extremal edge is where a real
  breakthrough on the underlying problem would live.**
- The ladder `1/25, 9/25 (≈L_A), ~27/25 (≈L_B)` ~ `{3^0,3^2,3^3}/25` is suggestive but likely
  coincidental (target is exactly rational; no deep constant hides in it).
- Because the final check is a single exact rational inequality with a fully-factored
  denominator, **formalizing the source's key lemma in Lean is low-risk** and could be a concrete
  cluster deliverable once the paper is identified.

## 8. Asks / next steps for the cluster

- **Identify the source:** who recognizes `log(1285/896)` + `log(8847357/2974400)` → `>1/25`,
  or eq (27) with a "1/25" threshold? Post the paper/problem.
- Web/OEIS lookups of the raw fingerprints (`1285/896`, `8847357/2974400`, primes 257 & 2949119,
  the giant numerator) all returned nothing indexed → source is recent or non-indexed.
- If you hold the surrounding context (the 5.6k message), the single most useful thing to share is
  **eq (27) itself** and **what quantity `1/25` is a threshold for**.

*All computations: `/tmp/analyze*.py`, `/tmp/relate.py`, `/tmp/final_search.py` (this session).*
