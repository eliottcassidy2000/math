---
source: claudebox-2026-06-02-S561
status: result (exact density of the LRC sieve-covered core + gap distribution; refines S557/HYP-2064)
tags:
  - lonely-runner
  - gap-bound
  - sieve
  - core-density
  - inclusion-exclusion
  - near-AP-slice
  - s557-followup
---

# The LRC sieve-covered core is not measure-zero — it has positive density; the *near-wall slice* inside it is what's thin

## Setup

oracle-S557 (HYP-2064) proved the repo's sieve gives the **full** lonely-runner conjecture
`g(v) := sup_t min_i ||v_i t|| >= 1/(k+1)` for every **non-sieve-covered** speed set, leaving
only the **sieve-covered** sets (a multiple of every `q in {2..k+1}` is present) as the locus
where Tao's general bound is the operative one. S557 described that core as "measure-zero-ish,"
with the genuinely hard part in an "even thinner near-AP slice."

This session makes both halves quantitative at `n = k+1`. (`g` is dilation-invariant,
`g(cv)=g(v)`, so primitive sets are WLOG.)

## Result 1 — the core has POSITIVE limiting density (exact)

The `R -> inf` density that a random primitive `k`-set is sieve-covered is closed-form by
inclusion-exclusion over `A_q = "no speed divisible by q"`:

```
P(covered) = 1 - sum_{∅≠S⊆{2..k+1}} (-1)^{|S|+1} rho(S)^k,
rho(S) = density of integers divisible by no q in S = sum_{T⊆S} (-1)^{|T|}/lcm(T).
```

| n=k+1 | exact P(sieve-covered) |
|------:|-----------------------:|
| 10    | 0.19753 |
| 12    | 0.15278 |
| 14    | **0.11203** |
| 16    | 0.09361 |

Monte-Carlo over growing speed ranges converges to the limit (n=14: 0.1016 → 0.1100 → 0.1109
at R=200/1000/5000), and the seed used by S557 reproduces its **2/40** headline exactly.
So at n=14 roughly **1 in 9** random primitive 13-sets is sieve-covered — a positive-density
set, **not** measure-zero. The density *rises* as n shrinks.

## Result 2 — generic core sets are FAR from the wall

Over 200 random sieve-covered primitive 13-sets at n=14, the gap distribution is
`min = 0.138, median = 0.201, max = 0.307` against the conjecture `1/14 = 0.0714`. The
**minimum** gap is `1.93×` the conjecture; **none** of the 200 lands within `1.15×` of `1/14`.
The wall `AP = {1,...,13}` is the one tight core set (`g = 1/14`); generic core sets miss it
by 2–4×.

## Reading

The "thin core" intuition is correct about the **extremizers** — the near-wall slice
(`g ≈ 1/(k+1)`) is a genuinely thin sub-locus that demands the rigid AP / dilated-AP structure —
but it **overstates the thinness of the sieve-covered set itself**. The honest picture at `n=14`:

- proven bulk (~89%): non-sieve-covered → `g ≥ 1/(k+1)` by the sieve (S557/THM-369);
- a **positive-density** (~11%) sieve-covered core, almost all of which still has `g ≈ 2–4 × 1/(k+1)`;
- a **thin near-AP slice inside the core** carrying the actual hardness (`g → 1/(k+1)`).

**Consequence for proof strategy.** A density/measure argument ("the open core is negligible")
cannot finish the problem — the core is not negligible. What is thin is cut by *structure*, not
*size*: the near-AP slice must be excluded by the AP/dilation machinery already in the repo
(S530 apex / co-observer, S556 local LP, S550 θ-deformed measure bound), applied **on** the
positive-density core rather than hoping the core itself vanishes.

## Follow-ups (small, handoff-sized)

1. **θ-deformed bound on the core (S557 Prediction 1).** For sieve-covered sets compute the
   largest `θ < 1/n` with the S550 energy `E(v) < (1-2θ)^k`, giving an explicit `g ≥ θ` off the
   AP. Does it cover the whole core minus the near-AP slice?
2. **AP-distance vs gap (S557 Prediction 2).** Define a dilated-AP distance and confirm small
   gap ⟺ near a dilated AP across the core (the 200-set data already shows the small-gap tail
   trends toward AP-like structure — quantify it).
3. **Density asymptotics.** Fit/derive `P(sieve-covered)` as `n → ∞`; the prime-power
   constraints `q ∈ {11,13,...}` dominate — does the core density → 0, and how fast?
