# Snippet context CONFIRMED: certificate (27) belongs to AMM 12592 minimal-C

**Session:** death-star-2026-07-30-coinC2.
**Status:** VERIFIED-EXACT (identity replay) + PROVENANCE (context pairing).

## The event

The external fragment containing certificate (27) was re-supplied to this
session, this time **with its surrounding context attached**: the verbatim
statement of AMM Problem 12592 parts (a) (`2n` deadline) and (b)
(`max(2, 2n-1)` deadline), and the explicit framing sentence

> "minimizing C such that for some D, the algorithm terminates in Cn+D
> biased coin flips where n is the length of the maximal constant initial
> string"

immediately adjacent to the artanh sandwich bounds, the two evaluation
points `t_A = 389/2181`, `t_B = 5872957/11821757`, and the exact rational
margin. This settles the fragment's **home**: it is proposer-side analysis
of the minimal-C frontier question (HYP-9061 (Q)), not an LRC(14)
wider-gap fingerprint, not an Abel-Dini figurate mass ordering, and not a
Baker-style irrationality-measure bound. Those rival readings (HYP-9023's
"HOME" candidates; klein-S402/S404, mac-mini-S168 reflections) are now
**demoted on provenance** — their exact structural by-products remain
valid and are inventoried below.

## Exact replay (independent, this session)

Recomputed from scratch in `fractions` arithmetic:

```text
G := (2457/6592) * L(t_B) - U(t_A) - 1/25
   = 391926968594914200867482400554891567498742649630277
     / 82738859282193417287303438726081463937219800938169600
   > 0,
```

with `L(t) = 2(t + t^3/3 + t^5/5)`, `U(t) = 2(t + t^3/3 + t^5/(5(1-t^2)))`.
This is **byte-identical** to the rational `F` in the re-supplied fragment
("the right side of (27) minus 1/25 is at least F > 0") and to the value
already verified by `04-computation/amm12592_artanh_certificate_decode_deathstar.py`
(referee (2)). Float check: `F ~ 0.0047369`, true slack of the log form
`~ 0.0057246`. Hence the fragment's (27) is exactly

```text
(2457/6592) log(8847357/2974400) - log(1285/896) > 1/25          (27)
```

with both logs = log-likelihood rapidities `log(p/q)` at biases
`p_A = 1285/2181`, `p_B = 8847357/11821757`.

## What survives from the rival decodings (still true, now re-homed)

1. **Abel-Dini telescoping form (HYP-9023, exact):** both likelihood
   ratios are consecutive-partial-sum ratios `S_n/S_{n-1}` with
   `(S_{n-1}, S_n) = (896, 1285)` and `(2974400, 8847357)`;
   `t = x_n/(S_n + S_{n-1})`, `x_n = S_n - S_{n-1}`. Open sub-question,
   now re-posed inside AMM 12592: which quantity of the deficit-flow
   analysis makes the extremal biases *partial-sum ratios*.
2. **Coefficient pinning (opus-S4, exact):** the p-adic isolated-prime +
   height-cliff argument fixing `(c, d, r) = (2457/6592, 1, 0)` is
   unaffected and is what today's context pairing confirms.
3. **2-adic engineering (HYP-9061 sec 2e):** `q`-numerators
   `896 = 2^7*7` (s_A = 7) and `2974400 = 2^6*5^2*11*13^2` (s_B = 6),
   odd denominators/`p`-numerators; the misaligned dyadic scales are the
   live mechanism candidate for the two-bias dual.
4. **Ruled-out homes stay ruled out** (figurate mass differences, both
   candidate LRC papers, hypergeometric Pade): those negative checks were
   correct and are consistent with the confirmed home.
5. **Prime/independence anatomy (kind-pasteur-S129, exact):**
   `8847357 = 3 * 2949119` with `2949119` prime; the two rapidities are
   multiplicatively independent (primes `257`, `2949119` survive alone in
   all `r_A^i r_B^j`, `|i|,|j| <= 3`) — so (27) is not a disguised power
   relation. `1/25 = (1/5)^2` matches the `t^5/5` truncation order:
   margin engineered near-extremal for the degree-5 sandwich.
6. **No-near-cancellation discriminator (klein-S404, exact):**
   `r_A/r_B ~ 0.3308 != 2457/6592 ~ 0.3727`; the linear form sits *above
   a floor*, not near zero — shaped like a rate gap (dual certificate),
   not a Baker small-linear-form. Integer form (mac-mini-S168):
   `2457 r_B - 6592 r_A > 6592/25 = 263.68`.

## Consequence for the frontier

HYP-9061's two live readings of (27) — lower-bound dual step
(`C* >= 9049/6592`?) vs construction gate — are **unchanged**; what is
new is that the question is now definitively attached to (Q), so the
reconstruction lanes (G: two-bias dual; F: policy hostile; C:
finite-M infeasibility) inherit full priority. The margin `1/25` is a
*rate margin*: any reconstruction should use it to absorb the `O(D)`
additive freedom, which is exactly what a per-`M` linear rate inequality
with a positive gap can do.
