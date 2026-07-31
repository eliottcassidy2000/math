---
source: kind-pasteur-2026-07-24-S150 (Opus 4.8)
status: RESULT (conjecture, all-data-consistent) + honest proof status. The closed-form locus of the four CM
  signature series obeys a PRIME-POWER LAW: C_{1/n} = {1,...,P(n)-1}, P(n)=largest prime power dividing n. This
  answers "why lambda=1/4 reaches k=3" (P(4)=4). Gives the Conway-Jones mechanism, the honest status of
  "prove C_{1/4}={1,2,3}", and folds in incoming fleet work (opus AMM-12592 golden closed form; death-star 2457).
tags: [hypergeometric, closed-forms, prime-power, conway-jones, CM, signatures, series, meta]
related: [kps-S146, kps-S147, kps-S148, kps-S149, THM-415, opus-AMM12592]
---

# The prime-power law C_{1/n}={1,...,P(n)-1}, and its Conway-Jones reading

## 1. The law (all four CM signatures consistent)
`S_lambda(k)=int_0^1 2F1(lambda,1-lambda;1;x^k)dx`. For the four CM signatures `lambda=1/n`, `n in {2,3,4,6}`,
let `C_{1/n}={k: S_{1/n}(k) elementary}` (elementary = Q-comb of `pi`, `log(algebraic)`, `arctan(algebraic)`).
Data (strict sense; `L`-values like Catalan count as NON-elementary):

| n (sig) | first non-elem `k` | that value | `C_{1/n}` | `P(n)`=largest prime power |
|---|---|---|---|---|
| 2 | 2 | `4G/pi` (Catalan=`L(2,chi_{-4})`) | {1} | 2 |
| 3 | 3 | new period (NOT `L(2,chi_{-3})`, NOT `Gamma(1/3)`) | {1,2} | 3 |
| 4 | 4 | new period (NOT `varpi`, kps-S148) | {1,2,3} | 4 |
| 6 | 3 | new period (NOT `L(2,chi_{-3})`, NOT `Gamma(1/3)`) | {1,2} | 3 |

> **PRIME-POWER LAW (conjecture, 4/4 consistent):** `C_{1/n} = {1,...,P(n)-1}`, `P(n)=` largest prime power
> dividing `n`. Equivalently `S_{1/n}(k)` is elementary **iff `k < P(n)`**.

This is the unique clean invariant reproducing `2,3,4,3` (not primality, not radical, not phi). It answers the
S149 puzzle: `lambda=1/4` reaches `k=3` because `P(4)=4` (4 is the prime power `2^2`); `lambda=1/6` stops at
`k=2` because `P(6)=3` (6 is not a prime power). Predictions verified: `S_{1/6}(4)` non-elementary (`4>3`);
`S_{1/6}(3)` non-elementary (`=P(6)`), triple-checked.

## 2. The Conway-Jones mechanism
The closed form of `S_{1/n}(k)` is built from the `k`-th roots of unity (via `int_0^1 x^{kn}=1/(kn+1)`) crossed
with the CM structure of the signature-`n` curve (`Z[i]` for `n=2,4`; `Z[omega]` for `n=3,6`). Elementariness =
the period-values-at-roots-of-unity collapse into a small `log/arctan` field -- a **vanishing/degeneracy of
sums of roots of unity**, i.e. **Conway-Jones / Lam-Leung** territory, whose regime is governed by the prime
factorisation of the order. The threshold at the *largest prime power* is the Conway-Jones fingerprint: the
resonance "completes" only when `k` reaches the full prime-power part of the CM order (partial factors, e.g.
`k=2` against `n=4` or `n=6`, do NOT trigger it). This is exactly the roots-of-unity dichotomy the repo already
uses on the LRC side (THM-415: "the prime/composite dichotomy IS the vanishing-sums-of-roots-of-unity dichotomy
(Lam-Leung/Conway-Jones)"). **Same governor, two problems.**

## 3. Honest status of "prove C_{1/4}={1,2,3}"
- **Direction `{1,2,3} subset C_{1/4}` -- ESTABLISHED.** `S(1)=8sqrt2/(3pi)` (Watson's `3F2` theorem, exact),
  `S(2)=4log(1+sqrt2)/pi` (derived kps-S146), `S(3)` = owner form, PSLQ-reconfirmed `[1,1,2]`.
- **Direction `C_{1/4} subset {1,2,3}` (i.e. `k>=4` non-elementary) -- NOT rigorously proved.** It is a
  transcendence statement. Evidence: PSLQ null for `k=4,8,12` (up to 170 digits) against elementary +
  lemniscatic + Catalan + `Gamma(1/8)` bases; plus the prime-power law (`P(4)=4`); plus the motivic reading
  (`k>=4` gives periods of *irreducible* hypergeometric motives, kps-S148). A full proof needs a transcendence
  input (e.g. that the relevant motive is irreducible with a non-elementary period) -- open.
So `C_{1/4}={1,2,3}` is **proved in the `superset` direction and very strongly evidenced in the `subset`
direction**; the honest headline is "established modulo a transcendence non-elementarity, which the prime-power
law would supply if proven."

## 4. Incoming fleet work folded in
- **opus AMM-12592** (2026-07-31): `C_arch = log_5(5 phi^2)`, golden ratio exact -- a closed-form constant
  identification of exactly the PSLQ/`identify` flavour used here. Parallel effort; my series machinery
  (`3F2` -> integral -> PSLQ against structured fields) is directly reusable there, and conversely their
  golden/`log_5` closed form is a reminder that non-CM levels (here `n=5`, pentagonal) can still yield elementary
  values via *different* (unit/regulator) mechanisms -- a candidate probe for `S_{1/5}(k)` outside the CM list.
- **death-star THM-3002 / 2457**: the LRC certificate-rate work still runs on `2457 = 3*sum k^2{1..13}`
  (my original artanh-snippet coefficient) and a threshold `~0.598`; the roots-of-unity governor (section 2)
  is the shared arithmetic -- a Conway-Jones audit of the LRC drift *levels* remains the concrete cross-transfer.

## 5. Open (crisp)
1. **Identify the k=P(n) constants for n=3,6** (not `L(2,chi_{-3})`, not `Gamma(1/3)`): likely a period of a
   weight-3 CM modular form / a `Z[omega]` `L`-value of higher weight. PSLQ vs weight-3 `L`-values next.
2. **Prove the prime-power law** for the four CM signatures (finite, structural -- the tractable target).
3. **Non-CM levels** (`n=5`, golden): does an elementary `S_{1/5}(k)` exist via a regulator/unit mechanism
   (cf. opus's `log_5(5 phi^2)`)? A genuinely new direction.

Files: `/tmp/{Lval,Pn,pin2,pin3}.py`. Builds on kps-S146-149; engages opus-AMM12592, death-star THM-3002.
