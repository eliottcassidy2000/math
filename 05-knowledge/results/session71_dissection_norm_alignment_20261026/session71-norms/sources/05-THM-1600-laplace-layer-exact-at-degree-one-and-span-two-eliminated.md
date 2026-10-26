---
id: THM-1600
title: "THE LAPLACE-GMC(1) LAYER IS AN IDENTITY AT DEGREE 1, AND THE ONE-SIDED CONJECTURE IS PROVED FOR CHARGE SPAN 2. (A) THM-1520's saddle lemma was an asymptotic sketch; at degree 1 it is EXACT: L((av+b)ᵐ) = m!·aᵐ·e_m(b/a) where e_m is the TRUNCATED EXPONENTIAL, because C(m,k)·k! = m!/(m−k)! turns the binomial expansion into a partial exponential sum. Hence L(pᵐ) ≠ 0 for all m past an EXPLICIT threshold, from |e_m(x) − eˣ| ≤ (|x|^{m+1}/(m+1)!)·e^{|x|} — no saddle point, no error analysis. For p = v−1 this gives L(pᵐ) = !m, THE DERANGEMENT NUMBERS 0,1,2,9,44,265, so nonvanishing is the classical fact that derangements exist and the S135 saddle limit 1/e IS the derangement asymptotic. (B) CHARGE SPAN 2 ELIMINATED EXACTLY. For charges {−1,0,+1} — the smallest case beyond THM-1540's two-charge theorem, and precisely the S137 trichotomy — write P = W·s(V) + q(V) + Z·r(V); then E[Pᵐ] = Σ_k [m!/(k!²(m−2k)!)]·L(v^k(rs)^k q^{m−2k}). With constant coefficients this eliminates in three lines: m=1 forces q=0, after which only m=2k survives with E[P^{2k}] = (2k)!(rs)^k/k!, forcing rs=0. So P is ONE-SIDED. Confirmed on 729 triples; extended computationally to degree-1 coefficients (15 624 triples, zero exceptions). (C) The requested form ∫₀^∞ CT_w(Pᵐ)e^{−s}ds ≠ 0 for some m is verified on every two-sided case, with one-sided cases identically zero. (D) THE TEN ERDŐS PROBLEMS ARE A CLEAN NEGATIVE — none involves polynomials, power sums, moments, constant terms, Laurent polynomials, truncated exponentials, or Mathieu-Zhao/DvdK structure. #475's 'partial sums' is a false friend. But 8 of the 10 carry the DECIDABLE badge — 'true for all sufficiently large n, small cases open' — which is EXACTLY the shape of the m≫0 quantifier that defines a Mathieu-Zhao subspace, and exactly the shape of the saddle lemma. The resonance is of PROOF SHAPE, not content."
status: >
  (A) PROVED — an exact algebraic identity (verified for several (a,b), m = 1..8) plus a
  standard tail bound giving an explicit threshold. This CLOSES the degree-1 case of
  HYP-8350 outright, replacing the Laplace estimate at that degree.
  (B) PROVED for constant coefficients (three-line elimination, brute-force confirmed on
  729 triples with m = 1..10). VERIFIED but NOT proved for degree-1 coefficients
  (15 624 triples, coefficients in [−2,2], m = 1..8, zero exceptions). A Gröbner
  elimination for degree ≥ 1 was attempted and did not complete.
  (C) VERIFIED on the cases tested; it is a restatement of (B) plus THM-1540(A).
  (D) SOURCED — all ten pages retrieved and read; the negative is reported as a negative.
source: mac-mini-2026-07-20-S140 (owner: "prove the one sided conjecture for bounded charge
  span by exact elimination and settle the laplace-GMC(1) layer ... show
  ∫₀^∞ CT_w(Pᵐ)e^{−s}ds ≠ 0 for some m", plus ten Erdős problem numbers)
depends_on:
  - THM-1540  # the nullcone conjecture; the two-charge theorem; the polar reduction
  - THM-1520  # the telescoping lemma and the saddle lemma this sharpens
script: 04-computation/laplace_layer_and_span_elimination_macmini_S140.py,
        05-knowledge/results/span2_elimination_quick_macmini_S140.out
---

# THM-1600 — the Laplace layer at degree 1, and charge span 2

## (A) The Laplace-GMC(1) layer is an identity at degree 1

THM-1520's saddle lemma was an asymptotic derivation without error bounds (HYP-8350). At
degree 1 it is not asymptotic at all. Since `C(m,k)·k! = m!/(m−k)!`, the binomial expansion of
`L((av+b)ᵐ)` reindexes into a **partial exponential sum**:

> **`L((av+b)ᵐ) = m!·aᵐ·e_m(b/a)`**, where `e_m(x) = Σ_{j≤m} xʲ/j!` is the truncated
> exponential.

Verified exactly for several `(a,b)`, `m = 1..8`. Three consequences, all rigorous:

- `e_m(x) → eˣ`, so S135's saddle limit `exp(a_{D−1}/(D·a_D))` is **exact** here, equal to
  `e^{b/a}`.
- `|e_m(x) − eˣ| ≤ (|x|^{m+1}/(m+1)!)·e^{|x|}`, so `e_m(x) ≠ 0` as soon as that tail drops
  below `e^{Re x}` — an **explicit threshold `m₀(x)`**, no saddle point required.
- **`p = v − 1` gives `L(pᵐ) = m!·e_m(−1) = !m`, the derangement numbers `0, 1, 2, 9, 44, 265`.**
  Nonvanishing for `m ≠ 1` is the classical fact that derangements exist, and `!m/m! → 1/e`
  **is** the saddle limit — the constant `1/e` from S135 was the derangement asymptotic all
  along.

> **The degree-1 case of HYP-8350 is closed.** Degrees `≥ 2` remain (see scope).

## (B) Charge span 2, eliminated exactly

The smallest case not covered by THM-1540's two-charge theorem is charge support
`{−1, 0, +1}` — which is precisely the trichotomy of the S137 reflection. Write

`P = W·s(V) + q(V) + Z·r(V)`, `V = ZW`.

Charge-0 needs equally many `+1` and `−1` factors, so

> `E[Pᵐ] = Σ_k [m!/(k!²(m−2k)!)]·L(v^k·(rs)^k·q^{m−2k})`.

**With constant `r, q, s` this eliminates in three lines:**

1. `m = 1` gives `E[P] = q`, so **`q = 0`**.
2. With `q = 0` only `m = 2k` survives: `E[P^{2k}] = (2k)!·(rs)^k/k!`.
3. Setting that to zero forces **`rs = 0`**, i.e. `r = 0` or `s = 0`.

So `P = Z·r` (charges `{+1}`) or `P = W·s` (charges `{−1}`) — **one-sided**. ∎

Brute-force confirmation: `r,q,s ∈ [−4,4]` (729 triples), `m = 1..10` — **zero** two-sided
nullcone elements. Extended to **degree-1 coefficients**: 15 624 triples with coefficients in
`[−2,2]`, `m = 1..8` — **zero** exceptions.

## (C) The requested integral form

THM-1540(A) gives `E[Pᵐ] = ∫₀^∞ CT_w(H_{√s}(w)ᵐ)·e^{−s}ds`. Tested:

| `(r,q,s)` | charges | `E[Pᵐ]`, `m=1..6` | nonzero? |
|---|---|---|---|
| `(1,0,1)` | `{−1,+1}` | `0, 2, 0, 12, 0, 120` | yes |
| `(3,0,−2)` | `{−1,+1}` | `0, −12, 0, 432, 0, −25920` | yes |
| `(1,1,1)` | `{−1,0,+1}` | `1, 3, 7, 25, 81, 331` | yes |
| `(2,−1,3)` | `{−1,0,+1}` | `−1, 13, −37, 505, −2281, 32581` | yes |
| `(1,0,0)` | `{+1}` | all `0` | **no** |
| `(0,0,1)` | `{−1}` | all `0` | **no** |

Every two-sided case is nonzero at some `m`; the one-sided cases are identically zero — they
**are** in the nullcone, as the conjecture says. Note the `{−1,+1}` rows are nonzero exactly at
even `m`, matching THM-1540(B)'s prediction (`B'+C' = 2`).

## (D) The ten Erdős problems: a clean negative, with one real resonance

All ten pages (#1016, 506, 742, 19, 580, 547, 460, 556, 475, 848) were retrieved and read.

**None of them involves** polynomials, power sums, moments, constant terms, Laurent
polynomials, truncated exponentials, zeros of partial sums, sum–product structure, or
Mathieu-Zhao / Duistermaat–van der Kallen material. The DvdK circle needs a polynomial
algebra, a linear functional, and a vanishing-on-all-powers hypothesis; **not one of those
three appears anywhere on the list.** Six are extremal graph theory, one combinatorial
geometry, three elementary number theory.

**#475 is a false friend.** It is literally phrased with "partial sums," but they are prefix
sums of a permuted finite subset of `𝔽_p` — a sequenceability / Hall–Paige problem, unrelated
to truncating a power series.

**But the list is not random, and the axis that organizes it is worth having.** Eight of ten
carry Bloom's **DECIDABLE** badge: *true for all sufficiently large `n`, small cases open.*
That is an unusual concentration — and it is **exactly the shape of the `m ≫ 0` quantifier
that defines a Mathieu-Zhao subspace**, and exactly the shape of (A) above ("nonzero past an
explicit threshold, small `m` by hand"). The resonance is of **proof shape, not content**, and
that is the honest way to state it.

## Honest scope

- **(A) closes degree 1 only.** Degrees `≥ 2` have no comparable identity; symbolic
  elimination clears `D ≤ 3` (S135) but a uniform-in-`D` argument is still missing, so
  HYP-8350 is *reduced*, not closed.
- **(B) proves the constant-coefficient case only.** The degree-1 coefficient case is
  computational evidence over a bounded box, not a proof; a Gröbner elimination was attempted
  and **did not complete**. Charge spans `≥ 3` are untouched.
- (C) is a restatement of (B) plus THM-1540(A) on the tested cases, not an independent result.
- **(D) is a negative and should be cited as one.** The "DECIDABLE" resonance is an
  observation about how the list was probably filtered, not a mathematical connection; no
  Erdős problem on it bears on GMC.
- The overall nullcone conjecture (THM-1540) remains **open**: proved now for one-sided
  support, exactly two charges, and span-2-with-constant-coefficients.

*Artifacts:* `04-computation/laplace_layer_and_span_elimination_macmini_S140.py`,
`05-knowledge/results/span2_elimination_quick_macmini_S140.out`.
