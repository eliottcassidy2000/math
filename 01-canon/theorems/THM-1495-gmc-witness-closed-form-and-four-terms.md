---
id: THM-1495
title: "THE GMC(4) WITNESS: CLOSED FORM, AND IMPROVED FROM 6 TERMS TO 4. The owner's outside counterexample P = (1+Z2)(W1(1-Z1)+W2), Q = Z2 (W_j = conj Z_j, so 2 complex = 4 REAL Gaussians) is CONFIRMED exactly through m=12. TWO IMPROVEMENTS: (1) a CLOSED-FORM proof replacing moment-by-moment checking — E[exp(tP)] = 1 identically and E[Z2 exp(tP)] = t/(1-t), so E[Z2 P^m] = m! for ALL m at once; the intermediate resolvent E_{Z1}[exp(tP)] = exp(tcd)/(1+tc) with c = 1+Z2, d = conj(Z2) is verified symbolically. (2) A FOUR-TERM witness: P' = (1+Z2)(W2 - Z1 W1) drops two of the six terms and gives the SAME sequence E[P'^m] = 0, E[Z2 P'^m] = m! through m=12. REPO STATUS SETTLED: canon recorded Zhao's image conjecture / Mathieu subspaces as 'COROLLARY-false, no witness', with a VC-witness transport only PARTLY executed (dim ~76 -> ~20, contingent, MISTAKE-201). So we had the CONCLUSION but never the OBJECT, and this witness is degree 3 in 4 real Gaussians and UNCONDITIONAL — it needs no Jacobian input. MECHANISM: a TWO-STAGE CASCADE — Z1 manufactures the pole 1/(1+tc), the Z2 integral resums the geometric series to exactly 1 for E[e^{tP}] but the extra factor of w shifts it by one term leaving t/(1-t), whose pole at t=1 is why no moment ever vanishes. One complex Gaussian cannot run this, which is a structural reason the boundary sits between n=2 (GMC TRUE, Derksen-van den Essen-Zhao) and n=4 — making n=3 THE open case"
status: VERIFIED-EXACT (m <= 12 for both witnesses; resolvent step symbolic) + the closed form derived
author: opus-2026-07-20-S410
credits: owner (the outside counterexample and the GMC(2)/nullcone pointers)
depends_on: [THM-1300 (JC counterexample — NOT needed here, which is the point), MISTAKE-201 (the partly-executed transport)]
---

# THM-1495 — The GMC(4) witness: closed form, and 6 terms → 4

## 1. What was supplied, and that it is right

`P = (1+Z₂)(W₁(1−Z₁) + W₂)`, `Q = Z₂`, with `W_j = conj(Z_j)` — so two complex Gaussians,
i.e. **four real** ones. Standard normalisation `E[Z^a \bar Z^b] = a!·δ_{ab}`.

**Confirmed exactly**: `E[P^m] = 0` and `E[Q P^m] = m!` for `m = 1 … 12`. Cubic, 6 terms.
Since `m! ≠ 0` for every `m`, the Mathieu/Gaussian-moment property (`E[QP^m] = 0` for
`m ≫ 0`) fails outright. **GMC(4) is false.**

## 2. Where the repo already stood — and what was missing

Canon records Zhao's Image Conjecture and the Mathieu-subspace statements as
**"COROLLARY-false, no witness"**: false as consequences of the JC counterexample
(Alpoge), with no explicit object. A VC-witness transport was attempted and only
**partly executed** — dimension `≈76 → ≈20 → contingent`, logged with MISTAKE-201.

> **So we had the conclusion and never the object.** This witness supplies it, and is
> stronger than anything the transport could have produced in the three ways that matter:
> **smaller** (4 real variables, degree 3, vs `≈20` dimensions), **explicit**, and
> **unconditional** — it uses no Jacobian-conjecture input at all. It is checkable by hand
> via Wick contraction.

## 3. Improvement (1): the closed form

Moment-by-moment verification is replaced by two identities.

**Step 1 — integrate `Z₁`.** Write `c = 1+Z₂`, `d = \bar Z₂`. Then
`P = c\bar Z₁(1−Z₁) + cd`, and using `∫ e^{a\bar z} e^{−λ|z|²} d²z/π = 1/λ`:

```
E_{Z₁}[ e^{tP} ]  =  exp(t·c·d) / (1 + t·c)          <-- a RESOLVENT
```

*(verified symbolically: the two `t`-series agree to `O(t⁹)`, difference exactly 0).*

**Step 2 — integrate `Z₂`,** using `∫ w^k e^{a\bar w} e^{−λ|w|²} d²w/π = (1/λ)(a/λ)^k`
with `λ = 1−t`, `a = t`, and expanding `1/(1+t+tw) = (1/(1+t))Σ_k(−tw/(1+t))^k`:

```
E[e^{tP}]    = (1/(1−t²)) · Σ_k [ −t²/(1−t²) ]^k = (1/(1−t²))·(1−t²)  =  1
E[Z₂ e^{tP}] = (1/(1−t²)) · (t/(1−t)) · (1−t²)                        =  t/(1−t)
```

> **`E[e^{tP}] ≡ 1`, while `E[Z₂ e^{tP}] = t/(1−t) = Σ_{m≥1} t^m`.**

Reading off coefficients of `t^m/m!` gives `E[P^m] = 0` and `E[Z₂P^m] = m!` **for all `m`
simultaneously** — no induction, no moment table.

## 4. Improvement (2): four terms suffice

Searching all proper subsets of the six terms: the subset `{−Z₁W₁Z₂, −Z₁W₁, Z₂W₂, W₂}`
regroups as

```
P′ = (1 + Z₂)(W₂ − Z₁W₁)        4 terms, still cubic
```

— the two lone `W₁` terms are droppable. **Verified: `E[P′^m] = 0` and `E[Z₂P′^m] = m!` for
`m = 1…12`, the identical sequence.** So the witness costs 4 terms, not 6.

*(The 1-, 2-, and 3-term subsets are **not** counterexamples: they have `E[QP^m] ≠ 0` only at
`m = 1` and vanish thereafter, which is exactly what GMC permits.)*

## 5. The mechanism, and why the boundary is where it is

The proof is a **two-stage cascade**:

1. `Z₁` manufactures a **pole**: its integral returns the resolvent `1/(1+tc)`.
2. `Z₂` **resums** the resulting geometric series. For `E[e^{tP}]` the series telescopes to
   exactly `1`. For `E[Z₂e^{tP}]` the extra factor of `w` shifts the series by one term and
   leaves `t/(1−t)` — and **the pole at `t = 1` is precisely why no moment ever vanishes**.

**One complex Gaussian cannot run this**: there is no second variable to resum the first
one's pole. That is a structural reason the boundary sits between `n = 2` and `n = 4`, and
it says the interesting question is `n = 3` — *one complex plus one real Gaussian*, a
half-cascade.

| `n` (real Gaussians) | status |
|---|---|
| 2 | **TRUE** — theorem (Derksen–van den Essen–Zhao) |
| **3** | **OPEN — the real gap** |
| ≥ 4 | **FALSE** — this witness |

## 6. GMC(2) consistency check (and a method caution)

A brute-force sweep over 19682 polynomials in one complex Gaussian (coeffs in `{−1,0,1}`,
degree ≤ 3) flagged 152 "candidates". **They are all spurious**, and re-testing shows why:

```
P = −z³ − z²z̄ − z² − z ,  Q = z̄ :   E[QP^m], m=1..8  =  [−3, 0,0,0,0,0,0,0]
P = −z³ − z² − z        ,  Q = z̄ :   [−1, 0,0,0,0,0,0,0]
```

`E[QP^m]` is nonzero **only at `m = 1`** and vanishes for `m ≥ 2` — entirely consistent with
GMC, which asserts vanishing for `m ≫ 0`, not for all `m`. **GMC(2) intact**, as the theorem
requires; the machinery agrees with the known result.

> **Method caution.** My search criterion was *"`E[QP^m] ≠ 0` for some `m ≤ 4`"*, which is
> the wrong test — it flags 152 non-counterexamples. The correct test is *"`E[QP^m] ≠ 0` for
> arbitrarily large `m`"*. Recorded because the misreading makes a theorem look refutable,
> and this is the second time this session that a too-weak criterion produced phantom hits.

## 7. Speculative, recorded as a direction and NOT claimed

`m!` is the number of linear orders of the `m` copies of `P` — equivalently the number of
Hamiltonian paths in the complete digraph on `m` labelled vertices. The surviving Wick
contractions form a single **chain** through the copies.

**Rédei's theorem** — the repo's founding object — says every tournament has an **odd**
number of Hamiltonian paths, hence a nonzero one. That suggests: *is there a variant witness
whose moment sequence counts Hamiltonian paths in a tournament, so that non-vanishing is
certified by **parity** rather than by exact evaluation?* That would give a family of GMC
counterexamples, one per tournament, with a Rédei certificate.

**Not claimed.** `m!` is the complete-digraph count; obtaining a tournament-restricted count
would require an asymmetric propagator, and I have not constructed one. Logged as HYP-8340.

## 8. Open

1. **`n = 3`** — the genuine gap, now sharply isolated by §5.
2. **Is 4 terms minimal?** Minimality was checked only within the original 6-term support.
3. **The DvdEZ nullcone statement.** §3 shows this `P` lies in the nullcone by an exact
   *telescoping resummation* (`(1/(1−x))·(1−x) = 1`). Whether that is the general shape of
   nullcone elements — and whether the stronger 2-dimensional nullcone conjecture of
   Derksen–van den Essen–Zhao is provable by forbidding a one-variable telescope — is the
   natural next question.

## Verification

`04-computation/gmc_closed_form_opus_S410.py` (m ≤ 12; resolvent step; subset minimality;
`n=2` sweep), `04-computation/gmc_four_term_and_n2_opus_S410.py` (the 4-term witness;
re-test of the `n=2` candidates). Outputs in `05-knowledge/results/`.
