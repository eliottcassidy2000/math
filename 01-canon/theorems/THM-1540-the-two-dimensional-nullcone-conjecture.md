---
id: THM-1540
title: "THE 2-DIMENSIONAL NULLCONE CONJECTURE — stated, two of its three cases PROVED, and GMC(2) shown to be a two-line corollary. CONJECTURE: N(E) = {P : all charges ≥ 1} ∪ {P : all charges ≤ −1}, i.e. the Gaussian nullcone at n = 2 is exactly the STRICTLY ONE-SIDED charge support. (A) THE EXACT POLAR REDUCTION: E[Pᵐ] = L(v ↦ CT_u(H_√v(u)ᵐ)) with H_r(u) := P(ru, r/u) — so n = 2 is literally DUISTERMAAT–VAN DER KALLEN's constant-term problem, Gaussian-averaged over the radius. Verified exactly. (B) THE TWO-CHARGE THEOREM, PROVED: if P has exactly two charges C > 0 > −B then charge-0 forces k_C = tB', k_{−B} = tC' uniquely, so E[Pᵐ] = 0 unless (B'+C') | m and otherwise equals C(m; tB', tC')·L(Fᵗ) with F = v^{CB'}q^{B'}s^{C'} ≠ 0 — nonzero for large t by the saddle lemma. So two-charge P is NEVER in the nullcone. (C) THE GENERAL REDUCTION: for ≥3 charges the leading term of E[Pᵐ] is governed by the TOP EDGE of the (charge, degree) Newton polygon where it crosses charge 0, and that edge polynomial is a ONE-VARIABLE LAURENT polynomial — so the remaining case follows from the 1-variable Laurent nullcone lemma plus domination. (D) EVIDENCE, NOW CONTROLLED: 886 800 two-sided P tested, zero in the nullcone — and unlike the S135 sweep this one HAS A WORKING POSITIVE CONTROL, since the conjecture predicts one-sided P are in the nullcone and all five controls pass. (E) GMC(2) FOLLOWS IN TWO LINES from the conjecture."
status: >
  (A) PROVED (direct substitution) and VERIFIED-EXACT on random P for m = 1..4.
  (B) PROVED. The uniqueness of the charge-0 splitting is elementary; the nonvanishing is
  THM-1520's saddle lemma, so (B) inherits that lemma's one gap (Laplace estimate without
  explicit error bounds — HYP-8350). Verified on 6 explicit P: E[Pᵐ] is nonzero exactly at
  the predicted multiples of B'+C'.
  (C) A REDUCTION, not a proof. The Newton-polygon argument is sketched and the required
  1-variable input is isolated and tested (10 320 two-sided Laurent h, none in the 1-variable
  nullcone), but domination of the leading term is NOT established.
  (D) BOUNDED but CONTROLLED — materially stronger than the uncontrolled S135 sweep.
  (E) PROVED from the conjecture.
  THE CONJECTURE ITSELF IS OPEN; only the ≥3-charge case of the hard inclusion remains.
source: mac-mini-2026-07-20-S136 (owner: "aim to finish GMC(2) by finishing the stronger
  2 dimensional nullcone conjecture")
depends_on:
  - THM-1520  # charge telescoping + the saddle lemma; the one-sided case
  - THM-1500  # the GMC master theorem; GMC(N) false for all N >= 3
related:
  - 07-reflections/the-telescoping-principle-macmini-S135.md
script: 04-computation/gmc2_nullcone_macmini_S136.py (+ .out)
---

# THM-1540 — the 2-dimensional nullcone conjecture

**One line.** Asking for the *nullcone* rather than for GMC(2) is the better question: it is
strictly stronger, GMC(2) drops out of it in two lines, and — decisively — it makes the
computational search **controllable**, which the direct GMC(2) search could never be.

## The conjecture

For `E` on `ℂ[Z,W]` (`W = Z̄`, `E[Z^aW^b] = a!δ_{ab}`), let
`N(E) := {P : E[Pᵐ] = 0 for all m ≥ 1}` and grade by **charge** `c = deg_Z − deg_W`.

> **NULLCONE CONJECTURE (n=2).**
> `N(E) = {P : all charges ≥ 1} ∪ {P : all charges ≤ −1}` — strictly one-sided support.

**Why it is the better target.** A direct GMC(2) search has *no possible positive control*: a
positive control would BE the counterexample being sought (the flaw I flagged in S135). The
nullcone conjecture makes a **positive prediction** — one-sided `P` *are* in the nullcone — so
the machinery can be validated before its negatives are trusted. That is the whole
methodological gain, and all five controls pass.

## (A) The exact polar reduction

With `H_r(u) := P(ru, r/u) = Σ p_{ab} r^{a+b} u^{a−b}`, a Laurent polynomial in `u` whose
**exponents are exactly the charges** and whose `r`-degree is the total degree:

> **`E[Pᵐ] = L( v ↦ CT_u(H_{√v}(u)ᵐ) )`**, `L(vᵏ) = k!`.

Verified exactly. So `CT_u` is the charge-0 projection and `L` is the Gaussian average over
the radius (`r² ~ Exp(1)`): **GMC(2) is Duistermaat–van der Kallen's constant-term problem,
Gaussian-averaged.** That names the right context — DvdK is the model theorem of this area.

## (B) The two-charge theorem — proved

Let `P = P_C + P_{−B}` with exactly two charges `C > 0 > −B`. Write `P_C = Z^C q(V)`,
`P_{−B} = W^B s(V)` with `V = ZW`; let `g = gcd(B,C)`, `B' = B/g`, `C' = C/g`.

Charge-0 requires `k_C·C = k_{−B}·B` with `k_C + k_{−B} = m`, forcing
`k_C = tB'`, `k_{−B} = tC'`, `m = t(B'+C')` — **unique when it exists, impossible otherwise**:

> `E[Pᵐ] = 0` unless `(B'+C') | m`; and at `m = t(B'+C')`,
> `E[Pᵐ] = C(m; tB', tC')·L(Fᵗ)` with `F = v^{CB'}q^{B'}s^{C'} ≠ 0`.

By THM-1520's saddle lemma, `L(Fᵗ) ≠ 0` for all large `t`. **So a two-charge `P` is never in
the nullcone.** Verified on six explicit `P` — nonzero exactly at the predicted multiples:

| `P` | `B,C` | `B'+C'` | `E[Pᵐ]`, `m = 1..8` |
|---|---|---|---|
| `W+Z` | 1,1 | 2 | `0,2,0,12,0,120,0,1680` |
| `W+Z²` | 1,2 | 3 | `0,0,6,0,0,360,0,0` |
| `W²+Z³` | 2,3 | 5 | `0,0,0,0,7200,0,0,0` |
| `−2W³+Z²W` | 3,1 | 4 | `0,0,0,−5760,0,0,0,53648179200` |

## (C) The general reduction — what the ≥3-charge case needs

`E[Pᵐ] = Σ_d [charge-0, degree-`d` coefficient of `Pᵐ`]·(d/2)!`, and the factorials weight the
**largest** degree overwhelmingly. Plot the support in `(charge, degree)` coordinates; the
maximum degree at charge 0 in `Pᵐ` is `m·D*` where `D*` is where the **upper hull** crosses
`c = 0`. So the leading coefficient comes from the **top edge** of the Newton polygon at
charge 0 — and the monomials on a single edge form a **one-variable Laurent polynomial**.

> **The ≥3-charge case reduces to: (i) the 1-variable Laurent nullcone lemma — `CT(hᵐ) = 0`
> for all `m` implies `h` is one-sided — plus (ii) domination of the leading term.**

If the top edge has only two lattice points, (i) is exactly (B). Input (i) tested: **10 320
two-sided Laurent `h` (exponents in `[−3,3]`, coefficients in `{±1,±2}`), none with
`CT(hᵐ) = 0` for `m = 1..10`.**

## (D) Evidence, now controlled

| support size | two-sided `P` tested | in nullcone |
|---|---|---|
| 2 | 400 | 0 |
| 3 | 10 000 | 0 |
| 4 | 118 800 | 0 |
| 5 | 886 800 | **0** |

with the positive control passing on all five one-sided test cases.

## (E) GMC(2) is a two-line corollary

Assume the conjecture and let `E[Pᵐ] = 0` for all `m`. Then `P` is strictly one-sided, say all
charges `≥ 1`; so every charge of `Pᵐ` is `≥ m`; so `QPᵐ` can have charge 0 only if
`charge(Q) ≤ −m`, i.e. only if `m ≤ deg_W(Q)`. Hence **`E[QPᵐ] = 0` for every
`m > deg_W(Q)`.** ∎

## Honest scope

- **GMC(2) is NOT finished.** What is finished: the one-sided case (THM-1520), the
  two-charge case (B), the reduction of the rest (C), and the corollary (E). The open part is
  precisely: *two-sided `P` with three or more distinct charges cannot lie in the nullcone.*
- **(B) inherits THM-1520's gap.** The saddle lemma is a Laplace estimate without explicit
  error bounds (HYP-8350). Until that is closed, (B) is conditional — though its combinatorial
  half (uniqueness of the charge-0 splitting) is unconditional and elementary.
- **(C) is a reduction sketch, not a proof.** Domination of the leading term is exactly what
  is not established; the same `O(1)`-ratio difficulty as in the 1-variable case appears here
  and needs the same saddle treatment.
- **(D) is bounded.** Coefficients in `{±1,±2}`, `deg_Z, deg_W ≤ 3`, support `≤ 5`. The control
  makes it meaningful, not conclusive.
- The 1-variable Laurent nullcone statement in (C) is presumably classical (it is the
  Newton-polytope criterion, and the `n`-variable Mathieu-subspace statement is
  Duistermaat–van der Kallen 1998). **No priority is claimed for it**; it is isolated here as
  the required input and tested, not proved.

*Artifacts:* `04-computation/gmc2_nullcone_macmini_S136.py` (+out).
