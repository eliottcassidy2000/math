---
id: THM-1720
title: "THE FULL GMC(2) NULLCONE IS ONE NULLSTELLENSATZ EMPTINESS TEST — unifying opus THM-1685 (angular) and klein THM-1700 (bottom-up radial), and CLOSING the charge-0 cancellation residual. With W = Z̄ and E[Z^A W^B] = A!·δ_{AB}, the coefficients c_i of P are free complex symbols and the MOMENT IDEAL I = ⟨E[P^m] : m ≥ 1⟩ ⊆ ℂ[c_1..c_k] is honest. The conjecture 'GMC(2) nullcone = charge-one-sided' is EXACTLY V(I) ∩ {genuinely two-sided} = ∅, a Nullstellensatz emptiness test: for every positive-charge coeff c_p and negative-charge coeff c_n, 1 ∈ sat(I, c_p·c_n) (Rabinowitsch). This is the SAME object as THM-1685's angular TNC test V(I_CT) ∩ (ℂ*)^{k−2} = ∅ — one decision procedure for both the angular (DvdK) and radial (cross-shell) layers of GMC(2). It CLOSES the exact residual klein THM-1700 left open — cancellation among charge-0 pairs: 6/6 two-sided patterns close by exact Gröbner over ℚ, INCLUDING the flagged cancellation cases {+1,−1} with two and with THREE charge-0 terms (1, ZW, Z²W²), {±2}+two charge-0, and {±2,±1}+two charge-0. Each is a rigorous per-pattern PROOF (a finite Nullstellensatz certificate, ≤ 8 moments), since the nullcone ⊆ V(⟨E[P^{1..M}]⟩) for every M. klein THM-1700's bottom-up descent IS the triangular Gröbner structure of I (the lowest straddle appears in the lowest-m generator). The one remaining step is UNIFORM: a moment-count bound M(pattern) — the exact analogue of TNC's HYP-8505 CT-level bound — would upgrade the per-pattern procedure to full GMC(2)."
status: >
  PROVED per pattern (each 'closed' is a finite Nullstellensatz certificate: 1 ∈ the ideal
  saturated by c_p·c_n from ≤ 8 moments, so V(I) ∩ {c_p ≠ 0, c_n ≠ 0} = ∅ for every pos-neg
  pair, hence V(I) ∩ two-sided = ∅). VERIFIED-EXACT for 6 charge patterns over ℚ with the
  RIGOROUS all-pairs saturation (not just top×bottom). The UNIFICATION with THM-1685 is a
  structural identity of the two decision procedures. NOT a uniform proof of GMC(2) — the
  moment-count bound is open (see HYP-8535), exactly as for TNC (HYP-8505).
source: mac-mini-2026-07-20-S148 (owner: follow the descent direction, mine past threads for
  connections, take the cutting edge as far as possible)
depends_on:
  - THM-1695  # complex radial closed (Cauchy transform) -- the charge-0 layer this rests on
  - THM-1700  # klein: cross-shell runs bottom-up; the residual this closes
  - THM-1685  # opus: TNC as a Nullstellensatz emptiness test -- the object this unifies with
related:
  - THM-1510  # klein EMP: the two-weight radial case
  - THM-1540  # the 2-D nullcone conjecture
  - HYP-8470  # cross-shell coupling
  - HYP-8505  # opus: uniform CT-level bound for TNC -- the parallel open step
---

# THM-1720 — the GMC(2) nullcone is one Nullstellensatz emptiness test

## The unification

Three recent threads converge:

- **mac-mini THM-1695:** the complex radial (charge-0) layer is closed (Cauchy transform).
- **klein THM-1700:** the cross-shell descent runs **bottom-up**; its residual is *"the general
  HYP-8470 (several straddling shells, whose charge-0 pairings could cancel at low order) is
  not closed."*
- **opus THM-1685:** TNC for a `k`-nomial charge pattern is a **Nullstellensatz emptiness
  test** `V(I_CT) ∩ (ℂ*)^{k−2} = ∅`, one Gröbner per pattern.

They are the same object. With `W = Z̄` and `E[Z^A W^B] = A!·δ_{AB}`, treat the coefficients
`c_1,…,c_k` of `P` as free complex symbols. The **moment ideal**

> `I = ⟨E[P^m] : m ≥ 1⟩ ⊆ ℂ[c_1,…,c_k]`

is honest, and the GMC(2) conjecture *"the nullcone is exactly the charge-one-sided
polynomials"* is precisely

> **`V(I) ∩ {genuinely two-sided} = ∅`**, i.e. for every positive-charge coeff `c_p` and
> negative-charge coeff `c_n`, **`1 ∈ sat(I, c_p·c_n)`** (Rabinowitsch).

This is *structurally identical* to THM-1685's angular test — the angular (DvdK/TNC) and radial
(cross-shell) layers of GMC(2) are **one decision procedure**, run on the moment ideal instead
of the constant-term ideal.

## Closing the charge-0 cancellation residual

THM-1700's open case was cancellation among charge-0 pairs. Testing the full moment ideal by
**exact Gröbner over ℚ**, with the rigorous **all-pairs** saturation (every `c_p·c_n`, not just
the top×bottom pair):

| pattern | charges | charge-0 terms | result |
|---|---|---|---|
| `{+1,−1}` + one charge-0 | `−1,0,1` | 1 | **closed** |
| `{+1,−1}` + **two** charge-0 (`1, ZW`) | `−1,0,1` | 2 | **closed** |
| klein witness `aZ³+bZ̄+cZ` | `−1,1,3` | 0 | **closed** |
| `{±2}` + two charge-0 | `−2,0,2` | 2 | **closed** |
| `{±2,±1}` + two charge-0 | `−2..2` | 2 | **closed** |
| `{+1,−1}` + **three** charge-0 (`1, ZW, Z²W²`) | `−1,0,1` | 3 | **closed** |

> **6/6, including the flagged cancellation cases with up to THREE charge-0 terms.** Charge-0
> cancellation does **not** create a two-sided nullcone member.

**Each is a rigorous per-pattern proof.** If `1 ∈ ⟨E[P^{1..M}]⟩ + ⟨1 − w·c_p c_n⟩` for finite
`M ≤ 8`, then `V(⟨E[P^{1..M}]⟩) ∩ {c_p≠0, c_n≠0} = ∅`; and the nullcone `⊆ V(⟨E[P^{1..M}]⟩)`
for every `M`, so the nullcone contains no two-sided `P` with `c_p, c_n ≠ 0`. Taking the union
over pos-neg pairs gives `nullcone ∩ two-sided = ∅` for that pattern. A finite certificate.

## The bottom-up descent is the triangular Gröbner structure

klein THM-1700's mechanism — `E[P²]` kills the lowest straddle first, higher moments force the
rest — is exactly the **triangular structure of a Gröbner basis** of `I`: the lowest straddle
product `c_p c_n` appears in the lowest-`m` generator, and elimination proceeds upward. The
"descent" and the "Gröbner reduction" are one process seen twice.

## What remains: the uniform moment bound

Every "closed" above is per-pattern. The one missing step is **uniform**: a bound `M(pattern)`
on how many moments suffice, as a function of the charge pattern (number of terms `k`, charge
span). This is the **exact analogue of TNC's HYP-8505** (opus's uniform CT-level bound). A
uniform bound on either side would likely transfer, and it would upgrade the per-pattern
procedure to **full GMC(2)** — the last real gap, given the complex radial layer (THM-1695) and
the span-2 base (THM-1600) are closed.

## Honest scope

- **Per-pattern rigor, not uniform.** Each of the 6 is a genuine proof for its charge pattern;
  GMC(2) as a whole needs the uniform bound (HYP-8535 / HYP-8505). Nothing here proves GMC(2).
- **Six patterns, chosen to hit the residual.** They are not a census; they were picked to
  stress charge-0 cancellation (the flagged worry), which they close.
- **The unification with THM-1685 is structural**, an identity of the two decision procedures,
  not a new algorithm — both are "saturate the vanishing-ideal, test for `1`."
- The moment ideal uses the free-complex-coefficient model (GMC(2) allows complex `P`); the
  Hermitian sub-case is a real slice and is a fortiori covered.
- **THM-1695 collision flagged:** two canon files share id THM-1695 — mine (complex-radial,
  first-pushed 16:06:38) and klein's emp/cross-shell (16:07:50). I hold precedence by 72 s;
  klein's should renumber. Recorded, not unilaterally edited.

*Artifacts:* `04-computation/gmc2_moment_ideal_nullstellensatz_macmini_S148.py` (+out).
*Credits:* opus THM-1685 (the Nullstellensatz framing), klein THM-1700 (bottom-up descent +
the residual), mac-mini THM-1695 (complex radial). This is their synthesis.
