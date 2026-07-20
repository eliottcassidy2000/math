---
id: THM-1640
title: "THE GAMMA-BRIDGE MECHANISM IS POSITIVITY, NOT DOMINATION — and the algebraic core is now rigorous. (1) THE TOP-r-COEFFICIENT IDENTITY, proved and exact-tested: for P(r,u) = Σ_q g_q(r)u^q with D = max_q deg_r g_q and leading symbol Λ(u) = Σ_{deg g_q=D} lc(g_q)u^q, [r^{mD}] CT_u(P^m) = CT_u(Λ^m) exactly — the rigorous form of 'the top coefficient is the toral quantity', pure algebra, no analysis. (2) DOMINATION IS FALSE: the top term's fraction of E_r[L_m] collapses (0.67, 0.36, 0.18, 0.11, 0.068 at m = 2,6,12,18,24), independently confirming death-star's MISTAKE-202 — and this RETRACTS the mechanism of my own klein-S351 Gamma bridge, which invoked exactly 'the top term dominates the r-average'. The bridge's conclusion survives; its stated reason does not. (3) POSITIVITY is the real mechanism on the radial {−1,+1} span: α = r|c|² ≥ 0 makes every E_r[α^j] > 0 termwise, for ALL c and no degree bound (matching opus THM-1535/1540 Hankel-PD) — immune to the collapse in (2) because it never compares terms. (4) THE PRECISE GAP: on the both-signs {−1,0,1} span the integrand L_m(r) is SIGN-INDEFINITE, so positivity is unavailable and domination is dead; GMC(2)'s remaining content is to show E_r[L_m] ≠ 0 for some m using neither."
status: >
  (1) PROVED (elementary — the r-degree of Σ_{Σq_i=0} Π g_{q_i} is mD exactly on the
      all-top-charge terms) + EXACT-TESTED with a control that first caught MY OWN wrong test
      predicate (it demanded deg = mD; the correct test is [r^{mD}] directly, since both
      sides vanish together when the top charges cannot sum to 0). Fixed; identity holds.
  (2) VERIFIED-EXACT and it is a RETRACTION. My klein-S351 Gamma bridge's mechanism is void;
      HYP-8395 §3 is superseded. TNC (constant-coeff) is unaffected — it is DvdK 1998
      (THM-1630) and needs no r-average at all.
  (3) PROVED for the radial {−1,+1} span, any degree (the conjugacy a = c-bar makes
      α = r|c|² ≥ 0). This is not new — it is opus/boxeph's positivity route — but it is the
      correct mechanism and is recorded as such against the dead domination one.
  (4) SCOPE STATEMENT + exhibited sign-indefiniteness. The both-signs case is OPEN; this
      file does not close it.
  GMC(2) REMAINS OPEN.
source: klein-2026-07-20-S363 (owner: work GMC(2) through TNC)
depends_on:
  - THM-1630  # TNC = Duistermaat–van der Kallen 1998 (the constant-coeff layer, cited)
  - THM-1510  # NC2 ⇒ GMC(2) (now Lean-formalized by death-star)
related:
  - THM-1535  # opus: Hankel positive-definiteness (the positivity mechanism, sign-coherent)
  - THM-1540  # opus: charge-0 vanishes over C
  - "death-star MISTAKE-202 / HYP-8426: the domination retraction this confirms"
  - "klein-S351 HYP-8395: the Gamma bridge whose MECHANISM this supersedes"
script: 04-computation/bridge_mechanism_klein_S363.py (+ .out)
---

# THM-1640 — the bridge mechanism is positivity, not domination

## Where this sits

The GMC(2) chain is `TNC ⟹ NC2 ⟹ GMC(2)`:

- **TNC** (constant-coefficient toral nullcone) = **Duistermaat–van der Kallen 1998**, Thm 2 +
  Rmk 3 (THM-1630). A citation. No `r`-average involved.
- **NC2 ⟹ GMC(2)**: Lean-formalized (death-star `mathieuZhao_of_charge_pos`, kernel-pure).
- **TNC ⟹ NC2** — the `r`-dependent lift, the "Gamma bridge" (klein-S351) — is the live gap,
  and its domination mechanism was refuted (death-star MISTAKE-202). This file confirms the
  refutation, makes the bridge's *algebraic* core rigorous, and identifies the *correct*
  analytic mechanism.

## 1. The top-`r`-coefficient identity (rigorous core)

Write `P(r,u) = Σ_q g_q(r) u^q`, `D = max_q deg_r g_q`, and the **leading symbol**
`Λ(u) = Σ_{q : deg g_q = D} lc(g_q) u^q`.

> **Identity.** `[r^{mD}] CT_u(P(r,u)^m) = CT_u(Λ(u)^m)`.

*Proof.* `CT_u(P^m) = Σ_{q_1+…+q_m=0} Π_i g_{q_i}(r)`. A product `Π_i g_{q_i}` has `r`-degree
`Σ_i deg g_{q_i} ≤ mD`, with equality iff every `q_i` is a top-degree charge; there its
`r^{mD}` coefficient is `Π_i lc(g_{q_i})`. Summing over the all-top-charge tuples with
`Σq_i = 0` gives `CT_u(Λ^m)`. ∎

This is the precise, analysis-free form of S351's "the top coefficient is the toral quantity".
**Note both sides vanish together** when the top charges cannot sum to `0` (e.g. `Λ = u^{−1}+u`,
`m` odd) — a fact my first *test* got wrong by demanding `deg = mD`; the correct test extracts
`[r^{mD}]` directly. (Control caught it; identity holds on every case.)

## 2. Domination is false — retracting my own S351 mechanism

For `P = Z̄ + 1 + Z` (`α = r`, `β = 1`, `E[P^m]` known nonzero), the top term's share of
`E_r[L_m] = Σ_k ([r^k]L_m) k!`:

| `m` | 2 | 6 | 12 | 18 | 24 |
|---|---|---|---|---|---|
| `|top|/Σ|terms|` | 0.67 | 0.36 | 0.18 | 0.11 | **0.068** |

The top term's fraction **collapses** — independently reproducing death-star's `0.068` at
`m = 24`. So *"the top term dominates the `r`-average"* is false, and that was exactly the
mechanism klein-S351 (HYP-8395 §3) invoked. **I retract the mechanism.** The bridge's
*conclusion* (`E[P^m] ≠ 0`) still holds here — but not for that reason. `k! = Γ(k+1)` growing
is real, yet the lower coefficients grow too, and no single term wins.

## 3. Positivity is the real mechanism (radial span)

On the `{−1,+1}` span (`β = 0`), the generating function collapses to
`E[P^{2j}] = C(2j,j)·E_r[α^j]`, `α = r·a·c`. For `P ∈ ℂ[Z,Z̄]` the charge `±1` coefficients are
conjugate, so `a·c = |c|²` and

> `α = r·|c|² ≥ 0` pointwise ⟹ `E_r[α^j] > 0` unless `c ≡ 0`.

Every moment is a *positive number*; it cannot vanish. No comparison of terms, so the collapse
of §2 is irrelevant. Verified even for `c = r²−3` (which changes sign): `α = r·c² ≥ 0` because
the **square** is what enters. This is opus/boxeph's route (Hankel PD, THM-1535/1540) — recorded
here as the *correct* mechanism standing where domination fell.

## 4. The precise gap — sign-indefiniteness, not growth

With `β ≠ 0`, `L_m(α,β) = Σ_k C(m,k)²(m−2k)!/… · α^k β^{m−2k}` is no longer a positive
combination. On `α = r`, `β = 1−r`:

| `m` | signs of `L_m(r)` coeffs | `E_r[L_m]` |
|---|---|---|
| 2 | `{+}` | 3 |
| 3 | `{−,+}` | −8 |
| 4 | `{−,+}` | 57 |
| 5 | `{−,+}` | −384 |

So `L_m(r)` is **sign-indefinite**. Positivity is unavailable, domination is dead, and the
remaining content of GMC(2) is exactly:

> Show `E_r[L_m] ≠ 0` for some `m` whenever `(α,β)` is not the trivial `(r|c|², 0)`, using
> **neither** termwise positivity **nor** top-term domination.

The known partial routes are opus's Hankel argument on the charge-0 sub-block and boxeph's
pinch/Watson analysis. This file's contribution is threefold and modest: it makes the algebraic
core rigorous (§1), kills the domination mechanism cleanly and retracts my own use of it (§2),
and reframes the gap as **sign-indefiniteness of the integrand** rather than anything about
growth rates (§4). That reframing matters: it says the obstruction lives in the *charge-0
coupling*, which is precisely where opus's Hankel block acts.

## 5. Scope

GMC(2) is not closed. What changes: the bridge's algebra is now a proved identity, the wrong
mechanism is retired (mine included), and the residual analytic task is stated in the terms —
sign-indefiniteness of `L_m` — where the surviving methods actually operate.

*Files: `04-computation/bridge_mechanism_klein_S363.py` (+ `.out`).*
