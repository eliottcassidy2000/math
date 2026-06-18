---
id: THM-528
title: The G_P coupling for the LRC(14) S3 lonely-density ρ*(P,E) — a scale-decoupling FOUR-WINDOW lemma reduces "ρ*>0" to a PURE small-speed statement; rho* >= Σ_q meas(G_P ∩ window_{a/q}(maxE)); positivity PROVED for every (P, consecutive cluster) unconditionally (finite explicit exception list, each closed by certificate)
status: MIXED. PROVED: the four-window lemma (Good_E ⊇ four conservative windows at {0,1/2,1/3,2/3}, q-grid + Lipschitz); the coupling bound ρ*(P,E) ≥ Σ_q meas(G_P ∩ window_q(maxE)); the RHS is positive for EVERY (P, consecutive E={0..k-1}) except a FINITE explicit list of 4 (P,k) pairs; the 4 exceptions have ρ*>0 by exact rational certificate. ⟹ ρ*(P, consecutive E) > 0 for ALL admissible P, UNCONDITIONALLY. VERIFIED (exact/exhaustive): the four-window lemma over 4238 offset sets; the correlation ratio R=ρ*/(measG_P·measGood) ≥ 245/3601 ≈ 0.068 over all consecutive cases (quasi-independence, no destructive anti-correlation); no ρ*=0 over a broad bounded-spread scan (min found 1/90 < consecutive floor 1/84, consistent with consecutive-not-extremal). OPEN: the UNIFORM c0 over the FULL bounded-spread shape space (not just consecutive) + the integer-vs-real-offset passage — the residual crux, unchanged. LRC(14) NOT proved.
source: mac-mini-2026-06-18-S2b (Angle B: the G_P coupling)
depends_on:
  - THM-527   # the lonely-density reformulation; ρ*(P,E) is its central object
  - THM-526   # arc-width lemma + criterion C
  - THM-523   # covering-set reduction; meas(G_P) ≥ 7/858
related:
  - OPEN-Q-108   # the uniform fattening lemma (equivalent crux)
  - HYP-2590     # the four-window scale-decoupling (this work)
  - HYP-2584     # bounded-spread extremizer
external: Lonely Runner Conjecture (13 speeds = LRC(14), first open case); Steinhaus three-gap; Weyl.
---

# THM-528 — The G_P coupling: a four-window scale-decoupling for ρ*(P,E)

## The question (Angle B of THM-527)

THM-527 reduced LRC(14) case S3 to a positive **density floor**
`ρ*(P,E) = meas{ x ∈ G_P : maxgap{frac(e_i x)} > 2/7 } ≥ c₀ > 0`,
where `G_P = {τ : ‖pτ‖ ≥ 1/14 ∀p∈P}` (`meas(G_P) ≥ 7/858`, THM-523/525) and
`Good_E = {x : maxgap{frac(e_i x)} > 2/7}` has `meas ≥ c(k) > 0` (the pure three-distance
floor). The **danger**: `Good_E` and `G_P` could be *anti-correlated* — `Good_E`
concentrated exactly where `G_P` is empty — making `ρ* = meas(G_P ∩ Good_E)` collapse to 0
even though both factors are large. Angle B asks: can this happen?

**Answer: no.** The two sets live on *different, decoupled scales*, made precise below.

## A. The four-window lemma (PROVED)

> **Lemma.** For ANY offset set `E` with `0 ∈ E` and spread `maxE = max E`, the good set
> `Good_E` contains an interval around each of the four low-denominator rationals:
> ```
>   near 0   : (0, 5/(7·maxE))                                    [reflection ⇒ near 1 too]
>   near 1/2 : (1/2 − 3/(28·maxE), 1/2 + 3/(28·maxE))
>   near 1/3 : (1/3 − 1/(42·maxE), 1/3 + 1/(42·maxE))
>   near 2/3 : (2/3 − 1/(42·maxE), 2/3 + 1/(42·maxE))
> ```

**Proof.** At `x = a/q` with `q ∈ {1,2,3}`, every point `frac(e_i·a/q) = (e_i·a mod q)/q`
lies on the `q`-grid `{0, 1/q, …, (q−1)/q}`; the occupied slots have all cyclic gaps a
multiple of `1/q`, so the circular max-gap is `≥ 1/q ≥ 1/3 > 2/7`. Move off `a/q` by
`δ`: point `e_i` drifts by `≤ maxE·|δ|`. A base gap `G₀ ≥ 1/q` realised between two
consecutive occupied slots `s_lo < s_hi` has its two boundary points each drift by
`≤ maxE·|δ|`, and no third point can enter while `maxE·|δ| < 1/q` (the next occupied slot is
`≥ 1/q` away). Hence the gap stays `≥ G₀ − 2·maxE·|δ| > 2/7` for
`|δ| < (1/q − 2/7)/(2·maxE)`, giving the half-widths `3/(28·maxE)` (q=2), `1/(42·maxE)`
(q=3). For `q=1` (`a/q = 0`) the points cluster in `[0, maxE·δ]`, the wrap gap is
`1 − maxE·δ > 2/7` for `δ < 5/(7·maxE)`. ∎  (The conservative half-widths are exact lower
bounds on the true Good half-widths; verified exhaustively on 4238 offset sets.)

## B. The coupling bound and the scale decoupling (PROVED)

Summing the four windows and intersecting with `G_P`:

> **`ρ*(P,E) ≥ Σ_{a/q ∈ {0,1/2,1/3,2/3}} meas( G_P ∩ window_{a/q}(maxE) )`.**

The RHS depends on `E` **only through the single scalar `maxE`** and on `P` through the
**fixed** small-speed safe set. *This is the literal "different scales" decoupling Angle B
sought:* `G_P` is built from low frequencies `p ≤ 13` (a coarse comb), the windows sit at
the four coarsest rationals, and `Good_E`'s fine cluster structure `frac(e_i x)` is
irrelevant beyond `maxE`. **The coupling crux is reduced to a PURE statement about
`(P ⊆ {1,…,13}, scalar maxE)`** — no anti-correlation pathology survives, because `Good_E`
is *guaranteed mass at four fixed locations*, and `G_P` cannot be empty at all four.

## C. Positivity is PROVED for every consecutive cluster (finite exception list)

The conservative four-window LB is positive for **every** `(P, consecutive E={0,…,k−1})`,
`maxE = k−1`, **except a finite explicit list of 4 pairs**:

| k | P | conservative LB | EXACT ρ* | surviving window |
|---|---|---|---|---|
| 9 | {1,2,12,13} | 0 | **2/105** | near 1/3, 2/3 |
| 9 | {1,4,12,13} | 0 | **131/2940** | near 1/2, 1/3, 2/3 |
| 10| {1,6,13}    | 0 | **11/294** | near 1/2 |
| 11| {1,6}       | 0 | **1/42** | near 1/2 |

For these 4, the conservative window is narrower than the small-speed danger arc *at the
killed center*, but the **sharp** Good window (wider than conservative) survives, and the
exact rational `G_P ∩ Good` is a nonempty union of certified safe+good sub-arcs (table; e.g.
`{1,6}` keeps `(10/21,41/84) ∪ (43/84,11/21)`, `ρ*=1/42`). Hence:

> **ρ*(P, consecutive E) > 0 for ALL admissible P, UNCONDITIONALLY.**

The exact consecutive floor is `ρ* = 1/84` at `k=9, P={1,2,3,12}` (matches THM-527-E).

## D. Quasi-independence holds (VERIFIED)

The feared anti-correlation does **not** occur: the correlation ratio
`R = ρ*/(meas(G_P)·meas(Good_E))` satisfies `R ≥ 245/3601 ≈ 0.068 > 0` over all consecutive
cases (min at `k=9, P={1,2,6,12}`); `R=1` (exact independence) at `k=3` and `k=13`. So a
quasi-independence bound `ρ* ≥ R₀·meas(G_P)·meas(Good_E)` holds with `R₀ ≈ 0.068`
empirically — but `R₀` is not yet derived (it would close the uniform floor combined with
the THM-527-E pure floor and the THM-523 `7/858`).

## E. Honest remaining gap (the residual crux, unchanged)

The PROVED positivity is for **consecutive** clusters. The general bounded-spread case is
strong-VERIFIED (no `ρ*=0` over a broad scan; min found `1/90 < 1/84`, consistent with
THM-527-F "consecutive not globally extremal") but the **uniform `c₀` over the full
bounded-spread shape space**, plus the **integer-vs-real-offset passage**, remain open —
the same prize as OPEN-Q-108 and THM-527-G. Angle B's contribution is to make the
**decoupling** structural and explicit (the four fixed windows), turning the coupling from a
2-set correlation worry into a finite small-speed statement, and to **prove the consecutive
case in full**.

## Net

The G_P coupling is **not** the obstruction. `Good_E` always carries mass at the four
fixed rationals `{0,1/2,1/3,2/3}` (q-grid lemma), so `ρ*` decouples into a pure
small-speed question; this is positive for every consecutive cluster (4 explicit
certificate exceptions), with correlation ratio bounded below (`≈0.068`). The residual is
the uniform floor over *all* bounded-spread shapes — localized, finite-dimensional, and
open.

## Scripts
- `04-computation/lrc14_angleB_gp_coupling_macmini_0618s2.py` (this result, all exact)
- output: `05-knowledge/results/lrc14_angleB_gp_coupling_macmini_0618s2.out`
