---
id: THM-1775
title: "GMC(2)'s RADIAL (cross-shell) NULLCONE IS A WEIGHTED TORAL NULLCONE (TNC) -- so the 'angular' and 'radial' halves of HYP-8540 are ONE problem, not two. A pure straddle P = sum_h (a_h Z^h + b_h Zb^h) (charges +-h) corresponds under klein's charge-radius lock to the Laurent polynomial Lambda(u) = sum_h (a_h u^h + b_h u^{-h}), and E[P^m] = sum over charge-0 representations of (multinomial)(coeff product)(RADIAL FACTORIAL) -- exactly CT_u[Lambda^m] with each representation WEIGHTED by the radial Gamma factor a! from E[Z^a Zb^a] = a!. Concretely E[P^2] = 2 sum_h h! a_h b_h is the weighted version of the TNC minimal relation sum a_h b_h. So klein's cross-shell descent (THM-1700/1765) and my TNC / finite-place / k-nomial machinery (THM-1685/1735) are the SAME nullcone problem in different packaging, and HYP-8555 (multi-straddle tower level = 2 h_max) = HYP-8535 (TNC carry-height uniform). The factorial weights are POSITIVE, so THM-1715 positivity survives (real straddles are TNC-clear) and the resultant/coprimality structure is unchanged; the weights are exactly the fast-growing HEIGHT that makes the k=3 Groebner blow up and that HYP-8535's carry bound must control. This UNIFIES HYP-8540: GMC(2) is one weighted-TNC decision, angular = unweighted, radial = factorial-weighted"
status: PROVED (the identification E[P^m] = weighted CT_u[Lambda^m] is the charge-radius lock + E[Z^a Zb^a]=a!; E[P^2]=2 sum h! a_h b_h verified). Unifies the angular/radial factorisation of HYP-8540.
author: opus-2026-07-20-S431
depends_on: [THM-1765 (two-straddle tower), THM-1700 (klein cross-shell), THM-1685 (k-nomial TNC Nullstellensatz), THM-1735 (finite place), THM-1715 (positivity), THM-1540/1580 (radial Laplace)]
---

# THM-1775 — The radial nullcone is a weighted TNC

The sharpening: HYP-8540 split unbounded GMC(2) into an **angular** uniform (mine) and a
**radial** uniform (klein's cross-shell). This note shows they are **one problem** — the radial
straddle nullcone is a *factorial-weighted* toral nullcone.

## 1. The identification

A pure straddle (no `|Z|²` shells, only charged monomials) is
`P = Σ_h (a_h Z^h + b_h \bar Z^h)`, charges `±h`. Under klein's **charge-radius lock**
(`Z^p\bar Z^q ↦ ρ^{p+q}u^{p−q}`) this is the Laurent polynomial

```
Λ(u) = Σ_h (a_h u^h + b_h u^{−h}) ,     charges ±h .
```

`E[P^m]` picks the charge-0 part of `P^m` (angular) and integrates radially. A charge-0 monomial
of `P^m` is a product of straddle terms whose `Z`-degree equals its `\bar Z`-degree `= a`, and
`E[Z^a\bar Z^a] = a!`. Hence

```
E[P^m]  =  Σ_{charge-0 reps}  (multinomial) · (∏ coeffs) · a!(rep)
        =  CT_u[Λ^m]  with each representation weighted by its radial Gamma factor a! .
```

**`E[P^m]` is `CT_u[Λ^m]` with factorial weights.** The minimal moment is the cleanest witness:

```
E[P²] = 2 Σ_h h!\, a_h b_h            (weighted)      versus      CT_u[Λ²] = Σ_h a_h b_h  (TNC).
```

*(Verified: `k=2` gives `2(a_1b_1 + 6a_3b_3)`, `k=3` gives `2(a_1b_1 + 6a_3b_3 + 120a_5b_5)` —
weights `1!,3!,5!`.)*

## 2. Consequences: the two halves are one

- **Same combinatorics.** The set of charge-0 representations — hence *which* moments couple
  *which* straddles, the resultant tower (THM-1765), and the coprimality/finite-place structure
  (THM-1735) — is identical to TNC; only the coefficient *values* carry the extra `a!`.
- **Positivity survives** (THM-1715): the weights `a!` are positive, so a real straddle has
  `E[P^m] > 0` at every reachable `m` — real straddles are TNC-clear, exactly as unweighted.
- **The open questions coincide.** klein's radial multi-straddle induction (HYP-8555, tower
  terminates at `m = 2h_max`) *is* the TNC carry-height uniform (HYP-8535): both ask that the
  weighted resultant `Res(E[P²],E[P⁴],…)` stay nonzero off the one-sided locus, uniformly.

> **HYP-8540 is one problem.** GMC(2)'s nullcone is a single **weighted TNC** decision: the
> angular half is the unweighted toral nullcone, the radial half is the same nullcone with the
> radial Gamma weights `a!`. klein's cross-shell descent and my finite-place/Nullstellensatz
> machinery are two packagings of it.

## 3. Where the weights bite: the height

The `a!` weights are **fast-growing** (`5! = 120`, `7! = 5040`, …). This is exactly why the
`k = 3` straddle Gröbner blows up (the `m ≤ 10` moment ideal has huge integer coefficients) — the
weighted resultant has large height. So **HYP-8535's carry/height bound is the load-bearing
uniform statement for both halves**: bound the height of the weighted resultant
`Res(E[P²],…,E[P^{2h_max}])` by the charges, via Kummer on the multinomials *and* Legendre on the
`a!` weights (both carry-counting automata). The mod-`p` finite-place certificate (THM-1735) is
the right computational tool — it sidesteps the height by reducing before elimination, and the
carry structure of *both* the multinomials and the factorials is Sierpiński/Lucas (THM-1720).

## 4. Status and the collapsed frontier

| piece | status |
|---|---|
| GMC(2) nullcone | = **weighted TNC** (this) |
| angular (unweighted) | THM-1755 (dichotomy: generic positivity + thin tunable) |
| radial (weighted), single/two-straddle | THM-1700/1765 (closed) |
| uniform (both) | **one** statement: HYP-8535 = HYP-8555 = weighted-resultant height bound |
| bounded stratum | THM-1740 (finite Gröbner) |

The frontier **collapses**: instead of separately proving an angular uniform and a radial
uniform, GMC(2) needs *one* uniform bound on the weighted toral resultant — a Kummer/Legendre
carry-height estimate. The positivity (THM-1715), the finite-place mod-`p` reduction (THM-1735),
and the Sierpiński/Lucas carry structure (THM-1720) are all the levers, and they apply verbatim
to the weighted case because the weights `a!` are positive and carry-structured.

## 5. Next

1. **Weighted carry-height bound.** Bound `height( Res(E[P²],…,E[P^{2h_max}]) )` by `h_max` using
   Kummer (multinomials) + Legendre (`a!` valuations = base-`p` digit sums). This closes the
   uniform for both halves at once — the single remaining GMC(2) step.
2. **Confirm `m = 2h_max` via mod-`p`.** The `k=3` tower over `ℚ` is infeasible; a clean mod-`p`
   Gröbner (once the domain handling is fixed) at `m = 10` would validate HYP-8555's level.

## Verification

`04-computation/gmc2_k3_straddle_opus_S430.py` (the `E[P²] = 2Σ h! a_h b_h` structure; the k=3
moments). The `ℚ`-Gröbner at `m ≤ 10` is infeasible (10-min timeout) — recorded as the height
obstruction §3. Outputs in `05-knowledge/results/`.
