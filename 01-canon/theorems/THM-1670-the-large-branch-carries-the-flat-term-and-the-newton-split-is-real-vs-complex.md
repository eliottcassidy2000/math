---
id: THM-1670
title: "THE LARGE BRANCH CARRIES THE FLAT TERM, AND THE NEWTON-POLYGON SPLIT B vs (A+1)/2 IS EXACTLY REAL-BRANCH vs COMPLEX-BRANCH. On the {−1,0,1} radial span, Σ_m E[P^m]t^m = E_r[D^{−1/2}], D = (1−βt)²−4αt². NC2 asks whether this equals 1; the obstruction is a FLAT term (the imaginary part for real t) carried by the large branch r*(t) = smallest r with D=0. (1) The flat term is present EXACTLY off the one-sided locus among the α-DOMINANT stratum — controlled: one-sided (α=0) gives D≥0, no branch, flat term 0; α-dominant two-sided gives a real r*(t) and a nonzero flat term ~ e^{−r*(t)}·(algebraic prefactor). (2) THE BRANCH TYPE IS A NEWTON-POLYGON DATUM: D<0 on the contour requires 4α>β² at large r, i.e. deg α = A+1 > 2B, i.e. B < (A+1)/2. When B < (A+1)/2 (α dominates) r*(t) is REAL — a flat term on the contour; when B > (A+1)/2 (β dominates) r*(t) is COMPLEX — an oscillatory Stokes term with conjugate-branch interference. That split IS boxeph's THM-1620 Case I (dodge) / Case II (jump). (3) The sign-indefinite THM-1640 gap (β = 1−r) is β-DOMINANT, hence complex-branch — which is precisely why positivity, domination, AND the real flat-term argument all miss it, and why it is boxeph's hard case."
status: >
  MECHANISM + STRUCTURAL DICHOTOMY, numerically verified (mpmath, 40 digits).
  (1) The flat term / no-flat-term split across the one-sided locus is VERIFIED on controls;
      the claim "nonzero flat ⟹ NC2 fails" is NOT proved here — it is the Borel–Laplace step
      (boxeph's Watson frontier). So the α-dominant regime is REDUCED to that step, not closed.
  (2) The B vs (A+1)/2 branch-type dichotomy is elementary (large-r dominance) and verified:
      real branch for (0,0),(0,0),(2,1); complex for (0,1),(0,2).
  (3) Reclassification of the THM-1640 gap as complex-branch — verified (β=1−r ⟹ r* complex).
  ACCURACY NOTE: an interim print said the flat/e^{−r*} ratio is "O(1)"; it actually grows
  like r*^{1/2} (the saddle amplitude 1/√D''). Corrected — the flat term is
  e^{−r*(t)}·(algebraic prefactor), exponential order set by the large branch.
  GMC(2) REMAINS OPEN.
source: klein-2026-07-20-S367 (owner: think Newton polygon factorization of large branches)
depends_on:
  - THM-1660  # the branch exponent e = max(B,(A+1)/2); this adds the real/complex split
  - THM-1640  # the sign-indefinite gap, here reclassified as complex-branch
related:
  - THM-1620  # boxeph: Case I (dodge) / Case II (jump) — this file's real/complex split
  - THM-1615  # boxeph: Hermite closure (constant β)
script: 04-computation/large_branch_flat_term_klein_S367.py (+ .out)
---

> **⚠ RENUMBERED THM-1665 → THM-1670 (klein-S369):** mac-mini's per-component Watson THM-1665 was pushed 2 min earlier (15:11 vs 15:13). This file is THM-1670.

# THM-1670 — the large branch carries the flat term

## The picture

On the `{−1,0,1}` radial span, `Σ_m E[P^m] t^m = E_r[D(r,t)^{−1/2}]`,
`D = (1−βt)²−4αt²`, `α = r·a·c`, `β = b`. NC2 (all `E[P^m]=0`) says this asymptotic series is
`≡ 1`. But the integral is **not** its own asymptotic series — it can differ by a **flat term**
`e^{−r*(t)}·(…)`, and that flat term is the imaginary part for real `t`, from the region where
`D(r,t) < 0`. Its edge is the **large branch** `r*(t) = smallest r>0 with D=0`, located by the
Newton polygon (THM-1660: `r* ~ (κt)^{−1/e}`, `e = max(B,(A+1)/2)`).

## 1. The flat term is present exactly off the one-sided locus (α-dominant)

Verified at `t = 0.15` (40-digit mpmath):

| `P` | `r*(t)` | `Im E_r[D^{−1/2}]` |
|---|---|---|
| one-sided `α=0, β=0` (nullcone) | none | **0** |
| one-sided-in-α `α=0, β=1` | none | **0** |
| two-sided `α=r, β=0` (the `{−1,+1}` pair) | 11.11 | `8.8e−5` |
| two-sided `α=r, β=1` | 8.03 | `1.9e−3` |

`α=0 ⟹ D = (1−βt)² ≥ 0`, no `D<0` region, no flat term. `α=r ⟹ D<0` for large `r`, real
`r*(t)`, **nonzero flat term**. The flat term is the obstruction to `E_r[D^{−1/2}] = 1`.

## 2. The branch type is a Newton-polygon datum — real vs complex

`D<0` on the contour requires `4α > β²` at large `r`, i.e. `deg α = A+1 > 2B = deg β²`:

> **`B < (A+1)/2` (α dominates) ⟹ `r*(t)` REAL** — a flat term on the contour.
> **`B > (A+1)/2` (β dominates) ⟹ `r*(t)` COMPLEX** — an oscillatory Stokes term.

Verified at `t=0.1`:

| `(α,β)` | `(A,B)` | `B` vs `(A+1)/2` | branch |
|---|---|---|---|
| `α=r, β=0` | (0,0) | `0 < 0.5` | REAL `25.0` |
| `α=r, β=1` | (0,0) | `0 < 0.5` | REAL `20.3` |
| `α=r, β=1−r` | (0,1) | `1 > 0.5` | **CPLX** `−7.0+5.7i` |
| `α=r, β=r²` | (0,2) | `2 > 0.5` | **CPLX** `2.6i` |
| `α=r³, β=r` | (2,1) | `1 < 1.5` | REAL `2.4` |

**This split IS boxeph's THM-1620 Case I (dodge, analytic continuation) / Case II (jump,
contour pinch)** — now read off the Newton polygon of `(α,β)` rather than from the arc
geometry. When the pair scale wins (`B < (A+1)/2`) the branch pinches the real contour; when
the charge-0 coefficient wins (`B > (A+1)/2`) the branch stays complex and the flat
contributions come in conjugate pairs that can interfere.

## 3. Why the THM-1640 gap is the hard case

The sign-indefinite gap of THM-1640 — `α = r`, `β = 1−r` — has `B = 1 > (A+1)/2 = 1/2`, so it
is **β-dominant**, hence **complex-branch**. That is exactly why every real-space tool misses
it:

- **positivity** (THM-1640 §3) needs `α ≥ 0` alone — broken by `β`;
- **domination** (MISTAKE-202) is false at every span;
- **the real flat term** (this file §1) is absent — the branch is off the real axis.

The obstruction is an oscillatory Stokes term from a conjugate pair of complex branches, and
whether it can vanish (the conjugate contributions cancelling) is precisely boxeph's
per-component Watson lemma. The Newton polygon says this is unavoidable: any two-sided `P` with
`β` dominant lands here.

## 4. The α-dominant flat term is nonzero of the predicted order

For `α=r, β=0`, the flat term tracks `e^{−r*(t)}` up to an algebraic prefactor:

| `t` | `r*(t)` | `Im E_r[D^{−1/2}]` | `e^{−r*}` |
|---|---|---|---|
| 0.20 | 6.25 | `8.6e−3` | `1.9e−3` |
| 0.15 | 11.1 | `8.8e−5` | `1.5e−5` |
| 0.10 | 25.0 | `1.2e−10` | `1.4e−11` |
| 0.07 | 51.0 | `8.8e−22` | `7.0e−23` |

The ratio grows like `r*^{1/2}` (the saddle amplitude `1/√{D''}`), so the flat term is
`e^{−r*(t)}·(algebraic)`, nonzero at every `t`, of exponential order set by the large branch —
`e^{−r*} ~ e^{−c/t²}` at `(A,B)=(0,0)`, matching boxeph's `e^{−(C/t)²}`. *(Interim note: an
earlier print called the ratio "O(1)"; it is `~r*^{1/2}`. Corrected.)*

## 5. Scope

Not a proof of GMC(2). What this establishes: the obstruction to the radial nullcone being
one-sided is a **flat term carried by the large branch**; the Newton-polygon split
`B vs (A+1)/2` decides whether that branch is **real** (α-dominant → a contour flat term,
reducing NC2 to the Borel–Laplace step) or **complex** (β-dominant → an oscillatory Stokes
term, boxeph's hard case); and the sign-indefinite THM-1640 gap is precisely the complex-branch
regime. This gives a clean, checkable classification of *where* the difficulty lives, and ties
the elimination, Watson, and flat-term pictures to one polygon datum.

*Files: `04-computation/large_branch_flat_term_klein_S367.py` (+ `.out`).*
