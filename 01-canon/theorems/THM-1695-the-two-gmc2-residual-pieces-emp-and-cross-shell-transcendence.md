---
id: THM-1695
title: "THE TWO GMC(2) RESIDUAL PIECES, BOTH RADIAL, WORKED. GMC(2) = HYP-8350 + HYP-8470 (the angular layer is DvdK-closed, THM-1645). PIECE 1 (HYP-8350, the radial one-variable Laplace nullcone) IS MY OWN EMP (THM-1510): the stated obstruction ker(L)≠0 — L(s−1)=1!−0!=0 — is exactly what EMP overcomes, since L(p^m)=0 for ALL m forces p=0 even though L has a kernel; verified that every kernel element fails at a higher moment, and that EMP holds under the Rayleigh weight 2ρe^{−ρ²} too. So Piece 1 is CLOSED for polynomial p. PIECE 2 (HYP-8470, the cross-shell coupling) worked on the family Λ_s = √s·u + (a/u+b+cu) (top shell one-sided, bottom straddles): closed EXACTLY over ℂ by elimination (the ideal ⟨E[P^m]⟩ contains (ac)^k, so every nullcone member has ac=0 = one-sided), all signs; and — the mechanism — closed at the SECOND MOMENT ALONE over ℚ by the TRANSCENDENCE of √π: E[P²] = (2ac+b²) + a√π, so vanishing over ℚ forces a=0 then b=0 (one-sided), because the half-integer Gamma weights Γ(k/2+1) put consecutive shells on ℚ-linearly-independent lines. The definite-sign sub-locus (a,b,c ≥ 0) gives positive EVEN moments."
status: >
  PIECE 1: CLOSED for polynomial p — it is EMP (THM-1510, proved by me via Laplace; boxeph
  re-proved via Hermite, THM-1615). This file's contribution is the IDENTIFICATION that
  HYP-8350's ker(L)≠0 obstruction is exactly EMP's content, plus the Rayleigh-weight extension
  (Gröbner, deg ≤ 2, variety = origin). If HYP-8350's object is a single polynomial p, done;
  if it is broader (algebraic p), that is flagged as the residual.
  PIECE 2: the specific cross-shell family Λ_s = √s·u + (a/u+b+cu) is CLOSED — over ℂ by exact
  Gröbner ((ac)^k in the ideal), over ℚ at m=2 by the √π transcendence. This is ONE shell
  family, not the general HYP-8470; the general case (higher/multiple straddling shells) is
  NOT closed here.
  SELF-CORRECTION: an interim claim "a,b,c ≥ 0 ⟹ E[P^m] > 0 for all m" is FALSE (odd moments
  vanish by parity at b=0); corrected to even moments. Control caught it. See §2.3.
  GMC(2) REMAINS OPEN.
source: klein-2026-07-20-S369 (owner: work the 2 GMC(2) residual pieces)
depends_on:
  - THM-1510  # EMP — this is Piece 1
  - THM-1645  # the polar bridge factoring GMC(2) into angular (DvdK) + radial (the two pieces)
related:
  - THM-1615  # boxeph: Hermite closure of the radial layer (same as EMP)
  - THM-1670  # klein: the large-branch flat term / Newton split (the other radial view)
  - "klein-S353 Vieta note: transcendence upgrades disjunction to rigidity — the §2.2 mechanism"
script: 04-computation/gmc2_two_residuals_klein_S369.py (+ .out)
---

> **⚠ RENUMBERED THM-1675 → THM-1695 (klein-S371), triple collision with opus THM-1675 (trinomial-gcd) and a third.**
>
> **⚠ SCOPE CORRECTION TO PIECE 2 (klein-S371): the family `Λ_s = √s·u + (a/u+b+cu)` is UNLOCKED — not a valid `Λ_s`.** A genuine `Λ_s(u) = P(√s·u, √s/u)` obeys the charge–radius LOCK: a monomial `Z^pZ̄^q ↦ ρ^{p+q}u^{p−q}`, so charge `p−q` and radius power `p+q` have the SAME PARITY. My `h=0` shell `(a/u+b+cu)` carries charges `±1` at `ρ^0` — impossible (charge ±1 needs odd `ρ`). So the `√π` mechanism is a fact about the ABSTRACT (unlocked) nullcone, NOT about GMC(2): for genuine locked `P`, charge balance forces CT_u[Λ_s^m] to have only EVEN `ρ`-powers, i.e. integer `s`-powers, so the radial functional is the pure exponential (integer moments, NO `√π`). The genuine cross-shell descent is done in THM-1700. Piece 1 (= EMP) is unaffected by this correction.

# THM-1695 — the two GMC(2) residual pieces

Via the polar bridge (THM-1645), `E[P^m] = ∫₀^∞ CT_u[Λ_s^m] e^{−s} ds`, and the **angular**
functional `CT_u` is DvdK-closed uniformly in `s`. So GMC(2) reduces to two **radial** pieces:
`HYP-8350` (the one-variable Laplace nullcone) and `HYP-8470` (the cross-shell coupling).

## 1. Piece 1 (HYP-8350) is EMP

The radial Laplace nullcone asks: `L(p^m) = 0` for all `m ≥ 1` ⟹ `p ≡ 0`, where
`L(g) = ∫₀^∞ g(s) e^{−s} ds` (so `L(s^k) = k!`). death-star flags the obstruction as
**`ker(L) ≠ 0`**: `L(s−1) = 1! − 0! = 0`, so integrated vanishing at one `m` does not force
pointwise vanishing.

> **This is exactly EMP (THM-1510).** EMP is the statement that the *full family*
> `L(p^m) = 0 ∀m` forces `p = 0` — precisely overcoming `ker(L) ≠ 0`.

Verified: every kernel/near-kernel element fails at a higher moment —
`p = s−1` gives `L(p^m) = 0, 1, 2, 9, 44, 265` (first nonzero at `m=2`). EMP was proved by me
(Laplace: `L(p^m) ~ c_d^m (dm)! e^{c_{d−1}/(c_d d)}`, amplitude nonzero) and re-proved by boxeph
(Hermite, THM-1615). And it survives the change of radial variable: under the **Rayleigh**
weight `2ρe^{−ρ²}` (the `|Z|` law, if the radial variable is `ρ = √s`), Gröbner over `ℚ` gives
variety = origin for `deg p ≤ 2`. **Piece 1 is closed for polynomial `p`.** If HYP-8350's object
is a single polynomial, done; a broader (algebraic) `p` is the only residual.

## 2. Piece 2 (HYP-8470) — the cross-shell coupling

The open step: the top shell `λ_D` is one-sided (DvdK-safe, its `CT` vanishes) but lower shells
straddle, and `L` mixes the cross-shell `CT_u[λ_D^{m−j} λ_{D'}^j]`. Worked on the family

```text
Λ_s = √s·u  (top shell, ONE-SIDED, charge +1)  +  (a/u + b + c·u)  (bottom shell, straddles).
E[P^m] = Σ_k C(m,k)·Γ(k/2+1)·[u^{−k}](a/u+b+cu)^{m−k},   top term CT_u[u^m] = 0.
```

### 2.1 Closed over ℂ by elimination

Treating `E[P^m]` (`m = 1..7`) as polynomials in `(a,b,c)`: the Gröbner basis contains
`(a·c)^k`, so the variety lies in `{ac = 0}` — one of the outer charges absent, i.e.
**one-sided**. No two-sided cross-shell nullcone member exists. Closed, all signs. (This is the
"normalise-or-radicalise" test of S361: the degenerate locus is `{ac=0}`, and radical
membership is the right predicate.)

### 2.2 The mechanism — transcendence of `√π`, closing it at `m = 2` over ℚ

```text
E[P²] = (2ac + b²) + a·√π.
```

`√π = 2Γ(3/2)` is **transcendental**, so over `ℚ` this vanishes iff *both* parts vanish:
`a = 0` (from the `√π` coefficient) and then `2ac + b² = b² = 0`, i.e. `a = b = 0`. But `a = 0`
leaves `Λ_s = (√s + c)u`, **pure charge +1 — one-sided**. So over `ℚ` the family is closed at the
**second moment alone**.

The reason is structural: `Γ(k/2+1)` is *rational* for `k` even and *rational·√π* for `k` odd,
so **consecutive shells sit on `ℚ`-linearly-independent lines** and cannot cancel unless the
outer-charge coefficient itself vanishes. This is the Vieta/transcendence pattern of klein-S353
(*one relation buys a disjunction; a family buys rigidity*) appearing exactly where HYP-8470
lives — the half-integer Gamma values are the family of relations that force rigidity. *(Over
`ℂ`, `m=2` alone is insufficient — e.g. `a=1,b=0,c=−√π/2` kills `E[P²]` — and the elimination of
§2.1 is what closes it.)*

### 2.3 The definite-sign sub-locus — and a correction

For `a,b,c ≥ 0` (`a,c` not both 0), `[u^{−k}](a/u+b+cu)^{m−k}` is a nonnegative lattice-path
count, so the **even** moments `E[P^{2j}] = ` (positive sum, with `k=0` term `2ac+b² > 0`) are
strictly positive — some moment is nonzero, all NC2 needs.

*Correction on record:* I first wrote "`a,b,c ≥ 0 ⟹ E[P^m] > 0` for **all** `m`"; the control
returned False, because at `b = 0` the **odd** moments vanish by parity (the recurring S345
artefact). The correct statement is even moments, above.

## 3. Scope

- **Piece 1:** closed for polynomial `p` (= EMP); the connection to HYP-8350's `ker(L)≠0` is the
  contribution.
- **Piece 2:** the one cross-shell family `√s·u + (a/u+b+cu)` is closed (ℂ: elimination; ℚ: `√π`
  transcendence at `m=2`). The **general** HYP-8470 — multiple straddling shells, higher top
  degree — is **not** closed here. What generalises is the mechanism: half-integer Gamma weights
  put shells on `ℚ`-independent lines, and this should be the lever for the general cross-shell
  descent.

GMC(2) remains open. Piece 1 is identified as already-proved (EMP); Piece 2 is closed on a
concrete family with a mechanism (Gamma-half-integer transcendence) that points at the general
case.

*Files: `04-computation/gmc2_two_residuals_klein_S369.py` (+ `.out`).*
