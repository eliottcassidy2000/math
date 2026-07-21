---
id: THM-1875
title: "THE TRANSITIVE TOURNAMENT'S SKEW CHARACTERISTIC FORM IS THE BINOMIAL REFLECTION OF THE NULLCONE MONOMIAL — char_S = ((x+1)ⁿ + (x−1)ⁿ)/2, spectrum i·cot((2k−1)π/2n), reviving the forgotten era-1 tangent-number (A000182) connection and closing the even/Weyl companion of the GIT nullcone vertex. The transitive tournament is the GIT nullcone on the ADJACENCY side (char_A = xⁿ, the maximally-unstable monomial, THM-1810). Its SKEW matrix S = sign(j−i) (the all-+1-above staircase skew matrix) has, for every n: (1) char_S(x) = ((x+1)ⁿ + (x−1)ⁿ)/2 = Σⱼ C(n,2j) x^{n−2j} — the EVEN part of (x+1)ⁿ, i.e. the binomial reflection of the nullcone monomial xⁿ — verified exactly n = 3..8; (2) spectrum {i·cot((2k−1)π/2n) : k = 1..n}, a pure COTANGENT / roots-of-unity-adjacent spectrum (the classical eigenvalues of the staircase skew matrix), with a 0 eigenvalue at odd n (k = (n+1)/2 gives cot(π/2) = 0); (3) char_S(1) = 2^{n−1}. This REVIVES a niche forgotten thread found by deep OEIS-layer archaeology: an era-1 session (~S33, in the session log only, never canonized) verified 'Σ(−1)ᵏ A(n,k) = ±T_{(n+1)/2} connecting W(0) of the transitive tournament to OEIS A000182 (tangent numbers)'. The tangent numbers are exactly the odd-part generating function tan(x) = Σ T_m x^{2m−1}/(2m−1)!, and cot is its reciprocal-shift — so the transitive tournament's skew spectrum (cot at odd multiples of π/2n) and its W-generating function (tangent numbers) are two faces of the SAME odd/even (sin/cos/tan) structure that runs through char_S being even (THM-1810 §5), the half-dictionary's ½ (THM-1555), and the regular-tournament Re = −½ line. The GIT nullcone vertex xⁿ (A-side) and its binomial-reflection companion ((x+1)ⁿ+(x−1)ⁿ)/2 (S-side) are the adjacency/skew pair of the single most degenerate tournament"
status: >
  (1) VERIFIED exactly n = 3..8 (char_S matches ((x+1)ⁿ+(x−1)ⁿ)/2 on the nose; e.g. n=8 gives
  x⁸+28x⁶+70x⁴+28x²+1, n=6 gives x⁶+15x⁴+15x²+1).  The closed form is the classical
  characteristic polynomial of the staircase skew matrix and is PROVABLE (the cotangent-product
  identity ∏(x − i·cot((2k−1)π/2n)) = ((x+1)ⁿ+(x−1)ⁿ)/2); stated as verified + cited, not reproved.
  (2) VERIFIED numerically n = 3..8 (eigenvalues match cot((2k−1)π/2n) to 1e−6).
  (3) Immediate: char_S(1) = ((2)ⁿ + 0)/2 = 2^{n−1}.
  The tangent-number connection to W(0) is the era-1 result RECOVERED verbatim from the session
  log (A000182); the reciprocal cot/tan link is the reading offered here, not a re-derivation of
  that identity (which concerns the succession-GF W = THM-293, also forgotten).
  This is a REVIVAL + a clean closed form, not a new open-problem advance.
source: kind-pasteur-2026-07-21-S128c137 (owner: find niche forgotten ideas via deep archaeology; collaboratively map past ideas)
depends_on:
  - THM-1810    # transitive = GIT nullcone char_A = xⁿ; char_S is the even companion
  - THM-1555    # the half-dictionary (the ½ / odd-even axis)
related: [THM-293, THM-1725]   # THM-293 = the forgotten succession-GF W; THM-1725 = roots-of-unity exponents
external:
  - "OEIS A000182 (tangent numbers); the staircase/Kac-type skew-matrix cotangent spectrum."
script: 04-computation/deep_archaeology_kps_S128c137.py, /tmp tanrev check (+ .out)
---

# THM-1875 — the transitive skew form is the binomial reflection, and the tangent thread revived

Deep OEIS-layer archaeology (`deep_archaeology_kps_S128c137.py`, mining A-numbers mentioned once)
surfaced a niche forgotten line, present only in an era-1 session log:

> "Tangent number connection verified: `Σ(−1)ᵏ A(n,k) = (−1)^{(n−1)/2} T_{(n+1)/2}`, connecting
> `W(0)` of the transitive tournament to OEIS **A000182**."

Chasing it into the current binary-form frame gives a clean closed form and a place on the map.

## The closed form

The transitive tournament's skew matrix is `S = sign(j − i)` (the all-`+1`-above staircase skew
matrix). For every `n`:

> **`char_S(x) = ((x+1)ⁿ + (x−1)ⁿ)/2 = Σⱼ C(n, 2j) x^{n−2j}`**,

the **even part of `(x+1)ⁿ`** — the *binomial reflection* of the nullcone monomial `xⁿ`. Verified
`n = 3..8` on the nose (`n=8`: `x⁸+28x⁶+70x⁴+28x²+1`; `n=6`: `x⁶+15x⁴+15x²+1`; `n=5`:
`x(x⁴+10x²+5)`). Its roots are

> **`spec(S) = { i·cot((2k−1)π/2n) : k = 1..n }`**

— a pure cotangent spectrum (the classical staircase-skew eigenvalues), with a `0` at odd `n`. And
`char_S(1) = 2^{n−1}`.

## Why this is the even companion of the nullcone vertex

THM-1810 established: the transitive tournament is the GIT nullcone on the **adjacency** side,
`char_A = xⁿ` (nilpotent, one root of multiplicity `n`), and `char_S` is its **even/Weyl
companion** — but the nullcone cannot live on the skew side (`tr S² = −n(n−1) ≠ 0`). This gives
that companion its **closed form**: while `char_A = xⁿ` is the pure monomial, `char_S =
((x+1)ⁿ+(x−1)ⁿ)/2` is its binomial reflection. The single most degenerate tournament has the
adjacency/skew pair

> `( xⁿ , ((x+1)ⁿ+(x−1)ⁿ)/2 )`.

## The odd/even structure, unified

The revival closes a loop the owner has pulled for many sessions:

- `char_S` is **even** (`spec ⊂ iℝ`), so it is a form in `x²` — the **cos/secant** side;
- the spectrum is **cotangent** at odd multiples of `π/2n` — roots-of-unity-adjacent;
- the forgotten `W(0)` generating function is the **tangent numbers** `A000182` (`tan(x) =
  Σ T_m x^{2m−1}/(2m−1)!`, the **odd/sine** side), and `cot = 1/tan`;
- the half-dictionary's `½` (THM-1555) is the shift onto the `Re = −½` line where the regular
  tournament's `char_A` spectrum sits.

So *"tangent numbers,"* *"cotangent spectrum,"* *"`char_S` even,"* and *"the `½`"* are one
odd/even (sin/cos/tan) axis — the Weyl involution `x ↦ −x` of the characteristic binary form —
and the transitive tournament is where it is cleanest, because there both `char_A` and `char_S`
are exactly known.

## Named next

- **Revive `W` itself (THM-293, the forgotten succession GF), and re-derive `Σ(−1)ᵏA(n,k) =
  ±T_{(n+1)/2}`** in the char_S / tangent frame — the tangent identity should be `char_S` or `W`
  evaluated at a root of unity.
- **The regular tournament's `char_S`.** Transitive gives `((x+1)ⁿ+(x−1)ⁿ)/2`; the Paley/regular
  `char_S` is the Gauss-sum object (opus THM-1810 §Q2). Between them is the full skew-spectrum
  map of the tournament, indexed by the odd/even axis.
- **Secant numbers on the even side.** `A000182` (tangent) is the odd companion; `A000364`
  (secant/Euler) should be the even companion — check whether it indexes `char_S` directly.
