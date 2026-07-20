---
id: THM-1650
title: "THE NEWTON POLYGON OF THE EFFECTIVE DvdK BOUND — u^M − tR(u) has a V-shaped Newton polygon (vertices (0,1)–(M,0)–(d,1), d=M+N) with a LEFT edge of length M (slope −1/M, the M SMALL branches u ~ t^{1/M}) and a RIGHT edge of length N (slope +1/N, the N LARGE branches u ~ t^{−1/N}); verified n=1..4 that small-branch count = M and large-branch count = N exactly. The two edges are exchanged by u ↦ 1/u, which IS the duality TNC(M,N) ⇔ TNC(N,M): the constant-term series is computable from EITHER edge by residues, with small-branch sum = −(large-branch sum) since Res_∞ = 0 for N ≥ 1 (verified). The Sturmfels/ESV effective bound (first nonzero CT(Λ^m) by index m+n) is exactly the TOTAL branch count M+N = both edge lengths summed; verified ≤ M+N and never exceeded across nine bidegrees, tight at the coprime ones. BUT MY OWN 'sparsest R is extremal' HEURISTIC IS REFUTED BY THE COMPUTATION: the two-monomial R = 1+u^d (the two polygon vertices) gives first-nonzero (M+N)/gcd(M,N), which is the extremal ONLY when gcd(M,N)=1; when gcd > 1 it is fast (m=2 at (3,3)) and interior monomials RAISE the first-nonzero index in 54/3875 sampled cases. So the extremal structure is gcd-graded, and the two-vertex picture is the coprime special case."
status: >
  Newton-polygon structure (V, two edges, branch counts M and N): PROVED (monomial support
  is (M,0) and the (j,1); the lower hull is the stated V) and VERIFIED numerically n = 1..4.
  The residue duality small = −large (Res_∞ = 0 for N ≥ 1): PROVED and verified.
  The effective bound ≤ M+N: VERIFIED across nine bidegrees in a bounded coefficient box,
  NOT proved (it is the OPEN Sturmfels/ESV conjecture, HYP-8460). Tightness holds at the
  coprime bidegrees tested; at (3,3) the box max was 4 < 6, so tightness is NOT universal
  and was not observed there.
  The 'sparsest-R-is-extremal' heuristic: REFUTED by the same computation (54/3875 cases
  raise the index above the two-monomial value). The gcd-grading is the correction.
source: mac-mini-2026-07-20-S144 (owner: "think newton polygon factorization of large branches")
depends_on:
  - THM-1630  # TNC = DvdK; the effective bound is the surviving open problem
related:
  - THM-1595  # min(M,N) indexing = the polygon duality
  - THM-1610  # the coefficient ladder that finds the ESV-extremal-approaching R
script: 04-computation/newton_polygon_branches_macmini_S144.py (+ .out)
---

# THM-1650 — the Newton polygon of the effective DvdK bound

## The polygon

For `Λ = u^{−M}R(u)` (`R(0)=1`, `deg R = d = M+N`), the constant-term generating function is
`F(t) = CT_u[u^M/(u^M − tR(u))] − 1`, and the denominator `Φ(u,t) = u^M − tR(u)` has `d` roots.
Its Newton polygon in `(u`-exponent`, t`-exponent`)` has support `(M,0)` (from `u^M`) and
`(j,1)` for `j = 0..d` (from `−t r_j u^j`). The lower hull is a **V** with vertices
`(0,1) – (M,0) – (d,1)`:

| edge | slope | length | branches |
|---|---|---|---|
| left `(0,1)–(M,0)` | `−1/M` | `M` | **`M` small** branches `u ~ t^{1/M}` |
| right `(M,0)–(d,1)` | `+1/N` | `N` | **`N` large** branches `u ~ t^{−1/N}` |

Verified `n = 1..4`: exactly `M` roots inside the unit disk and `N` outside, with the predicted
valuations (e.g. `(2,3)`: small `~ 10^{−3} = t^{1/2}`, large `~ 10^{2} = t^{−1/3}`).

## The duality is the two edges

`CT_u[u^M/(u^M − tR)]` is the sum of residues of `u^{M−1}/(u^M − tR)` inside `|u|<1` — the
**small** branches (`u=0` is a zero of order `M` of the numerator, not a pole). By the residue
theorem this equals `−(`residues outside `+` at `∞)` — the **large** branches — and `Res_∞ = 0`
for `N ≥ 1` (the integrand `~ u^{−N−1}` at `∞`). Verified: small `+` large `= 0` at
`(2,2),(2,3),(3,3)`.

> **So the constant-term series is computable from EITHER edge, and `u ↦ 1/u` (which sends
> `R ↦ R*`, `M ↔ N`, small ↔ large) is exactly `TNC(M,N) ⇔ TNC(N,M)`** — the duality THM-1595
> already used to index by `min(M,N)`.

## The effective bound is the total branch count

The surviving open problem (HYP-8460, Sturmfels; Erman–Smith–Várilly-Alvarado arXiv:0908.2609):
the first nonzero `CT(Λ^m)` occurs by index `m+n`. In this language `m+n = M+N = d = ` **the
total branch count = both edge lengths summed.** Searching `R` that delays the first nonzero
longest:

| `(M,N)` | `d` | max first-nonzero (in box) | `≤ d`? | `= d`? |
|---|---|---|---|---|
| (1,1) | 2 | 2 | yes | yes |
| (2,2) | 4 | 4 | yes | yes |
| (2,3) | 5 | 5 | yes | yes |
| (3,2) | 5 | 5 | yes | yes |
| (1,4) | 5 | 5 | yes | yes |
| **(3,3)** | **6** | **4** | yes | **no** |

Never exceeded; tight at the coprime bidegrees; **not tight at `(3,3)`** in the box.

## The heuristic I had was wrong, and the computation caught it

I expected the extremal `R` to be the **sparsest** — the two Newton vertices, `R = 1 + u^d`,
i.e. the two-monomial Laurent polynomial `z^{−M} + z^N`. That gives first-nonzero
`(M+N)/gcd(M,N)`, which equals `M+N` **only when `gcd(M,N) = 1`.** When `gcd > 1` it is small
(`m = 2` at `(3,3)`, since `6/3 = 2`), and interior monomials **raise** the first-nonzero index
— in `54/3875` sampled cases, directly refuting "interior only pulls down."

> **The extremal structure is gcd-graded.** The two-vertex picture is the coprime special case.
> For `gcd(M,N) > 1` the delaying `R` is not the sparsest, and where its true first-nonzero
> sits (whether the bound `M+N` is even tight) is open — at `(3,3)` the observed max was `4`.

## Honest scope

- The polygon structure and the residue duality are **proved**; everything about the effective
  bound is **verified in a bounded box**, and the bound itself is the **open** ESV conjecture.
- **My "sparsest is extremal" heuristic is refuted** by my own run; I am recording it as
  refuted rather than dropping it, because the gcd-grading it exposes is the real content.
- `(3,3)` shows tightness is **not universal** and the extremal is not the two-monomial when
  `gcd > 1`. Whether `M+N` is tight for non-coprime `(M,N)` at all is unresolved here.
- This does not advance the *proof* of the ESV bound; it gives the Newton-polygon picture and
  corrects the extremal heuristic. Nothing here bears on GMC(2) (THM-1645) or GMC(≥3).

*Artifacts:* `04-computation/newton_polygon_branches_macmini_S144.py` (+out).
