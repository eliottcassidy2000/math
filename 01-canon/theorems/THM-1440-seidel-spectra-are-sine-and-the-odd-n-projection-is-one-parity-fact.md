---
id: THM-1440
title: "THE TWO 'ODD n ONLY' FACTS ARE ONE PARITY FACT, AND IT IS THE PARITY OF sin AND cos. (A) S = A − Aᵀ is skew, so p(−x) = (−1)ⁿ p(x): the characteristic polynomial is an EVEN function at even n (like cos) and an ODD function at odd n (like sin), which is exactly why odd n forces a ZERO eigenvalue — sin(0) = 0 and cos(0) ≠ 0. Rank is always even. Verified n = 3..7. (B) THE SKEW-SEIDEL SPECTRUM IS A COMPLETE INVARIANT OF THE SWITCHING CLASS FOR n ≤ 6 AND FIRST FAILS AT n = 7: distinct spectra 1, 2, 2, 6, 11 against A049313 = 1, 2, 2, 6, 12, i.e. exactly ONE cospectral pair of switching classes at n = 7. Switching acts by S ↦ DSD, a similarity, so unlike H / min-FAS / cyclic-triangle count the spectrum survives it. (C) The repo's odd-n-only projection T_cycle = (I+L(K_n))T mod 2 is explained twice over: the image degree at w is (n−1)·d_w mod 2, vanishing identically iff n is ODD; and the BICYCLE space Cut ∩ Cycle of K_n is 0 iff n is odd (dimension n−2 when n is even, verified n = 3..9), which is exactly when Cut ⊕ Cycle is a direct sum. (D) Circulant tournaments make the sine LITERAL: μ_j = 2i·Σ_{k∈C} sin(2πjk/n), a pure sine sum, while the cosine part is PINNED at Re(z_j) = −1/2 for every j ≠ 0 — a critical line. Paley gives {0} ∪ {±i√q} exactly. (E) Cauchy interlacing under vertex deletion holds on all 16/60/336/3192 deletions at n = 4..7, reproducing the interlacing of the zeros of sin and cos."
status: >
  (A) PROVED (three lines from Sᵀ = −S) and VERIFIED-EXACT n = 3..7 over all iso classes.
  (B) Switching-invariance PROVED (D = D⁻¹ makes S ↦ DSD a similarity) and verified
  numerically on every class at every n; the COMPLETENESS boundary is VERIFIED-EXACT by
  full enumeration against A049313 — complete at n ≤ 6, exactly one cospectral pair at n = 7.
  (C) Both mechanisms PROVED (the degree computation is three lines; the bicycle argument is
  two) and VERIFIED n = 3..9 / n ≤ 8.
  (D) PROVED (circulant diagonalization + C ⊔ −C = Z_n∖{0}) and verified to 1e-14 over all
  circulant tournaments at n = 5,7,9,11; Paley verified q = 3,7,11,19,23.
  (E) PROVED (Cauchy interlacing for the Hermitian iS) and verified on every deletion n = 4..7.
  The sin/cos reading is a FRAMING built on the repo's existing odd/even-function axis, not a
  new theorem — see Honest scope.
source: mac-mini-2026-07-20-S130 (owner: "think seidel spectra and even graph bijections that
  are odd n only, and how these two interlaced can be seen in sin and cosine")
depends_on:
  - THM-474   # the Gauge Theorem: tilings ARE switching classes of tournaments
related:
  - THM-499   # H is finer than the ADJACENCY spectrum, breaks at n=6
  - THM-500   # odd-cycle count non-spectral from n=7 -- SAME boundary, DIFFERENT matrix
  - THM-1415  # switching is the canonical star quotient (+ the S129 A049313 addendum)
  - THM-1430  # graph switching lands exactly on E_n (the undirected companion)
  - THM-1420  # no F_2-linear invariants -- the sharp contrast drawn in §F
  - 07-reflections/the-determinant-lens-sgn-vs-chi-and-the-three-geometries.md
script: 04-computation/seidel_spectra_sine_cosine_macmini_S130.py,
        04-computation/seidel_cospectral_pair_n7_macmini_S130.py (+ .outs)
---

# THM-1440 — Seidel spectra are sine, and the odd-`n` projection is the same parity fact

**One line.** `sin` is odd and vanishes at `0`; `cos` is even and does not. A tournament's
skew-Seidel characteristic polynomial is sin-like at odd `n` and cos-like at even `n` — and
*that single parity* is simultaneously why odd `n` forces a zero eigenvalue and why the
even-graph projection exists only at odd `n`.

## (A) The characteristic polynomial has the parity of sin/cos

Let `S_ij = +1` if `i→j`, `−1` if `j→i`, so `Sᵀ = −S`. Then

`p(−x) = det(−xI − S) = (−1)ⁿ det(xI + S) = (−1)ⁿ det(xI − Sᵀ) = (−1)ⁿ p(x).`

> **`p` is an EVEN polynomial at even `n` and an ODD polynomial at odd `n`.**
> Hence at odd `n` it has the factor `x`: a **forced zero eigenvalue**.

Concretely `p(x) = ∏(x² + λ_k²)` at even `n` and `x·∏(x² + λ_k²)` at odd `n` — the same
"extra linear factor" that distinguishes `sin z = z∏(1 − z²/k²π²)` from
`cos z = ∏(1 − z²/(k−½)²π²)`. Verified over all iso classes at `n = 3..7`: every eigenvalue
purely imaginary, the parity exact, `0` present at every odd `n`, and the rank always even
(skew matrices have even rank, so odd `n` forces corank `≥ 1`).

## (B) The spectrum is switching-invariant, and complete exactly up to `n = 6`

Switching at `X` sends `S ↦ DSD` with `D = diag(±1)`. Since `D = D⁻¹` this is a **similarity**,
so the spectrum cannot move. This is a sharp contrast with the S129 findings: `H`, min-FAS and
the cyclic-triangle count are **none of them** switching-invariant, but the spectrum is.

| `n` | iso classes | distinct skew-Seidel spectra | switching classes (A049313) | complete? |
|---|---|---|---|---|
| 3 | 2 | 1 | 1 | yes |
| 4 | 4 | 2 | 2 | yes |
| 5 | 12 | 2 | 2 | yes |
| 6 | 56 | 6 | 6 | yes |
| 7 | 456 | **11** | **12** | **NO — one cospectral pair** |

> **The skew-Seidel spectrum determines the switching class (equivalently the oriented
> two-graph) for `n ≤ 6`, and first fails at `n = 7` with exactly one cospectral pair.**

### The pair, exhibited

Independent recomputation confirms `456` iso classes, **`12` switching classes** (matching
`A049313(7)` — an independent verification of the OEIS term) and **`11` spectra**. The two
cospectral switching classes each contain **exactly 64 iso classes**, and share

> `spec(iS) = {0, ±√7, ±√(7−4√2), ±√(7+4√2)}`,
> i.e. `p(x) = x(x² + 7)(x⁴ + 14x² + 17)`,

verified to the integer: coefficients `[1, 0, 21, 0, 115, 0, 119, 0]`. Note the coefficients
are supported on degrees `7, 5, 3, 1` — **all odd**, which is (A) at work. Two remarks:

- The splitting field is `ℚ(√2)`: the quartic `x⁴ + 14x² + 17` is irreducible, and `√2` is
  the repo's hypotenuse/leg ratio. `Σλ² = 42 = n(n−1)` as forced by `tr S² = −n(n−1)`.
- `√7` appears as an eigenvalue, which is the Paley-heptagon value — but the pair is **not**
  Paley: the Paley tournament on 7 vertices has the degenerate spectrum `{0, ±i√7}` with
  multiplicity 3, whereas here all three magnitudes are distinct. Recorded as an observation,
  **not** a claimed connection.
- Representatives have score sequences `[0,1,3,4,4,4,5]` and `[0,2,2,4,4,4,5]` with `4` and
  `5` cyclic triangles — but neither statistic separates the classes, since (per S129) neither
  is switching-invariant. The `H`-value sets overlap heavily (both contain
  `29, 33, 73, 89, 91, 95, 109, 115`), for the same reason.

**A boundary worth flagging.** THM-499/500 place the *adjacency*-spectrum boundaries at
`n = 6` (`H` non-spectral) and `n = 7` (odd-cycle count non-spectral). This is a **different
matrix** — `S = 2A + I − J`, and the `S`-spectrum is switching-invariant while the
`A`-spectrum is not — yet the completeness boundary lands at the *same* `n = 7`. Whether that
is one mechanism or two is **open**; nothing here establishes a connection.

## (C) Why the even-graph projection is odd-`n` only — two proofs

**C1, directly.** For the repo's map `T_cycle = (I + L(K_n))T mod 2`, the image at edge
`e = {u,v}` is `T_e + d_u + d_v` where `d_v = Σ_{f at v} T_f`. Summing over edges at `w`:

`Σ_{e at w} [(I+L)T]_e = (n−1)·d_w + Σ_x d_x = (n−1)·d_w (mod 2)`,

since `Σ_x d_x = Σ_f 2T_f = 0`. **So the image is an even graph iff `n` is odd.** Verified
empirically at `n = 3..8`.

**C2, structurally.** A cut `δ(S)` has all degrees even iff `|S|` and `|S^c|` are both even —
impossible when `n` is odd. So the **bicycle space** `Cut ∩ Cycle` of `K_n` is `0` exactly
when `n` is odd, which is exactly when `Cut ⊕ Cycle` is a genuine direct sum and the
projection onto the cycle space is well defined:

| `n` | 3 | 4 | 5 | 6 | 7 | 8 | 9 |
|---|---|---|---|---|---|---|---|
| bicycle dim | 0 | 2 | 0 | 4 | 0 | 6 | 0 |
| predicted (`0` odd, `n−2` even) | 0 | 2 | 0 | 4 | 0 | 6 | 0 |

This also completes THM-1405, which noted the bicycle space as a caveat to its "gauge bits /
holonomy bits" splitting: **the splitting is canonical exactly at odd `n`.**

## (D) Circulants: the sine is literal, and the cosine is pinned

For a circulant tournament with connection set `C` (so `C ⊔ (−C) = Z_n∖{0}`),

> `μ_j = Σ_{k∈C}(ω^{jk} − ω^{−jk}) = **2i·Σ_{k∈C} sin(2πjk/n)**` — a pure **sine** sum.

And because `C ⊔ (−C)` exhausts the nonzero residues, `z_j := Σ_{k∈C} ω^{jk}` satisfies
`z_j + z̄_j = −1`, i.e.

> **`Re(z_j) = −1/2` for every `j ≠ 0`** — the cosine part is pinned to a critical line, and
> *all* the content of a circulant tournament sits in its sine part.

Verified to `< 2·10⁻¹⁴` over all `4/8/16/32` circulant tournaments at `n = 5,7,9,11`. The
Paley tournaments are the degenerate extreme: spectrum exactly `{0} ∪ {±i√q}`, verified
`q = 3,7,11,19,23`.

## (E) Interlacing is the interlacing of the zeros of sin and cos

`iS` is Hermitian, so Cauchy interlacing holds for principal submatrices — i.e. under **vertex
deletion**. Verified on all `16 / 60 / 336 / 3192` deletions at `n = 4..7` (strict in
`8/16`, `44/60`, `220/336`, `2342/3192`).

The spectra alternate exactly as sin/cos zeros do:

| | spectrum | zeros of |
|---|---|---|
| odd `n` | symmetric about `0`, **contains** `0` | `sin`: `0, ±π, ±2π, …` |
| even `n` | symmetric about `0`, **omits** `0` | `cos`: `±π/2, ±3π/2, …` |

Deleting a vertex flips the parity, and the two spectra must interlace — which is precisely
the classical statement that the zeros of `sin` and `cos` interlace.

## (F) The contrast with THM-1420

S129 proved there are **no `F₂`-linear isomorphism invariants at all**. This file exhibits a
genuine `ℝ`-spectral invariant that is not merely isomorphism-invariant but *switching*-
invariant. There is no contradiction — and the pair draws a clean three-tier picture:

| invariant | iso-invariant | switching-invariant |
|---|---|---|
| any `F₂`-linear functional | **no** (THM-1420: none exist) | — |
| `H`, min-FAS, cyclic triangles | yes | **no** (S129) |
| skew-Seidel spectrum | yes | **yes** |

So the tiers are strictly nested, and the spectrum is the coarsest of the three — which is
also why it can be complete for switching classes at small `n` while `H` is not even defined
on them.

## Honest scope

- **The sin/cos reading is a framing, not a new theorem.** The repo's determinant-lens
  reflection already records the load-bearing axis ("graph = even (symmetric); tournament =
  odd (skew)"). What is added is that the *same* parity controls (A), (C) and (E), and that
  (D) makes the sine literal rather than metaphorical.
- (A) and (D) are classical facts about skew matrices, circulants and Paley tournaments,
  re-derived and verified here; they are not claimed as new.
- **(B) is the new content**, together with (C)'s bicycle identification. The `n = 7`
  cospectral pair is exhibited in the companion script.
- The `n = 7` coincidence with THM-499/500's adjacency boundary is **flagged, not explained**.
- Completeness is checked only to `n = 7`; `A049313(8) = 79` is the obvious next test and is
  not done here.
- No claim is made that the "critical line" `Re(z_j) = −1/2` in (D) has any connection to any
  other object bearing that name. It is an elementary consequence of `C ⊔ (−C) = Z_n∖{0}`.

*Artifacts:* `04-computation/seidel_spectra_sine_cosine_macmini_S130.py`,
`seidel_cospectral_pair_n7_macmini_S130.py` (+outs).
*Credits:* THM-474 (tilings are switching classes) for the framework; Babai–Cameron for the
switching-class/oriented-two-graph dictionary; THM-1430 for the undirected companion; the
determinant-lens reflection for the odd/even-function axis this builds on.
