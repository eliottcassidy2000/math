---
id: THM-1920
title: "THE SPECTRAL INSERTION-RESPONSE: how a = vertex-insertion acts on the SKEW characteristic polynomial for GENERAL tournaments, with kps's transitive a/b Chebyshev-Pell frame (THM-1880) as the special case. The combinatorial a of THM-1900 (add a vertex u beating exactly the subset P) acts on char_S by an exact bordered recursion, verified for ALL tournaments n<=4 and all 2^m patterns: char_S(T+u_P)(x) = x*char_S(T)(x) + B_P(x), where B_P(x) = s_P^T adj(xI-S) s_P (s_P the +-1 sign vector of P). THREE STRUCTURAL LAWS. (1) THE CORRECTION SPLITS: B_P = char_S'(x) + (signed off-diagonal cofactor sum), because the DIAGONAL cofactor sum is exactly char_S'(x) = sum_i char_{S-i}(x) (the derivative = the vertex-deletion sum) -- so the pattern-INDEPENDENT core of every insertion is the derivative of the characteristic polynomial, and P enters only through the quadratic form s_P^T (off-diagonal) s_P. (2) SOURCE = SINK: the all-+1 and all--1 patterns give the SAME B (s_i s_j identical), so char_S(T+source) = char_S(T+sink); hence the COMPLEMENT b = T -> T^op is char_S-INVARIANT, which is exactly why the b-quotient merged metagraph G_n/Z2 is the natural SPECTRAL object (the skew spectrum is complement-symmetric). (3) INTERLACING: the m+1 eigenvalues of T+u_P interlace the m eigenvalues of T on the imaginary axis (Cauchy interlacing on the Hermitian iS), verified 2000 random n=5 -- inserting a vertex threads one new imaginary eigenvalue between each consecutive pair. THE TRANSITIVE TOWER RECOVERS kps THM-1880: building the transitive tournament by n all-source insertions (a^n) gives char_S = E_n = ((x+1)^n+(x-1)^n)/2 via the source-recursion E_n = x*E_{n-1} + B, so kps's ALGEBRAIC a(x)=x+1 is this COMBINATORIAL a (THM-1900's insertion) restricted to the transitive tower. TRIANGULAR NUMBERS ARE UNIVERSAL: the x^{n-2} coefficient of char_S is C(n,2)=T_{n-1}=#arcs for EVERY tournament (= (1/2)(-tr S^2)), extending kps THM-1880(5) beyond the transitive case; under a, #arcs grows T_{n-1} -> T_n. So a=insertion (combinatorial, spectrum interlaces) and b=complement (spectrum-invariant) ARE the owner's a,b at the OBJECT level, dual to kps's polynomial-level a(x)=x+1, b(x)=x/2"
status: PROVED (bordered recursion exact by Schur complement, verified all (T,P) n<=4; diagonal=char_S' is Jacobi's formula; source=sink exact; interlacing = Cauchy on iS, numerically confirmed 2000 random n=5; transitive tower = ((x+1)^n+(x-1)^n)/2 with the source-recursion, exact n<=6). Complementary to kps THM-1880 (transitive only); this is the a-action on ALL tournaments' spectra + the combinatorial<->algebraic a bridge.
author: opus-2026-07-20-S440
extends: THM-1900 (insertion-response calculus: H,c3 under insertion -> now char_S/the spectrum)
recovers: THM-1880 (kps: transitive a/b Chebyshev-Pell, = the source-insertion special case)
depends_on: [THM-1900 (combinatorial a = vertex-insertion), THM-1880 (kps a/b frame), THM-1875-transitive-skew-char (kps: ((x+1)^n+(x-1)^n)/2), THM-1810 (transitive = GIT nullcone), THM-474 (skew d(T)=det(I+S)/2^{n-1}), THM-012b-insertion-decomposition (era-0 H(T)-H(T-v) parent), THM-1440-seidel-spectra-are-sine (Cauchy interlacing under DELETION -- this is the insertion dual), THM-1560-the-halving-dictionary (priority owner of b: x->(1+x)/2=(J-I+S)/2), THM-1830-unstable-non-transitive (char_S MULTIPLICATIVE under order-join = the down-set B_P case)]
cite_by_filename: true  # duplicate THM numbers exist (1875, 1810, 1830, 1440, ...); cite files not numbers
---

# THM-1920 — The spectral insertion-response

Owner (S440): think of triangular numbers and tournaments as `a(x)=x+1`, `b(x)=x/2` composed
recursively; think functionally/trigonometrically. kps THM-1875/1880 built this frame for the
**transitive** tournament (`E_n = ((x+1)^n+(x-1)^n)/2`, Pell `E_n²−O_n²=(x²−1)^n`, cotangent
spectra). This is the **complementary object-level law**: the combinatorial `a` of [THM-1900](THM-1900-the-insertion-response-calculus.md)
(insert a vertex `u` beating exactly `P`) acting on the skew characteristic polynomial of **every**
tournament — recovering kps's algebraic `a(x)=x+1` as the transitive special case.

## The bordered recursion (exact, all tournaments)

Skew matrix `S` (`S_{ij}=+1` if `i→j`). Insert `u` beating `P` (sign vector `s_P`, `s_j=+1 ⟺ j∈P`).
By Schur complement on the bordered skew `[[0, s^T],[−s, S]]`:

> **`char_S(T + u_P)(x) = x·char_S(T)(x) + B_P(x)`, where `B_P = s_P^T \operatorname{adj}(xI−S)\, s_P`.**

Verified exactly for all tournaments `n≤4` and all `2^m` patterns.

## Three structural laws

**(1) The correction splits: `B_P = char_S′ + signed-off-diagonal(P)`.** The diagonal cofactor
sum is `Σ_i \operatorname{adj}(xI−S)_{ii} = Σ_i char_{S−i}(x) = char_S′(x)` (Jacobi's formula: the
derivative is the vertex-deletion sum). Since `s_i² = 1` always, the diagonal part of `B_P` is
`char_S′` for **every** `P` — the **pattern-independent core of any insertion is the derivative of
the characteristic polynomial**. `P` enters only through the off-diagonal quadratic form
`Σ_{i≠j} s_i s_j \operatorname{adj}(xI−S)_{ij}`. (Verified `74/74`.)

**(2) Source = sink ⟹ the complement `b` is char_S-invariant.** All-`+1` and all-`−1` give the
same `B` (`s_i s_j` unchanged), so `char_S(T+\text{source}) = char_S(T+\text{sink})`. Consequently
`T ↦ T^{op}` (reverse all arcs = the complement `b`) preserves `char_S`: **the skew spectrum is
complement-symmetric**, which is exactly why the `b`-quotient **merged metagraph `G_n/ℤ₂`** is the
natural *spectral* object (S211+ threads), and why `b(x)=x/2` = `(A000568+SC)/2` is the right count.

**(3) Interlacing.** The `m+1` eigenvalues of `T+u_P` **interlace** the `m` eigenvalues of `T` on
the imaginary axis (Cauchy interlacing applied to the Hermitian `iS`) — one new imaginary
eigenvalue threads between each consecutive old pair. Verified on 2000 random `n=5` insertions.
This is the spectral face of the insertion-response: `a` cannot move an eigenvalue past its
neighbours. It is the **insertion dual** of `THM-1440-seidel-spectra-are-sine` (which proves Cauchy
interlacing under vertex **deletion**) — `a` and `a⁻¹` interlace in opposite directions.

## The transitive tower = kps THM-1880

Building the transitive tournament by `n` all-**source** insertions (the combinatorial `a^n`)
gives, step by step, `char_S = E_n = ((x+1)^n+(x−1)^n)/2` via the source-recursion
`E_n = x·E_{n−1} + B` (verified `n≤6`). So kps's **algebraic** `a(x)=x+1` is precisely **this
combinatorial `a`** (THM-1900's insertion) restricted to the transitive tower, and `E_n±O_n =
(x±1)^n` is the eigenpolynomial of that tower. The two `a`'s — insert-a-vertex (object) and
multiply-by-`(x+1)` (polynomial) — are one functor.

## Triangular numbers, universally

For **every** tournament the `x^{n−2}` coefficient of `char_S` is
`C(n,2) = T_{n−1} = #\text{arcs} = \tfrac12(−\operatorname{tr}S²)` (extends kps THM-1880(5), which
noted it for the transitive case, to all tournaments — it is `e₂` of the skew spectrum, always the
arc-count). Under `a` (insertion) the arc-count grows `T_{n−1} → T_n` (add `n` arcs), the
combinatorial shadow of `char_S ↦ x·char_S + B`.

## The dictionary (object level ↔ kps's polynomial level)

| owner's generator | object level (this THM) | polynomial level (kps THM-1880) |
|---|---|---|
| `a(x)=x+1` | insert a vertex `n→n+1` (THM-1900); `char_S ↦ x·char_S + B_P`; spectrum interlaces | multiply the eigenpoly by `x+1`: `E_n±O_n=(x±1)^n` |
| `b(x)=x/2` | complement `T↦T^{op}` — **char_S-invariant** — quotient to `G_n/ℤ₂` | the symmetriser `b(a^n±ā^n)` and the `½` half-dictionary (THM-1555) |
| triangular `T_{n−1}` | `#arcs`, universally the `x^{n−2}` coeff of `char_S` | the `C(n,2j)` binomial coefficients of `E_n` |

## Open

1. **The off-diagonal quadratic form (Q1 — partially resolved).** For `P` a **down-set** (THM-1900
   H-neutral), the insertion is an **order-join at an SCC boundary**, so `char_S` becomes
   **multiplicative** — `char_S(T+u_P) = char_S` of the two order-join factors (this is
   `THM-1830-unstable-non-transitive`, the `B_P`-structured case). So the down-set signature is
   *char_S-factorisation*, not H-neutrality's combinatorial form. The tempting guess "down-set
   insertion pins a zero eigenvalue" is **REFUTED** (`spectral_downset_probe_opus_S440.py`): a
   `0` eigenvalue appears iff `n` is **odd** (skew-parity), for down-sets and non-down-sets alike
   (`0/28` at `n=4`, `216/216` and `808/808` at `n=5`) — it is nothing to do with `P`. The general
   off-diagonal form (non-down-set `P`) remains the open cell.
2. **The regular pole.** The transitive (nullcone) tournament is the Chebyshev-Pell object; its dual,
   the **regular/Paley** tournament, has skew spectrum `±i√p` (THM-1810) — a *single* repeated
   eigenvalue, **not** cotangent. So the a/b Chebyshev frame is special to the nullcone vertex; what
   is the regular pole's functional frame? **(Partially answered concurrently by kps-S128c139:** the
   scalar `var(λ²)` of the squared skew spectrum is a one-number GIT-instability measure — **maximal
   at the transitive vertex, exactly `0` at Paley** (all `λ²=p`) — and the deformation family
   `b((x+c)^n+(x−c)^n)` interpolates transitive `↔` Paley. So the two poles are the `var(λ²)`-extreme
   and `var(λ²)=0` ends; under insertion `a`, how does `var(λ²)` move? — the natural next cell of the
   spectral insertion-response matrix.)**

## Verification

`04-computation/spectral_insertion_response_opus_S440.py` (+ `.out`) — the bordered recursion for
all `(T,P)` `n≤4`; `B_P` diagonal `= char_S′` (`74/74`); source `=` sink; the transitive tower
`= ((x+1)^n+(x−1)^n)/2`; interlacing on 2000 random `n=5`.
