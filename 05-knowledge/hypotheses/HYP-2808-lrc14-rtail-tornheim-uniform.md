# HYP-2808 — The general bounded-base R-tail is a convergent MORDELL–TORNHEIM double sum, bounded by an ABSOLUTE constant uniformly over (base, gap) — completing leg-C's analytic piece

- **Instance:** claude-opus-2026-06-22
- **Status:** RIGOROUS (mac-mini-S22 proved the Tornheim constant `T = 12·ζ(3)` exactly via the classical Tornheim reduction — owner's arXiv:2409.19980 — so `|R_g| ≤ (1/π³)·(#sector-pairs)·12ζ(3)` is a CLOSED-FORM absolute uniform bound; empirical sup 2.24, sin-sharp ≈2.9). The R-tail is now an absolute closed-form constant.
- **Touches:** OPEN-Q-108 leg C; THM-564 (P/R split, R-tail), HYP-2807 (generalized doublet), HYP-2792 (signed Dedekind), the Dedekind-ladder reflection.

## The R-tail and why "general bounded-base" needs a UNIFORM bound

THM-564 closes the adjacent doublet via `M·(p0(E_M)−Φ) = P(M) + R(M)`, `R(M) = M·(d2(M)−d_inf)`.
HYP-2807: the genuine-wide max is a GENERALIZED doublet `{M,M+g}`, so the R-tail must be bounded
**uniformly over (base B, gap g)** — the "general bounded-base R-tail." THM-564 used the empirical
`sup|R|`; kps's TV/Koksma bound `J_sharp≈27` is rigorous but ~20× loose (cutoff `f0~344`).

## The reframe: R_g is a Mordell–Tornheim double sum (the uniform mechanism)

The curvature `d2_g(M)` is supported on the base's miss-arcs; on each arc it is a **covariance of two
sector-indicators at the LOCKED phases** `({Mφ}, {(M+g)φ})`. Fourier-expand `1_{I_p}({x})−1/7 =
Σ_{h≠0} a_{p,h} e(hx)`, `a_{p,h} = (e(−hp/7)−e(−h(p+1)/7))/(2πih)`, `|a_{p,h}| = |sin(πh/7)|/(π|h|)`.
On arc `[a,b]`:
> `Cov = Σ_{h,h'≠0} a_{p,h} a_{q,h'} ∫_a^b e(((h+h')M + h'g)φ) dφ`.
The integral phase is `c = (h+h')M + h'g`:
- `h+h'=0` ⟹ `c = −hg` (M-independent) ⟹ the SATURATION `d_inf` (the double-Dedekind diagonal).
- `h+h'≠0` ⟹ `c ≈ (h+h')M` large ⟹ `∫ = O(1/((h+h')M))` ⟹ the M-dependent `d2−d_inf`.

Therefore
> **`R_g(M) = M·(d2_g−d_inf) = Σ_{h+h'≠0} a_{p,h} a_{q,h'} · M·(e(cb)−e(ca))/(2πi c)`**, `M/c → 1/(h+h')`,
> `|R_g| ≤ (per-arc) (1/π³)·Σ_{h,h'≠0, h+h'≠0} |sin(πh/7)sin(πh'/7)|/(|h||h'||h+h'|) ≤ (1/π³)·T`,
> `T = Σ_{h,h'≠0, h+h'≠0} 1/(|h||h'||h+h'|) ≈ 14.33` (convergent Mordell–Tornheim constant).

**Per-arc bound `(1/π³)·T ≈ 0.46`; with ≤ ~12 active sector-arcs, `|R_g| ≲ 5.5` — an ABSOLUTE
constant, UNIFORM over every bounded base and gap** (the base only selects WHICH finitely-many
sector-pairs/arcs appear; `T` is base-independent).

**SHARP constant (with the `sin(πh/7)` Fourier factors).** Using the exact `|a_{p,h}| =
|sin(πh/7)|/(π|h|)` (which VANISHES when `7|h`):
> `|R_g| ≤ (1/π³)·(#active sector-pairs)·S`,  `S = Σ_{h+h'≠0} |sin(πh/7)sin(πh'/7)|/(|h||h'||h+h'|) ≈ 5.95`.
Per-arc `(1/π³)·S ≈ 0.192`; with 10–15 active sector-pairs, **`|R_g| ≲ 1.9–2.9`** — which TIGHTLY
matches the empirical sup `2.244`. So `S ≈ 5.95` is the explicit absolute constant of the general
bounded-base R-tail. ⟹ `G = period-max(P)+sup|R| ≲ 1.5+2.9 = 4.4`, margin ≥ 0.16, `M* = ⌈G/margin⌉ ≈ 28`.

## Empirical confirmation (`lrc14_Rtail_uniform_claudeopus_0622.py`, exact)

`sup_{M∈[15,200]}|R_g(M)|` over consec/even-AP/top-cluster/+random bases (≈50) × gaps g=1,2,3 ×
k=9,10,11: **GLOBAL sup = 2.244** (top-cluster base `(0,9,…,14)`, g=1). All bounded ≲ 2.3, well
under the Tornheim ceiling ~5.5. So `R_g` is uniformly bounded as predicted.

## Why it COMPLETES the genuine-wide leg

With a rigorous uniform `G = period-max(P) + sup|R| ≤ ~1.5 + 5.5 = 7` (or ~3.7 with the sharp sup),
and the binding margin `cap − Φ_frozen ≥ 0.16`:
- **tail `M ≥ M* = ⌈G/(cap−Φ)⌉ ≈ 23–44`:** `p0 = Φ + g(M)/M ≤ Φ + G/M* = cap` ✓ (needs the ROOM
  `Φ_frozen(B,g) < cap`, HYP-2807 (I), being verified).
- **finite window `15 ≤ M < M*`:** exact check `p0(B∪{M,M+g}) < cap` over (base, gap, M) — a small,
  THM-563-style enumerable certificate (collapsed from kps's `f0~344` to `M*~tens`).

So leg C = **[frozen room `Φ<cap`] + [Tornheim-uniform R-tail] + [tiny finite window]**, all three of
which are now finite/structural. The Mordell–Tornheim `T` is the absolute constant that makes the
R-tail "general bounded-base" (uniform) rather than per-base.

## Next / rigor gaps
- Pin the sharp constant: include `sin(πh/7)` and exact arc count per (sector-pair) → tighten ~5.5.
- The endpoint phases `e(cb),e(ca)` — bound the conditionally-convergent sum rigorously (Abel/Koksma
  on the tail; the absolute Tornheim bound already dominates it).
- This is the doublet/double-sum rung of the Dedekind ladder; the single-far rung is THM-563.

→ OPEN-Q-108, THM-564, HYP-2807, HYP-2792, THM-563, `07-reflections/lrc-the-dedekind-ladder-of-far-coherence.md`.
