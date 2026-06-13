# THM-485 — The two temperatures of the hard-core transfer operator: Zeckendorf, the tournament H, and Viswanath

**Status:** PROVED (deterministic parts) + VERIFIED (Lyapunov/Monte-Carlo) — claudebox-2026-06-11-S5.
**Provenance:** user dispatch (Viswanath's constant + Zeckendorf + the transfer-spectrum thread).
**Builds on:** THM-326 (H = I(Ω,2)), THM-337 (base-path order-3 recurrence, root 3.38298),
the S720/721 temperature ladder, fibonacci_tiling_galois.out (F_n = I(P_{n−2},1), the
Fibonacci-on-paths identity), tribonacci_base_tournament.out (Zeckendorf in tribonacci base).
**External (cited):** Viswanath 2000 (Math. Comp. 69, the constant 1.13198824…, Lyapunov of
random [[0,1],[±1,1]] products via the Stern–Brocot stationary measure); Embree–Trefethen 1999
(critical β* ≈ 0.70258 for f_n = f_{n−1} ± β f_{n−2}); Lekkerkerker (Zeckendorf summand average).

## Setup

The **hard-core / independence transfer operator** of the path is
I(P_n, x) = I(P_{n−1}, x) + x·I(P_{n−2}, x) (independence polynomial of P_n at activity/fugacity
x). The admissible configurations are 0/1 strings with no two adjacent 1s — the **golden-mean
shift of finite type**, which is exactly the **Zeckendorf numeration** rule (no two consecutive
Fibonacci summands). The same operator is the repo's tournament partition function in the
path case: H = I(Ω, 2) (THM-326). It carries TWO independent temperature axes.

## Statement

**(1) Deterministic fugacity temperature (annealed).** I(P_n, x) grows like
λ(x) = (1 + √(1+4x))/2 — the dominant eigenvalue of [[1, x],[1, 0]]. Three landmarks on one
operator (verified exact, I_40/I_39):
- x = 1: λ = φ = 1.61803… — **Zeckendorf / golden-mean SFT**, topological entropy log φ
  (admissible length-n strings = F_{n+2}; Lekkerkerker average summand density 1/(φ²+1)).
- x = 2: λ = 2 — the repo's **H = I(Ω,2)**, Jacobsthal J_n = (2^{n+1}+(−1)^n)/3 = 1,3,5,11,21,…
- general x: the S720 "fugacity temperature." x = 1 (Zeckendorf) and x = 2 (tournament H) are
  two readings of ONE transfer matrix at two activities.

**(2) Quenched sign-disorder temperature.** Randomize the recurrence SIGNS:
f_n = ±f_{n−1} ± x·f_{n−2} (i.i.d. fair signs). The growth rate drops from the eigenvalue λ(x)
to the top **Lyapunov exponent** of the random matrix cocycle. At x = 1 this is exactly
**Viswanath's constant 1.13198824…** (verified: Monte-Carlo 1.1321, 4 i.i.d. ×1.5M-step runs).
Disorder strictly LOWERS growth: φ = 1.618 (ordered) → 1.132 (disordered). NEW — the disordered
("Viswanath") constants of the repo's transfer families (verified Lyapunov, leading
geometry/coefficient frozen, signs randomized):

| family | deterministic root | disordered (Lyapunov) |
|--------|--------------------|------------------------|
| Fibonacci (Zeckendorf) | φ = 1.61803 | **1.1320 = Viswanath** |
| tribonacci (up-staircase) | 1.83929 | **1.2228** |
| base-path (THM-337) | 3.38298 | **2.9786** |

**(3) Disorder-induced phase transition (new).** The DETERMINISTIC operator grows for every
x > 0 (λ(x) > 1). The QUENCHED operator does NOT: its Lyapunov exponent crosses zero at the
**Embree–Trefethen activity β* ≈ 0.70258** — the random hard-core partition function DECAYS for
x < β* and grows for x > β* (verified: growth 0.92 at x=0.5, ≈1.000 at x=0.70258, 1.13 at x=1).
So disorder CREATES a phase transition the ordered system lacks; Zeckendorf (x=1) and the
tournament H (x=2) both sit in the grown phase. β* = the critical activity of the
**quenched-random hard-core lattice gas**.

## Proof / verification

(1) is the Perron root of [[1,x],[1,0]], char t² = t + x ⟹ t = (1+√(1+4x))/2; exact ratio
check I_40/I_39 matches to 1e−6 at x = 1,2,3. The Zeckendorf=golden-mean-SFT identity is the
no-11 count F_{n+2} (standard; = independent sets of P_n, the repo's fibonacci_tiling_galois).
(2) is the Furstenberg–Kesten Lyapunov exponent of the i.i.d. matrix cocycle; Viswanath proved
the x=1 value via the Stern–Brocot stationary measure. Monte-Carlo with renormalization
reproduces 1.1320 (Fibonacci/Viswanath, agreeing to 4 digits) and gives the tribonacci/base-path
constants. (3) is Embree–Trefethen's β* (decay↔growth threshold); the new content is the
identification of x with the hard-core fugacity and of the grown phase with Zeckendorf + the
tournament H. Script transfer_temperatures_viswanath_cbx5.py (+ .out).

## Remarks

1. **Two orthogonal temperatures.** The S720 ladder Λ(a) = x³−3x²+ax−1 varies a COEFFICIENT
   deterministically (the annealed/fugacity axis, frozen leading 3 = geometry). Randomizing
   SIGNS is a SECOND, quenched axis. The eigenvalue (annealed) ≥ Lyapunov (quenched) always;
   the gap is the disorder/"glass" cost. This is the transfer-operator face of the S637 glass
   transition: quenched disorder = the spin-glass/Lyapunov regime, deterministic = the crystal.
2. **Viswanath = the disordered golden ratio = the disordered Zeckendorf entropy.** The golden
   ratio is the entropy of the Zeckendorf SFT; Viswanath's constant is what that entropy becomes
   under quenched sign-disorder. The two famous Fibonacci constants (φ, Viswanath) are the
   ordered and disordered endpoints of ONE transfer operator — the same operator whose x=2
   reading is the tournament independence polynomial H.
3. **Open (toward HYP-2416):** is there a closed form for the base-path disordered constant
   ≈ 2.9786? Viswanath's own constant has no known closed form (only the Stern–Brocot stationary
   measure), so a closed form here is unlikely; the honest targets are rigorous interval-arithmetic
   digits and the leading-order disorder correction (3 − λ_disordered) as a function of the frozen
   geometry. (No claimed match to the 0.014/Cl₂(π/3) surplus — the magnitudes differ.)
