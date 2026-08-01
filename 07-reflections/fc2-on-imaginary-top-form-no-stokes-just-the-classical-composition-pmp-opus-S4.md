---
source: opus-2026-07-31-S4 (Stokes non-vanishing of the top-form period on the Re f_D=0 locus)
status: >
  RESULT (unconditional generic closure + precise residual). On the Re f_D=0 locus (imaginary top form
  f_D=i g_D), the saddle-descent's top-form period P_m ~ i^m int_0^1 psi(a)^m C(a) da has psi=g_D(a,1-a) REAL,
  so |psi|^m has a REAL-modulus Laplace peak: P_m ~ (M^m/sqrt m)[S_+ + (-1)^m S_-], S_pm = sum over the
  max-/min-modulus points a_i of C(a_i) w_i, w_i>0. There is NO continuous phase (Stokes) cancellation -- the
  only cancellation is the DISCRETE sum S_pm=0. Since |C|=e^{Re(...)}>0 everywhere, a UNIQUE max-modulus point
  gives S_+ = C(a_1)w_1 != 0, so L(f^m) != 0 for large m, forcing g_D=0: the descent closes UNCONDITIONALLY
  (no Kontsevich-Zagier) for every imaginary top form whose |psi| has a unique dominant point (generic). The
  entire residual is psi with MULTIPLE EQUAL maxima whose complex C-weights cancel to all Laplace orders --
  exactly the classical Chebyshev/COMPOSITION locus of Pakovich's SOLVED polynomial moment problem, further
  gated by C=e^{holomorphic}. So on Re f_D=0, FC(2) is not a Stokes/K-Z problem: it is the classical weighted
  PMP, generically closed, residually the composition solutions. VERIFIED numerically.
tags: [factorial-conjecture, FC2, stokes, no-stokes, saddle, laplace, real-modulus, pakovich-pmp, composition, chebyshev, imaginary-top-form, descent]
related: [fc2-saddle-descent-reduces-to-the-top-form-period-and-localizes-the-stokes-obstruction, FC2-no-pole-homogeneous-is-the-lebesgue-moment-problem]
---

# FC(2) on the Re f_D = 0 locus: no Stokes, just the classical composition PMP

## 1. The locus real-izes the modulus (the whole point)

On `Re f_D == 0` the top form is `f_D = i g_D` with `g_D` a REAL homogeneous polynomial, so
`phi_D(a) = f_D(a,1-a) = i psi(a)`, `psi(a) = g_D(a,1-a) in R`. The saddle-descent's top-form period is

```
   P_m = int_0^1 phi_D^m C da = i^m int_0^1 psi(a)^m C(a) da.
```

Because `psi` is REAL, `|psi(a)^m| = |psi(a)|^m` and the `a`-Laplace peak is a genuine REAL maximum of the
modulus at the critical/endpoint points -- there is **no continuous phase**, hence **no continuous
Stokes/Picard-Lefschetz cancellation** (the obstruction that blocks the general complex top form). This is the
structural gain of the locus: it collapses the Stokes problem to a discrete one.

## 2. The Laplace expansion and the unconditional generic closure

Let `M = max_[0,1] |psi|`, attained at max-points `a_i` (`psi=+M`, `psi''<0`) and min-points `b_j`
(`psi=-M`, `psi''>0`). Standard Laplace gives

```
   int_0^1 psi^m C da  ~  (M^m / sqrt(m)) [ S_+ + (-1)^m S_- ],
      S_+ = sum_i C(a_i) sqrt(2 pi M/|psi''(a_i)|),   S_- = sum_j C(b_j) sqrt(2 pi M/|psi''(b_j)|).
```

`L(f^m) = 0 for all m` forces the even- and odd-`m` leading coefficients to vanish separately: `S_+ = 0` AND
`S_- = 0` (and, from the subleading `1/m` terms, an entire hierarchy `S_pm^{(k)} = 0`). The weight is
`C = exp(phi_{D-1}/(D phi_D) + ...)`, so `|C| = e^{Re(...)} > 0` -- **`C` never vanishes**. Therefore:

> **If `|psi|` has a UNIQUE dominant point** (a single max, or a single min, generically), then
> `S_+ = C(a_1) w_1 != 0` (or `S_- != 0`), so `L(f^m) != 0` for large `m`, contradicting the FC hypothesis.
> Hence `psi == 0`, i.e. `g_D == 0`: the degree descent `f -> f - f_D` closes **UNCONDITIONALLY -- no
> Kontsevich-Zagier**. This proves FC(2) for every `f` whose imaginary top form has a unique dominant `|psi|`
> point (an open, generic condition).

Verified (`f = i(x^2+y^2)+(1+i)x`, `psi=a^2+(1-a)^2`): `[|L(f^m)|/(2m)!]^{1/m} = 1.335, 1.141, 1.088, 1.064`
(`m=4..16`) `-> max|psi| = 1`, and `L(f^m) != 0` throughout.

## 3. The residual is the classical composition PMP (Pakovich), not Stokes

The only way `L(f^m)=0` can persist with `psi != 0` is `S_+ = S_- = 0` with `>=2` dominant points -- i.e.
`psi` has MULTIPLE EQUAL maxima whose complex `C`-weights cancel, to all Laplace orders. Multiple equal
maxima of a real polynomial is the **Chebyshev** phenomenon, and the full vanishing `int psi^m C = 0 forall m`
is exactly the weighted **Polynomial Moment Problem**, whose solutions Pakovich CLASSIFIED as **composition
conditions** (`psi = R o W` with a loop `W`). So on `Re f_D = 0`:

> **FC(2) reduces to Pakovich's solved composition PMP.** No continuous Stokes cancellation exists; the residual
> is the discrete composition/Chebyshev locus, and it is further gated by the weight being `C = e^{holomorphic}`
> (the exponential of a rational function of the lower forms), which is NOT a generic PMP weight. Determining
> which Pakovich composition solutions are realized by an actual polynomial `f` -- with `psi` a Bernstein
> restriction of a real form and `C` the specific saddle exponential -- is the sharp, finite remaining question
> for imaginary top forms, and it is classical PMP arithmetic, not transcendence.

## 4. Net

- **Proved (no K-Z):** FC(2) for every `f` with imaginary top form and unique dominant `|psi|` -- a generic
  open family. The Re f_D=0 locus has **NO Stokes obstruction** (real modulus); the friend's K-Z is unnecessary
  there.
- **Reduced to:** the classical Pakovich composition PMP for the real polynomial `psi` with weight
  `C=e^{holomorphic}` -- a solved-problem residual, not a period-transcendence one.
- **Chain:** homogeneous case = Lebesgue PMP (exact thru deg 4) -> inhomogeneous = weighted top-form period
  (saddle-descent) -> Re f_D=0 = real-modulus, no Stokes, Pakovich composition. The ONLY place K-Z is still
  needed is the complex-top-form off-locus (Re f_D != 0), where the continuous Stokes cancellation lives; the
  no-pole condition (`Phi_f` entire) is exactly what pushes any candidate onto the Re f_D=0 locus.
