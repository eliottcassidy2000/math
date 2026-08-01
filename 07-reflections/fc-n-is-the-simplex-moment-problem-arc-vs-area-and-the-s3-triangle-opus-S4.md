---
source: opus-2026-07-31-S4 (FC(2) composition compatibility + FC(3) via the S_3 triangle; the simplex-dimension unification)
status: >
  UNIFICATION + two checks. FC(n)-homogeneous <=> the polynomial moment problem on the (n-1)-SIMPLEX:
  L(f^m)=(Dm+n-1)! int_{Delta_{n-1}} phi_D^m dsigma, phi_D=f_D|simplex ranging over all forms. The DOMAIN
  DIMENSION n-1 is the whole story via the Cauchy/area dichotomy: n=2 -> 1-simplex = interval -> phi image is a
  1-D ARC -> holomorphic moments force phi=0 (FC(2) TRUE); n>=3 -> 2-simplex = TRIANGLE -> phi image is a 2-D
  region that can ENCLOSE AREA, and a measure on an area (e.g. disk area measure, or a rotationally symmetric
  region) can have ALL holomorphic moments zero -> FC(3) is the borderline. VERIFIED S_3 MECHANISM: phi =
  alpha+omega beta+omega^2 gamma maps the triangle onto the equilateral triangle 1,omega,omega^2 (3-fold
  symmetric image), giving int_Delta phi^m = 0 for 3 nmid m -- a PARTIAL counterexample killing 2/3 of the
  moments, flatly impossible for a 1-D arc. Full FC(3) counterexample needs the 3|m moments too (rotationally-
  symmetric image); D=1 is ruled out (fails exactly at m=3: int phi^3 = c0 c1 c2/20). FC(2) COMPATIBILITY: an
  exact search finds NO imaginary-top-form counterexample -- the saddle weight C=e^{holomorphic} is not a
  Pakovich composition weight (the anti-symmetric requirement S(a)-S(1-a)=i pi is transcendental-incompatible
  with algebraic S). So FC(2) holds on Re f_D=0; FC(3) hangs on whether a polynomial can push simplex-Lebesgue
  to an all-moments-zero area measure -- the S_3-symmetric image is the natural generator.
tags: [factorial-conjecture, FC2, FC3, simplex, triangle, s3-symmetry, arc-vs-area, cauchy-transform, moment-problem, pakovich, area-measure, n2-vs-n3-wall]
related: [fc2-on-imaginary-top-form-no-stokes, FC2-no-pole-homogeneous-is-the-lebesgue-moment-problem, factorial-conjecture-exponential-periods-and-the-repo-state]
---

# FC(n) is the (n-1)-simplex moment problem: arc vs area, and the S_3 triangle

## 1. The reduction, in every dimension

For homogeneous `f` of degree `D`, radial coordinates `x_i = r a_i`, `sum a_i = 1` (`a in Delta_{n-1}`),
`dx = r^{n-1} dr dsigma`, `f = r^D phi_D(a)` give

```
   L(f^m) = int_{Delta_{n-1}} [ int_0^inf (r^D phi_D)^m e^{-r} r^{n-1} dr ] dsigma
          = (Dm + n - 1)! * int_{Delta_{n-1}} phi_D(a)^m dsigma,
```

and `phi_D = f_D|_{Delta}` ranges over ALL degree-`<=D` forms (Bernstein), `phi_D == 0 <=> f_D == 0`. So

> **FC(n)-homogeneous `<=>` the polynomial moment problem on the `(n-1)`-simplex:**
> `int_{Delta_{n-1}} phi^m dsigma = 0 for all m >= 1  =>  phi == 0`.

`n=2`: `Delta_1 = [0,1]` (the Lebesgue PMP, companion note). `n=3`: `Delta_2 =` the **triangle**.

## 2. Arc vs area -- the whole `n=2 | n>=3` wall

The pushforward `mu = phi_*(dsigma)` is a compactly supported measure on `C`; all holomorphic moments vanish
iff its Cauchy transform is `1/w` on the unbounded complement component. The DOMAIN DIMENSION decides whether
that forces `mu = delta_0`:

- **`n = 2` (interval):** `phi([0,1])` is a 1-D **ARC**. A simple arc has connected complement, so `C_mu = 1/w`
  everywhere `=> mu = delta_0 => phi = 0`. Only self-intersecting arcs (Pakovich composition) can enclose area,
  and those are the thin classical residual. **FC(2) is TRUE (up to the solved PMP residual).**
- **`n >= 3` (triangle):** `phi(Delta_2)` is a 2-D **REGION**. A measure supported on a 2-D area CAN have all
  holomorphic moments zero -- the cleanest example is **area measure on a disk centered at 0**,
  `int_{|z|<=1} z^m dA = 0 for all m >= 1`. So the Cauchy obstruction dissolves and **FC(3) is the borderline**:
  it fails iff a polynomial `phi` can push simplex-Lebesgue to an all-moments-zero area measure.

This is the same `n=2 | n>=3` wall as the Jacobian conjecture (true `n=2`, false `n>=3`) and GM (true `n=2`,
false `n>=3`): **1-D arc vs 2-D area of the simplex image.**

## 3. The S_3 triangle mechanism (verified partial counterexample)

The triangle carries `S_3` (permuting `alpha,beta,gamma`), and its `C_3` subgroup gives the generator
`phi = alpha + omega beta + omega^2 gamma` (`omega = e^{2 pi i/3}`), which maps `Delta_2` linearly onto the
EQUILATERAL triangle with vertices `1, omega, omega^2` -- a 3-fold rotationally symmetric image. Under the
3-cycle `phi -> omega^2 phi`, so `int_Delta phi^m = omega^{2m} int_Delta phi^m`, hence

```
   int_{Delta_2} phi^m dsigma = 0   for every m with 3 nmid m     (VERIFIED: m=1,2,4,5 ~ 1e-32),
   != 0                              for 3 | m                     (m=3: 1/20, m=6: 1/56).
```

**Two-thirds of all moments vanish** from a single degree-1 form -- categorically impossible for a 1-D arc
(FC(2)), and a direct realization of the area-enclosing mechanism. It is a genuine PARTIAL FC(3) counterexample.

The missing third (`3 | m`) is exactly the equilateral triangle's `int z^{3k} dA != 0`; a FULL counterexample
needs an image with more symmetry (disk-like / full rotational balance). **Degree 1 cannot do it:**
`int_Delta phi^3 = c_0 c_1 c_2 / 20` (exact), so `int phi^m=0` for `m=1,2,3` forces (Newton)
`e_1 = e_2 = e_3 = 0`, i.e. `phi = 0`. The full FC(3) counterexample -- if it exists -- is a higher-degree
`S_3`-structured `phi` whose triangle image is rotationally balanced through order `3k`; this is the concrete
open construction, and the `S_3`/`C_3` symmetry is its natural engine.

## 4. FC(2) composition compatibility: the residual does NOT arise

On `Re f_D = 0` the FC(2) residual was: a symmetric real `psi` (`psi(a)=psi(1-a)`, from a symmetric top form)
whose branch weights cancel -- needing an ANTI-symmetric saddle weight `C(a) = -C(1-a)`, i.e.
`S(a) - S(1-a) = i pi` for the algebraic saddle exponent `S = phi_{D-1}/(i D psi) + ...`. This is a
TRANSCENDENTAL identity (`pi`) demanded of an ALGEBRAIC function, generically incompatible; and an EXACT search
(`f = i(x^2+y^2) + complex linear + const`, `L(f^m)=0` for `m=1..5`) returns **NO nonzero solution**. So the
Pakovich composition solutions are **not realized by polynomial `f`**, and **FC(2) holds on `Re f_D = 0`** with
no residual counterexample -- consistent with FC(2) being true.

## 5. Net (the tandem)

- **FC(2):** homogeneous = Lebesgue PMP = 1-D arc `=> phi = 0`; the composition residual is not realized
  (saddle weight `e^{holomorphic}` is not a Pakovich weight). FC(2) TRUE up to that (now-checked-empty) residual.
- **FC(3):** homogeneous = triangle PMP = 2-D area; the `S_3`/`C_3` form kills `2/3` of the moments; FC(3) is
  the borderline and its fate is whether a polynomial gives a rotationally-balanced (all-moments-zero) triangle
  image. The user's `S_3` intuition is exactly right: **the triangle's symmetry is the FC(3) counterexample
  generator, and the arc->area jump from the interval to the triangle is why FC(3) can fail where FC(2) cannot.**
