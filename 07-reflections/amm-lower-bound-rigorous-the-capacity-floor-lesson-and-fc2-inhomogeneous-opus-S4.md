---
source: opus-2026-07-31-S4 (rigorous AMM lower bound; the capacity-floor lesson; FC(2) inhomogeneous coupling)
status: >
  THEOREM (rigorous lower bound) + method-lesson + FC(2) inhomogeneous step. (1) C*_block >= log_5(5 phi^2)
  is now a THEOREM, not a numerical evaluation: THM-3009's (ARCH) is necessary, its supply is MONOTONE in the
  profile a_k (extremal profile maximizes it), Stirling/Laplace reduces (ARCH) to the entropy comparison (T)
  with uniform control, and my exact Lagrangian gives the (T) threshold = 2 ln(phi)/ln 5 EXACTLY; below it the
  demand exponent H(delta) strictly exceeds the supply exponent so (ARCH) fails for large m. With the cross-
  shell Gale argument, C*_general >= log_5(5 phi^2). (2) LESSON (transferable): a sharp lower/impossibility
  floor = threshold of a SCALE-INVARIANT necessary CAPACITY inequality, obtained by (a) a necessary counting
  inequality, (b) its Laplace/entropy asymptotic (a sup of an entropy functional), (c) a closed-form Lagrangian
  threshold = a generating-function EDGE value. (3) FC(2) INHOMOGENEOUS: the no-pole condition forces a growth
  constraint on the top form f_D (Re f_D <= 0 on the quadrant for entireness), and the saddle/Laplace
  asymptotic of L(f^m) is governed by f_D -- the SAME capacity/Laplace lesson -- giving a degree-descent to the
  homogeneous (=Lebesgue PMP) case, modulo the Stokes structure of the exponential integral.
tags: [amm12592, lower-bound, rigorous, capacity, entropy, laplace, golden, lesson, transfer, factorial-conjecture, FC2, inhomogeneous, no-pole, saddle-point]
related: [THM-3009, cross-shell-floor-is-golden, FC2-no-pole-homogeneous-is-the-lebesgue-moment-problem, amm12592-two-ray-threshold-is-golden]
---

# The AMM lower bound, rigorously; the capacity-floor lesson; and FC(2)'s inhomogeneous coupling

## 1. THEOREM: `C*_block >= log_5(5 phi^2)` (rigorous, upgrading THM-3009's numerics)

THM-3009 states the block floor `~1.598` as "a numerical evaluation of (T)." It is in fact a THEOREM. Four
pieces, three prior + one mine:

1. **(ARCH) is necessary** [THM-3009 sec 2]. Any balanced-block fair extractor on shell `[m,2m)` solves the
   deviation system `sum_k E_k(u)(1+u)^{L_k} = eps u^{m-1}`, `|[u^i]E_k| <= binom(a_k,i)`; evaluating at
   `u=-1` gives, for every `d`,
   ```
      binom(m-1,d)  <=  sum_k binom(a_k, d-L_k) 2^{a_k-d+L_k}.        (ARCH)
   ```
2. **The supply is MONOTONE in the profile.** With `L_k=m-1-k-a_k`, the summand simplifies to
   `binom(a_k, m-1-k-d) 2^{m-1-k-d}` -- the `2`-power is independent of `a_k` and the binomial is
   nondecreasing in `a_k`. So the **extremal profile** `a_k=min(m-1-k, gamma(m+k))` maximizes the supply; if
   (ARCH) fails there, it fails for every profile. (This closes the "which profile?" gap: the constructor's
   best choice is the extremal one.)
3. **Stirling/Laplace reduces (ARCH) to (T)** [standard, verified R-stable]. At `k=xm, d=delta m`,
   `binom(m-1,delta m) = 2^{m H(delta)+o(m)}` and each summand `= 2^{m[alpha H(r/alpha)+(alpha-r)]+o(m)}`
   (`r=delta-ell`); the `<= m+1` terms give `sum = 2^{m*max_x[...]+o(m)}`. Verified: the finite exponent
   `(1/m)log_2(sum)` climbs `0.925, 0.954, 0.958, 0.959` (`m=64..4096`) to `max_x[...] = H(1/phi)` at
   `gamma=golden`.
4. **The (T) threshold is EXACTLY `2 ln(phi)/ln 5`** [mine, symbolic Lagrangian, `amm12592-two-ray-threshold`].
   So for `gamma < golden` there is a `delta` (near `1/phi`) with `H(delta) > max_x[...]` by a POSITIVE
   exponent gap; then `binom(m-1,delta m)` exceeds the supply for large `m` (exponential beats the `m+1`
   polynomial factor), so (ARCH) fails and no balanced-block scheme with that `gamma` exists.

Hence `gamma*_block >= 2 ln(phi)/ln 5` and **`C*_block = 1+gamma*_block >= log_5(5 phi^2) = 1.59799...`,
rigorously.** With the forward-routing Gale + scale-invariance reduction (`cross-shell-floor-is-golden`),
`C*_general >= log_5(5 phi^2)` as well. (The matching upper bound -- a construction achieving `gamma=golden`
for all `R` -- remains open; death-star is at `C=8/5` through `n<=127`.)

## 2. THE LESSON (transferable): a sharp floor is the edge of a scale-invariant capacity inequality

The proof has a shape that recurs and transfers. To prove a sharp impossibility floor:

> **(a)** Find a NECESSARY COUNTING inequality `demand(scale) <= supply(scale)` that any solution must satisfy
> (here: forced parity deficit `<=` box capacity). **(b)** Take its LAPLACE/ENTROPY asymptotic: both sides are
> exponential in the scale, `demand = 2^{m D(delta)}`, `supply = 2^{m sup_x S(x,delta)}`, so the condition is
> the ENTROPY COMPARISON `D(delta) <= sup_x S(x,delta)`. **(c)** Solve the threshold in CLOSED FORM by
> Lagrange multipliers; it is the EDGE VALUE of a generating function (here golden `= C(-1)`, the Catalan GF at
> its branch point). Below the threshold the exponent gap is positive, killing the inequality for large scale.

This is exactly the same machine as: THM-3009 `(ARCH)`; my two-ray comparison; the FMM/multipole
truncation-plus-geometric-tail certificate; and (dually) the S(k) motivic-dimension gate. The transferable
core is **"scale-invariance turns a hard exact problem into an entropy optimization whose Lagrangian threshold
is a generating-function edge value."** The universal `n=2 | n>=3` walls (JC, GM) and level bounds (`C_{1/4}`)
are the same phenomenon: capacity holds until an entropy threshold, then fails.

## 3. FC(2), the inhomogeneous coupling: the same capacity/Laplace lesson

Homogeneous FC(2) is the Lebesgue PMP (`FC2-no-pole-...`). The inhomogeneous coupling yields to the lesson.
Write `f = sum_{d<=D} f_d`, `f_D != 0` the top form. Two capacity/Laplace facts:

- **No-pole constrains the top form.** `Phi_f(t) = int_{[0,inf)^2} e^{t f} e^{-x-y}` has
  `|integrand| = e^{t Re f - x - y}`. For `Phi_f` to be ENTIRE (the FC hypothesis `L(f^m)=0 forall m` forces
  `Phi_f == 1`, hence entire), the exponent `t Re f(x,y) - (x+y)` must stay bounded above for every real `t`,
  as `(x,y) -> infinity`. Since `Re f ~ Re f_D = r^D Re phi_D(a)` (`r=x+y`, `a=x/r`), for `D>=2` this forces
  `Re phi_D(a) <= 0` for all `a in [0,1]` (`t->+inf`) AND `>= 0` (`t->-inf`), i.e. **`Re f_D == 0` on the
  quadrant** -- the top form is purely imaginary-valued on `[0,inf)^2`. (For `D=1`, linear, direct.)
- **The Laplace asymptotic of `L(f^m)` is governed by `f_D`.** `L(f^m)=E[f^m]` is, for large `m`, a saddle
  integral `int e^{m log f - x - y}`; the saddle `m grad f / f = (1,1)` moves to large `(x,y)` where
  `f ~ f_D`, so the leading `m -> infinity` asymptotic of `L(f^m)` is a homogeneous-`f_D` quantity -- exactly
  the Lebesgue-PMP moment of `phi_D`. `L(f^m)=0 forall m` therefore forces the **`f_D` PMP moments to vanish**,
  and by sec 2 of `FC2-no-pole` that gives `phi_D == 0`, i.e. `f_D == 0`: a **degree descent** `f -> f - f_D`.
  Iterating kills all homogeneous parts, so `f == 0`.

**Status of the descent.** The descent is rigorous MODULO the Stokes/Picard-Lefschetz structure of the
exponential integral `Phi_f(t)` -- for complex `f` the saddle contributions can cancel across Stokes rays, and
the "leading asymptotic is the `f_D` PMP moment" step needs the exponential-period monodromy control that is
the Kontsevich-Zagier content. So the inhomogeneous coupling is: **no-pole `=>` `Re f_D == 0` on the quadrant
(rigorous) + saddle/K-Z asymptotic `=>` `f_D` PMP-vanishes `=>` `f_D=0` (modulo Stokes) + descent.** This is
precisely where the friend's Kontsevich-Zagier conjecture enters: it is the tool that legitimizes the
saddle-asymptotic-vanishing step, turning the descent into a full proof of FC(2). The purely imaginary top
form `Re f_D == 0` is a sharp, checkable necessary condition that a fresh FC(2) attack can start from.

## 4. One-line transfers

- **AMM:** floor = entropy-threshold = Catalan edge (golden). RIGOROUS lower bound now.
- **FC(2):** no-pole = entire `Phi_f` => `Re f_D == 0` + saddle-descent to the Lebesgue PMP (homogeneous solved).
- **LRC:** covering floor = capacity of a graded shared-atom family (klein) -- same "demand vs box" shape.
- **S(k):** elementary iff below the motivic-dimension / quadratic-transformation threshold -- the dual gate.
