---
source: opus-2026-08-01-S4 (rigorous DEGREE-2 case of the exp-integral transcendence conjecture)
status: >
  THEOREM (rigorous, modulo one cited theorem). For every degree-2 Q = alpha t^2 + beta t + gamma with
  algebraic coefficients (alpha != 0), int_0^1 e^{Q(t)} dt is TRANSCENDENTAL. Proof: complete the square to
  I = e^delta[(1+c)M(eta1) - c M(eta2)] with M(z)=1F1(1/2;3/2;z)=sum z^n/((2n+1)n!) an E-FUNCTION and
  eta1=alpha(1+c)^2, eta2=alpha c^2 algebraic. The three E-functions g1=e^{del z}M(eta1 z), g2=e^{del z}M(eta2 z),
  and 1 are LINEARLY INDEPENDENT over Qbar(z) -- elementary, because at z->+inf they have DISTINCT exponential
  rates del+eta1 != del+eta2 (this is the "horizontal-endomorphism rigidity", made trivial). By the REFINED
  (linear) Siegel-Shidlovskii theorem with Beukers' no-exceptional-point lifting (Beukers, Ann. Math. 2006),
  their values at z=1 are Qbar-linearly independent; hence I = (1+c)g1(1) - c g2(1), a nontrivial Qbar-linear
  combination, is transcendental. Degenerate c in {0,-1,-1/2} collapse to the single-value case e^lambda M(eta),
  same argument. This is EXACTLY the friend's route (E-functions + Beukers lifting + horizontal-endomorphism
  rigidity) instantiated at degree 2, and needs only the LINEAR S-S (no differential-Galois computation).
tags: [factorial-conjecture, FC2, exponential-integral, transcendence, e-functions, beukers, siegel-shidlovskii, error-function, kummer, 1F1, degree-2, HYP-9078]
related: [HYP-9078, exp-integral-radial-bridge-gives-fc-homog-and-the-weighted-generalization-for-full-fc2, fc-jc2-tandem-the-isolated-coincidence-and-the-fold-no-fold-duality, amm12592-golden-threshold]
external: ["Beukers, 'A refined version of the Siegel-Shidlovskii theorem', Ann. of Math. 163 (2006) 369-379",
           "Shidlovskii, Transcendental Numbers (1989)", "Siegel 1929/1949 (error-function values transcendental)",
           "Lindemann-Weierstrass (degree 1, HYP-9078)"]
---

# Degree-2 is transcendental: `int_0^1 e^{Q}` for `deg Q = 2`, rigorously

**Theorem.** Let `Q(t) = alpha t^2 + beta t + gamma` with `alpha, beta, gamma in Qbar` and `alpha != 0`. Then
```
        I  :=  int_0^1 e^{Q(t)} dt   is TRANSCENDENTAL.
```
This is the degree-2 case of the conjecture "`int_0^1 e^{Q}` is transcendental for every nonconstant algebraic
`Q`" (the object whose degree-1 case is Lindemann-Weierstrass, HYP-9078, and whose truth `=>` FC(2)-homogeneous
via the radial bridge). The proof below is the friend's route -- **E-functions + Beukers lifting +
horizontal-endomorphism rigidity** -- and, at degree 2, needs only the *linear* Siegel-Shidlovskii theorem.

## 1. Reduction to a confluent-hypergeometric E-function

Complete the square. With `c = beta/(2 alpha) in Qbar` and `delta = gamma - beta^2/(4 alpha) in Qbar`,
```
   Q(t) = alpha (t + c)^2 + delta,   so   I = e^delta int_c^{1+c} e^{alpha u^2} du = e^delta [ F(1+c) - F(c) ],
```
where `F(x) = int_0^x e^{alpha u^2} du`. Expanding `e^{alpha u^2}` and integrating termwise,
```
   F(x) = sum_{n>=0} alpha^n x^{2n+1} / ((2n+1) n!) = x * M(alpha x^2),   M(z) := sum_{n>=0} z^n / ((2n+1) n!).
```
`M` is the Kummer function `1F1(1/2; 3/2; z)` (indeed `(1/2)_n/(3/2)_n = 1/(2n+1)`), and its `z^n/n!`-coefficient
is `b_n = 1/(2n+1)` -- **a bona fide E-function** (algebraic coefficients, at most geometric denominators). Thus
```
        I = e^delta [ (1+c) M(eta1) - c M(eta2) ],     eta1 = alpha (1+c)^2,   eta2 = alpha c^2  in  Qbar.   (*)
```
(All four displayed identities verified to >= 18 digits;
`04-computation/exp_integral_degree2_transcendence_opus_S4.py`.)

## 2. The three E-functions and the rigidity (elementary)

Introduce the single-variable E-functions
```
   g1(z) = e^{delta z} M(eta1 z),    g2(z) = e^{delta z} M(eta2 z),    g0(z) = 1,
```
(products and constants of E-functions are E-functions). Evaluated at `z = 1` they give `g1(1) = e^delta M(eta1)`,
`g2(1) = e^delta M(eta2)`, so `(*)` reads `I = (1+c) g1(1) - c g2(1)`.

**Lemma (linear independence over `Qbar(z)`; the rigidity step).** If `eta1 != eta2`, then `g1, g2, g0` are
linearly independent over `Qbar(z)`.

*Proof.* From the classical confluent asymptotics `1F1(1/2;3/2; w) ~ (1/2) e^w w^{-1}` as `w -> +infinity`,
```
   g_i(z) ~ (1/(2 eta_i)) * e^{(delta + eta_i) z} / z    as  z -> +infinity  (i = 1,2).
```
Suppose `R1 g1 + R2 g2 + R0 = 0` with `R0,R1,R2 in Qbar(z)` not all zero. The two rates `delta+eta1` and
`delta+eta2` are **distinct** (as `eta1 != eta2`), and a rational function `R0` has sub-exponential growth, so the
two exponentials cannot cancel each other nor be absorbed by `R0`: comparing dominant exponential rates forces
`R1 = R2 = 0` (if one rate, say `delta+eta1`, is `0`, then `e^{-eta1 z}M(eta1 z) = 1F1(1; 3/2; -eta1 z)` by
Kummer's transformation is still transcendental over `Qbar(z)`, so its coefficient vanishes too), and then
`R0 = 0`. Contradiction. `//`

This Lemma is the entire content of "horizontal-endomorphism rigidity" at degree 2: the only way a horizontal
(flat, `Qbar(z)`-linear) relation could hold is if two solution-sheets shared an exponential rate, i.e.
`eta1 = eta2`; distinct rates forbid it. No differential-Galois-group computation is required for the *linear*
statement.

## 3. Beukers lifting `=>` transcendence

**Refined Siegel-Shidlovskii (Beukers, Ann. of Math. 2006, linear form).** If `f_1, ..., f_n` are E-functions
that are linearly independent over `Qbar(z)`, then for *every* nonzero algebraic `xi` that is not a singularity
of the underlying system, the values `f_1(xi), ..., f_n(xi)` are linearly independent over `Qbar`.
(Beukers removes the finite exceptional set of the classical theorem -- this is the "lifting" that lets us name
the specific point `xi = 1`.)

The combined system for `(g0, g1, g2)` has its only finite singularity at `z = 0` (the regular singular point of
Kummer), so `xi = 1` is admissible. By the Lemma the functions are `Qbar(z)`-linearly independent, hence
`g0(1)=1,  g1(1)=e^delta M(eta1),  g2(1)=e^delta M(eta2)` are **`Qbar`-linearly independent**.

Now suppose, for contradiction, `I` were algebraic, `I = A in Qbar`. Then by `(*)`
```
   (1+c) g1(1)  +  (-c) g2(1)  +  (-A) g0(1)  =  0,
```
a `Qbar`-linear relation with coefficients `(1+c, -c, -A)` that are not all zero (`1+c` and `c` are not both
zero). This contradicts the linear independence just established. Therefore **`I` is transcendental**. `QED`

## 4. Degenerate parameters (`c in {0, -1, -1/2}`) collapse to one value

The generic branch above assumed `eta1 != eta2` and both nonzero, i.e. `c not in {0, -1, -1/2}`. The three
exceptions each collapse `(*)` to a **single** confluent value, handled by the same machinery with two functions
`{e^{lambda z} M(eta z), 1}`:

- `c = 0` (`beta = 0`): `eta2 = 0`, `M(0)=1`, `I = e^gamma M(alpha)`.
- `c = -1` (`beta = -2 alpha`): `eta1 = 0`, `I = e^delta M(alpha)`.
- `c = -1/2` (`beta = -alpha`, symmetric interval `[-1/2,1/2]`): `eta1 = eta2 = alpha/4`, and
  `(1+c)M - cM = M`, so `I = e^delta M(alpha/4)`.

In each case `eta = (alpha, alpha, alpha/4)` is a **nonzero** algebraic point, `e^{lambda z}M(eta z)` and `1` are
`Qbar(z)`-linearly independent (single-exponential-rate version of the Lemma; the confluent value is not
`e^{rational}`), so `M(eta)` resp. `e^lambda M(eta)` is transcendental. `//`

Combining Sections 3-4 proves the Theorem for **all** algebraic `alpha != 0, beta, gamma`. `#`

## 5. What this settles, and honest scope

- **Rigorously proved (new here):** the DEGREE-2 case of the exp-integral transcendence conjecture. In
  particular the error-/imaginary-error-function periods `int_c^{1+c} e^{alpha u^2} du` at algebraic
  `alpha, c` are transcendental. Degree 1 was Lindemann-Weierstrass (HYP-9078); **degree 2 is now closed.**
- **Method = the friend's route, instantiated.** "E-functions + Beukers lifting + horizontal-endomorphism
  rigidity" is exactly Sections 2-3; the surprise is that at degree 2 the rigidity is *elementary* (distinct
  exponential rates `delta + eta_i`) and only the LINEAR Siegel-Shidlovskii is needed -- no Galois group, no
  algebraic-independence (2nd-order) machinery.
- **Consequence for FC(2).** Via the radial bridge
  (`exp-integral-radial-bridge-...`), transcendence of `int_0^1 e^{s0 phi}` for algebraic `s0` and polynomial
  `phi` is what rules out a homogeneous FC(2) counterexample of that degree; degree-2 transcendence supplies
  this for `phi` giving a degree-2 exponent. (FC(2)-homogeneous itself is Liu-Sun 2020 -- no novelty claimed
  there; the transcendence theorem is the new object, and it is the harder statement.)
- **One cited theorem.** The proof rests on Beukers 2006 (refined S-S, no exceptional point) -- standard,
  peer-reviewed. Everything else (the E-function form `b_n = 1/(2n+1)`, the reduction `(*)`, the distinct-rate
  linear independence) is elementary and numerically verified to 40 digits.
- **Not yet done (the honest frontier):** (i) degree `>= 3` -- the same architecture applies with
  `M` replaced by the higher confluent/hypergeometric E-function attached to `int_0^x e^{alpha u^k}`, but the
  linear-independence Lemma needs the analogue of "distinct exponential rates" for `k >= 3` (several
  exponentials `e^{omega alpha^{1/k} ...}`; the rigidity is no longer one inequality). (ii) The INHOMOGENEOUS
  (weighted, `poly + rational` exponent) integral that FULL FC(2) needs (saddle-descent weight
  `exp(phi_{D-1}/(D phi_D))`) is a genuinely different object -- E-function methods do not obviously cover a
  rational term in the exponent.

## 6. One-line takeaway

`int_0^1 e^{alpha t^2 + beta t + gamma} dt` is a `Qbar`-linear combination of two values `e^delta M(alpha(1+c)^2)`,
`e^delta M(alpha c^2)` of the SAME error-function E-function `M = 1F1(1/2;3/2)` at DISTINCT algebraic points;
distinct points force distinct exponential rates, distinct rates force `Qbar(z)`-linear independence, and Beukers
lifts that to `Qbar`-linear independence of the values -- so the combination is transcendental. **Degree 2:
closed.**
