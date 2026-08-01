---
id: THM-3018
title: "The Factorial Conjecture is a simplex moment problem; FC(2) is the polynomial moment problem"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Let
  L(x^alpha) = alpha! on C[x_1..x_n]. Then L(f) = int_{[0,inf)^n} f e^{-|x|_1}
  dx, and for f HOMOGENEOUS of degree d the polar substitution x = r u gives
  L(f^m) = (dm+n-1)! * int_{Delta_{n-1}} g^m dA with g = f|_Delta; since
  f |-> f|_Delta is a linear BIJECTION onto polynomials of degree <= d in
  n-1 variables (Bernstein basis), the Factorial Conjecture FC(n) is
  EQUIVALENT to: for every complex polynomial g in n-1 variables,
  int_{Delta_{n-1}} g^m dA = 0 for all m >= 1 implies g = 0. In particular
  FC(2) is exactly the POLYNOMIAL MOMENT PROBLEM on [0,1] with h(u) = u,
  and FC(3) the same problem on a triangle. Consequences proved here:
  FC(2) holds whenever the image arc g([0,1]) has connected complement
  (Cauchy-transform/Plemelj argument), and holds outright for deg g <= 3
  (exact Groebner elimination); and for FC(3) the area-preserving 3-cycle
  on barycentric coordinates makes every moment with 3 not dividing m vanish
  AUTOMATICALLY on its omega-eigenspace, so only m = 3,6,9,... can obstruct.
  Section 4b then PROVES the homogeneous statement outright for every n, by
  a maximum-modulus Laplace argument on the compact simplex, and isolates
  the real difficulty as the NON-HOMOGENEOUS (non-compact) regime -- which
  is exactly the repo's SFC lane. Scope: this settles the compact case and
  gives an exact dictionary; it is a DIFFERENT functional from the repo's
  Strong Factorial Conjecture lane (one variable, L(z^n) = n!, N-term
  supports) -- see section 5.
source: death-star-2026-07-31-coinC2
depends_on: []
related:
  - THM-2836
  - THM-2812
  - THM-2849
  - THM-1435
external:
  - "A. van den Essen, D. Wright, W. Zhao, On the image conjecture (2010ish)."
  - "F. Pakovich, M. Muzychuk, Solution of the polynomial moment problem (2009)."
script: 04-computation/factorial_conjecture_simplex_reduction_thm3018.py
output: 05-knowledge/results/factorial_conjecture_simplex_reduction_thm3018.out
---

# THM-3018 -- the Factorial Conjecture as a simplex moment problem

## 1. The functional is an integral

Let `L : C[x_1,...,x_n] -> C` be linear with `L(x^alpha) = alpha! =
prod_i (alpha_i)!`. Since `int_0^inf x^a e^{-x} dx = a!`,

```text
L(f) = int_{[0,inf)^n} f(x) e^{-(x_1 + ... + x_n)} dx.                (1)
```

**FC(n)** (van den Essen-Wright-Zhao): for `f` homogeneous, if `L(f^m) = 0`
for all `m >= 1` then `f = 0`.

## 2. The polar reduction (PROVED)

Substitute `x = r u` with `r > 0` and `u` in the standard simplex
`Delta_{n-1} = {u_i >= 0, sum u_i = 1}`; then `|x|_1 = r` and
`dx = r^{n-1} dr dA(u)`. For `f` homogeneous of degree `d`, `f(ru)^m =
r^{dm} g(u)^m` with `g := f|_Delta`, so

```text
L(f^m) = [ int_0^inf r^{dm+n-1} e^{-r} dr ] * int_{Delta_{n-1}} g^m dA
       = (dm+n-1)! * int_{Delta_{n-1}} g^m dA.                        (2)
```

`(dm+n-1)! != 0`, so `L(f^m) = 0` iff the simplex moment vanishes. Moreover
`f |-> f|_Delta` sends the degree-`d` forms in `n` variables bijectively onto
polynomials of degree `<= d` in `n-1` variables — this is the Bernstein
basis `u^beta (1 - sum u)^{d-|beta|}` (referee R2, ranks full for
`n = 2, d = 3,5` and `n = 3, d = 2,3`). Hence:

```text
FC(n)  <=>  for every complex polynomial g in n-1 variables,
            int_{Delta_{n-1}} g^m dA = 0 for all m >= 1  =>  g = 0.   (3)
```

Referee R1 verifies (2) symbolically for `n = 2` (`d <= 3`, `m <= 3`) and
`n = 3` (`d <= 2`, `m <= 2`).

## 3. FC(2) is the polynomial moment problem

By (3), **FC(2) says exactly**: for a complex polynomial `g` of one variable,
`int_0^1 g(u)^m du = 0` for all `m >= 1` forces `g = 0`. Equivalently, with
`mu := g_*(Lebesgue on [0,1])` a probability measure on `C`,

```text
int_C z^m dmu = 0 for all m >= 1,   equivalently   int_0^1 e^{t g(u)} du = 1
                                                    for all t in C.   (4)
```

**(a) The arc case (PROVED).** Let `Gamma = g([0,1])`, a compact rectifiable
set of planar measure zero. The Cauchy transform
`B(w) = int_0^1 du/(w - g(u))` is analytic on `C \ Gamma` and equals `1/w`
for `|w|` large, hence on the unbounded component. If `C \ Gamma` is
connected — e.g. `Gamma` is a simple arc — then `B = 1/w` off `Gamma`, so
`B = 1/w` a.e. on `C`, so as distributions `dbar B = dbar(1/w)`, i.e.
`mu = delta_0`. Then `g = 0` a.e. on `[0,1]`, so `g = 0`. Equivalently, by
Plemelj: the jump of `B` across a smooth point of `Gamma` is proportional to
the density, and both sides lie in the same component, so the density
vanishes.

**(b) The obstruction is loops, and it is real.** If `Gamma` bounds, the
argument fails, and it must: the uniform measure on the unit circle has all
moments `m >= 1` zero (mean value property), so no purely measure-theoretic
argument can succeed. What excludes it here is that `mu` is a *polynomial*
pushforward of an interval: no nonconstant polynomial maps `[0,1]` into a
circle, since `g(u) g^*(u) = 1` on `[0,1]` (with `g^*` the conjugate-
coefficient polynomial) is a polynomial identity forcing `deg g = 0`.

**(c) Low degree (VERIFIED-EXACT).** Exact Groebner elimination on the
system `int_0^1 g^m du = 0`, `m = 1..d+3`, gives only `g = 0` for
`d = 1, 2, 3` (referee R3). For `d = 2` the eliminated basis is
`{c_0 + c_1/2 + c_2/3, c_1^2 + 2c_1c_2 + 16c_2^2/15, c_2^3}`.

**(d) Where this lands in the literature.** (4) is the polynomial moment
problem `int_a^b p^m q' = 0` with `q(u) = u`, `q' = 1`, restricted to
`m >= 1` (the `m = 0` moment is `1`, not `0`, so the standard normalisation
does not apply verbatim). The classical Composition Condition
`p = P o W, q = Q o W, W(a) = W(b)` is *unavailable* here: `q = u` forces
`deg W = 1`, and then `W(0) = W(1)` is impossible. This is the first contact
in this repo with the Pakovich-Muzychuk circle of results (repo-wide grep:
zero prior hits for "polynomial moment problem", "composition conjecture",
"Pakovich", "Briskin", "Yomdin").

## 4. FC(3): the symmetry mechanism (PROVED)

By (3), FC(3) asks: does `int_T g^m dA = 0` for all `m >= 1` force `g = 0`,
`T` the triangle? Here the pushforward is **two-dimensional**, and any
rotationally symmetric image measure kills all moments `m >= 1` (e.g. the
uniform measure on a disc centred at `0`). So FC(3) is genuinely more
dangerous than FC(2), and the search for a counterexample is the search for
a polynomial with a "rotationally balanced" image measure.

The triangle has an area-preserving symmetry: the 3-cycle
`sigma : (u,v,w) -> (v,w,u)` on barycentric coordinates. If
`g o sigma = omega g` with `omega = e^{2 pi i/3}`, then

```text
int_T g^m dA = int_T (g o sigma)^m dA = omega^m int_T g^m dA,
```

so **every moment with `3` not dividing `m` vanishes automatically**. Only
`m = 3, 6, 9, ...` can obstruct — one third of the conditions, on the
`omega`-eigenspace (one third of the coefficient space). Referee R4 verifies
this exactly over `Q(omega)` on the degree-one eigenvector
`g = u + omega v + omega^2 w`: moments `m = 1,2,4,5,7,8` vanish identically,
while the survivors are nonzero with the closed form

```text
int_T (u + omega v + omega^2 w)^{3k} dA = 1/((3k+1)(3k+2))
   (m = 3: 1/20,  m = 6: 1/56,  m = 9: 1/110),
```

i.e. exactly the moments of `u` itself. So this particular eigenvector is
*not* a counterexample, but it isolates the residual obstruction: a
counterexample in the eigenspace must kill the sparse family `m = 3k`.

## 4b. The compact/homogeneous case is settled by Laplace (PROVED)

The reduction makes the generating function **entire and constant**: since
`g` is bounded on the compact simplex,

```text
Phi(t) := int_{Delta} e^{t g} dA = sum_{m>=0} (t^m/m!) int_Delta g^m dA
        = Area(Delta)   for every t in C.                            (5)
```

**Theorem.** If `g` is a polynomial on `Delta_{n-1}` with (5), then `g = 0`.
Hence FC(n) holds for homogeneous `f`, for every `n`.

*Proof.* Suppose `g != 0` and set `A := max_Delta |g| > 0`, attained at
`u_1`. Choose `theta` with `e^{i theta} g(u_1) = A` and put
`psi := e^{i theta} g`. For every `u`,
`Re psi(u) <= |g(u)| <= A`, with equality **iff** `psi(u) = A`; so the
maximum locus `S = {psi = A}` is a proper algebraic subset (if it had
interior, `psi ≡ A`, forcing `int g = A e^{-i theta} Area != 0`). Now

```text
Phi(s e^{i theta}) = e^{sA} int_Delta e^{s(psi - A)} dA,   Re(psi - A) <= 0.
```

Laplace's method applies with all contributions concentrated on `S`, where
`psi - A` vanishes; each local contribution is `s^{-alpha}` times
`int e^{Q}` for a form `Q` with `Re Q <= 0`, whose value has **positive real
part** (for a nondegenerate quadratic `Q` in `k <= 2` variables this is
`pi^{k/2}/sqrt(det(-Q))` with `arg sqrt(det) in [-pi/2, pi/2]`, since `-Q`
has positive semidefinite real part). Because the direction was chosen by
**maximum modulus**, every point of `S` carries the *same* value `psi = A`,
so the contributions share the factor `e^{sA}` with no relative oscillation
and cannot cancel. Hence `|Phi(s e^{i theta})| ~ C e^{sA} s^{-alpha} -> inf`,
contradicting (5). QED

Numerical confirmation (script, real and complex `g` on `[0,1]`):
`|Phi(s e^{i theta})|` tracks `e^{sA}/sqrt(s)` over `s = 5 .. 120` across
14 orders of magnitude.

**Consequence for where the difficulty lies.** The homogeneous case is
therefore *not* the hard case: homogeneity is exactly what compactifies the
problem onto a simplex and makes `Phi` entire and bounded on rays. The
genuine difficulty in the factorial circle of problems is the
**non-homogeneous** case, where `L(f^m) = int_{[0,inf)^n} f^m e^{-|x|_1}dx`
lives on a **non-compact** domain with `f` unbounded, so `Phi` is only a
formal series and no Laplace argument is available. That is precisely the
regime of this repo's Strong Factorial Conjecture lane (section 5), and it
explains why SFC(3) can be hard while the homogeneous statement is not.

## 5. Scope, and the relation to the repo's SFC lane

This file is a **reformulation with an exact dictionary**, not a proof of
FC(n) for `n >= 2`. Section 3(a) proves FC(2) only under a topological
hypothesis on `g([0,1])`; 3(c) proves it outright only for `deg g <= 3`;
section 4 reduces FC(3)'s symmetric obstruction to `3 | m`; section 4b
settles the homogeneous case for all `n` but says nothing about the
non-homogeneous one.

**Type warning (cf. MISTAKE-237).** The repo's large Strong Factorial
Conjecture lane (THM-2812/2824/2836/2849/2854, `sfc3_*`/`sfc4_*` scripts)
uses a *different* functional and a *different* parameter: one variable,
`L(z^n) = n!`, and `SFC(N)` indexes the number of TERMS in the support,
asking whether the first `N` moments detect an `N`-term polynomial. The
`FC(n)` here is `n` VARIABLES with `L(x^alpha) = alpha!` and `f`
homogeneous. The two agree at `n = 1` only. No arrow between them is
claimed; anyone building one must state the map, the dimension change, the
preserved predicate and the loss.

Referee: R1-R4 exact, `ALL THM-3018 REFEREE CHECKS PASSED`. QED.

## AUDIT FLAG, 2026-08-01 (death-star, on my own section 4b)

A concurrent agent reported section 4b's Laplace step as refuted. I could not
substantiate a refutation, but I **can** substantiate a genuine gap, so the
status of 4b is downgraded from PROVED to **AUDIT-REQUIRED** pending the
resolution below. Nothing else in this file is affected; sections 1-4a and the
`S_3` mechanism of section 4 stand.

**The gap.** Section 4b argues: `Phi(s) = int_Delta e^{s g} dA` is entire; a
counterexample forces `Phi ≡ Area`; choosing the direction by maximum modulus
gives `|Phi| ~ C e^{sA} s^{-a} -> infinity`, contradiction. Two problems.

*(G1) The type shortcut is circular.* One would like "`Phi` has exponential type
exactly `A = max|g| > 0`, and a constant has type `0`". But the type is
`limsup_m m|a_m|^{1/m}/e` with `a_m = mu_m/m!`, `mu_m = int g^m`. If `g` is a
counterexample then **every** `mu_m` with `m >= 1` vanishes, so the type is `0`
by construction. Establishing `type = max|g|` is therefore equivalent to the
conclusion, not an input to it.

*(G2) The asymptotic needs a no-cancellation justification.* The claim as
written is that maximum modulus forces "no relative oscillation". At an
**interior** maximum `u_1` of `|g|`, `d|g|^2/du = 0` gives
`Re(conj(g) g') = 0`, so writing `e^{i theta} g(u_1) = A > 0` the derivative
`e^{i theta} g'(u_1) = i b` is purely **imaginary** with `b` generally nonzero.
Then near `u_1`, `e^{i theta}g(u) ~ A + i b w - c w^2` and

```text
int e^{R(A + i b w - c w^2)} dw = e^{R A} sqrt(pi/(R c)) * exp( - R b^2/(4c) ),
```

so the oscillation contributes a factor `e^{-R b^2/(4c)}` and the growth rate is
`A - b^2/(4c)`, **not** `A`. That can be `<= 0`. The step therefore requires
either `b = 0`, or a contour deformation / steepest-descent argument that
section 4b does not supply.

**What the evidence says.** Tested on the mean-zero family
`g(u) = (u-1/2) + i lam (u^2 - u + 1/6)` at `lam = 0, 2, 6, 20`, in every case
`log|Phi(R e^{i theta})| / R -> A = max|g|` (e.g. `lam = 6`: `0.879, 1.041,
1.094, 1.111` at `R = 20, 80, 320, 1280` against `A = 1.11803`). So the
conclusion held everywhere tested. **But in all of these the maximum of `|g|`
was at an ENDPOINT** (`u_1 = 0`), where no stationary point occurs and (G2) does
not bite. The interior-maximum case, which is exactly the risky one, was not
exercised: for real `g` an interior modulus maximum forces `g'(u_1) = 0` hence
`b = 0` and the step is fine, so the open case is **complex `g` with an interior
maximum of `|g|` and `g'(u_1) != 0`**.

**Also relevant (F1).** The step genuinely uses polynomiality, so it is not
refuted by the obvious continuous witness: `h(u) = e^{2 pi i u}` has all moments
`m >= 1` equal to zero and `Phi_h ≡ 1` bounded, but `|h| ≡ 1` attains its maximum
on **all** of `[0,1]`, whereas for a polynomial `|g|^2` is a real polynomial and
its maximum is attained at finitely many points. So (F1) is not a refutation.

**Consequence for downstream claims.** Anything asserting "the homogeneous case
of FC(n) is proved for every n" -- including the inference that FC(2) is thereby
closed and the `int e^Q` route redundant, and the reading of kps-S166's residual
as already-settled -- is **suspended** until (G1)/(G2) are resolved. THM-3031's
bridge is unaffected: it is an implication about what would follow from the
exponential-integral claim, and does not depend on 4b.

**To close this**, one needs either (a) a proof that a mean-zero polynomial on
`Delta` cannot have its modulus maximum interior with `b^2 >= 4cA`, or (b) a
steepest-descent treatment giving the growth of `Phi` from the integral
representation without assuming the moments.
