---
id: THM-3018
title: "The homogeneous Factorial-Conjecture subclass as a simplex moment problem"
status: >
  CORRECTED / PARTIALLY PROVED + VERIFIED-EXACT + AUDIT-REQUIRED.  For the
  factorial functional on all polynomials in n variables, the exponential
  integral identity is exact.  On the homogeneous degree-d subclass only,
  polar substitution gives L(f^m)=(dm+n-1)! times the simplex moment of
  g=f|_Delta, and restriction is a Bernstein-basis bijection.  Hence the
  homogeneous FC(n) subclass, not full FC(n), is equivalent to a polynomial
  moment problem on Delta_(n-1).  The connected-complement argument, exact
  degree-at-most-three elimination, and the triangle 3-cycle selection rule
  survive with that homogeneous scope.  Section 4b's proposed Laplace
  closure remains AUDIT-REQUIRED because maximum-modulus localization does
  not control oscillatory cancellation at an interior saddle.  MISTAKE-350
  also corrects the former slot-count use of SFC(N): the univariate N-slot
  lane is a restriction of ambient SFC(1), not SFC(N).
source: death-star-2026-07-31-coinC2
depends_on: []
related:
  - THM-2836
  - THM-2812
  - THM-2849
  - THM-1435
external:
  - "E. Edo and A. van den Essen, The Strong Factorial Conjecture, arXiv:1304.3956v2, Definitions 2.1 and 2.7, Conjectures 2.4 and 2.8."
  - "A. van den Essen, D. Wright, W. Zhao, On the image conjecture (2010ish)."
  - "F. Pakovich, M. Muzychuk, Solution of the polynomial moment problem (2009)."
script: 04-computation/factorial_conjecture_simplex_reduction_thm3018.py
output: 05-knowledge/results/factorial_conjecture_simplex_reduction_thm3018.out
---

# THM-3018 -- the homogeneous Factorial-Conjecture subclass as a simplex moment problem

## 1. The functional is an integral

Let `L : C[x_1,...,x_n] -> C` be linear with `L(x^alpha) = alpha! =
prod_i (alpha_i)!`. Since `int_0^inf x^a e^{-x} dx = a!`,

```text
L(f) = int_{[0,inf)^n} f(x) e^{-(x_1 + ... + x_n)} dx.                (1)
```

**FC(n)** quantifies over every `f in C[x_1,...,x_n]`: if `L(f^m)=0` for all
`m>=1`, then `f=0`.  Homogeneity is not part of the original definition.
Write **HFC(n)** locally for its restriction to homogeneous polynomials.
The legal FC(1) input `1+x` is the minimal witness that this distinction is
load-bearing: it has no single polar degree.

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
HFC(n) <=>  for every complex polynomial g in n-1 variables,
            int_{Delta_{n-1}} g^m dA = 0 for all m >= 1  =>  g = 0.   (3)
```

Referee R1 verifies (2) symbolically for `n = 2` (`d <= 3`, `m <= 3`) and
`n = 3` (`d <= 2`, `m <= 2`).

## 3. The homogeneous FC(2) subclass is the polynomial moment problem

By (3), **HFC(2) says exactly**: for a complex polynomial `g` of one variable,
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

## 4. The homogeneous FC(3) subclass: the symmetry mechanism (PROVED)

By (3), HFC(3) asks: does `int_T g^m dA = 0` for all `m >= 1` force `g = 0`,
`T` the triangle? Here the pushforward is **two-dimensional**, and any
rotationally symmetric image measure kills all moments `m >= 1` (e.g. the
uniform measure on a disc centred at `0`). So HFC(3) is genuinely more
dangerous than HFC(2), and the search for a counterexample is the search for
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

## 4b. Proposed compact/homogeneous Laplace closure (AUDIT-REQUIRED / UNPROVED)

The argument in this subsection is retained as a proof candidate only. The
audit below identifies a circular exponential-type shortcut and an unhandled
interior-phase cancellation. Nothing here belongs to the proved dependency
graph.

The reduction makes the generating function **entire and constant**: since
`g` is bounded on the compact simplex,

```text
Phi(t) := int_{Delta} e^{t g} dA = sum_{m>=0} (t^m/m!) int_Delta g^m dA
        = Area(Delta)   for every t in C.                            (5)
```

**Candidate statement.** If `g` is a polynomial on `Delta_{n-1}` with (5),
then `g=0`.  This would prove HFC(n) for every `n`, but the proof below has
the oscillatory-saddle gap recorded at the end of this file and is not canon.

*Unproved proof candidate.* Suppose `g != 0` and set
`A := max_Delta |g| > 0`, attained at
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

**No consequence is licensed from this historical paragraph.** Homogeneity compactifies the integral onto
a simplex, but compactness does not itself prevent complex saddle
cancellation.  If the gap below is repaired, the homogeneous subclass would
close and the remaining full-FC difficulty would be nonhomogeneous on the
noncompact orthant.  At present both that homogeneous closure and full FC
remain open beyond the proved fragments above.

## 5. Scope, and the relation to the repo's SFC lane

This file gives a **homogeneous-subclass reformulation**, not a proof of full
`FC(n)` for any `n>=2`.  Section 3(a) proves the HFC(2) implication under a
topological hypothesis on `g([0,1])`; section 3(c) verifies it through degree
three; section 4 gives the HFC(3) cyclic selection rule; section 4b is an
unproved candidate.

**Type warning (MISTAKE-350).** The functional is not different from the
repo's univariate lane: setting `n=1` gives exactly `L(z^j)=j!`.  What differs
is the restriction.  Original `FC(n)` and `SFC(n)` index **ambient variable
dimension** and quantify over all polynomials.  For `SFC(n)`, the separate
quantity `N(f)` is the support size and the window length.  Therefore the
legacy `sfc3_*`/`sfc4_*` artifacts are three-/four-slot restrictions of
`SFC(1)`, not ambient `SFC(3)`/`SFC(4)`.  Their exact computations survive;
their old names are provenance only.

Referees R1--R4 verify the displayed polar identities, finite eliminations,
and cyclic selection rule.  They do not verify section 4b.

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
of FC(n) is proved for every n" -- including the inference that HFC(2) is thereby
closed and the `int e^Q` route redundant, and the reading of kps-S166's residual
as already-settled -- is **suspended** until (G1)/(G2) are resolved. THM-3031's
bridge is unaffected: it is an implication about what would follow from the
exponential-integral claim, and does not depend on 4b.

**To close this**, one needs either (a) a proof that a mean-zero polynomial on
`Delta` cannot have its modulus maximum interior with `b^2 >= 4cA`, or (b) a
steepest-descent treatment giving the growth of `Phi` from the integral
representation without assuming the moments.
