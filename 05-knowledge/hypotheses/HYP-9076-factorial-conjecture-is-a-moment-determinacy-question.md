---
id: HYP-9076
title: "The Factorial Conjecture is a moment-determinacy question, and the difficulty is indexed by DEGREE"
status: >
  FRAMING + CENSUS. L is the moment functional of the PRODUCT EXPONENTIAL
  measure, so FC(n) says: if every moment (m>=1) of the pushforward
  nu = f_* mu vanishes then f = 0 -- exactly a determinacy question. The
  Krein/Carleman criterion then locates the difficulty by DEGREE, not
  dimension: a Stieltjes density decaying like exp(-t^rho) is determinate iff
  rho >= 1/2, and the pushforward of e^{-|x|} under a degree-d polynomial has
  tail exp(-t^{1/d}), so the argument closes for deg f <= 2 in EVERY number of
  variables and fails from cubics on. Census: zero survivors among all small
  integer f of degree <= 2 in 2 and 3 variables with the first 3-4 moments
  vanishing -- but see section 6: ALL of those censuses are VACUOUS, because
  L is positive definite and so L(f^2) > 0 for every real f != 0. The real
  case of FC is trivial; the conjecture lives entirely over C, and the
  structural reduction is to pairs (g,h) of real polynomials that are
  L-mean-zero, ORTHOGONAL and EQUINORMAL in the Gram form <u,v> = L(uv),
  with the first genuine content at m = 3.
source: opus-2026-07-31-amm12592-writeup
related:
  - THM-2922  # first-window SFC(4), Macaulay-Newton atlas
  - THM-2891
  - HYP-9075  # the shell-imbalance module, same moment-vanishing shape
external:
  - "van den Essen, Wright, Zhao, On the image conjecture (Factorial Conjecture)."
script: 04-computation/fc_factorial_conjecture_census.py
output: 05-knowledge/results/fc_factorial_conjecture_census.out
---

# HYP-9076 -- FC as determinacy

## 1. The framing

With `L(x^alpha) = alpha!`,

```text
L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx,
```

so `L` is the moment functional of the product exponential measure `mu`, and
`L(f^m)` is the `m`-th moment of the pushforward `nu = f_* mu`. FC(n) reads:

```text
all moments (m >= 1) of nu vanish   =>   f = 0.
```

If `nu` is DETERMINATE this is immediate: vanishing moments plus total mass
force `nu = delta_0`, hence `f = 0` `mu`-a.e., hence `f = 0` as a polynomial.
**So FC is exactly a moment-determinacy statement.**

## 2. The difficulty is indexed by degree, not dimension

Krein/Carleman: a Stieltjes density decaying like `exp(-t^rho)` is determinate
iff `rho >= 1/2`. The pushforward of `e^{-|x|}` under a degree-`d` polynomial
has tail `exp(-t^{1/d})`, so

```text
deg f <= 2 :  rho = 1/d >= 1/2  ->  DETERMINATE  ->  FC closes
deg f >= 3 :  rho < 1/2         ->  INDETERMINATE -> the argument fails
```

in EVERY number of variables. This is worth stating because FC is indexed by
`n` in the literature and by supports/slots in the repo's SFC work, while the
determinacy obstruction is indexed by `d`. Caveat: for complex `f` the
pushforward is a complex measure and the criterion needs the real/imaginary
split; the clean statement is for real `f`.

## 3. Census

Small integer `f` with `L(f^m) = 0` for `m = 1..M`:

```text
n=2, deg<=2, coeffs in [-2,2], M=3 :  0 survivors / 3124
n=2, deg<=2, coeffs in [-3,3], M=4 :  0 survivors / 16806
n=3, deg<=1, coeffs in [-3,3], M=4 :  0 survivors / 342
n=3, deg<=2, coeffs in [-1,1], M=3 :  0 survivors / 19682
```

Nothing survives even the first three moments in degree `<= 2`, consistent
with (and much stronger than) the determinacy argument, which only says the
infinite family of moments suffices.

## 4. Degree three: still nothing

Pushed into the regime where the determinacy argument stops working:

```text
homogeneous cubic, n=2, coeffs in [-6,6], m=1..5 :  0 / 28560
homogeneous cubic, n=2, coeffs in [-9,9], m=1..4 :  0 / 130320
two-term cubics, coeffs in [-40,40], only L(f)=L(f^2)=0 : 0
n=2, deg <= 3 (all 9 monomials), coeffs in [-1,1], m=1..4 : 0 / 19682
n=2, deg <= 3 (all 9 monomials), coeffs in [-2,2], m=1..3 : 0 / 1953124
n=3, deg <= 3 : NOT COMPLETED -- 19 monomials, ~1.2e9 vectors, cut by timeout
```

The `n = 2` inhomogeneous sweep over `1953124` vectors is the substantial one:
allowing every monomial of total degree `<= 3` and coefficients up to `+-2`
still kills everything by the third moment.

The two-term case is the sharpest: `L(f) = 0` already pins the coefficient
ratio to a fixed rational `-c!d!/(a!b!)`, and the quadratic `L(f^2)` then has
no solution at all in that range. So degree three is if anything MORE rigid
in these windows, not less.

**Logical care.** This does not conflict with section 2. The determinacy
argument failing for `d >= 3` says we lose a PROOF of FC there, not that
counterexamples appear. A census finding nothing is exactly what one expects
if FC is true and merely harder to prove past quadratics. Larger degree, wider
coefficients and deeper windows all remain untested, as does the repo's own
regime (THM-2922's four-slot supports of diameter five).

## 5. Relation to HYP-9075

The shell-imbalance module has the same shape -- a signed object whose moments
against a PRODUCT measure vanish -- with Bernoulli in place of exponential and
the family indexed by `p` rather than by powers. The difference is that our
Bernstein basis makes the determinacy step trivial (`D_m = 0` iff the
composition vector vanishes), so all our difficulty is in the integer box,
whereas FC has no box and all its difficulty is in the determinacy. The two
problems are complementary in exactly the way HYP-9075 sec 8 found for
THM-2922's machinery.


## 6. CORRECTION: the censuses were vacuous, and the real reduction

`L(g) = int_{[0,inf)^n} g(x) e^{-(x_1+...+x_n)} dx` is a POSITIVE functional:
`L(g) > 0` whenever `g >= 0`, `g != 0`. Hence for **real** `f`,

```text
L(f^2) = int f(x)^2 e^{-|x|} dx > 0    unless f = 0.
```

Verified directly: the Gram matrix `<u,v> = L(uv) = sum (alpha+beta)!` is
positive definite (least eigenvalue `0.0253` at `n=2, deg<=3`; `0.0150` at
`n=3, deg<=3`).

**Every census in sections 3 and 4 imposed `m = 1..M` with `M >= 2` over
INTEGER, hence real, coefficients. Their emptiness was therefore forced by
positivity alone and carries NO information about FC.** They are marked
vacuous. This includes the degree-3 runs and the `1953124`-vector sweep that
were reported as substantial; they were not.

### The reduction

Write `f = g + i h` with `g, h` real. With `<u,v> := L(uv)`:

```text
L(f)   = 0  <=>  L(g) = 0  and  L(h) = 0
L(f^2) = 0  <=>  <g,g> = <h,h>   and   <g,h> = 0
L(f^3) = 0  <=>  L(g^3) = 3L(g h^2)   and   3L(g^2 h) = L(h^3)
```

So FC is a question of EUCLIDEAN GEOMETRY in `(R[x]_{<=d}, <,>_L)`: seek `g`
orthogonal to `h`, of equal norm, both `L`-mean zero, satisfying in addition
the cubic identities and their analogues at every `m`. The `m <= 2` level is
just "an orthogonal equinormal mean-zero pair", which is trivially non-empty
-- an explicit pair in `n = 3, deg <= 3` is produced by the referee, giving a
nonzero `f` with `L(f) = L(f^2) = 0`.

**Consequence for the n = 3, degree 3 problem.** The obstacle was never the
size of the box. Searching real coefficient vectors answers a trivial
question; the meaningful search space is the variety of orthogonal equinormal
mean-zero PAIRS, on which `m = 1, 2` are already satisfied by construction and
the first real conditions are the two cubic identities. That is the structural
reduction, and it also explains why every census so far returned zero: they
were all measuring positivity, not FC.

Referee: `04-computation/fc_real_case_is_trivial_and_the_complex_reduction.py`.


## 7. The m = 3 search: solvable, and no finite window can ever obstruct

Searching the variety of section 6 directly (least-squares on the real and
imaginary parts of `L(f^m)`, normalised by `||f||_G = 1` to kill the scale
symmetry and exclude `f = 0`), for `n = 2`, `deg <= 3`:

```text
M = 2 :  max|residual| = 7.2e-16    ||f||^2 = 1.0000    SOLVED
M = 3 :  max|residual| = 1.1e-14    ||f||^2 = 1.0000    SOLVED
M = 4 :  max|residual| = 1.9e-12    ||f||^2 = 1.0000    SOLVED
```

**So the `m = 3` conditions ARE satisfiable on the orthogonal / equinormal /
mean-zero variety** -- there is a nonzero complex `f` in two variables of
degree `<= 3` with `L(f) = L(f^2) = L(f^3) = 0`, found to machine precision.

### The dimension count that explains it

Each `L(f^m) = 0` is ONE COMPLEX equation, i.e. two real ones, so `m = 1..M`
imposes `2M` real conditions on the `2N` real parameters of `f`
(`N` = number of monomials). The symmetries `f -> lambda e^{i theta} f`
remove two more, so solutions are expected whenever

```text
2M < 2N - 2,     i.e.     M < N - 1.
```

For `n = 2, deg <= 3` that is `N = 10`, so every window `M <= 8` should be
solvable; for `n = 3, deg <= 3`, `N = 20` and every `M <= 18`.

**Consequence: no finite window can ever obstruct FC.** Every finite
truncation is under-determined, so a census over any window will find
solutions once it is allowed to search `C` rather than `R`. This closes the
loop on section 6: the earlier censuses returned zero not because finite
windows are restrictive but because they were confined to the reals, where
positivity forbids everything at `m = 2`. Over `C` the opposite is true --
finite windows are permissive, and FC is genuinely a statement about the
infinite family.

Referee: `04-computation/fc_orthogonal_equinormal_pair_search.py`.
