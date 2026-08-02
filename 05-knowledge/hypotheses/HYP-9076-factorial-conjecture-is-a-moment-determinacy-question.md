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
  - THM-2922  # first-window SlotSFC_1(4), Macaulay-Newton atlas
  - THM-2891
  - HYP-9075  # the shell-imbalance module, same moment-vanishing shape
external:
  - "van den Essen, Wright, Zhao, On the image conjecture (Factorial Conjecture)."
script: 04-computation/fc_factorial_conjecture_census.py
output: 05-knowledge/results/fc_factorial_conjecture_census.out
---

# HYP-9076 -- FC as determinacy

**Notation correction (MISTAKE-350).** `FC(n)` below uses the original
ambient-variable index.  Every `SlotSFC_1(N)` reference means an
`N`-monomial restriction inside ambient `SFC(1)`, not `SFC(N)`.

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
2M < 2N - 2,     i.e.     M < N - 1.        <-- WRONG, see section 8
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


## 8. CORRECTION to the count, and what it explains

The count in section 7 subtracted two for the symmetries `f -> lambda
e^{i theta} f`. That is wrong: the normalisation `||f|| = 1` already consumes
the scale freedom, and the phase symmetry does not remove an equation -- it
only says the solution set, when nonempty, contains circles. The correct
count is `2M` moment equations plus `1` normalisation against `2N` real
unknowns, so the window is solvable exactly when

```text
2M + 1 <= 2N,     i.e.     M <= N - 1,      failing first at M = N.
```

Tested in ONE variable, where FC is known true and the computation is cheap
(residual `~1e-12` or better when solvable, `~1e-2` when not -- a clean
transition, not a numerical grey zone):

```text
deg<=2, N=3 :  solvable M = 1,2      fails at M = 3 = N
deg<=3, N=4 :  solvable M = 1,2,3    fails at M = 4 = N
deg<=4, N=5 :  solvable M = 1,...,4  fails at M = 5 = N
```

Three independent confirmations of `M <= N-1`, and three refutations of the
`M < N-1` of section 7.

### This explains the repo's "first window"

The consequence is sharper than "no finite window obstructs". For a polynomial
with a FIXED support of `N` monomials, only about `N` moments can be killed at
all: `L(f^m) = 0` for `m = 1..N` already forces `f = 0` generically. **So FC
restricted to a fixed support is a FINITE problem, of window length exactly
the support size.**

That is precisely the shape of the repo's SFC results. THM-2922 proves
first-window `SlotSFC_1(4)` for FOUR-slot supports using exactly the four moments
`L(H), L(H^2), L(H^3), L(H^4)` -- `N = 4` slots, `M = 4` moments. The
dimension count above says `M = N` is exactly where the system first becomes
over-determined, so the "first window" is not a convenient truncation: it is
the threshold, and THM-2922 is proving the statement at the first `M` where it
could possibly be true.

Referee: `04-computation/fc_first_window_threshold.py`.


## 9. The threshold at N = 2, PROVED (not generic) -- and why it stops there

**Theorem.** Let `f = c_1 x^alpha + c_2 x^beta` with `alpha != beta`
multi-indices in any number of variables and `(c_1,c_2) != 0` complex. Then
`L(f)` and `L(f^2)` cannot both vanish.

*Proof.* Put `A = alpha!`, `B = beta!`, `C = (alpha+beta)!`, `D = (2alpha)!`,
`E = (2beta)!` (multi-index factorials). Then `L(f) = c_1 A + c_2 B` and
`L(f^2) = c_1^2 D + 2 c_1 c_2 C + c_2^2 E`. If `L(f) = 0` then
`c_2 = -c_1 A/B` and `c_1 != 0`, so

```text
L(f^2) = (c_1^2 / B^2) * Q,     Q := D B^2 - 2 A B C + A^2 E.
```

Two classical facts finish it. Coordinatewise **log-convexity of Gamma** gives
`((alpha_i+beta_i)!)^2 <= (2alpha_i)!(2beta_i)!` with equality iff
`alpha_i = beta_i`; multiplying, `C^2 <= D E`, strict unless `alpha = beta`.
**AM-GM** gives `D B^2 + A^2 E >= 2 A B sqrt(D E)`. Hence

```text
Q >= 2 A B ( sqrt(D E) - C ) > 0     when alpha != beta,
```

so `L(f^2) != 0`. QED

Verified on `29112` multi-index pairs in 1, 2 and 3 variables: no violation,
and `Q = 0` exactly in the excluded case `alpha = beta`.

### Why the argument does not extend past N = 2

After using `L(f) = 0` to eliminate one coefficient, `N` slots leave `N-1`
complex unknowns. At `N = 2` that is ONE unknown, and `L(f^2)` becomes
`c_1^2` times a positive real number -- so it vanishes only at `c_1 = 0`. From
`N = 3` on there are at least two complex unknowns, and **a complex quadratic
form in two or more variables always has nontrivial zeros** (a conic in
`P^1` has points). So `L(f^2) = 0` alone is satisfiable for every `N >= 3`,
the positivity of `L` gives nothing further, and the higher moments
`L(f^3), ..., L(f^N)` must be brought in and eliminated.

That is exactly why the repo's program is case-by-case with Macaulay
resultants rather than a single soft argument: the transition at `N = 3` is
real, not a gap in effort. **A general rigorous proof of the `M = N` threshold
is therefore not available here, and should not be expected from positivity
alone; `N = 2` is the base case and the elimination-theoretic work
(THM-2824/2917 at three slots, THM-2849/2922 at four) is what the rest costs.**

Referee: `04-computation/fc_two_slot_threshold_proof.py`.


## 10. The proof's engine is a Newton ratio -- the THM-3010 connection

The discriminant in section 9 is not a new object. klein's THM-3010 studies
the **Newton ratio** `R_k = h_k^2 / (h_(k-1) h_(k+1))` of an integer sequence,
log-concave exactly where `R_k > 1`. The quantity controlling the two-slot
proof is

```text
C^2 / (D E) = (alpha+beta)!^2 / ((2alpha)! (2beta)!),
```

which is precisely a Newton ratio of the factorial sequence taken at the
midpoint of the segment `[2alpha, 2beta]`. With `h_j = j!` one gets
`R_n = (n!)^2/((n-1)!(n+1)!) = n/(n+1) < 1`, i.e. **`j!` is log-CONVEX**, and
that single fact is what makes `Q > 0` and kills the two-slot case.

So two independently developed threads use the same engine: THM-3010 computes
exact rational Newton ratios `R_k = 1 - Q_(a,b)(k)/D_(a,b)(k)` for BALLOT
columns, and the factorial-moment program needs Newton ratios of the FACTORIAL
sequence. A grep confirms log-convexity appears nowhere else in canon --
THM-3010 and this entry are the only two occurrences -- so the tool is
essentially unused in the repo despite being decisive here.

**Concrete transfer to try.** THM-3010 gives `R_k` in closed rational form for
`binom(2k,k)`, Catalan, `binom(2k,k-1)` and `binom(2k+1,k-1)`. Supports whose
slots are built from those columns would then have *exactly computable*
two-slot discriminants, and the `Q > 0` step would become a rational-function
positivity check in `k` rather than an inequality -- which is the same shape as
the Gregory--Newton certificates THM-2922 already uses. That is a concrete,
untried bridge between the two programs.


## 11. CORRECTION: my dimension count is MISTAKE-246, and the lower half is THM-2173

An extensive repo search settles the status of the `M = N` threshold, and it
is not what section 8 claimed.

### The counting argument is an already-quarantined error

MISTAKE-246 (audit of THM-1790) records exactly this reasoning and rejects it:

> "`d` equations in `d` affine variables may cut out the empty set or only a
> forbidden point; equation counting alone supplies no solution."

Section 8 argued "solvable iff `2M+1 <= 2N`" by counting equations against
parameters. **The forward direction is true but not for that reason, and the
reverse direction does not follow at all.** Equation counting cannot show that
`M = N` forces `f = 0`; that inference is precisely the retracted step.

### The lower half is already a theorem, by Krull height

**THM-2173** (sparse projective factorial moment floor) proves, for every
`t >= 2` and every `t`-dimensional envelope `V ⊂ C[s]`, that there is a
nonzero `H in V` with `L(H^j) = 0` for `1 <= j <= t-1`. The proof is Krull's
height theorem: if the `t-1` homogeneous moment equations had only the origin
as common zero, their radical would have height `t`, which is impossible.

So the "solvable for `M <= N-1`" half of my numerics is **already proved,
rigorously and for all `N`** -- my least-squares runs rediscovered THM-2173.
Its own scope note is apt: *"a dimension obstruction, not a factorial
peculiarity."*

### The upper half is the open program

There is **no** general "`N` slots need `N` moments" theorem in the repo, and
no hypothesis file conjectures one; the phrase occurs only in my own script's
docstring. The upper half is strictly case-by-case: `N = 2` external
(Edo--van den Essen) and section 9 here; `N = 3` arbitrary support (THM-2824);
`N = 4` only in pieces (THM-2908/2940 consecutive, THM-2921 diameter 4,
THM-2922 diameter 5, THM-2849 bounded box). Every `N = 4` theorem explicitly
disclaims SFC(5), and no SFC(5) file exists.

**Net status of the request.** The threshold `M = N` is NOT proved rigorously
in general, cannot be by counting, and the honest decomposition is:
lower half `M <= N-1` = THM-2173 (proved, Krull height); upper half `M = N`
= the open SFC program, with `N = 2` closed by section 9.

## 12. What the search says about section 9 and 10

* The two-slot proof of section 9 is, per the survey, **the strongest
  uncanonized asset in this lane** -- rigorous, all dimensions, arbitrary
  distinct multi-indices -- and it is cited by nothing. THM-3019 and THM-3020
  re-derive two-slot results by unrelated resultant and ODE methods and would
  be shortened by it. It should be promoted to canon.
* Its inequality `C^2 <= D E` is a **Turan-type inequality** for the factorial
  sequence, and the survey confirms Turan inequalities are **entirely absent**
  from the repo -- every `Turan` hit is extremal graph theory or
  Erdos--Turan discrepancy. That is a naming/vocabulary gap worth closing,
  since the Gamma-Turan literature is exactly the right toolbox.
* The THM-3010 Newton-ratio identification of section 10 stands, and the
  survey confirms log-convexity appears in only two other places repo-wide.

### A free certificate at two slots, in the THM-2922 style

For a TRANSLATED two-slot support `alpha = n+a`, `beta = n+b` in one variable,
put `d = b-a >= 1` and `N_0 = 2n+2a`. Then the discriminant ratio is

```text
rho(n) = (2n+a+b)!^2 / ((2n+2a)!(2n+2b)!)
       = prod_(i=1)^(d) (N_0+i)/(N_0+d+i)   <  1,
```

each factor being `< 1`. Verified exactly for several `(a,b)` and `n`.
Positivity therefore holds for ALL translations at once with no resultant and
no case analysis -- the same kind of statement THM-2922 needs a 197-term
Gregory--Newton certificate for at four slots, free at two.
