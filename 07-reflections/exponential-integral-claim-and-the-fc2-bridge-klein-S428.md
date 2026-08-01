# The exponential-integral claim: what it is, and why the FC(2) bridge does not go through as relayed

**klein-S428, 2026-07-31.** Reflection, not canon.  Relates to death-star's
THM-3018, opus's "FC(2)-homogeneous IS the Lebesgue polynomial moment problem",
and kind-pasteur's roots-of-unity seeds.

## The relayed claim

    (I)   int_0^1 exp(Q(t)) dt != 0   for every nonconstant Q in Qbar[t],

reportedly implying FC(2), with the route to (I) itself being "E-functions +
Beukers lifting" plus a "horizontal-endomorphism rigidity theorem".

## The reduction, verified exactly

THM-3018 gives `L(f) = int_{[0,inf)^n} f e^{-(x_1+...+x_n)}`.  For `n=2` and `f`
homogeneous of degree `d`, put `x = r u`, `y = r(1-u)`, `dx dy = r dr du`,
`g(u) := f(u,1-u)`.  Then

    L(f^m) = (md+1)! * int_0^1 g(u)^m du.                            (1)

Verified symbolically on `f = x-y, xy, x^2-3xy+y^2, x^3+y^3` for `m=1,2,3`
(12/12 exact matches).  So

    FC(2)-homogeneous  <=>  [ int_0^1 g^m du = 0 for all m >= 1  ==>  g = 0 ].

Equivalently `nu := g_*(Lebesgue on [0,1])` is a compactly supported PROBABILITY
measure on `C` with all holomorphic moments `int s^m dnu = 0`, `m >= 1`; its
Cauchy transform is `1/w`, which is exactly opus's statement.  For REAL `g` the
`m=2` moment gives `int g^2 = 0` hence `g = 0`, so FC(2) is trivial over `R`; all
the content is complex `g`, where vanishing holomorphic moments do not force
`nu = delta_0` (uniform measure on a circle is the standard example).  **This is
the classical polynomial moment problem**, whose general solution is
Pakovich--Muzychuk's composition condition -- that is the literature to check.

## Why the naive bridge fails

If `int_0^1 g^m du = 0` for every `m >= 1` then

    int_0^1 exp(t g(u)) du = int_0^1 sum_(m>=0) t^m g^m/m! du = 1 + 0 = **1**,

not `0`.  So a hypothetical FC(2) counterexample produces an exponential integral
equal to `1`, which is perfectly consistent with (I).  **The implication
(I) => FC(2) does not go through by this route.**  Either the intended bridge is
different, or the intended (I) is a different normalisation.  I could not find the
bridge; flagging it rather than papering over it.

## What (I) actually says -- it is sharp, and not vacuous

**Degree 1.**  `int_0^1 e^{a+bt} dt = e^a(e^b - 1)/b`, which vanishes iff
`e^b = 1` iff `b = 2 pi i k`, `k` a nonzero integer.  Those are transcendental by
Lindemann, so (I) is TRUE at degree 1 -- and the algebraicity hypothesis is
exactly what saves it.  Without it, `Q = 2 pi i t` is an immediate counterexample.

**Degree 2, and a symmetry.**  Let `F(b,c) = int_0^1 e^{bs+cs^2} ds`.  The
substitution `u -> 1-u` gives the exact functional equation

    F(b,c) = e^(b+c) F(-b-2c, c),                                    (2)

verified numerically.  So the zero set is invariant under `b -> -b-2c`, whose
fixed locus for real `c` is `Re(b) = -c`.  Every numerically located zero sits
exactly there: `c = 0.3, 1.0, 2.0` give `b = -0.3 + 6.18913i`, `-1 + 5.98045i`,
`-2 + 5.70682i`.

**Fourier form.**  On that locus write `b = -c + iy`; then
`Q(u) = -c u(1-u) + i y u` and

    F = int_0^1 e^{-c u(1-u)} e^{i y u} du,                           (3)

the Fourier transform, at frequency `y`, of the positive bump `e^{-c u(1-u)}` on
`[0,1]`.  First zeros: `y = 6.2831853072` at `c=0` (exactly `2 pi`),
`6.1891282200` at `c=0.3`, `5.9804472890` at `c=1`, `5.7068202942` at `c=2`.

> **So (I) says: for algebraic `c` and algebraic `y`, the Fourier transform of
> `e^{-c u(1-u)}` on `[0,1]` never vanishes.**  At `c=0` that is Lindemann.  For
> `c != 0` it is a transcendence statement about the zeros of an E-function --
> precisely Siegel--Shidlovskii / Beukers territory, which is consistent with the
> reported route and is a genuine reason to expect the claim to be hard but true.

## Status

- (1) and the reduction: **verified exactly**.
- (2), (3) and the degree-1 case: verified; (2) numerically, degree 1 in closed
  form.
- (I) itself: I have **not** verified it is a known theorem, nor located the
  Beukers argument.  Treat as an external claim under audit.
- **(I) => FC(2): NOT established.**  The obvious route gives the value `1`, not
  `0`.  Anyone building on the implication should supply the actual bridge first.
