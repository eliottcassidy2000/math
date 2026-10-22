# Reciprocal-residue inverse clocks: distinguish the greedy rule

Status: analytic formulas PROVED below; later numerical indices are VERIFIED, not certified integers unless separately enclosed by the interval packet. The original supplied prompt does not define a unique greedy convention.

Put A(x)=(2/pi)x^{-1/x}sin(pi/x) for x>=3. It is positive, strictly decreasing to0. Its logarithmic derivative is

A'(x)/A(x)=[log x−1−pi cot(pi/x)]/x².

For x>=3, cos u>=1−u²/2 and sin u<=u give pi cot(pi/x)>=x−pi²/(2x)>log x−1. The last difference is positive at3 and increasing. Also pi cot(pi/x)<x and log x>=1, whence

0<−A'(x)<=2/x².

## Three distinct algorithms

For a residual r in (0,A(3)), let x(r) be the unique real inverse and m=floor x(r).

- Over-residue: choose m, so A(m)>=r.
- Under-residue: choose m+1, so A(m+1)<=r (with the obvious exact-integer exception).
- Nearest residue: choose whichever of m,m+1 minimizes |r−A(M)|.

Nearest in A-value is NOT always rounding x to the nearest integer, because A is nonlinear. A tie requires an explicit convention. Signed reconstruction tracks the sign of r−A(M) at every step; replacing a signed residual with its absolute value without recording that sign changes the represented constant.

For r=A(8)−C_MRB, the inverse is approximately

x(r)=105719003.855602567339775184066853029452834470099887...

Thus 105719004 is the under-residue and the nearest choice, not the over-residue choice. Numerically:

A(105719003)−r = +1.5310726756958271062715295939e−16,
A(105719004)−r = −2.5839445781801954275830241249e−17.

This distinction is material when describing how the later 'greedy' corrections are generated.

## Inverse expansion

As r decreases to0, put X=2/r and L=log X. Direct expansion gives

A(x)=2/x−2log(x)/x²+[(log x)²−pi²/3]/x³+O((log x)³/x⁴),

X=x+log x+[(log x)²/2+pi²/6]/x+O((log x)³/x²),

x(r)=X−L+[L−L²/2−pi²/6]/X+O(L³/X²).

These are inverse asymptotics, not rounded-integer certificates. Near an integer, the displayed remainder must be controlled before flooring/ceiling.

## Quadratic debt and integer clocks

The nearest residual r'=min(r−A(m+1),A(m)−r) obeys

0<=r' <= [A(m)−A(m+1)]/2 <=1/m².

Consequently r' <= r²/4 [1+O(r log(1/r))]. A zero residual terminates the process; otherwise the small residual creates a quadratically larger next inverse scale. This is a deterministic analytic bound with a rounding coordinate, not a source-supply theorem of the kind required for Collatz coverage.

More explicitly, let theta=x−floor x. Uniformly as x tends to infinity,

r' = (2/x²) min(theta,1−theta)+O(log(x)/x³).

The leading approximation is additive. It must not be divided by min(theta,1−theta) near an integer without a separate lower bound. Extremely large next indices come from this fractional-part clock, not from the fourth nome power alone.

For a coarse all-step growth bound, c=2 exp(−1/e) sin(pi/3)/(pi/3)>0 gives A(x)>=c/x for x>=3. Thus x(r')>=c/r'>=c m²; the next selected integer is at least c m²−1. Tracking the chosen M rather than floor x gives the safe inequality M_next>=c(M−1)²−1. These bounds apply when r'>0.

## Connection and boundary

The inverse series and the mirror tree function both require branch and discrete-index data. This is a valid shared mechanism, not identification of the functions: mirror inversion solves log r/r=z, while residue inversion additionally carries a sine factor. Integer clocks explain the correction cascade once the original defect is known; they do not explain why that defect has its original size.

## Independent integer certificate upgrade

`analytic-audit/AUDIT.md` and its interval script independently certify the first four nearest-amplitude indices by strict neighboring midpoint inequalities:

8,
105719004,
77401040908104410,
24845630402060641702669818431121505.

All four integer decisions, rather than only the leading defect digits, are certified. The companion numerical script is retained as a separate numerical reproduction. This cascade computes consequences of the initial defect; it still does not derive that defect from CM data.
