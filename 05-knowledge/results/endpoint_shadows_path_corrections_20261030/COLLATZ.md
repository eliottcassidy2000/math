# Negative-cycle shadows: exact endpoint phases and original-source tail bounds

2026-10-30. **PROVED elementary, self-audited; FINITE-EXACT controls.**
This bounds a specified family of actual prefixes without assuming that
their endpoints reach1. Universal ROOT coverage is **OPEN**.

## Inheritance and preserved quantities

The [October4 negative-shadow obstruction](../collatz_negative_cycle_shadow_20261004.md)
already proves the exact repeated(1,2) budget floor((v2(n+5)-1)/3), and
rules out bounded-depth checkpoint banks on specified mixed cells. The
[signed rational-anchor compiler](../rational_anchor_returns_20261004.md)
already retains denominator and exact-word guards. We use these mechanisms,
not a weaker finite-residue potential. In particular
[THM-4507, finite-valuation obstruction](../../../01-canon/theorems/THM-4507-finite-valuation-collatz-potential-obstruction.md)
already excludes arbitrary nonlinear corrections to logarithmic height
from finitely many polynomial valuations.

The [current signed-cutoff packet](../source_cutoffs_balanced_lifts_20261029/CUTOFFS.md)
uses positive odd sources n>=3 with their ORIGINAL weights

    w(n)=4/((n+1)(n+3)),        sum_(odd n>=3) w(n)=1/2.

We keep these weights at transported endpoints. We never replace w(n) by
w(U^t(n)) or by uniform mass on residue classes. The new calculation records
the reduced rational anchor, its finite residue phase, and endpoint tail
moments together. The fixed-source upper-cutoff rule in the incoming
[native survival note](../session106_collatz_transfer_20261029/recurrence/NATIVE_UPPER_CUTOFF.md)
is complementary: it compresses first-descent inequalities; the present
rule measures the supplied shadow cylinder independently of descent.

## C1. Exact source and endpoint lattices

Fix a nonempty positive valuation word u=(a1,...,ad), with A=sum a_i,
P=3^d, Q=2^A, and affine map F_u(x)=(Px+B)/Q. Assume P>Q. Write its
negative fixed point in lowest terms

    c=-b/D=-B/(P-Q),        b,D positive odd, gcd(b,D)=1.

Every such word has this fixed point with its exact 2-adic valuation
itinerary. To see exactness, iterate its inverse letters (2^a x-1)/3
backwards from the odd fixed endpoint: every intermediate value is odd
2-adically. The unique full-word source cylinder is therefore c modulo2Q.
The same argument for u repeated j times gives, with M=2Q^j,

    S_j={r_j+M t:t>=0},     r_j=-b D^(-1) mod M,  1<=r_j<M.       (C1)

These are actual positive odd sources. None is1: the map at1 with its
actual letters remains1, whereas the displayed expanding affine map sends
every positive source strictly above itself after a full block.

If M>b, there is a unique h_j in {1,...,D} with

    h_j M=b mod D,
    r_j=(h_j M-b)/D,
    U^(dj)(r_j+Mt)=(2h_j P^j-b)/D + 2P^j t.                    (C2)

For D=1 take h_j=1. The numerator is divisible by D because P=Q mod D.
The phase h_j is periodic in j, with period ord_D(Q) when D>1. Thus
the independent least-endpoint bound is

    U^(dj)(n)>=(2P^j-b)/D,              n in S_j.               (C3)

There is no density argument in (C2): it retains the actual source address.
More generally Dn+b is divisible by2Q^j exactly when n follows u^j, so
the number of repeated full blocks at a fixed positive source is exactly
floor((v2(Dn+b)-1)/A). This is the general anchor version of the inherited
consumed-precision budget, not a ROOT certificate after the mismatch.

## C2. Quadratic decay of original-source mass

Write W_j=sum_(n in S_j) w(n). Once M>=2b,

    4/((M+1)(M+3)) <= W_j <= (16D²+8D)/M².                    (C4)

For the lower bound keep the first term and use r_j<M. For the upper bound,
r_j>=(M-b)/D>=M/(2D), and the decreasing-function sum estimate gives

    W_j <= 4/(r_j+1)² + 4/[M(r_j+1)].                         (C5)

The second term follows by integrating4/(r_j+Mx+1)² over x>=0.
This proves (C4). The decay is quadratic in M for a fixed anchor, although
the odd-relative Haar density of S_j is only2/M. The additional power is
caused by its least positive representative escaping proportionally to M.

The exact phase constants are also visible. Along any subsequence on which
h_j=h, dominated convergence gives

    M² W_j -> 4 sum_(t>=0) (t+h/D)^(-2).                       (C6)

The right-hand side is4 times a Hurwitz-zeta value. This is an ordinary
convergent series identity, not a claim about irrationality of its values.
For integer negative anchors D=1 it becomes4 zeta(2). Distinct phases of
a nonintegral anchor can have different limits: discarding h loses the
leading constant, even though it leaves the coarse exponent unchanged.

## C3. Transported endpoint tails, and the first-moment boundary

Set alpha=h_j/D and L=2P^j. In the source parametrization of (C2),

    n_t+b/D=M(t+alpha),       y_t+b/D=L(t+alpha).

For M>=2b, n_t>=M(t+alpha)/2, and hence
w(n_t)<=16/[M²(t+alpha)²]. Let nu_j give mass w(n_t) to the actual
endpoint y_t. Then for every Y>0,

    Y nu_j{y:y+b/D>=Y}
      <=16(D+1)L/M²
       =8(D+1)(P/Q²)^j.                                      (C7)

Indeed the first retained index t0 has v=t0+alpha>=max(alpha,Y/L),
and its reciprocal-square tail is at most v^-2+v^-1. Multiply by Y,
use Y<=Lv and v^-1<=D. No endpoint multiplicity or weight is resampled.

More generally for 0<=beta<1,

    sum_(n in S_j) w(n)(U^(dj)(n)+b/D)^beta
      <=4*2^beta [D^(2-beta)+D^(1-beta)/(1-beta)]
                   * (P^beta/Q²)^j.                          (C8)

Integrate(t+alpha)^(beta-2) to obtain this inequality. The base of this
exponential is less than1: Q>=2^d, so P^beta/Q²<=(3^beta/4)^d.
Along each fixed phase the sharper scaled limit is

    (M²/L^beta) sum w(n)(y+b/D)^beta
        ->4 sum_(t>=0)(t+alpha)^(beta-2).                      (C9)

For beta>=1 that moment is infinite at EVERY finite j, since its summands
are asymptotic to a positive constant times t^(beta-2). In particular an
untruncated mean endpoint is not a suitable observable for these original
weights. A fractional moment, a weak-L1 tail, or a source cutoff is needed.
Growing individual endpoints and shrinking source-weighted packet mass
are compatible; neither statement by itself establishes ROOT entry.

## C4. Controls, limits, and the next independent target

| word | negative anchor | P/Q | phase period |
|---|---:|---:|---:|
|(1)|-1|3/2|1|
|(1,2)|-5|9/8|1|
|(1,1,2)|-19/11|27/16|5|
|(1,1,1,2)|-65/49|81/32|21|
|(1,1,1,2,1,1,4)|-17|2187/2048|1|

For -19/11 the phases start3,5,1,9,4 and repeat. For -65/49, (C2)'s
simple phase range starts once M>b; the first j need not satisfy that guard.
The appearance of5 and21 here is the multiplicative order of the actual
denominator clock, not a numerical identification with a tournament class.

Two hostile controls delimit the result. First, the contracting word(2)
has positive anchor1; its cylinders retain a root atom of mass1/2 if the
terminal is included. That atom does not vanish as the modulus grows.
The proof of (C4) uses the NEGATIVE anchor and cannot be carried across this
boundary unchanged. Second, do reduce B/(P-Q): replacing a word by its own
powers keeps its anchor but makes the unreduced denominator arbitrarily
large. A bound in that unreduced coordinate would charge artificial debt.

The fixed-word bounds are independent of convergence, but are not uniform
over a changing alphabet of anchors. D, b, phase, and the threshold M>=2b
all change. We have not proved that a positive orbit must return to one
of these charts, or that its post-mismatch words have summable costs.

A precise next target is a source-preserving stopping-prefix cover of the
remaining killed-ROOT residual, with a summable phase-sensitive budget of
the right-hand sides of (C5) or (C8). Any claimed cover must retain every
marked surviving source and count its mass only once. Union bounds suffice
for overlapping packets, but do not create a small budget. A finite catalogue
of repeated words solves only its own subfamily. The earlier signed-dual
obstruction also remains: optimizing an uninformative endpoint bound of1
does not manufacture a new terminal certificate.

## C5. The changing-word obstruction has an arithmetic bill

There is an explicit source-dependent reason the fixed-word decay cannot
simply be applied to longer and longer DIFFERENT words at one source n.
For any actual expanding prefix, with reduced anchor -b/D and M=2Q,

    Dn+b is a positive multiple of M;
    max(D,b)>=(2Q)/(n+1).                                     (C10)

Thus its reduced arithmetic height must grow at least proportionally to
the consumed binary modulus at a fixed source. This is just the exact
cylinder guard, followed by Dn+b<=max(D,b)(n+1). It neither assumes a
divergent source exists nor rules one out. If the real anchor is bounded
in magnitude by C, it specifically forces D>=2Q/(n+C).

On the guard M>=2b used in (C4), the same identity gives D>=M/(2n).
Consequently the upper-bound expression (16D²+8D)/M² cannot even drop
below4/n² as long as that packet still contains n. The fixed-source floor
w(n) is therefore not silently erased by the quadratic decay formula.
When b is too large, it is instead the entrance threshold M>=2b that fails.

This is the precise rerouting obligation: changing the controller must
pay reduced denominator growth or real-anchor growth. Cancelling a common
factor in B and P-Q avoids fictitious debt from repeated copies of one word,
but cannot remove the actual bill (C10). A successful ROOT argument must
use a further native mechanism to discharge it; norm cancellation or a
Haar-density estimate cannot by itself do so.
