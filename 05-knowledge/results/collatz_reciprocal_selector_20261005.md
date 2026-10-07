# Reversing the shell selector: a stable dyadic product with an integration boundary

2026-10-05. **PROVED:** coefficient reversal gives a stable entire product,
an exact reciprocal-coordinate shift law, and the trace-square recursion.
All negative shell moments of the actual leaf law diverge. **FINITE-EXACT:**
the finite identities and hostile controls below. **OPEN:** an independently
positive source measurement and universal Collatz coverage.

[Checker](../../04-computation/experiments/collatz_reciprocal_selector_20261005.py)
and [saved output](collatz_reciprocal_selector_20261005.out).

## 1. Inheritance and transfer contract

The nearest proved mechanism is the arithmetic-shell minorant in
[refinement energy duals](collatz_refinement_energy_dual_20261005.md), sections
2-5. Its uniform coefficient norm is below512, but the degree grows.
The canonical hostile is a zero-target shell law; the corrected near miss is
exchanging an infinite expansion with an expectation without paying its
integrability cost. The underused sidecar is the polynomial's degree and the
choice of expansion point. The companion
[denominator budget](collatz_shell_denominator_budget_20261005.md) independently
keeps the integer lattice lost by a small real norm.

Portfolio: **anchor** source-specific positive measurements; **niche** reverse
the finite selector instead of increasing its degree again; **wildcard**
reciprocal traces and dyadic products. Board: source / degree / norm /
arithmetic denominator / integration / normalization. The connection to
odd-zeta work is methodological: small linear forms need coefficient control
and nonvanishing separately. No zeta irrationality result is asserted here.

| Map | Preserved | Lost if omitted | Needed sidecar | Cheap test |
|---|---|---|---|---|
| Q_N(h) to z^N Q_N(1/z) | coefficients, norm, finite shell zeros | which end carries high degree | N and evaluation coordinate | independent factor multiplication |
| z to reciprocal pair alpha,1/alpha | z=(alpha+1/alpha+2)/4 | branch choice, harmless for symmetric product | nonzero alpha | exact rational factorization |
| finite product to entire limit | compact convergence and all shell zeros | uniform control at unbounded inverse coordinates | domain and integrability | explicit rooted ray |

## 2. Forward coefficients escape; reverse coefficients converge

Put x_d=4*2^d/(2^d+1)^2, d>=1, and

    Q_N(h)=product_(d=1)^N (h-x_d)/(1-x_d),
    R_N(z)=z^N Q_N(1/z)=product_(d=1)^N (1-x_d*z)/(1-x_d).

Here N may be any nonnegative integer. Only the odd-N Q_N are the inherited
minorants. Both coefficient l1 norms are below512, including even N by
comparison with the next odd degree.

**Proposition R1 — PROVED.** For every fixed k, the coefficient of h^k in Q_N
tends to zero as N tends to infinity, although Q_N(1)=1 always. More explicitly,

    |[h^k]Q_N| <= 512 binom(N,k) product_(d=1)^(N-k) x_d
               <= 512 binom(N,k) 2^(-(N-k)(N-k-3)/2).       (1)

Indeed the coefficient is the product of 1/(1-x_d) times the elementary
symmetric function of order N-k in the positive, decreasing x_d. Bound each
term by the product of the largest N-k entries. The prefactor is below512;
also x_d<=2^(2-d). For fixed k the quadratic exponent dominates the binomial.
Thus bounded l1 norm does not make evaluation at1 continuous for coordinatewise
coefficient convergence. Dropping N changes the represented function.

**Proposition R2 — PROVED.** R_N converges in coefficient l1 and uniformly
on every compact subset of the complex plane to the entire function

    R(z)=product_(d>=1) (1-x_d*z)/(1-x_d).                   (2)

It satisfies R(1)=1, and its zeros are exactly 1/x_d, each simple. For N>=5,

    ||R-R_N||_coefficient-l1 <= (65536/15)*2^(-N).           (3)

For d>=6, x_d<=1/16, and the l1 norm of the difference of the dth factor from1
is 2x_d/(1-x_d)<=(128/15)2^(-d). Telescope the products, using the norm bound512,
then sum the geometric tail. Compact convergence also follows directly from
sum x_d<infinity. Away from the displayed zeros the infinite numerator product
is nonzero; its denominator product is positive. This proves the zero claim.
Bound(3) controls the closed unit disk, not all inverse shell points.

## 3. The reciprocal trace gives an exact shift equation

Let q=1/2, and define (a;q)_N=product_(r=0)^(N-1)(1-a*q^r).
For alpha!=0 set z=(alpha+alpha^(-1)+2)/4. Direct factorization gives

    R_N(z)= (alpha*q;q)_N (alpha^(-1)*q;q)_N / (q;q)_N^2.   (4)

No modularity theorem is being imported. This is an identity between finite
products. If F(alpha)=R((alpha+alpha^(-1)+2)/4), shifting the index in the
convergent product proves

    (2*alpha-1) F(2*alpha)=2*alpha*(1-alpha) F(alpha),
    F(alpha)=F(alpha^(-1)),   F(1)=1.                       (5)

The finite version retains a terminal factor:

    (2*alpha-1)(1-alpha/2^N) F_N(2*alpha)
      =2*alpha*(1-alpha)(1-1/(2^(N+1)*alpha)) F_N(alpha).    (6)

These are cross-multiplied identities and remain valid at their vanishing
factors. Dividing by such a factor would lose information. For instance,
at alpha=1/2 equation(5) becomes0=0 and cannot recover F(1)=1. The zero
function satisfies the homogeneous shift equation too. Normalization must
therefore be retained; recurrence alone is not a positivity proof.

Write J=alpha+alpha^(-1)=4z-2. Squaring alpha gives

    J -> J^2-2,     z -> (2z-1)^2.

At the reciprocal shell coordinates this is exactly

    x_d -> x_(2d)=x_d^2/(2-x_d)^2.                         (7)

This is the precise common structure with the earlier quadratic trace
recursion: reciprocal multiplication becomes a polynomial operation on a
quotient coordinate. It does not identify dyadic shell distance with a
Collatz itinerary. A golden-coordinate example is z=9/4, where J=7 and the
reciprocal pair is phi^4,phi^(-4), since phi^4+phi^(-4)=7. This is an algebraic
specialization of the same map, not evidence for a source floor.

## 4. The pointwise exact selector cannot be integrated term by term

For a fixed source index m>=0, let h_m(j)=x_|j-m| with x_0=1. Let
p_j=lambda(6j+3)=W(6j+3) be the inherited normalized Collatz leaf law.
Universal positivity is not assumed. Equations(2) and the zero set prove

    R(1/h_m(j)) = 1 if j=m, and 0 otherwise.                (8)

Thus its pointwise expectation is exactly p_m. That fact alone is a
representation, not an independently computable positive measurement.

**Proposition R3 — PROVED.** For every fixed m>=0 and real s>0,

    sum_(j>=0) p_j h_m(j)^(-s)=infinity.                    (9)

Use the explicit rooted leaf ray

    n_t=(2^(6t)-1)/3,   j_t=(n_t-3)/6,    t>=1.

Each n_t is an odd multiple of3 and takes exactly one odd step to1, with
valuation6t. The established weight formula has L=0,K=3t-1, so

    p_(j_t)=2/((3t)(3t+1)).                                (10)

For large t, j_t>m and h_m(j_t)<=4*2^(m-j_t). Consequently the j_t summand
in(9) is at least

    [2/((3t)(3t+1))] * 2^(s*(j_t-m-2)),

which tends to infinity because j_t grows exponentially in t. In particular,
these nonnegative summands do not tend to zero. No assumption about any
unresolved orbit enters this argument.

The power series R(z)=sum c_k z^k has c_k nonzero for every k: its signs
alternate, and each absolute coefficient is a positive elementary symmetric
sum in the x_d, times a positive normalizer. Hence every k>=1 term would
require the divergent moment in(9). The formal expression

    E R(1/h) = sum_k c_k E h^(-k)

is invalid. It would combine infinite signed quantities; absolute integration
fails already in degree1. The zeros in(8) arise from cancellation performed
pointwise before integration. Coefficient l1 convergence near z=0 does not
license that exchange at the unbounded points1/h.

This obstruction is specific to inverse-power readout. It does not rule out
bounded rational functions, a different expansion, or independently proved
global signed inequalities. The ordinary bounded h^k moments remain valid.

## 5. A bounded carrier survives the pole cancellation

There is a constructive alternative to the divergent inverse-power readout.
Let F_m(z)=sum_(k>=0) H_k z^k be the ordinary-moment generating function,
and put T_m(z)=(1-z)F_m(z), initially in the unit disk. Its coefficients are

    t_0=1,     t_k=H_k-H_(k-1)<=0 for k>=1.

Since H_k decreases to p_m, telescoping gives

    ||T_m||_coefficient-l1=2-p_m<=2,
    sum_(k>M)|t_k|=H_M-p_m<=(8/9)^M.                     (11)

The companion denominator note proves that

    E_m(z)=R(z)T_m(z)

extends to an entire function: (1-z) removes the target pole, and R removes
every off-target pole. Its coefficient l1 norm is below1024 and
E_m(1)=p_m. Thus canceling the poles yields a stable ordinary-moment carrier,
without integrating any inverse moments.

Let T_[N] be the degree-N truncation of T_m. Equations(3),(11) give the
explicit finite approximation, of degree2N,

    ||E_m-R_N T_[N]||_coefficient-l1
      <= (131072/15)*2^(-N)+512*(8/9)^N,   N>=5.           (12)

This follows by adding (R-R_N)T_m and R_N(T_m-T_[N]) and using the l1
product inequality. It preserves the source but does not introduce a new
positive lower bound: at z=1 the approximant is exactly H_N. The elementary
one-moment interval H_N-(8/9)^N<=p_m<=H_N is already sharper there. Its
purpose is to identify a legal analytic carrier after the raw moments'
non-D-finite boundary, not to claim an improved source certificate.

## 6. Consequence and next target

We now have an exact compact recursive representation of all arithmetic-shell
zeros, with a paid approximation norm on a stated domain. The same operation
exposes its own failure boundary: normalization at a singular shift and
integration at inverse shell coordinates are indispensable extra data.

The productive next target is a bounded readout with either an independently
controlled tail or a source-sensitive global inequality. The companion
[spatial sidecar](collatz_specialization_sidecar_20261005.md) takes another
route: keep a fixed bounded moment degree and vary its source address.
Neither representation theorem supplies an all-source positive lower bound.

## 7. Reproduction and scope

    python 04-computation/experiments/collatz_reciprocal_selector_20261005.py
    python -O 04-computation/experiments/collatz_reciprocal_selector_20261005.py

There are2332 exact checks: independent factor multiplications and coefficient
reversal at degrees0..31; fixed-coordinate bounds; all finite reciprocal
shell zeros in that range; rational tail bounds; shift identities at degrees
0..18 over nine signed/singular rational reciprocal coordinates; doubled
shells through79; explicit rooted-ray weights through t20; and eight type
hostiles. No floats or assert-dependent checks are used. No infinite theorem
is inferred from these finite controls; the arguments above prove those
claims. Source21 is retained as the small ray boundary where the displayed
logarithmic growth control is not yet positive.
