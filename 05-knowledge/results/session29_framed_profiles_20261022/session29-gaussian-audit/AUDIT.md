# Independent audit: rank-88 scalar critical-point certificate

Verdict: ACCEPTED, conditional on the stated Session28 harmonic-theta/scalar reduction. This is an audit of the universal scalar function, not of the existence, perfection, or angular spectrum of any specific lattice. No Gaussian-energy local/global optimality claim is certified.

## Frozen inputs and independent outputs

Producer script SHA256: d55e4e5527bd7e68b01cc806b12b1d3f9db6b4ce878792372da20a4268479915.

Producer report SHA256 at review: 6f6cd32eeee5d865b382b7ce9aafd51ca0c81bfa290d90aa0dadcb5ff30fa570.

Producer JSON SHA256 at review: 59aa8f1c31a2edbba374b2e73aa680d4ccb17a1a0f2bbbd71606c6514e733621.

`producer_frozen.py` is the exact audited source. Its independent full replay is recorded in `replay.out` and `producer_frozen.json`. `audit.py` and `audit.json` supply additional independent arithmetic and algebraic controls. The producer JSON includes a runtime, so byte identity of independent output is not expected.

## What was checked independently

1. Reconstructed Delta through degree84 using q times the Euler product of (1-q^n)^24, instead of the producer's (E4^3-E6^2)/1728 identity. Reconstructed F from this Delta; every coefficient matches.
2. Verified all stationary numerator derivative arrays through ordinary x-derivative order3 using the Leibniz formula directly, not the producer recurrence.
3. Tested 1,698 exact rational interval combinations, including negative and mixed-sign multiplication, endpoint rounding, and negative-interval reciprocals. Reviewed floor/ceiling formulas for outwardness. Zero-crossing division is rejected.
4. Independently enclosed exp(-x) using 180th/181st alternating-series partial sums at x=0,1,6,7,14,16. Every enclosure lies within the producer interval. Reviewed positive Taylor remainder, the [0,16] guard, the alternating Machin bounds, and outward rational division. The routine uses integer endpoints only for proof decisions.
5. Replayed every adaptive sign-coverage interval and every root bisection from the frozen script. Counts: 1024 cells near1, 2755 descending, 440 in the minimum window, 866 ascending. The intervals cover all [1,2.2] without gaps, hence in particular all [1,2].

## Tail derivation and domain audit

At R=1/100, sigma3(n)<=n^4<=16^(n-1) and sigma5(n)<=n^6<=64^(n-1) justify the producer Cauchy majorant. Similarly sigma1(n)<=n^2<=4^(n-1) gives the weighted E2 coefficient norm at most5/4. Multiplication by E2 therefore gives the claimed component majorants

A_n=n^3+5n^2, B_n=45(n^2+5n), C_n=45n,

in units M R^-n, where n=k+4. On x<=16, put p0=256A_n+16B_n+C_n, p1=32A_n+B_n, p2=2A_n. The order-d derivative, d<=3, is bounded by

n^d p0 + d n^(d-1) p1 + binomial(d,2) n^(d-2) p2,

with absent terms omitted. Dividing by n^6 gives positive combinations of nonpositive powers of n, maximized at n=85. The largest of the four bounds there is 34544583/122825, less than282 and far less than the stated common constant20000. The separate constants2000 n^3 and3000 n^4 also hold. Successive common majorant terms have ratio at most (1/5)(86/85)^6<1/4. The geometric tail is therefore genuinely below10^-20, not a numerical coefficient extrapolation. The same common bound dominates the smaller numerator tails (x>=2*pi>6).

All finite evaluations satisfy x<=4.4*pi<16 and q<1/500. The centered mean-value interval uses a whole-cell derivative enclosure and an outward-rounded midpoint displacement, so it remains valid even when the original rational cell endpoints are not exactly dyadic.

For t>=2.2, each monomial x^j exp(-kx), k>=1 and 0<=j<=2, is decreasing: its derivative has sign j-kx<0. Consequently the absolute sum of all nonconstant terms is bounded by its value at x0=4.4*pi. The omitted tail is bounded at x0, where the finite-domain tail estimate applies; no claim of the same bounded-x estimate at arbitrarily large x is needed. The replay verifies the strict upper integer bound20,193,127,386, below the positive constant30,329,640,000. Thus P stays positive on the entire infinite tail, not merely to the last tested temperature.

The product remainder for J is also valid: the omitted negative logarithm is bounded by192 q^41/(1-q), and exp(-u)>=1-u. Outward interval operations preserve positivity before division.

## Critical-point and global argument

Direct differentiation of U=q^-4 DF and J=q^-4 Delta^4 gives

r_x=P/(990 x^2 J).

The denominator is positive. Session28 reciprocal symmetry r(t)=r(1/t) implies exactly P(1)=0; an approximate interval evaluation of that endpoint is not used to assert zero. P_x<0 on [1,1.02] then forces P<0 immediately to the right of1. The descending coverage reaches1.09; increasing P on [1.09,1.12], with the endpoint signs already certified by neighboring closed coverages, yields exactly one simple zero t0. P>0 after1.12 holds through the infinite tail as above. Hence t0 is the sole minimum on [1,infinity), and reciprocal symmetry supplies the unique other minimum. There are exactly three stationary temperatures on (0,infinity): t0, 1, and1/t0. No endpoint at0 is part of the domain.

Certified outward decimal consequences:

1.098604471215707 < t0 < 1.098604471215708;
11135.8763212 < r(t0) < 11135.8763213.

The root's two floating display endpoints coincide at binary64 precision; this does not collapse the distinct exact rational interval and is not a certification defect. Exact endpoints are retained in JSON.

The branch-sign diagram follows from this monotonicity, the inherited r(1) and limiting value2723200, and the positivity of b. The angular eigenvalues remain separate missing inputs. In particular, scalar positivity does not prove a positive-definite Gaussian Hessian.

## Scope and reproducibility

Run `python3 producer_frozen.py` and `python3 audit.py` in this directory. The first reproduces the entire finite/infinite-domain certificate. The second adds independent construction and arithmetic controls. No shell-vector enumeration, unpublished lattice computation, or novelty claim is part of this audit.
