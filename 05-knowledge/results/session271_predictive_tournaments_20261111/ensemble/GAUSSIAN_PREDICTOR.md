# A Gaussian prediction route, and why the choice of sampling measure matters

**PROVED elementary variance calculations.** This is a concrete connection to
[THM-4845's inherited Gaussian gcd field](../../../../01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md),
not a new derivation of that representation or of Cameron–Martin tilting.
It gives a useful counterindication for an otherwise attractive sixth
prediction strategy: do not sample the unshifted field to predict enormous
counts. No expensive Monte Carlo run is needed to discover that failure.

## 1. Exact identity, disastrous naive predictor

Use centered Gaussian X_j on odd j with
Cov(X_i,X_j)=log(2) gcd(i,j), and W_j=2^(-1/2) exp(X_j).
The inherited theorem identifies A000568 with the expected odd cycle-index
coefficient. The identity-cycle term is

    Y_n = W_1^n/n! = 2^(-n/2)exp(nX_1)/n!,
    E Y_n = b_n = 2^binom(n,2)/n!.

The normalized variance is exactly

    Var(Y_n)/(E Y_n)^2 = 2^(n^2)-1.                        (1)

Thus an M-sample ordinary Monte Carlo estimator of this term has relative
root-mean-square error sqrt((2^(n²)-1)/M). Reaching relative RMS ε needs
M >= (2^(n²)-1)/ε² for this estimator. This is a statement about this estimator's
RMS, not a universal lower bound on randomized counting or a guarantee for
any realized sample. At n=1000 the variance factor alone is 2^1000000.
A mathematically elegant stochastic representation is not automatically an
efficient stochastic predictor.

## 2. Identity tilt removes the dominant variance exactly

Exponentially tilt the Gaussian measure by nX_1. Under the new law Q_n,

    E_Q X_j = n log(2) for every odd j,
    Cov_Q(X_i,X_j)=log(2) gcd(i,j),
    dP/dQ_n = exp(-nX_1+n²log(2)/2).

The likelihood-weighted identity contribution is identically b_n, with zero
variance. This is the standard Gaussian exponential tilt applied to the
actual dominant Burnside sector, not a heuristic recentering.

For a core multiset μ of nonunit cycles, let m=sum j, r=number of cycles,
d=m-r, and z_μ its centralizer factor. After normalizing by b_n, its weighted
random contribution is

    V_μ = (n)_m/z_μ · 2^(d/2) exp(sum_(j in μ) X_j - mX_1). (2)

Under Q_n the exponent in (2) has mean -dn log2 and variance

    v_μ = log2 [sum_(i,j in μ) gcd(i,j) - 2mr + m²].         (3)

Consequently E_Q V_μ is exactly the deterministic cycle correction w_μ from
the forecast lane, and

    Var_Q(V_μ)/(E_Q V_μ)^2 = exp(v_μ)-1.                    (4)

At fixed support this variance factor no longer depends on n. For μ=(3),
v_μ=6log2, so the factor is exactly63. The relative RMS for estimating the
small correction itself is sqrt(63/M); its absolute RMS as a correction to
the normalized full count is δ_3 sqrt(63/M). This has a very different scale
from (1). Direct analytic evaluation still wins for the present problem.

The covariance of two log contributions is also explicit:

    c_(μ,ν)=log2 [sum_(i in μ,j in ν)gcd(i,j)
                    - m_ν r_μ - m_μ r_ν + m_μ m_ν].       (5)

Their normalized multiplicative correlation is exp(c_(μ,ν)). Two singleton
3-cycle cores have c=6log2, hence correlation factor64, the same enhancement
that appears between two disjoint 3-cycles in the deterministic correction
gas (up to falling-factorial depletion). The Gaussian and cycle-gas views
retain the same gcd interaction instead of pretending rare symmetries are
independent Poisson events.

## 3. What this opens, and what remains missing

A finite-support importance sampler is now mathematically specified and has
an exact variance matrix from (3)–(5). It might become useful after inserting
an additional observable whose Gaussian expectation lacks a cheap formula.
For unweighted A000568 all these moments are explicit, so deterministic
support summation is preferable. A uniform variance bound for sampling the
entire growing-support cycle index is NOT proved here; fixed-support variance
control must not be extrapolated to all sectors.

Connection contract:

- Source: the inherited Gaussian gcd-field cycle-index representation.
- Target: a normalized identity-tilted prediction ensemble.
- Map: Cameron–Martin shift along the unit-cycle field X_1.
- Preserved: each coefficient expectation and all retained sector means.
- Changed: sampling law and variance, not the counted isomorphism classes.
- Sidecar: likelihood ratio, support cutoff, covariance matrix and tail bound.
- Decisive test: analytic variance (1) versus (4), before any simulation.

This is a precise bridge between Gaussian large-deviation tilting, rare
symmetry interactions and predictive enumeration. It is not a proposed new
constant-time exact enumeration algorithm.
