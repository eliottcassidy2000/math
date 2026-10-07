# The first-moment field already determines every source weight

2026-10-05. **PROVED, using the cited classical identities:** the degree-one
shell observation operator on the half-line is boundedly invertible on l2,
with an explicit rational coercivity constant and finite-readout error bounds.
**FINITE-EXACT:** rational matrix, tail, and parameter controls. **OPEN:**
independently positive measurements at every Collatz source. This is an
information-recovery theorem, not universal source positivity.

[Checker](../../04-computation/experiments/collatz_spatial_fourier_inverse_20261005.py)
and [output](collatz_spatial_fourier_inverse_20261005.out).

## 1. Inheritance and the question that changes the representation

The [spatial sidecar](collatz_specialization_sidecar_20261005.md) reconstructs
the law from a fixed degree k>=7 field by a simple l-infinity Neumann bound.
Seven is the threshold for that diagonal-dominance argument. Does lower
degree actually lose information, or merely defeat that proof?

The closest proved mechanism is finite spatial inversion. The hostile is
p=delta_0: its smoothed field is positive everywhere, but every other atom
is zero. The corrected near miss is confusing a norm sufficient condition
with a necessary invertibility condition. The underused coordinate is the
Fourier frequency of the full source-indexed field.

Portfolio: **anchor** source measurement bounds; **niche** degree-one spatial
identifiability; **wildcard** the beta-integral formula connecting a discrete
kernel to a positive continuous Fourier transform. Board: source address /
degree / frequency / coercivity / finite tail / sign. No cosmetic tournament
is introduced: the kernel is symmetric and has no intrinsic orientation.

## 2. Continuous positivity supplies a discrete spectral gap

Write

    a=(ln2)/2,
    f(x)=sech^2(a*x),
    h_d=f(d)=4*2^d/(1+2^d)^2 for integer d>=0.

Use the Fourier convention fhat(xi)=integral_R f(x) exp(-i*xi*x) dx.
With s=xi/(2a) and u=exp(2ax), substitution gives

    fhat(xi)=(2/a) integral_0^infinity u^(-is)/(1+u)^2 du
            =(2/a) Gamma(1-is) Gamma(1+is)
            =pi*xi/[a^2 sinh(pi*xi/(2a))].                 (1)

The value at zero is its continuous limit2/a. The external inputs are
[Euler's beta integral](https://dlmf.nist.gov/5.12.E3) and the
[gamma recurrence and reflection formula](https://dlmf.nist.gov/5.5).
The integral has real parameter parts1, so those beta hypotheses hold.

In particular fhat is positive everywhere. It decreases with |xi|: the
derivative sign of u/sinh(u) follows from u*cosh(u)-sinh(u)>0 for u>0,
whose derivative is u*sinh(u)>0. Because f is smooth and decays exponentially,
[Poisson summation](https://dlmf.nist.gov/1.8.E14), with the convention above,
gives the absolutely convergent symbol

    S(theta)=sum_(d in Z) h_|d| exp(-id*theta)
            =sum_(r in Z) fhat(theta+2*pi*r),   |theta|<=pi.

There is no extra2*pi factor with this convention. Retaining just r=0 proves

    S(theta)>=c_star:=4*pi^2/[(ln2)^2 sinh(pi^2/ln2)]>0.    (2)

For a convenient rational constant put c=1/200000. The elementary bounds
3<pi<22/7, ln2>56/81, ln2<1 and e<3 give

    pi^2/ln2<15,
    4*pi^2/(ln2)^2>36,
    c_star>72/3^15>1/200000=c.                            (3)

The log bound follows from the first two positive terms of
ln2=2 sum_(j>=0) (1/3)^(2j+1)/(2j+1). No floating estimate is used for c.
Also sum_(d in Z) h_|d|<1+8 sum_(d>=1)2^(-d)=9.

Let A=(h_|i-j|)_(i,j>=0). Extend a finitely supported vector on the half-line
by zero to all integers. Parseval and(2)-(3), followed by l2 completion, show

    c ||v||_2^2 <= <v,Av> <= 9 ||v||_2^2.                 (4)

Thus A is bounded, self-adjoint, and boundedly invertible on l2(N_0), with
||A^(-1)||<=200000. One can also see surjectivity directly from the operator
norm convergence of (1/9)sum_(r>=0)(I-A/9)^r.

For any probability law p on N_0, its **first-moment field** is H=A p,

    H(m)=sum_(j>=0) p_j h_|m-j|.

Since p is in l1 and therefore l2, (4) proves that the whole field uniquely
determines p. More quantitatively, any two such laws obey

    ||p-p'||_2 <= 200000 ||H-H'||_2.                       (5)

This is a bound on the whole source-addressed field. A single H(m) neither
determines nor forces a positive p_m. The constant is deliberately coarse;
no claim of optimal conditioning or priority is made.

## 3. A finite spatial readout with an explicit error bill

Let A_R keep only entries with |i-j|<=R, and set

    delta_R=8*2^(-R),   gamma_R=c-delta_R>0,
    B_R=I-A_R/9,       rho_R=1-gamma_R/9.

R>=22 suffices (21 also gives a positive gamma). The omitted row and column
sums are bounded by delta_R, so ||A-A_R||_2<=delta_R. Moreover

    gamma_R I<=A_R<=9I,       0<=B_R<=rho_R I.              (6)

For an integer K>=0 define

    V_(R,K)=(1/9) sum_(r=0)^K B_R^r.

Its mth row uses only addresses j>=0 with |j-m|<=RK, and all its coefficients
are rational and computable. Substituting H=A p gives the exact identity

    V_(R,K)H-p=-B_R^(K+1)p+V_(R,K)(A-A_R)p.

Since ||p||_2<=1 and ||V_(R,K)||_2<=1/gamma_R,

    |(V_(R,K)H)_m-p_m|
       <=rho_R^(K+1)+delta_R/gamma_R.                     (7)

If independent simultaneous intervals enclose these finitely many H(j),
apply the exact signed row to the intervals and add the error(7). A uniform
absolute measurement error u costs at most u times the absolute row sum;
the latter is at most sqrt(2RK+1)/gamma_R by Cauchy-Schwarz, or the looser
rational bound (2RK+1)/gamma_R. Intersect successive justified intervals
to retain an earlier positive bound through refinement.

The checker includes a terminating parameter compiler for tolerance eps in
(0,1]. Choose delta_R<=c*eps/8, and K so that

    1/[1+(K+1)(1-rho_R)]<=eps/4.

The Bernoulli inequality bounds rho_R^(K+1) by that expression; meanwhile
delta_R/gamma_R<=eps/7. Choosing per-measurement error at most
eps*gamma_R/[4(2RK+1)] pays the remaining noise. All choices use rational
arithmetic and integer comparisons. This demonstrates finite computability,
not practical efficiency: the coarse plan for eps=1/16 uses radius28 and
depth114079969. Those enormous stencils were **not** expanded or executed.
The degree8/12 sidecar is the much better computational method currently
available in this comparison.

## 4. What the Fourier bridge preserves and what remains missing

| Source | Map | Preserved predicate | Necessary retained data |
|---|---|---|---|
| Continuous sech-squared kernel | Sampling and Fourier periodization | Strict positive transform gives a uniform discrete quadratic-form bound | Transform normalization and all frequency aliases |
| Half-line probability law | p to the whole H=A p field | Full law is recoverable in l2 | Every source address, not one scalar |
| Infinite inverse | Finite band and polynomial inverse | Source-specific interval with a proved remainder | Band, depth, and simultaneous measurement error |

The beta integral and the Collatz weight construction both use beta-function
arithmetic, but the measures are different: here the integral computes a
Fourier kernel. It does not assign Collatz ROOT membership. Similarly,
evenness makes this a cosine symbol; it does not erase source addresses.

For p=delta_0 every H(m)=h_m is positive while p_m=0 for m>0. Consequently
neither positivity of S nor positivity of the smoothed field proves the
target atom positive. In the Collatz application, independently bounding a
source-specific signed readout remains the decisive unproved step.

## 5. Validation

    python 04-computation/experiments/collatz_spatial_fourier_inverse_20261005.py
    python -O 04-computation/experiments/collatz_spatial_fourier_inverse_20261005.py

There are705 exact checks: LDL factorizations of both spectral inequalities
for principal matrices of sizes1..12; all243 signed vectors in{-1,0,1}^5;
rational finite-tail and Bernoulli controls for radii22..44; six tolerance
plans through2^(-32); the everywhere-positive-field hostile; and eight
malformed inputs. Normal, optimized, and saved output agree. The infinite
spectral theorem follows from the proof and cited identities, not from
finite eigenvalue sampling.
