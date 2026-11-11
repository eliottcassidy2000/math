# Common-mode stripping and constant-gcd sector convolution

Status: PROVED identities by the elementary derivation below; FINITE-EXACT
implementation controls. No theorem number reserved. No fastest-known or
full-A000568 polynomial-time claim. These are weighted Burnside sectors, not
counts of a specified subclass of unlabelled tournaments.

## Inheritance and live concepts

Inputs are THM-514 (`01-canon/theorems/THM-514-the-a000568-unit-tail-core-kernel-and-prime-tail-ladder.md`),
THM-4845 (`01-canon/theorems/THM-4845-a000568-as-a-gaussian-average-gcd-covariance-cycle-index-exponential-and-the-meromorphic-normalized-series-with-poles-at-powers-of-four.md`),
and THM-5097 (`01-canon/theorems/THM-5097-factorial-gauge-integral-cycle-response-recurrence.md`).
The first already solves unit tails and prime divisibility responses; the
second already supplies Gaussian divisor covariance; the third supplies
integral factorial transfer. The contribution here is their explicit
common-mode gauge, a budget-aware component factorization, and a restricted
sector implementation. No literature novelty claim is made.

Concept board: Gaussian common modes; arithmetic gcd components; factorial
integrality; coefficient truncation; pantograph equations. Anchor: accelerate
exact Burnside evaluation. Niche: composite constant-gcd sectors. Wildcard:
Gaussian independence becomes deterministic polynomial convolution.
Hostile: lengths 9 and15 have gcd3, not1. Least-used sidecar: the number of
cycles K, which carries the removed common-mode energy. Corrected near miss:
a Gaussian expression alone does not evaluate its expectation or prove a
finite-order ordinary recurrence.

## 1. Exact sector identity, including a boundary response

Let S be a finite set of positive odd lengths, g a positive odd common divisor
of every length, and N the mass horizon. Let r_l be nonnegative integral
charges from an already fixed outside multiplicity vector. Define

    B_n(S;r) = sum_{sum l m_l=n} 2^{E(m)+sum r_l m_l} / prod(l^{m_l}m_l!),
    E(m) = (sum_{l,h} gcd(l,h)m_l m_h - K)/2,
    K = sum_l m_l.

The exponent decomposes exactly as

    E(m) = g binom(K,2)
           + 1/2 sum_{l,h}(gcd(l,h)-g)m_l m_h + (g-1)K/2.

Put an edge between distinct l,h in S if gcd(l,h)>g AND l+h<=N.
For every connected component C form the truncated bivariate series

    F_C(x,y) = sum_m x^{sum l m_l} y^{sum m_l}
      * 2^{H_C(m)+sum r_l m_l} / prod(l^{m_l}m_l!),
    H_C(m) = 1/2 sum_{l,h in C}(gcd(l,h)-g)m_l m_h +(g-1)sum m_l/2.

All series are interpreted modulo x^{N+1}. Then

    B_n(S;r) = sum_K 2^{g binom(K,2)} [x^n y^K] prod_C F_C(x,y), n<=N.

Proof: expand the product and use the displayed quadratic identity. Any
omitted interaction with l+h>N multiplies a monomial of degree>N, so vanishes
in the truncation. All other nonzero residual interactions lie in one
component. This proves the formula coefficientwise. It preserves every
Burnside summand and charge, losing no data at a fixed boundary. An aggregated
outside sector still needs its response r; outside mass alone is insufficient.

Without deleting budget-infeasible edges, this is also the Gaussian common
mode X_l = sqrt(g log2) Z + Y_l, where
Cov(Y_l,Y_h)=log2*(gcd(l,h)-g). This residual matrix is positive semidefinite:
write l=g a and use g*(gcd(a,b)-1), the divisor covariance with divisor1
removed. Different graph components are independent. After budget deletion,
only the formal truncated identity is asserted, not a new Gaussian covariance.

## 2. Constant-gcd support is a polynomial-state kernel

If every feasible distinct pair has gcd g, every component is a singleton:

    F_l(x,y) = sum_{m>=0} 2^{((l-g)m^2+(g-1)m)/2+r_l m}
                           x^{lm} y^m / (l^m m!).

This includes composite supports g times pairwise coprime odd numbers, for
example {3,9,15,21,33,39,...}. Primality is not required. At horizon N even
non-coprime pairs are allowed if their total mass exceeds N.

Define the integral coefficient C[d,K] = d![x^d y^K] product F_l.
Adding length l with multiplicity m gives

    C_new[d+lm,K+m] += C[d,K]
      * (d+lm)!/(d! l^m m!)
      * 2^{((l-g)m^2+(g-1)m)/2+r_l m}.

The factorial quotient is integral: it chooses the new lm labels and then a
permutation consisting of m cycles of length l. All exponents are nonnegative
integers. Finally

    n! B_n = sum_K C[n,K] 2^{g binom(K,2)}.

There are at most (N+1)^2 states and O(|S| N^3) transfers, a loose polynomial
arithmetic bound, even when the support grows with N. Integer bit sizes are
polynomial for charges of polynomial magnitude. This is not a polynomial
algorithm for unrestricted A000568: its residual gcd graph has a large
connected component whose evaluation remains hard.

For a singleton let f_l(z)=sum_m 2^{((l-g)m^2+(g-1)m)/2+r_l m} z^m/(l^m m!).
Its exact formal pantograph equation is

    f_l'(z) = (2^{(l-1)/2+r_l}/l) f_l(2^{l-g}z), f_l(0)=1.

This is a dilation equation, not an ordinary D-finite/holonomic claim. For
l>g its coefficients grow too fast for a positive analytic radius. The final
common-mode transform weights coefficient y^K by 2^{g binom(K,2)}.

## 3. The computational gain is delayed energy, not fewer states

The direct constant-gcd response DP on the same (d,K) states uses exponent

    (l-1)m/2 + l binom(m,2) + g K m + r_l m.

Subtracting g binom(K,2) from the accumulated energy changes this exactly to
the exponent in section2. Thus the two implementations have identical states
and transfer counts; the gauged version carries smaller integers and restores
the large power of2 at the end. `gauge=False` in the implementation is this
same-state conventional response control, not naive partition enumeration.

The method also applies without constant-gcd support, for g=1 on all odd
lengths. Keep K and responses sum_h(gcd(l,h)-1)m_h. The universal transfer
has exponent (l-1)m^2/2 + m response_l. `common_mode.py` integrates this with
the inherited odd-factorial gauge: subtract v2(m!) and restore binom(K,2)
when no future length is feasible. Before the unit tail, shifts remain
nonnegative because l>=3; at the unit tail the restored binomial dominates
v2(m!). The low packed coordinate stores K, not the zero residual response
of length1. Restoration when a response is discarded is essential.

## 4. Finite controls and measured scope

Files `benchmark.jsonl`, `gauge_benchmark.jsonl`, and
`common_mode_benchmark.jsonl` are exact-run evidence. Timings are single-run
wall times on the shared environment, not robust hardware comparisons.

Constant-gcd support S=3*{1,3,5,7,11,13,17,19,23,29,31,37,41,43}:

- N=300: gauge0.0796s, same-state response0.1226s, 1.54x;
  76,509 transfers and2,507 peak states.
- N=600: gauge2.1177s, response4.7881s, 2.26x;
  528,343 transfers and9,908 peak states.
- N=900: gauge15.4027s, response34.4880s, 2.24x;
  1,694,066 transfers and22,195 peak states.
- Independent original-gcd multiplicity recursion through N=390 compares
  all7,573,727 feasible vectors, with exact prefix equality. That direct
  recursion takes63.44s versus0.251s for convolution; this252x comparison is
  against a naive baseline and is not the principal speed claim.

Charges and budget controls: support{3,9,15,21}, N90, r9=7,r15=4 agrees with
original-gcd recursion; support{9,15}, g1,N20 is permitted and agrees because
no two distinct lengths coexist; the same unsupported g1 atN24 is rejected.
Small positive controls include unit and prime supports and repeated cycles.

Universal A000568 gauge agrees with `local_scale.py` at all coefficients for
horizons0,1,2,3,4,5,10,20,30,50,100,150,200. At200 it uses2.675s versus2.947s,
identical2,994,814 transfers, a modest1.10x observation. At150 it is slower.
Do not describe this as a decisive general enumeration improvement. The
strong measured improvement is in constant-gcd composite sectors.

Reproduce from the repository root:

    python3 05-knowledge/results/session270_chanpak_enumeration_20261111/terminal/coprime_kernel.py
    python3 05-knowledge/results/session270_chanpak_enumeration_20261111/terminal/benchmark_gauge.py
    python3 05-knowledge/results/session270_chanpak_enumeration_20261111/terminal/common_mode.py

Next useful target: detect residual gcd components dynamically at a response
boundary and choose component convolution only when its cost is below transfer.
No full-counter improvement from that adaptive strategy has yet been measured.
