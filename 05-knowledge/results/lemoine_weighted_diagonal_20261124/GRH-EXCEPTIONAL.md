# Weighted two-prime circle method: a square-root exceptional set

Assume RH and Dirichlet GRH (the interfaces in Glazer's final link).
Fix 0<η<1/8 and X large. Choose smooth nonnegative f,g compactly
supported in (0,1) and positive on [1/4,3/4]. Write

    S_f(α) = ∑_p (log p) f(p/X)e(αp),
    V_f(β) = ∫ f(t/X)e(βt) dt,
    R(N) = ∫ S_f(α) S_g(2α)e(−Nα) dα.

For odd X≤N≤2X its geometric main term is
J(N)=X∫ f(N/X−2t)g(t)dt ≍X, by taking a fixed-width interval
around t=N/(3X). All weights here are on primes, not prime powers.

## Major arcs: narrow width, many denominators, numerator energy

Put Q=floor(X^(1/2−η)) and T=X^(η/4). Take arcs around a/d,
(a,d)=1, d≤Q, all with radius T/X on the unit circle. They are
disjoint since Q²T=o(X). Write c_d(k)=∑_(r,d)=1 e(kr/d).
For k=1,2 the principal approximation is

    M_{f,k}(a/d+β) = c_d(ka)/φ(d) · V_f(kβ).

In a reduced residue class r, subtract V_f(kβ)/φ(d) from the
corresponding twisted prime sum and call the residual A_r(β).
The multiplicative character Fourier expansion, Parseval on the
units mod d, Dirichlet GRH's explicit formula and partial
summation give

    ∑_r |A_r(β)|²
     = (1/φ(d)) ∑_χ |B_χ(kβ)|²
     ≪ X (1+X|β|)² log^C X.

Here B_χ is the character-twisted prime sum minus the principal
integral for the principal character. Replacing Λ by log-prime
weights costs O(√X log²X), and removing primes dividing d
from the principal term costs O(log d). The residual primes
dividing d in S_f itself contribute O(log d) on each arc.

For *fixed β*, sum squared additive errors over all a mod d
before restricting to reduced numerators. The map a→ka has
fibers of size ≤k, so additive orthogonality bounds the sum by
kd∑_r |A_r(β)|² (plus the harmless nonunit term just noted).
Consequently

    ∥S_f(k·) − M_{f,k}∥²_{2,major} ≪ Q² T³ log^C X = o(X).

This choice of widths is useful: charging each numerator
separately before summing would lose the desired cutoff.
Also ∥M_{f,k}∥²_{2,major} ≪ X log^C X: because |c_d(ka)|≤2
for a reduced numerator and k=1,2, the sum over d is bounded
by a constant times ∑_{d≤Q}1/φ(d), and ∫|V_f(kβ)|²dβ≪X.
Cauchy on each cross term yields total major replacement error o(X).

The main integral factors into a truncated real integral
and truncated arithmetic series. The real tail is O_A(X T^(-A))
by smoothness. Summing over reduced a leaves c_d(N) times
c_d(1)c_d(2)/φ(d)² (both numerator factors are independent of a).
The factor is zero for nonsquarefree d, since c_d(1)=μ(d).
For odd primes the local factor is 1+c_r(N)/(r−1)²;
at 2 it is 1−c_2(N)=2. Thus for odd N

    S(N) = 2 ∏_{r>2} (1−1/(r−1)²)
               ∏_{r|N, r>2} (r−1)/(r−2).

This equals the usual weight at 2N and is bounded below by
the positive twice-twin-prime constant. Uniform truncation
is elementary: multiply absolute summands by d^(1/2).
The resulting Euler product is O(exp(O(√log X))) for |N|≤2X,
since the only non-summable factors come from primes r|N,
with increments O(r^(-1/2)). Split at log X and use
∑_{r|N}log r≤log(2X). So the arithmetic tail beyond Q is
≪Q^(-1/2)X^o(1). The real truncation, even multiplied by the
arithmetic absolute sum X^o(1), is o(X). Uniformly over odd N
in the window the major coefficient is S(N)J(N)+o(X).

## Minor arcs: borrow the additive estimate, not a third prime

Use [Vaughan, Distribution of primes in arithmetic progression,
§7 lemma](https://personal.science.psu.edu/rcv4/personal/Publications/vrprapg.pdf).
In the present smooth form its additive GRH bound is

    |S_f(a/d+β) − μ(d)V_f(β)/φ(d)|
      ≪ (√(Xd) + X√(d|β|)) log^C X

for the near-rational arguments below. Apply the lemma to the
sharp sums for each endpoint, then partial summation with f.
The integral can replace Vaughan's discrete principal sum with
error O(1+X|β|), and the prime-power removal is lower order.
(The elementary GRH estimate used above allows 1+X|β|;
the sharper square root here prevents a loss at small d
outside the narrow arcs.)

By Dirichlet approximation, take d≤ceil(X/Q) and
|β|≤Q/(dX); both errors are ≪X/√Q times logs for Q≤√X.
For d>Q the principal factor is at most X/φ(d),
≪X log^C X/Q. For d≤Q and α on the minor arcs,
|β|>T/X, and smoothness makes that principal term
O_A(X T^(-A)). Thus sup_minor |S_f|≪X log^C X/√Q.
Bessel in N and Parseval in α (doubling preserves Haar measure)
give

    ∑_N |∫_minor S_f(α)S_g(2α)e(−Nα)dα|²
     ≪ (X²/Q) log^C X · ∥S_g∥²_2 ≪ (X³/Q) log^C X.

Every exceptional odd N forces a minor coefficient of magnitude
≫X. Hence the number of exceptions in [X,2X] is
O((X/Q)log^C X), or O_ε(X^(1/2+ε)) after taking η<ε and
absorbing logs. Dyadic summation gives the same order to X.
In fact the same conclusion applies to the set whose weighted
representation count is less than a sufficiently small
fixed multiple of X. This does not make the exceptional set empty.
