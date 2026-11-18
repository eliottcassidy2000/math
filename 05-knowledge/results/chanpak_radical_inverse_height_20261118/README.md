# A sharp address-height law for the clone family

A refinement of the [fixed-core height selector](../chanpak_sun_height_selector_20261118/README.md).
This counts **addresses**, not distinct Hamiltonian totals. It gives
an upper bound on attainable target density, not an equality for
the spectrum.

## 1. Arc cost becomes a binary radical inverse

Use the earlier H(k,u)=W_k u+h_k(2−s(u))−5, where u is odd,
W_k=36·3^k, h_k=12·2^k−1, R_k=W_k−h_k, and s(u)=u−F(u).
Write the bits of u as b_i, least significant first, and set
ρ(u)=Σ b_i 2^(−i−1). The append-zero formula for F shows that
F(u)−1 is the weighted sum of 1-before-0 bit pairs, with an
extra unit for each adjacent such pair. For odd u, subtract
this from u−1. The sum of weights at all lower positions
under a 1 at position j is 2^j−1. Since b_0=1, the remaining
population count minus one counts adjacent pairs with high bit
1; remove the 10 bonuses to leave only adjacent 11 pairs.
Consequently

    s(u) = Σ_(i>j) b_i b_j 2^(i−j−1) + Σ_i b_i b_(i+1).

The first sum is the upper triangle of uρ. Within a run of
l consecutive ones, the complementary triangle weighs
l−1+2^(−l). Subtract that run's l−1 adjacent pairs, and its
contribution is positive. Contributions between runs are
nonnegative. Also the weight over j≥i for each fixed i is at
most 1. Thus for r=bit length(u),

    0 < E(u):=uρ(u)−s(u) ≤ r,
    H(k,u) = (W_k−h_kρ(u))u + 2h_k−5 + h_k E(u).

In particular the radical-inverse formula gives a one-sided
numerical height certificate with an explicit logarithmic error,
not an average guess about the word.

## 2. Exact leading coefficients and a uniform error bound

For X≥2, write A_k(X) for the number of positive odd u with
H(k,u)≤X. The following strengthened count holds uniformly:

    Σ_(k≥0) A_k(X) = X Σ_(k≥0)c_k + O(√X+log X),

where

    c_k = ∫_(1/2)^1 dt/(W_k−h_k t)
        = log((W_k−h_k/2)/(W_k−h_k))/h_k.

Proof: restrict first to k with 24·3^k≤X; there are O(log X)
such k, and no addresses at larger k. Split odd u by its low
b bits with B=2^b≥2. Each of the B/2 classes is one integer
arithmetic progression of step B. Its truncated radical
inverse p takes the values 1/2,1/2+1/B,...,1−1/B in some
order, and p≤ρ(u)<p+1/B.

Write g(t)=W_k−h_k t≥R_k>0 on that interval. Because the
remaining term in H is positive and bounded by
O(h_k(1+log X)) for u≤X, the number in each class is between
the counts below (X−O(h_k(1+log X)))/g(p) and below
X/g(p+1/B), with the lower cutoff clipped at zero. The
arithmetic-progression rounding errors total O(B).
The variation between the two Riemann sums is
O(X h_k/(B R_k²)) since 1/g has derivative at most
h_k/R_k². The total height remainder contributes
O(h_k(1+log X)/R_k). Therefore
A_k(X)=c_k X+O(B+Xh_k/(BR_k²)+h_k(1+log X)/R_k).
Choose B dyadic comparable to max(2,√(Xh_k)/R_k).
The √X and log terms now sum geometrically because
R_k≥24·3^k and h_k≤12·2^k. At omitted k,
c_k≤1/(2R_k), so the main-term tail times X is O(1).
This proves the uniform estimate. No injectivity in u or
between clone depths was assumed.

Thus at most 2Σ c_k = 0.052301466... of odd absolute targets
can be covered asymptotically. For k=0, if representatives of
all odd classes modulo an even M have count at most T_M,
then A_0(T_M)≥M/2. Consequently
liminf T_M/M≥1/(2c_0)=27.6589200... . This improves the
earlier 25 lower bound. It is a necessary condition, not a
claim that the bound is attainable. The series tails are
geometric; the decimals are rounded.

## 3. Typical height of the earlier greedy lift

For the earlier high-bit algorithm modulo 2^L, any decision
sequence ε_b∈{0,1}, b=1,...,L−1, has a different output residue:
the first differing decision flips bit b and later decisions
cannot change it. All sequences occur for uniform odd targets.
The tag width therefore has exactly the distribution

    r = 1 + Σ_(b=1..L−1)(b+1) ε_b,

with independent fair bits. Its mean is
L(L+1)/4+1/2 and its variance is
(Σ_(j=2..L)j²)/4 = Θ(L³). Thus r/L² tends to 1/4 in probability,
and the same holds for log₂ H/L² by the fixed-core linear
bounds on H in u. The quasipolynomial overshoot is typical
for this algorithm, not just an adversarial worst case.
It leaves open different compilers and different cores.

Run `node verify.js` for the exact binary identities and
small-height recurrence controls.
