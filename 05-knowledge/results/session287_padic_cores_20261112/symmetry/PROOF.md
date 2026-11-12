# Two different symmetry-breaking scales in ternary tournament cores

**Status: PROVED IN NOTE; pending independent session review.** Exact laws
below concern the explicit balanced ternary family and uniformly independent
attachment bits, not uniformly random tournament isomorphism classes. The
Hamilton-count divisibility input is inherited. Chan–Pak remains open.

## 0. Inheritance, objects, and the new question

Write `T_0=K1`, `T_e=C3[T_(e-1),T_(e-1),T_(e-1)]`, `n=3^e`, and
`I_e=(3^e-1)/2`. Its full automorphism group is the cyclic wreath group
`W_e`, of order `3^I_e`. The local rotations are C3, not S3.

[Session284 wildcard §§3–4](../../session284_structural_births_20261112/wildcard/PROOF.md)
already proves the complete q-color orbit polynomial, every root-stabilizer
exponent, and the endpoint-stabilizer floor. We use these as inherited
mechanisms, not new results. Its earlier modular-centers source already
proves the two free binary-color orbits. The free action of automorphisms
on ordered Hamilton paths is also inherited (LEM-003).

New here is a probabilistic asymptotic classification of the **visible
3-primary symmetry** and its separation from the **protected path-count
valuation**. The minimal hostile to identifying the two is already the
one-apex extension of the nine-tower: the root-fixing group may be trivial
while every path count is still divisible by3. We make this discrepancy
quantitative at every scale.

Color each of the n leaves independently and uniformly from q colors,
q>=2. Let `S_e` be the base-3 logarithm of its color-preserving automorphism
group. For q=2^k this is exactly the group fixing k named newborn vertices
pointwise, when each old-new arc is an independent fair bit. Arcs among
newborns may be fixed arbitrarily. Keep the actual colors/attachment bits;
a color frequency histogram is not the sufficient recursive state.

**Connection contract.** A tournament core plus a named interface maps to
an iid colored ternary tree. Cyclic child isomorphism preserves the exact
pointwise stabilizer. Root-forgetting loses the orbit of the newborn tuple,
which is explicitly bounded below. Endpoint occupation, rather than this
color symmetry, controls the retained Hamilton-path divisor. No conclusion
about arbitrary odd-target realization follows from a symmetry distribution.

## 1. Exact arithmetic distribution, at every size

Let `O_e(z)` count color-orbits with weight `z^S`. The inherited recurrence is

    O_0(z)=q,
    O_e(z) = [O_(e-1)(z)^3 +(3z-1) O_(e-1)(z^3)]/3.             (1)

Let `L_e(z)` count **labelled** colorings with the same weight. Orbit–stabilizer
converts this into

    L_e(z)=3^I_e O_e(z/3),
    L_e(z)=L_(e-1)(z)^3 +(z-1)3^(2I_(e-1))L_(e-1)(z^3/9).     (2)

Although the second expression contains denominators in its argument, the
result has integer nonnegative coefficients. The first expression is a
safe integer implementation. Normalization is `L_e(1)=q^(3^e)`.

For the local mechanism, let `J_e` indicate that the three colored child
subtrees belong to the same W_(e-1) orbit. Then, exactly,

    S_e=S_(e-1)^(1)+S_(e-1)^(2)+S_(e-1)^(3)+J_e.              (3)

The three child variables are independent, but J_e is not independent of
them. Treating all internal-node symmetry bits as independent is false.
If child orbit t has probability p_t and stabilizer exponent a_t, put

    alpha_e=sum_t p_t^3,
    beta_e=sum_t a_t p_t^3.

With `M_e=E S_e`, `V_e=Var S_e`, (3) gives

    M_e=3M_(e-1)+alpha_e,
    V_e=3V_(e-1)+delta_e,
    delta_e=alpha_e(1-alpha_e)+6(beta_e-M_(e-1)alpha_e).         (4)

Thus the variance's correlation term is retained, not discarded. Exact
alpha and beta are obtained from the same orbit coefficients: each child
orbit of exponent a has probability `3^(I_(e-1)-a)/q^(3^(e-1))`.

Set `rho=3/q^2<1`. Every child orbit probability is at most
`3^((m-1)/2)/q^m`, where m=3^(e-1). Consequently

    alpha_e <= (1/3) rho^(3^(e-1)),
    |delta_e| <= (1+6 I_(e-1)) alpha_e < 3^e alpha_e.           (5)

The first bound follows from `sum p_t^3 <= (max p_t)^2 sum p_t`.
The second uses `0<=a_t,M_(e-1)<=I_(e-1)`.

## 2. Law of large numbers and central limit theorem

Define the computable constants

    mu_q=sum_(j>=1) alpha_j/3^j,
    sigma_q^2=sum_(j>=1) delta_j/3^j.                         (6)

Both converge, and sigma_q^2>0 for every integer q>=2. Explicit errors after
depth r are

    0 <= mu_q-M_r/3^r
       <= rho^(3^r) / [3^(r+2)(1-rho^(2*3^r)/3)],            (7)
    |sigma_q^2-V_r/3^r|
       <= rho^(3^r) / [3(1-rho^(2*3^r))].                   (8)

These follow by bounding the successive terms with a geometric series;
their exponents actually grow faster. They are certified intervals, not
fits to decimal data.

For positivity, when q>=3 the first variance term is
`q^-2(1-q^-2)/3`. The absolute remaining tail is at most
`rho^3/[3(1-rho^6)]`. The first is strictly larger: equivalently
`q^4(1-q^-2)(1-729/q^12)>27`, already true at q=3 and increasing.
For q=2, exact depth3 arithmetic gives `V_3/27 > 56/1000`, whereas
`rho^27/[3(1-rho^54)] < 1/1000`. Thus the limit is positive.

As e tends to infinity,

    S_e/3^e -> mu_q in L2 and probability,
    (S_e-mu_q*3^e)/sqrt(sigma_q^2*3^e) -> N(0,1).             (9)

Here is a proof of the second statement retaining all dependence. Choose
r=r(e) with e^2 <=3^r<3e^2. Partition the tree into n/3^r independent
bottom subtrees. The expected number of symmetry tolls above those blocks
is at most `n*sum_(j>r) alpha_j/3^j`, which tends to zero by (7).
So with probability tending to1, S_e is exactly the sum of these independent
block exponents. Their variances divided by3^r tend to sigma_q^2; their
centered absolute values are bounded by3^r. The total normalized third
absolute moments are O(3^r/sqrt(n)), tending to zero (use
`|X|^3 <=3^r X^2` for centered block X). The elementary characteristic-function
Taylor proof of the triangular-array central limit theorem applies. Formula
(7) also makes the difference between the block mean and mu_q*n tend to0.
This proves (9), not merely an independent-node heuristic.

Certified leading values (the verifier supplies outward decimal brackets):

| q | mu_q | sigma_q^2 |
|---:|---:|---:|
| 2 | 0.09583543106990665... | 0.05658834209466877... |
| 4 | 0.021063913562701205... | 0.019685054809346778... |
| 8 | 0.005212095048705426... | 0.005130401788629269... |

Thus repeated symmetry cores have a complete arithmetic stochastic law,
not just a list of the first few sizes. It concerns a specified model,
not the typical structure under uniform unlabelled enumeration.

## 3. The logarithmic attachment threshold and Poisson law

Now q may vary with e. Bottom-node tolls are independent Bernoulli(q^-2),
and there are n/3 of them. Let B_e be their sum. The remaining tolls are
nonnegative, and (5) gives the explicit coupling bound

    P(S_e != B_e) <= n/q^6 / [1-(3/q^2)^6/3].                 (10)

Indeed the first omitted level contributes at most n/q^6; each subsequent
upper-bound term has ratio at most rho^6/3. This is an upper bound and may
exceed1 at small parameters; it is not an equality.

Consequently the total variation distance to a Poisson variable of mean
`lambda_e=n/(3q^2)` is at most

    n/q^6 / [1-(3/q^2)^6/3] + 2n/(3q^4).                   (11)

The second term is the usual elementary Bernoulli-to-Poisson coupling
bound, obtained by coupling each Bernoulli(p) with Poisson(p) and summing
at most2p^2 error per coordinate. If `n/q^2 -> lambda in (0,infinity)`,
then `S_e -> Poisson(lambda/3)` and the probability of trivial pointwise
stabilizer tends to `exp(-lambda/3)`. If q^2/n tends to infinity, it tends
to1. If q^2/n tends to0, it tends to0, since already P(B_e=0)<=exp(-n/(3q^2)).

For k actual newborn vertices q=2^k. The critical window is therefore

    k = (1/2) log_2 n + O(1).                              (12)

Integer k can cause oscillation in lambda_e; the Poisson limit is stated
along subsequences on which n/4^k converges, not as a nonexistent unique
limit for every rounding convention. This is the **visible symmetry
threshold**, not a path-count divisibility threshold.

## 4. Protected p-adic mass survives after the visible group disappears

Let U extend T_e by k>=1 named vertices. Inherited Session284 gives

    E_e(2k)=sum_(d=0)^(e-1) max(3^d-2k,0),
    v_3 H(U) >= E_e(2k).                                   (13)

This holds for every attachment, including those with no pointwise
symmetry; its proof groups paths by the at most2k old endpoints adjacent
to newborns, then acts with the pointwise stabilizer of those old endpoints.
It does not assert the same group acts on the whole extension.

Let A=Aut(U), and let S_e be the colored pointwise stabilizer exponent.
The orbit of the ordered newborn tuple has at most `(n+k)_k` elements.
Orbit–stabilizer gives the precise correction

    v_3 |A| = S_e + v_3 |A·(newborn tuple)|,
    S_e <=v_3|A|<=S_e+floor(log_3((n+k)_k)).                 (14)

In particular, with Q(U)=H(U)/|A|, which is an integer by the free path action,

    v_3 Q(U) >= max(0,E_e(2k)-S_e-floor(log_3((n+k)_k))).     (15)

No claim that A itself is a3-group is needed: fixing the tuple embeds its
stabilizer in W_e, while the tuple-orbit factor can have other primes.

For fixed k, random attachment therefore gives

    v_3|Aut(U)|/n -> mu_(2^k),
    liminf v_3 Q(U)/n >= 1/2-mu_(2^k) in probability.         (16)

The liminf assertion means: for every epsilon>0 the probability of the
corresponding lower bound with minus epsilon tends to1. For one new vertex,
the protected quotient exponent has density at least
`0.4041645689300933...`, while visible symmetry has density about0.09583543.
These are lower bounds for the quotient valuation, not an exact H formula.

More strikingly, choose k=O(log n) with4^k/n tending to infinity. Then the
pointwise stabilizer is trivial with probability tending to1, yet

    v_3 Q(U) >= n/2-O((log n)^2)                            (17)

with probability tending to1. Deterministically, (13) has loss at most
2ke from I_e, and (14) loses at most k log_3(n+k). Thus (17) follows without
an unjustified equality between count valuation and group valuation.

This separates **color-distinguishing budget** from **path-contact budget**.
Logarithmically many random boundary bits can distinguish the whole tree,
but each Hamilton path touches only2k old boundary vertices; almost all
ancestral rotations still act on its fixed-endpoint class. The surviving
3-adic valuation is a hidden path-fiber symmetry, not a visible automorphism
of U. The first obstruction to a 3-adic unit target remains even after
ordinary symmetry-breaking succeeds.

## 5. Scope and constructive consequence

This provides a precise large-scale classification of one emerging core
family and a sharp probabilistic symmetry-erasure threshold. It also gives
a negative design test: attaching o(n/log n) vertices to a large balanced
ternary core cannot generally supply3-adic unit counts; for logarithmic k,
(13) alone forces a vast divisor for every attachment. To obtain arbitrary
odd targets a constructor must change the underlying core, use enough path
contacts, or apply a lawful reset that breaks the protected endpoint action.
The Session286 reset is such a different operation in its own family; no
interoperability with these ternary cores is asserted here.

The stochastic law, Poisson threshold, and visible/protected valuation
separation are the new conclusions. The orbit recurrence and endpoint
floor are inherited. No universal classification of all prime tournament
cores, all H valuations, or Chan–Pak coverage is claimed.
