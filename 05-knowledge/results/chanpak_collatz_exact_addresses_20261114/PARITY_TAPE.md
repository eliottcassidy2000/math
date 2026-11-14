# Complete binary addresses with an exactly sparse absolute spectrum

This is a scoped theorem about one fixed physical grammar, not a
restriction on all one-apex tournaments or the unrestricted spectrum.

Take a transitive spine, writing a 1 for each vertex beaten by its
apex and 0 otherwise. Start with word 10. At indices j=1,2,... in
increasing order, either skip j or append the block B_j=1^j 0. No
other appends or changes of source are allowed. If P is the binary
address of the spine and H the Hamilton-path count, the inherited
terminal-zero transfer is

    B_j: (P,H) → (2^(j+1)(P+1)−2, H+2^j(P+1)),

starting with (P,H)=(2,3). For the original recurrence see
[the earlier modular-escape proof](https://github.com/eliottcassidy2000/math/blob/138bab67a4a6c0ddab82269e295a5be27a5fa877/05-knowledge/results/chanpak_least_missing_jacobian_20261114/MODULAR_ESCAPE.md).
The new point here is the cost of the strictly increasing itinerary,
not the individual append formula.

## Prefix isometry

For two choices, let j be their first difference. They start that
index at the same state. Because P is even, exactly one of the
choices adds 2^j times an odd number. Every later increment vanishes
modulo 2^(j+1). Thus the two scalar outputs differ at precisely bit j,
whether they have finite support or are taken as 2-adic limits.

In particular all 2^(a−1) residue classes of odd targets modulo 2^a
correspond bijectively to the first a−1 choices. The compatible limit
is a homeomorphism from binary itineraries onto the odd 2-adic
integers and carries fair independent choice bits to normalized Haar
measure. This matches the finite-prefix mechanism of the
[Bernstein–Lagarias Collatz coding](https://doi.org/10.4153/CJM-1996-060-x);
it is not a dynamical conjugacy between these tournament operations
and the Collatz map.

The integer tape compiler starts at H=3, compares bit j of the
current H with bit j of the target N, and appends B_j just when it
must change. It recognises the *absolute* image too: stop on equality
or overshoot. If not already equal before processing the first
j with 2^j>N, then the congruence and positivity below that step
preclude indefinitely remaining below N. More explicitly, after
setting bits through j with 2^(j+1)≥N, a positive count in [1,N]
congruent to N modulo 2^(j+1) must equal N. The earlier overshoot
cannot be repaired because every nontrivial append is strictly
increasing.

## Exact budget and typical inflation

In this grammar, the shortest spine that attains a given residue
modulo 2^a has length

    m = 2 + sum_{j=1}^{a−1} (j+1) b_j.

The indices below a are forced, and higher ones are unnecessary.
If the target residue is uniformly random, then the b_j are
independent fair bits by the isometry. The mean is
3/2+a(a+1)/4 and the variance is
(a(a+1)(2a+1)/6−1)/4. Worst length is 1+a(a+1)/2,
and m/a^2 tends to 1/4 in probability.

Each literal word starts with 1 and ends in 0. For m≥3 the 1st and
last vertices contribute at least 2^(m−2) to H−1, since any even
subset of the m−2 interior vertices forms one directed odd cycle
with the apex and those endpoints; multiply its count by 2.
For the upper bound, the contribution from a 1 before a 0 at
distance d is 2 when d=1 and 2^(d−1) otherwise; this is the
same odd-cycle subset count. Sum the geometric bounds on contributions
over all possible right endpoints, and add at most m adjacent
corrections, to get H≤2^m+m. Therefore log₂H=m+O(1).
The ordinary magnitude of a typical least-cost modular output is
2^((1/4+o(1))a^2), despite the modulus being 2^a.

## Exact counting asymptotic

Let A(X) count the distinct absolute outputs of finite tapes at
most X, and L=log₂X. The prefix isometry makes all finite distinct
itineraries have distinct absolute counts. By the two-sided
bounds just proved, A(X) is sandwiched by the numbers of
partitions into distinct parts each at least 2 with total at
most L−O(1) and L+O(1), respectively. The parts are exactly
the block lengths j+1.

The classical [Hardy–Ramanujan asymptotic](https://ramanujan.sirinudi.org/Volumes/published/ram36.html)
for the number q(n) of partitions into distinct parts is

    q(n) ~ [1/(4 3^(1/4))] n^(−3/4) exp(π√(n/3)).

Its cumulative count Q(n) is asymptotic to
3^(1/4)/(2π) n^(−1/4) exp(π√(n/3)), by integral comparison.
If F(n) is the cumulative count with part 1 forbidden, then
Q(n)=F(n)+F(n−1), while F(n)−F(n−1)≤q(n)=o(Q(n)).
Thus F(n)~Q(n)/2. A bounded shift of n makes no difference to
this equivalent; consequently

    A(X) ~ [3^(1/4)/(4π)] L^(−1/4) exp(π√(L/3)).

The resulting separation is sharp: complete binary congruence
coverage in this physical grammar gives only a subpolynomial
number of exact positive integers. An inductive Chan–Pak solver
cannot merely complete this 2-adic address; it needs other appends
or changing sources supplying an independent magnitude control.
