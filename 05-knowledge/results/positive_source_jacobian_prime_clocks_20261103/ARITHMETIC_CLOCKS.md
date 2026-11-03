# 290, 227, and the difference between density and first arrival

Status: **PROVED elementary clock and CRT identities / FINITE-EXACT census**.
The supplied arithmetic table is correct when the clock counts all primes
coprime to m, starting from the smallest prime. Its probabilistic explanation
needs a different status.

## 1. What the seven rows actually count

For even m>=4 let p_m be the least prime congruent to1 or-1 modulo m and
tau_m the number of eligible primes <=p_m. Since p_m>=m-1,

    tau_m=pi(p_m)-omega(m).                             (1)

The distinct prime factors of even m are all <=m/2 and are precisely the
excluded primes. For m=2 the two target classes coincide; it is outside
both the two-class probability model and the following table.

| m | phi(m)/2 | p_m | tau_m | tau_m/(phi(m)/2) |
|---:|---:|---:|---:|---:|
|262|65|263|54|54/65|
|278|69|277|57|19/23|
|292|72|293|60|5/6|
|120|16|239|49|49/16|
|144|24|431|81|27/8|
|208|48|1039|173|173/48|
|290|56|1451|227|227/56|

For290 the entire exclusion certificate before1451 is

    289=17^2,        291=3*97,
    579=3*193,       581=7*83,
    869=11*79,       871=13*67,
    1159=19*61,      1161=27*43,
    1449=9*7*23,     1451 prime.

These are all candidates km+/-1 in increasing order, through5m+1.
There are230 primes through1451, of which2,5,29 are excluded, giving227.
This finite factorization explains the delay exactly; a distributional
heuristic is not needed for this example.

The complete census over even4<=m<=20000 shows that290 is a record at its
location, but is not an ultimate normalized record. Later records include

| m | p_m | tau_m | normalized ratio |
|---:|---:|---:|---:|
|1688|15193|1772|443/105|
|1820|12739|1516|379/72|
|4106|53377|5441|5441/1026|
|9970|149551|13810|6905/996|

The last ratio is about6.933, the maximum in this stated finite universe.
No all-size maximum or unboundedness conclusion is inferred. The full
record sequence, sieve bound2000000 and reproduction code are stored.

## 2. Which expectation is meant?

Independent draws uniformly distributed over the phi(m) units have success
probability2/phi(m); their geometric waiting time has mean phi(m)/2.
That is a valid model calculation. In the actual prime ordering, however,
all candidates below m-1 have zero chance to hit: the starting phase is
already incompatible with an iid constant-hazard model. It is not a theorem for primes examined
in increasing order. For a fixed modulus, the prime number theorem in
progressions supplies asymptotic equal proportions; it supplies neither
independence nor a random distribution for the first prime. Dirichlet's
infinitude theorem alone is weaker still.

[Granville--Martin's primary survey](https://arxiv.org/html/math/0408319v1)
distinguishes these distribution questions. The
[Rubinstein--Sarnak paper](https://doi.org/10.1080/10586458.1994.10504289)
studies prime races and, under its stated hypotheses, logarithmic densities.
It is not a formula for this first-arrival statistic. The supplied
[arXiv:1608.08647](https://arxiv.org/abs/1608.08647) concerns Gaussian
primes; it is not direct evidence for the ordinary-prime claim made here.

For290 both1 and-1 are quadratic residues. The square roots of-1 are
17,133,157,273; roots of1 are1,59,231,289. Thus a quadratic-residue race is
a sensible separate object to examine. This observation does not show
that Chebyshev bias caused the finite delay. In120,144,208, -1 is not a
square; the union{1,-1} does not even have the same residue-type profile.

## 3. A precise phase/variance theorem

Take a cyclic binary word of length N with K>=1 marked hits, and let g_i
be its successive positive hit gaps, sum g_i=N. Two clocks give

    mean next-return gap, starting at a uniform hit = N/K;
    mean inclusive wait, starting at a uniform position
       = (sum_i g_i^2+N)/(2N)
       = (mean(g)+1)/2 + Var(g)/(2 mean(g)).             (2)

Proof: each gap of length g contributes waits1,2,...,g exactly once.
For N=12,K=3, gaps(4,4,4) give wait5/2; gaps(1,1,10) give19/4.
Both next-return means are4. Density fixes the return mean but not the
uniform-phase waiting mean. This is the finite exact form of the missing
clock coordinate, not an assertion that primes form a periodic word.

The connection to the tournament calculation is variance after changing
the sampling rule: in (2) it is variance of gaps, while in TREE_RECURRENCE
it is variance of conditional Hamiltonian means across retained sources.
These are two different instances of an explicit second-order correction.

## 4. Arbitrarily long composite candidate prefixes by CRT

For any integer L>=1 choose2L distinct odd primes q_(k,e)>L, one for each
1<=k<=L and sign e in{-1,+1}. Impose

    m=0 mod2,  m=-e/k mod q_(k,e).                     (3)

All inverses exist, so CRT gives an infinite arithmetic progression of m.
Take m large enough that km+e>q_(k,e). Then every km+e for k<=L is composite.
Consequently p_m> Lm+1. Dirichlet ensures some later hit exists, since the
two target residues are units. This produces certified arbitrarily long
prefixes in the multiple clock, with every excluded candidate supplied
its divisor. `verify.py` gives a concrete L=6 certificate.

The construction does NOT prove that tau_m/(phi(m)/2) is unbounded: CRT
also enlarges m and changes the conversion between candidate and prime
clocks. A covering certificate for a finite initial segment is not a
cover of every subsequent k, and is not a Fermat-compositeness theorem.

## 5. Recovering older 227 and 1451 without equating their roles

The least-used relevant sidecar is
[TWO_RUN, R4](../overnight_rank_return_20261029/TWO_RUN.md). Its actual
Paley7 transporter determinant has a constituent norm

    N_-=2^6*7*1451,  1451=38^2+7.

The equality here is exact, so the arithmetic entry p=1451 has two named
sources: a quadratic norm from a21-element transporter representation,
and the first prime in the two progressions modulo290. They meet on

    5m-a^2=6,  (a,m)=(38,290).                         (4)

More generally any even a=+/-2 mod5 gives even m=(a^2+6)/5 and
a^2+7=5m+1. Whether this candidate is prime and is the FIRST signed hit
requires additional arithmetic; (4) guarantees neither. This is a useful
searchable family with a stated failure boundary, not a claimed spectral
explanation of prime onset.

[The old source splice](../collatz_atom_route_memory_20261003.md) has
151->227->341->1 under accelerated3x+1, valuations(1,1,10), and the family
n_k=(4^(9k+6)-19)/27 with n_(k+1)=262144n_k+184471.
Here227 is a source value, not a prime-index clock. The separate
[reversal study](../collatz_five_mirrors_20260929.md) has the3x+13 cycles

    (227,347,527,797,601), valuations(1,1,1,2,3);
    (259,395,599,905,341), valuations(1,1,1,3,2).

Reversal up to rotation exchanges the valuation words. Both lengths are5,
valuation sums8, cycle denominators2^8-3^5=13, and element sums2499.
`sidecars.py` replays these exact sources and clocks. Numerical227 alone
transports no prime-arrival theorem into Collatz.
