# Chan–Pak as positivity of an arithmetic renewal kernel

**Status:** exact structural reformulation and independent small controls, not a new proof of Chan–Pak. The SCC multiplicity laws are inherited. This note isolates what a number-theoretic attack must establish and why local information cannot do it.

## Provenance and priority

The class series, logarithmic strong-core cutoff, and fixed-count multiplicity polynomials already occur in `05-knowledge/results/h_spectrum_multiplicities_20261021.md`, which credits Session24 and `spectrum_basins_20261020e/spectrum.md`. The unweighted converse sequence formula is also in `converse_tournament_count_channels_20261021/session22-converse-channels/README.md`. The omissions7,21 and component product reduction are recorded in THM-1370. No novelty for those formulas is claimed. We use the classical Camion theorem (every nontrivial strong tournament has a Hamiltonian directed cycle) only to give a deliberately elementary, nonoptimal cutoff.

## 1. The arithmetic object is an ordered factorization system

Let a_n(N) count tournament isomorphism classes of order n and Hamilton count N. Let s_n(N) count strong classes, with the singleton counted strong. Use a_0(1)=1. A tournament's strong components have a unique total order. Every Hamilton path visits them in that order, so its count is the product of their counts. Consequently

    a_n(N)=sum_{k=1}^n sum_{d|N} s_k(d) a_{n-k}(N/d).       (1)

This is a renewal recurrence with **Dirichlet convolution** in the count variable. It is not an Euler product: component order matters, so A→B and B→A need not be isomorphic. Formal symbols [d][e]=[de] give A(x)=1/(1−S(x)); every coefficient of x is a finite sum, so no analytic convergence assumption is hidden here.

The number-theoretic analogy is representability by an ordered multiplicative system whose atoms are actual strong tournaments. Ordinary prime factorization predicts necessary arithmetic shapes, but it cannot manufacture the atoms.

## 2. A finite cutoff, including negative answers

If a nontrivial strong tournament has order r, breaking one Hamiltonian cycle at each of its r arcs gives r distinct Hamilton paths. Thus h≥r≥3.

Given a witness for N>1, remove every singleton strong component. This does not change N. Write the remaining orders/counts as r_i,h_i. Then

    sum_i r_i <= sum_i h_i <= product_i h_i=N.

The second inequality holds for one factor with equality and for two or more integers≥3 by induction. Therefore **every realizable N>1 has a witness of order at most N**. N=1 has the singleton witness. This bound is intentionally much weaker than the inherited Busch/logarithmic cutoff.

Thus an explicit algorithm on a prescribed positive odd N is:
1. Enumerate all labeled tournaments of orders1..N.
2. Compute H by subset dynamic programming.
3. Return the first witness if one occurs; otherwise return the complete finite absence certificate.

It terminates for every N, whether or not Chan–Pak is true. It can be made least-order by enumerating orders in increasing order. It is spectacularly inefficient, roughly exponential in N² before the path DP, and does not solve the universal assertion that its negative outputs are exactly7 and21. In particular, undecidability of **individual spectrum membership** is not a plausible escape route here. No claim about independence of the universal conjecture from an axiomatic theory follows.

## 3. Exact eventual laws and the prime/7p lanes

For fixed N let b_{m,k}(N) count ordered lists of k nontrivial strong classes, total order m, with product of counts N. There are finitely many: m≤N and k≤floor(log_3 N). Inserting n−m singleton components into the k+1 gaps gives

    a_n(N)=sum_{m,k} b_{m,k}(N) binom(n−m+k,k),             (2)

where a term is0 if n<m. Hence a_n(N) is eventually a polynomial of degree at most floor(log_3 N). Its degree and leading coefficient depend on actual strong inventories, not merely on the prime factors of N. This is the inherited fixed-count multiplicity law in arithmetic form.

Now take an odd prime p, or7p with p prime other than7. A realization cannot have two nontrivial components: in the prime case factorization forbids it; in the7p case one factor would be7, already impossible. The same argument applies to49,63,343 because every nontrivial factorization contains7 or21. This statement includes the zero inventories at7 and21.

Let S_N(x)=sum_r s_r(N)x^r. For these target lanes,

    A_N(x)=sum_n a_n(N)x^n = S_N(x)/(1−x)²,                (3)
    a_n(N)=sum_{r≤n} (n−r+1)s_r(N).

Since s_r(N)=0 for r>N, for every n≥N,

    a_n(N)=alpha_N(n+1)−beta_N,
    alpha_N=sum_r s_r(N),  beta_N=sum_r r s_r(N).          (4)

Thus realization of N is equivalent to **alpha_N>0**, or to positivity of the double-pole coefficient

    lim_{x→1−}(1−x)² A_N(x)=S_N(1).

This is an exact alternate formulation, not an analytic proof: proving positivity still requires physical strong sources. It does identify why asymptotic counts in n can be deceptive. For each fixed prime target there are only two possible eventual behaviors: identically zero, or affine with positive integer slope. There is no slowly emerging third regime.

The standard multiplicative reduction of Chan–Pak asks for prime counts p≠7, sevenfold primes7p with p≠3,7, and the exceptional generating counts49,63,343. The kernel above applies to each. To see sufficiency, factor any odd target as7^a m with gcd(m,7)=1. Exponents a≥2 are generated by2 and3 using49 and343, except a=1; for a=1 pair7 with a prime factor of m other than3, or use63 when m is a power of3 of exponent≥2. Only7 and21 remain. This factorization step is inherited, and does not imply positivity of any missing kernel.

## 4. Converse carries square divisors

Let c_n(N) count self-converse classes and u_n(N) strong self-converse classes. A self-converse SCC word consists of an arbitrary left sequence, its forced reversed converse on the right, and either no central component or one self-converse strong component. Therefore

    c_n(N)= [n even, N square] a_{n/2}(sqrt N)
            + sum_{r≥1, n−r even} sum_{d²|N}
                a_{(n−r)/2}(d) u_r(N/d²).              (5)

This is a square-divisor convolution, not ordinary prime splitting. A tournament has H=1 exactly when it is transitive: a directed triangle has three paths, and iterative vertex insertion injects its paths into those of any larger tournament containing it. Thus a_j(1)=1. For squarefree N>1, d=1, hence

    c_n(N)=sum_{r≤n, r≡n mod2} u_r(N).

It is eventually constant on each parity of n. A prime/7p ordinary fibre can grow linearly while its self-converse count stays bounded. Conversely, absence of self-converse witnesses does not imply absence of ordinary witnesses. This distinguishes count arithmetic from the converse geometry relevant to the merged metagraph.

## 5. Why the analogy is not merely local–global rhetoric

For any odd prime p other than7, the abstract set

    M_p = positive odd integers except {7,21,p}

is multiplicatively closed. Indeed a product of members cannot equal p or7, and a nontrivial factorization of21 includes7. It has density one within the odds, and meets every parity-admissible residue class modulo every positive modulus, because it omits only finitely many odds.

So **multiplicative closure + density one + complete congruence coverage** is logically compatible with an additional prime hole. These properties cannot prove Chan–Pak without a genuinely tournament-specific source theorem. This is a countermodel to an inference, not a proposed tournament spectrum.

A useful number-theoretic attack would prove alpha_N>0 through a positive lower bound for an actual source representation function, or build a descent preserving endpoint/slot legality. Replacing that function by unrestricted norm representations, Goldbach pairs, or local residues discards precisely the physical condition whose positivity is at issue. The goal is an analogue of a representation theorem with a quantified compatibility loss, not merely a stronger numerical census.

## Reproduction

`python arithmetic_renewal_replay.py` independently enumerates every labeled tournament through order5, quotients by all permutations, computes strong connectivity and converse classes, and counts Hamilton paths by subset DP. It verifies(1),(5), and the exact linear kernel on the listed prime/sevenfold target controls. Totals are1,1,1,2,4,12 including order0. These checks do not replace the proofs or establish any new global coverage.
