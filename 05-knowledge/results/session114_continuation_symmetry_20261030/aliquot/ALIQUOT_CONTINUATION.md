# A complete dilation quotient and its infinite one-step split

2026-10-30. **PROVED scoped**, with reproducible finite controls. Local labels A1–A4 refer only to this note. No literature-priority claim. No aliquot termination claim.

Write s(n)=sigma(n)−n for positive integers n, and let o(n)=s(n) mod 2. The legal dilation operations are multiplication by arbitrary named primes; their composites are multiplication by arbitrary positive integers. Define

- e(n)=n mod 2;
- k(n)=the squarefree kernel of the odd part of n: the product of odd primes whose valuations in n are odd;
- Q(n)=(e(n),k(n)).

The inherited Session 113-adjacent result gave the collision 25 versus 49 when an aliquot step was added to the binary dilation language. The present result classifies that entire quotient and proves that one aliquot step splits even one quotient class into infinitely many classes. This concerns actual positive sources, not formal linear mixtures.

## A1. Exact binary dilation equivalence

For all positive n,m,

    o(nm) = 1[k(n)=k(m)] XOR (e(n)e(m)).                 (1)

Indeed sigma(nm) is odd exactly when every odd-prime exponent in nm is even, equivalently k(n)=k(m). Subtracting nm modulo 2 gives (1).

Consequently, two positive sources have identical binary observations after EVERY dilation if and only if their Q-values agree.

Sufficiency is immediate. For necessity, suppose k(n) differs from k(n'). Choose the even multiplier m=2k(n). The parity term vanishes, and (1) distinguishes the sources. If the kernels agree but e(n) differs from e(n'), choose an odd prime q outside their common kernel. Then the indicator is zero and the observations distinguish their parities.

The same classification holds even if ONLY odd multipliers are allowed. Different kernels with equal source parity are separated by m=k(n); opposite source parity is separated by an odd prime outside both kernels. Such a prime exists by Euclid's theorem.

This is a deterministic quotient, not merely an invariant. Under multiplication by a prime q:

- q=2 sends (e,k) to (0,k);
- q odd sends (e,k) to (e,kq/gcd(k,q)^2).

The observation is 1[k=1] XOR e. Therefore Q is the coarsest exact deterministic quotient for all future dilation observations. This is a countable state set; no bounded-size or finite-precision representation is asserted.

## A2. An infinite split caused by one aliquot step

For each odd prime p, take the actual positive source

    n_p = 2^p,       y_p = s(n_p)=2^p−1,
    K_p = k(y_p).

Every n_p has Q(n_p)=(0,1), so all these sources are indistinguishable under EVERY dilation. Nevertheless:

1. y_p is odd and congruent to 7 modulo 8, so it is not a square and K_p>1.
2. For distinct primes p,q,

       gcd(y_p,y_q)=2^gcd(p,q)−1=1.

   The elementary identity follows by the Euclidean algorithm on exponents. Thus K_p and K_q are coprime and, since both exceed 1, distinct.
3. After the single added aliquot operation, the legal odd multiplier K_p distinguishes n_p from EVERY n_q with q≠p:

       o(K_p s(n_p))=0,        o(K_p s(n_q))=1.

   This is because multiplying y_p by its squarefree kernel makes a square, while y_q K_p is a nonsquare odd integer.

Thus the words “apply s, multiply by the named prime factors of K_p, then read o” distinguish infinitely many sources inside ONE original Q-class. All the sources are reachable from 2 by multiplication by 2 alone. Every displayed aliquot successor is positive, so the terminal value s(1)=0 is irrelevant. If a total transition system is desired, adjoin an absorbing terminal state 0 and define its observation separately; the obstruction is unchanged.

## A3. No uniformly bounded sidecar repairs this quotient

Suppose a deterministic carrier R(n)=(Q(n),c(n)) preserves the binary observer and is closed under prime dilations and an actual aliquot step. Then c(n_p) must be pairwise distinct. Otherwise deterministic updates would give the same observation on the distinguishing word from A2, a contradiction.

Therefore the fiber of any such carrier over Q=(0,1) must be infinite. In particular no finite-valued sidecar, and no uniform fixed number of additional bits, repairs closure. On any N of these sources, at least N distinct sidecar values and at least ceil(log2 N) bits are necessary.

This does NOT exclude an unbounded integer-valued sidecar: retaining n itself is an exact sufficient state. It does not show that the quotient for the enlarged language is equality of integers. It does not infer divergent aliquot trajectories.

**Alphabet boundary.** The infinite statement permits all prime-dilation letters. Its distinguishing multipliers depend on p. For a chosen finite collection of N sources, the finite alphabet consisting of 2 and all prime factors of their K_p already gives the N-way split. No infinite lower bound for one fixed finite alphabet is claimed. In fact, with r fixed odd dilation primes, binary dilation responses after a positive aliquot step have only finitely many possibilities; see A4.

## A4. Exact fixed-alphabet boundary

Let P be a finite set of r odd primes, and allow only their nonnegative powers as dilation multipliers. Let S(n) be the subset of P dividing k(n), and let b(n) indicate whether k(n) has any prime factor outside P.

If b(n)=1, every response is the constant e(n). If b(n)=0, the response to exponent-parity vector t in {0,1}^r is

    1[t=S(n)] XOR e(n).                                   (2)

For r≥2 these give exactly 2+2^(r+1) distinct response functions over all actual positive sources: two constants and two polarities of each singleton indicator. All are realized (use an outside prime for constants; use product(S) or twice that product for the singleton functions). Different singleton supports have different indicators, and when the cube has at least four points no singleton indicator is the complement of another. For r=1 the count is four, because the two singleton indicators are complements; for r=0 it is two.

Hence for a fixed finite P, adding one aliquot step followed only by P-dilations partitions ANY source class into at most these many response classes. A2 requires an expanding alphabet. This finite upper bound prevents misreading the infinite split as a fixed-alphabet impossibility.

## Inheritance, mechanism, and next question

Closest proved mechanism: binary prime-dilation responses and the 25/49 aliquot collision in [D7–D8 of the prime-dilation packet](../../prime_dilation_response_20261030/PROOFS.md). The missing coordinate is the odd squarefree kernel of the actual successor, not the source's existing dilation phase. The canonical hostile is now an entire family of even powers of 2. The retained sidecar must be unbounded on that one class when all primes remain legal.

The exact map Q preserves all dilation futures and loses the magnitudes of prime exponents and the even valuation. Adding s destroys that preservation because new factor support appears. The concrete decisive test is the successor kernel, exposed by a multiplier drawn from its prime factors.

An open next question is to classify the equivalence relation for the full mixed language, or to prove useful lower bounds for a FIXED finite multiplication alphabet with repeated aliquot steps. Neither is settled here.

## Reproduction

Run `python3 check_aliquot.py` and `python3 -O check_aliquot.py` in this directory. The checker uses explicit exception gates, direct divisor summation for observation controls, exhaustive small-source dilation tests, finite-alphabet function counts, and seven pairwise distinguishing Mersenne-source witnesses. It is evidence for the proofs above, not an exhaustive test of infinite universes.
