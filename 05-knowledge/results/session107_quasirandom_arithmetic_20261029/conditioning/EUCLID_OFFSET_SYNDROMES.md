# Euclid-product offsets create exact character correlations

**PROVED elementary identities and scoped consequences; FINITE-EXACT controls.**
The Sylvester and Fermat identities are classical and inherited. This note
assembles their factor-sign constraints as a precise boundary for the recovered
prime-tournament quasirandomness theorem. No literature-priority or infinitude
of Sylvester/Fermat prime terms is claimed.

## 1. A common character mechanism

Let p_1,...,p_r be distinct odd primes, none dividing the nonzero integer c.
Suppose N>1 satisfies N=c modulo every p_i, and write N=product q_j^e_j.
Define its factor-sign column v(q_j) in F2^r by

    (-1)^v(q_j)_i = (q_j/p_i).

All q_j are coprime to every p_i. Multiplicativity of the Legendre symbol gives

    sum_j (e_j mod2) v(q_j) = v(c).                 (1)

This is an exact affine syndrome, not independence of factor columns. If N is
squarefree with t prime factors and c is a square modulo every p_i, the t
columns sum to zero, so their rank is at most t−1. For general exponents only
odd multiplicities contribute; squarefreeness is not silently assumed.
If the offset syndrome is nonzero, appending its column gives a homogeneous
relation, but the t factor columns need not themselves be dependent.

This identity requires only congruences, not primality of N or distribution
of its factors. It provides a clean distinction between random-looking
individual character observations and exact joint arithmetic dependence.

## 2. Sylvester: offset one, and a transitive prime-term tournament

Set s_0=2 and s_(j+1)=s_j^2−s_j+1. Induction gives

    s_j=1+product_(i<j)s_i       (j>=1).

Therefore distinct terms are pairwise coprime, and every prime factor p of an
earlier term sees each later term as1 modulo p. Equation(1) has zero syndrome.
This is the inherited Phi6/Egyptian-fraction recursion, not a new sequence law.

Since s_1=3 and x=3 modulo4 implies x²−x+1=3 modulo4, all s_j with j>=1
are3 modulo4. Select any finite set of odd Sylvester terms (indices at least1) which are
themselves prime, in increasing index order p_1,...,p_t. For i<j,

    (p_j/p_i)=1,          (p_i/p_j)=-1.

The first equality follows from p_j=1 modulo p_i; the second from quadratic
reciprocity since both primes are3 modulo4. Hence the Legendre tournament
p→q iff(p/q)=1 on these prime terms is **transitive**, oriented from later
terms to earlier ones. Every triple is transitive, rather than the cyclic
triple density1/4 of a quasirandom tournament.

This is unconditional for each finite selection. Saying it forms an infinite
non-quasirandom family would require infinitely many prime Sylvester terms,
which is not assumed. The explicitly verified terms are3,7,43,3263443.
Do not replace "prime terms" by "all prime divisors of all terms": composite
terms have only the joint constraint(1), not individually constant columns.

There is no conflict with global prime-prefix quasirandomness. The prime-term
selection is extremely sparse: s_(j+1)−1=s_j(s_j−1)>=(s_j−1)^2, with s_1−1=2,
so s_j−1>=2^(2^(j−1)). Thus only O(log log X) Sylvester terms are at mostX.
The recovered large-sieve theorem guarantees robustness for sets larger than
X^(1/2+eta), not for this selection. This sparseness argument concerns terms,
not all prime factors of arbitrarily late terms.

## 3. The 13/139 split: equal columns in one gauge, opposite in another

The first composite term is s_4=1807=13·139. Against previous odd prime roots
P=(3,7,43), its two incoming sign columns are both

    ((13/p))_(p in P) = ((139/p))_(p in P) = (+1,-1,+1).

They have binary rank1, exactly as (1) predicts. However the fixed-root
features used by a prime-indexed tournament usually have the reversed
arguments (p/q). Here13=1 modulo4 and139=3 modulo4, so reciprocity gives

    ((p/13))_(p in P) = (+1,-1,+1),
    ((p/139))_(p in P) = (-1,+1,-1).

Thus the columns are opposite in this second gauge. The distinction is
structural, not a numerical discrepancy. Furthermore13 is outside the
prime3-mod4 tournament domain: the raw Legendre relation on arbitrary primes
is not antisymmetric. These are rectangular character observations, not a
license to insert13 as a tournament vertex with the same rule.

The later composite s_6=10650056950807=547·607·1033·31051 likewise has a zero
factor-sign syndrome against every earlier odd prime factor. The verifier
checks the full displayed factorizations and their primality directly.

## 4. Fermat: offset two changes the syndrome

Let F_j=2^(2^j)+1. The classical product identity gives

    F_j=2+product_(i<j)F_i       (j>=1).

For every prime p dividing an earlier F_i, later F_j is2 modulo p. Therefore
its factor-sign syndrome is v(2), not zero in general.

For p|F_i, the order of2 modulo p is exactly2^(i+1), since
2^(2^i)=-1 modulo p. For i>=2 this forces p=1 modulo8, and (2/p)=1.
The earlier primes3 and5 instead have(2/p)=-1. Thus the two exceptional
root rows have negative syndrome while subsequent rows have positive syndrome.
For example F_5=641·6700417 has syndrome(-1,-1,+1,+1,+1) against
(3,5,17,257,65537).

The known Fermat prime vertices17,257,65537 form a symmetric all-residue
clique: both(p/q) and(q/p) are+1 for distinct vertices. This is a graph,
not the Legendre tournament on primes3 modulo4. Including5 would spoil
the all-residue statement because every later F_j is the nonresidue2 modulo5;
that exception is retained explicitly. No claim that the known Fermat primes
are all Fermat primes is used.

## 5. How this helps interpret quasirandomness

The global character matrix has cancellation at scale, while an arithmetic
construction can force its selected columns into a code or a transitive pattern.
These are compatible facts:

* squareclass relations constrain the attainable sign space exactly;
* fixed-modulus prime equidistribution mixes within that space;
* large-sieve cancellation controls many varying prime pairs;
* recursive Euclid offsets condition on relations among several prime factors.

The character and group notes supply the quantitative statements in their
respective domains. They are not interchangeable randomness assumptions.
Our factor-syndrome theorem itself is deterministic and needs no analytic input.

Closest prior source actually read: HYP-3724,
`HYP-3724-phi6-is-the-sylvester-map-egyptian-fraction-tower.md`.
The historical covering-min interpretation is retracted in HYP-3730;
we use only the sequence identity, not its old extremality claim.
`number_tournament_bridges_20261014.md` labels the separate equality
139=3^7−2^11 as numerology; no causal link to its Sylvester occurrence is proved
here. The current global theorem is P1 of
[tournament_open_problem_routes at96b41ec](https://github.com/eliottcassidy2000/math/blob/96b41ec57202ea9cc89dc301abd7e8053fbce179/05-knowledge/results/tournament_open_problem_routes_20261029/PROOFS.md).

Reproduce with `python3 check_conditioning.py` and
`python3 -O check_conditioning.py`. Standard library only, exact Euler-symbol
calculations, trial-division primality and explicit non-removable checks.
