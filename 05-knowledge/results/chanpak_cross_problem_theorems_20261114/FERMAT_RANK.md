# A physical bank of cyclic sources with unbounded square-class rank

The [prime-wrap obstruction](PRIME_WRAP.md) requires unbounded
square-class rank of the source ratios Q/h for a universal prime or
7p strategy using cyclic wraps. Here is an explicit, unconditional
way to obtain that rank, without assuming there are infinitely many
Fermat primes.

For m≥1 take the strong tournament D_m=C3[b,a,TT_m] with two
singleton portals. Every directed cycle goes through both portals,
then uses an increasing nonempty subset of TT_m. It is odd exactly
when that subset has odd size. In particular two cycles never appear
together in a disjoint packing. The exact profile is

    h_m = 1+2^m,
    Q_m = 2^(m+2)+2(3^m−1).

For h use the [odd-cycle formula](https://arxiv.org/abs/2307.05569):
the 2^(m−1) odd subsets contribute weight two each. For Q use the
same formula on each of the two labeled path slots. The empty
cycle packing lets each of m+2 vertices choose a slot. For a
fixed odd subset of size s the corresponding cycle goes into one
of two slots, gives weight two in that slot, and the other m−s
vertices choose freely. Thus

    Q_m = 2^(m+2) + 4 sum over odd s of binom(m,s)·2^(m−s)
        = 2^(m+2) + 2(3^m−1).

Now set m=2^j, j≥1. Since h_m≡2 mod 3 and 3∤h_m, some prime
q≡2 mod 3 occurs to odd valuation in h_m. The order of 2 modulo
q is exactly 2^(j+1): its 2^j-th power is −1 and the possible
orders divide the next power of 2. In particular q≡1 mod 4 and
q≥2^(j+1)+1. Thus q≡5 mod 12 and 3 is a quadratic nonresidue
modulo q, by quadratic reciprocity. Reduction gives

    Q_m ≡ 2(3^m−3) mod q.

If this were zero, 3^(m−1)=1 mod q, so 3 would have odd order
and would be a square modulo q. Contradiction. Therefore the ratio
Q_m/h_m has odd valuation at q, with q absent from its numerator.

Choose the indices j successively large enough that 2^(j+1)+1 exceeds
every prime divisor of all earlier selected profiles. Each new ratio
has an odd valuation at a prime absent from all preceding ones. The
resulting triangular valuation matrix proves that the selected ratios
are linearly independent in the group of positive rational square classes.
This remains true if the class of 2 is adjoined: the witnesses q are odd.

This is a structural use of Fermat numbers even if every sufficiently
large Fermat number is composite. Their prime factors still supply new
torsion tests. What remains unproved is the *dual* condition: given
a prescribed prime target p, find a source whose ratio has the required
Legendre sign **at p** and whose actual positive values h and Q satisfy
the exact target equation. Unbounded square-class rank supplies fresh
directions; it does not automatically span the right directions at
each prescribed prime, or set the absolute height.

## A density-one consequence without a density-one count theorem

Allow arbitrary finite ordinal words built from this bank and
singletons as sources of the outer cyclic wrap. For a bank of the
first k selected ratios, the class of 2 gives an additional independent
quadratic character. An odd prime outside the finite profile support
that divides no outer output must see all k+1 classes as squares:
the exact subgroup criterion says every generator has odd order,
which implies it is a square even for primes 1 mod 4.
By the prime number theorem in arithmetic progressions (or the
multiquadratic splitting law), the condition has density 2^(-k-1)
among the primes. Thus **a finite bank can make the upper density of
missed prime divisors arbitrarily small**, even though each finite
bank still misses a whole prime progression. Taking the countable
bank makes the density of missed prime divisors zero. This is
unconditional, and is a density result only for divisibility.
It does not imply a single additional prime is realized *as* a count.

The negative and positive results suggest a concrete next test:
try to make the subset of ratios with a prescribed negative Legendre
sign carry a bounded-height inverse for that target. A character sieve
can choose a candidate family; only a separate height argument can
turn it into exact equality.
