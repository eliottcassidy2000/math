# Recovering the 3-smooth decoder and sharpening its exact complexity

The November 13 [twenty-lane inventory and decoder](https://github.com/eliottcassidy2000/math/blob/0423b7f9b4a4f8126dd74d05e2f760bced205e0a/05-knowledge/results/two_three_problems_census_20261113/README.md)
already proves that every positive integer is a sum of distinct
2^a3^b terms forming a divisibility antichain. This addresses the old
request to prove the **arithmetic** lemma. It does not permit choosing
these summands as Hamiltonian frames.

For even n divide by 2 and then double the answer. For odd n select
the largest 3^b≤n, remove it, and halve the remaining even part.
Every smaller doubled term is less than 2·3^b and even, so cannot be
divisible by 3^b; being even it cannot divide 3^b either. Continue
until zero. Across output terms, a strictly increases and b strictly
decreases.

## The greedy choice has an exact structural meaning

In *any* divisibility-antichain representation there is at most one
odd term, because powers of 3 are comparable. So for odd n there
is exactly one such term; for even n there is none. More generally,
order an antichain by strictly increasing powers of 2; the powers
of 3 then strictly decrease, and the smallest power of 2 is exactly
v₂(n), since all other summands are divisible by a higher power.
This shows the largest first 3-exponent allowed by size is the
lexicographically largest attainable one (the existence proof
realizes it); after removing it and halving the same argument
applies recursively. For a given descending 3-exponent list the
2-exponents are unique if realizable: repeatedly read the valuation
of the positive remainder and peel off its forced lowest-valuation
term. Thus the greedy decoder selects a canonical antichain among
possibly many, not just a period of base representations.

## The exact worst-case size of this canonical representation

Let f(n) be the number of summands emitted by this decoder.
For every integer L≥1,

    min{n : f(n)≥L} = 3^L−2^L.

For L=1 the minimum is 1. A minimum n for L>1 must be odd,
since even n has the same f as n/2. Write n=3^k+2m with
0<m<3^k, as selected by the decoder. The smallest possible
m is m_*=3^(L−1)−2^(L−1). It lies in
[3^(L−2),3^(L−1)), so the least k with 3^k>m_* is L−1.
For any larger m both m and this least k are nondecreasing.
Hence n is minimized at 3^(L−1)+2m_* =3^L−2^L, which
indeed gives the prescribed greedy depth. This is about worst
case for the greedy normal form, not about minimum representation
length over all antichains: 19=16+3 is shorter than the greedy
19=9+6+4.

## The physical side of the decoder remains the hardest step

The [two-hub three-slot forest inventory](https://github.com/eliottcassidy2000/math/blob/962d1028aeb5a2d1448fdaed5e8c12eba71cf963/05-knowledge/results/session335_three_slot_forests_pythagorean_addresses_20261114/inventory/TWO_THREE_OPEN_PROBLEMS.md)
and [frame-lineage waves](https://github.com/eliottcassidy2000/math/blob/0423b7f9b4a4f8126dd74d05e2f760bced205e0a/05-knowledge/results/chanpak_waves_20261113/README.md)
show why these coefficients are relevant: an admissible frame assigns
vertices with two or three allowed runs, giving the multiplicative
factors of 2 and 3. But the same endpoints and actual tournament
must supply the *entire* compatible frame family. For example the
canonical arithmetic recipes for 7 and 21 are 3+4 and 9+12.
They do not contradict the proved tournament exclusions.
The new exact bound specifies how many pieces this particular
target-first decoder demands; it is not a theorem that those pieces
are physically selectable or that all count realizations need
that many frames.
