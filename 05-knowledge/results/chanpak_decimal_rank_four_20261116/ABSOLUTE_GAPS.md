# Absolute gaps persist after joint binary/decimal steering

The [mixed-power controller](MIXED_POWERS.md) uses only the 139 physical return words whose determinant is nonzero mod 5. The exact replay finds **111 distinct integer actions** among them; our separate short decimal controller uses the same invertible-at-5 repertoire. Thus the power-sparsity obstruction extends to the actual arbitrary-precision policies, not just a simpler clock. It remains an obstruction for this architecture, not for span-3 as a whole.

Let B=13,460,883 be the [certified minimum](https://github.com/eliottcassidy2000/math/blob/3a2423f0bb11b8a63dc81b7b3616a881c101b3f6/05-knowledge/results/chanpak_eight_braid_20261116/RETURN_GROWTH.md) of H on any isolated ordered triple of old return blocks. The same lower bound holds when restricted to the invertible repertoire. For m complete blocks, independent odd-cycle packs on disjoint groups of three give H≥25B^floor(m/3). Stopping inside the next block cannot decrease H. At a boundary the exact semantic vector depends only on the composed integer action on the fixed rank-four lattice; a partial next block has only a constant-size menu of words, each evaluated from that same semantic state. Thus with m complete blocks there are at most a constant times 111^m values even allowing every partial stop. A finite collection of fixed seeds incurs only another constant.

Consequently the repertoire has at most O(X^α) distinct H≤X, where
```
α = 3 log(111)/log B = 0.8606965398...
```
(For comparison the original 168-word upper bound has exponent 0.9364369436....) This has two concrete consequences.

1. In every interval [X,2X] at large X there are Ω(X^(1−α)) consecutive odd targets all absent from the repertoire: distribute the roughly X/2 odd sites among the O(X^α) gaps cut by the attained ones.
2. For any sufficiently large even modulus M=2^s5^d, every odd residue is achievable, but the **worst smallest witness** among the M/2 residues is Ω(M^(1/α)). Indeed one needs M/2 distinct integers to witness them all. Here 1/α≈1.16185, so increasing precision cannot keep every absolute witness within a constant factor of its modulus.

Even closing under ordinal sums does not restore positive density. Let S be the set of attained values >1. For any σ>α the sum of h^(−σ) over h∈S converges by dyadic summation. Since ordinal sum multiplies exact Hamiltonian counts, the counting Dirichlet series for its product closure is dominated by the convergent product ∏_(h∈S)(1−h^(−σ))^(−1): the product counts multisets even when different multisets give the same value. The number of closure values ≤X is O_σ(X^σ), still power-sparse for α<σ<1. A move outside the return architecture or a genuinely different absolute remainder mechanism is necessary. These conclusions say nothing about arbitrary tournaments.
