# The Mersenne seam reads an arbitrary two-cover profile

The [Mersenne seam](MERSENNE.md) does not need its block B to be transitive. Write h=H(B), let P=P₂(B) count ordered pairs of possibly empty paths partitioning B, and let R=P−2h. Instead of 1 and 2^m−2, the one-run and two-run allocations contribute h and R. The same at-most-two-runs proof therefore gives

    a = 129h+64R,
    b = 190h+96R,
    H(assembled seam) = a+4b = 7(64P−h).

The terminal-left Hamiltonian frame still survives: replace the old slot on it by any Hamiltonian path of B. The assembled tournament is strong by the same closing arc. Thus this seam is a uniform **physical** reading of an arbitrary pair (H,P₂), not just a recurrence on transitive blocks.

Appending an ordinal singleton to B leaves h unchanged and doubles P: it can be attached to the end of either slot, counting the empty slot as a legitimate choice. Consequently q=64P−h changes by the exact recursion q↦2q+h; h=1 is exactly the Mersenne recursion. For a nontransitive example, B with outgoing masks (22,28,24,17,0) has (h,P)=(5,64) and gives q=4091, prime, with a literal strong count 28637=7·4091. This is not on the Mersenne ray. The replay also checks q=104059 from the profile (133,1628) and checks the unscaled quotients are prime by trial division.

For a **prescribed** q, the seam would require P=(q+h)/64. At block order n≥2, P≡0 (mod 4) by complement pairing and oddness, so h≡−q (mod 256). Splitting paths gives P≥(n+1)h and hence h≤q/(64n+63). These are exact constraints on which block to choose. No theorem here says that this line, or its union with the [two aliquot readers](HETEROGENEOUS.md), meets a physical pair for every prime.

Replay: `python replay_mersenne_modules.py`.
