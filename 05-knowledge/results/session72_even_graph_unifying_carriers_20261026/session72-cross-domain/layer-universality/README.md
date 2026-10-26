# Every finite set family occurs as a cyclic phase-layer design problem

**Status: PROVED elementary corollary of CA2–CA3.** Suggested by the parent after the core was frozen; this supplement changes none of its bytes. This is representation universality, **not a computational-hardness theorem**.

Let `S_1,...,S_t` be any finite indexed family of subsets of `[r]`, with r>=1. Empty or repeated subsets are allowed. Choose distinct odd primes `p_1,...,p_t` and put

    d_i=product_(j:i in S_j) p_j,
    m=product_j p_j, if t>0; m=3, if t=0.

Empty products are1. Every d_i divides the odd integer m, and m>=3. For a graph E on `[r]`, retain cyclic phase differences on its edges as in CA2.

**Corollary LU.** This observation separates actual joint cyclic phase states modulo common translation if and only if every nonempty induced graph `E[S_j]` is connected, with singletons counted connected.

**Proof.** Each d_i is squarefree, and the prime-power layer for p_j is exactly S_j at exponent1 and empty at every larger exponent. No other prime occurs. CA2's complete list of nontrivial layer conditions is therefore precisely the displayed family of induced-connectivity constraints. Empty sets impose no condition and duplicate subsets repeat an existing one. QED.

This has a genuine native tournament realization. For each i, take the binary length-d_i word `10^(d_i-1)` (the word1 when d_i=1) and repeat it to length m. Its least period is d_i: its ones occur exactly at the multiples of d_i. Apply CA3's r-bit block-colour carousel construction. It gives4rm vertices, with the simultaneous phase quotient and sparse-observation criterion exactly those in LU.

Thus general cyclic phase-query design is not limited to trees, nested/laminar supports, or one unweighted connectivity condition. Any finite family of overlapping set constraints can arise. This statement does not establish NP-hardness, a minimal-edge algorithm, a polynomial-size reduction to explicitly listed tournaments, or a generic graph-completion theorem. The arithmetic periods are succinct, while m and the native vertex count can be very large relative to the set-family input. Complexity transport would need its own size proof.

## Reproducible controls

The checker pins the exact parent CA2–CA3 theorem bytes. For r1,2,3 it tests **every set family** formed from the2^r subsets, including empty and full families:276 families total. For every graph on those r vertices it compares all actual prime-power layer constraints with the supplied set-family constraints:2,084 graph/family cases. It separately tests duplicate supports. This is a bounded control of the representation, not a computation of every corresponding enormous tournament or phase space.

Run `python3 -B check_layers.py` and `python3 -O -B check_layers.py`; output must byte-match the stored receipts. `--inject-failure` must fail in both modes. Only the standard library is used.
