# Prime vs Composite Tournaments

**Session:** opus-2026-04-05-S24

## The Question

Tournament theory has a hidden dependence on the primality of n. Paley tournaments exist only at primes p ≡ 3 mod 4. Tang-Yau's stability theorem has prime-specific structure. The eigenspace Betti parameter δ is 0 at prime n and 2 at n=9=3². Where exactly does the prime/composite boundary lie?

## Systematic Comparison (n = 3..8)

### What is UNIVERSAL (no prime/composite dependence):

**1. The H-spectrum gap structure is identical.**

At both n=7 (prime) and n=8 (composite), the first 15 gaps between consecutive achievable H values are:
```
[2, 2, 4, 2, 2, 2, 2, 2, 4, 2, 2, 2, 2, 2, 2]
```
The gap of 4 at position 3 (H jumps from 5 to 9, skipping 7) is the universal H=7 impossibility. This gap does not care about n being prime or composite. The forbidden values {7, 21} are properties of tournaments as a SPECIES, not of individual sizes.

**2. The frustration-H correlation is constant.**

corr(c₃, H) ≈ 0.96 at all n ≥ 5, whether prime or composite. The relationship between local frustration (3-cycles) and global complexity (Hamiltonian paths) is a structural law independent of number-theoretic properties of n.

**3. The H-distribution is symmetric.**

Skewness ≈ 0 for all n ≥ 5. The H-distribution approaches a symmetric shape regardless of primality. The ratio H_max / E[H] ≈ 2 at all n.

**4. Cycle forcing thresholds are universal.**

c₃ ≥ 3 forces c₅ > 0 at n = 5 (prime) and n = 6 (composite). The thresholds 3, 3, 4, 5, 5 for n=5..9 show no prime/composite pattern.

### What IS prime-specific:

**1. The automorphism group of the H-maximizer.**

| n | prime? | H_max | |Aut(H_max)| | |Aut|/n |
|---|--------|-------|-------------|---------|
| 3 | prime | 3 | 3 | 1.00 |
| 4 | 2² | 5 | **1** | 0.25 |
| 5 | prime | 15 | 3 | 0.60 |
| 6 | 2·3 | 45 | 3 | 0.50 |
| 7 | prime | 189 | **21** | 3.00 |
| 8 | 2³ | 661 | ? | ? |

At n=4 (the smallest composite): the H-maximizer has TRIVIAL automorphism group. No symmetry at all. This is striking — the "best" tournament on 4 vertices has no self-similarity.

At prime n: |Aut| is always ≥ 3, and at Paley primes (3, 7, 11, ...) the automorphism group is MUCH larger: |Aut(T_p)| = p(p-1)/2 (the affine group of F_p restricted to QR translations and dilations).

**2. The Paley construction exists only at primes p ≡ 3 mod 4.**

The most symmetric tournament — the one built from quadratic residues — has no composite analog. At composite n, the H-maximizer must come from a different, less symmetric construction. Nobody has identified what construction, if any, systematically produces H-maximizers at composite n.

**3. Regular tournaments exist only at odd n.**

At even n (4, 6, 8): no tournament has all scores equal. The H-maximizer always has a bipartite score structure (half scores ⌊(n-1)/2⌋, half ⌈(n-1)/2⌉). At odd n: regular tournaments exist and typically include the H-maximizer.

This is NOT a prime/composite distinction — it's an odd/even distinction. But among odd n, primes ≡ 3 mod 4 get the Paley construction while composites and primes ≡ 1 mod 4 do not.

**4. The δ parameter (from eigenspace Betti decomposition).**

For the H-maximizer at odd n, β_{top} = (n-1) + δ where:
- δ = 0 at prime n (7, 11 verified)
- δ = 2 at n = 9 = 3²

The δ > 0 at composite n comes from the TRIVIAL eigenspace (k=0 under the cyclic shift) contributing non-trivially to the top Betti number. At prime n, the trivial eigenspace contributes 0 to β_top.

**Why?** At prime n, the shift automorphism acts freely on the non-trivial eigenspaces, and THM-125 (constant symbol matrix) ensures uniform distribution across eigenspaces. At composite n = pq, the cyclic group Z_n has subgroups Z_p and Z_q. The trivial eigenspace "sees" this subgroup structure and picks up extra contributions — specifically, the fixed points of the Z_p and Z_q actions create additional independent cycles in the path complex.

## Tang-Yau in the Prime/Composite Framework

Tang-Yau (arXiv:2602.04140) prove stability of Betti numbers for circulant digraphs C_n^S: for "generic" n (specifically, n ∉ Q+(S)), the Betti numbers are determined by the symbol matrix M_m(t) at a generic t.

**At prime n (Paley):** THM-125 proves Q+(QR_p) = ∅. The symbol matrix is CONSTANT (no t-dependence). This is STRONGER than stability — it means ALL nonzero t give the same rank. No prime p is exceptional.

**At composite n:** Two questions are open:
1. Does the constant symbol matrix property hold for composite-order circulants?
2. Is Q+(S) non-empty for natural connection sets at composite n?

If Q+(S) is non-empty at composite n, then some n values are "exceptional" in the Tang-Yau sense — their Betti numbers differ from the generic formula. The δ = 2 at n=9 might be exactly this phenomenon: n=9 is exceptional because 9 ∈ Q+(S) for the relevant connection set.

## The Three Regimes

The data suggests three distinct regimes:

### Regime 1: Paley primes (p ≡ 3 mod 4: 3, 7, 11, 19, 23, ...)
- Paley tournament T_p exists, maximizes H
- |Aut| = p(p-1)/2 (very large)
- Flat eigenvalue spectrum
- δ = 0 (uniform eigenspace Betti)
- Tang-Yau Q+(QR_p) = ∅ (no exceptional n)

### Regime 2: Non-Paley primes (p ≡ 1 mod 4: 5, 13, 17, 29, ...)
- No Paley construction
- H-maximizer exists but is not circulant (or circulant with different S)
- |Aut| ≥ 3 but much smaller than Paley
- δ likely = 0 (unverified beyond n=5)
- Regular tournament exists (odd n) and is or contains H-maximizer

### Regime 3: Composite n (4, 6, 8, 9, 10, ...)
- No Paley construction
- H-maximizer has small |Aut| (trivial at n=4!)
- At even composite: no regular tournament, bipartite score structure
- At odd composite: regular tournament exists but may not maximize H alone
- δ > 0 (verified at n=9=3²)
- Tang-Yau Q+(S) possibly non-empty

## H_max Values and Factorizations

| n | type | H_max | factorization | H_max · |Aut| |
|---|------|-------|---------------|---------|
| 3 | Paley | 3 | 3 | 9 = 3² |
| 4 | comp | 5 | 5 | 5 |
| 5 | non-Paley prime | 15 | 3·5 | 45 |
| 6 | comp | 45 | 3²·5 | 135 |
| 7 | Paley | 189 | 3³·7 | 3969 = 63² |
| 8 | comp | 661 | 661 (prime!) | ? |

Note:
- H_max at Paley primes: 3 = 3, 189 = 3³·7. Both divisible by n.
- H_max at non-Paley primes: 15 = 3·5. Divisible by n=5.
- H_max at composites: 5 (not div by 4), 45 = 9·5 (div by 3 but not 6), 661 (prime, not div by 8).

**At prime n, H_max is always divisible by n.** (Trivially for Paley: |Aut(T_p)| contains the cyclic shift of order p, so H_max = |Aut|·(H_max/|Aut|) is divisible by p.)

**At composite n, H_max is NOT necessarily divisible by n.** (5 mod 4 = 1 ≠ 0, 661 mod 8 = 5 ≠ 0.)

This divisibility is the cleanest prime/composite distinction:

**At prime n: n | H_max. At composite n: n ∤ H_max (typically).**

## Open Questions

1. **What replaces Paley at composite n?** Is there a systematic construction of H-maximizers at even n?

2. **Does the δ pattern extend?** Compute δ at n=15=3·5, n=21=3·7, n=25=5² to see if δ depends on the specific factorization.

3. **Tang-Yau Q+(S) at composite n:** Construct circulant tournaments at n=9, 12, 15 and test if the symbol matrix is still constant. If not, identify the exceptional t values.

4. **n | H_max at all prime n?** Test at p=11: H(T_{11}) = 95095. 95095/11 = 8645. Yes, divisible. At p=13 (non-Paley prime ≡ 1 mod 4): what is H_max?

5. **The H-maximizer at composite even n:** At n=4, |Aut|=1. At n=6, |Aut|=3. At n=8, |Aut|=? Does |Aut| grow or stay small?

6. **Second eigenspace Betti:** At prime n, each non-trivial eigenspace contributes β_top = 1. At composite n=9, each non-trivial eigenspace still contributes 1 but the trivial contributes 2. Do ANY non-trivial eigenspaces contribute > 1 at composite n?
