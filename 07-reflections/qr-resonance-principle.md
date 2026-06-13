# The QR Resonance Principle

**opus-2026-04-05-S26**

## The Meta-Theorem

Every Burnside counting problem with a non-square base exhibits quadratic residue resonance.

Specifically: let G be a finite group acting on a set X of size m. The orbit count under b-colorings is N = (1/|G|) sum_g b^{orb(g)}. For any prime p | |G| with gcd(b,p) = 1:

1. The labeled count b^m mod p² is determined by the Legendre symbol (b/p)
2. The p-cycle contribution encodes (b/p) through the Euler criterion
3. When the p-cycle uniquely minimizes v_p, the orbit count's p-adic valuation is determined by a single number-theoretic quantity

This is proved for tournaments (THM-305/306/307) and general labeled counts (THM-308). The principle is a consequence of a single identity: **b^{(p-1)/2} ≡ (b/p) mod p** (Euler's criterion).

## Why Tournaments Are Special

Among all Burnside counting problems with base 2, ONLY tournaments have a clean v_2 formula. The reason is the **odd-cycle restriction**: tournaments count only permutations with all-odd cycles. This restriction ensures:

1. The n-cycle (which has the smallest orbit count) is always present
2. No even-cycle permutation can interfere at the minimum v_2 level
3. The n-cycle's v_2 contribution (n-1)/2 uniquely dominates

For **graphs** (no odd-cycle restriction): even-cycle permutations contribute and pollute the v_2 picture. The graph v_2 sequence 0,1,2,0,1,2,2,1,2,4,3,4... has no simple formula.

For **digraphs** (base 4 = 2²): 4 is always a QR, so (4/p) = 1 for all p. There is no QR discrimination — the labeled digraph count is always ≡ 1 mod p². Digraphs have **trivial** QR resonance.

The hierarchy:
- **Trivial resonance**: perfect square bases (b = 4, 9, 16, ...)
- **Nontrivial but messy**: graphs (base 2, all permutations)
- **Clean formula**: tournaments (base 2, odd-cycle restriction)

Tournaments sit at a sweet spot: the parity constraint on cycle types creates a clean algebraic structure that lets the QR resonance propagate all the way to an exact formula.

## Cross-Field Manifestations

### Statistical Mechanics
The Burnside sum IS a partition function Z(β=0) for the symmetry-reduced Ising model on K_n. The v_2 theorem constrains the "combinatorial free energy" at all temperatures. The antiferromagnetic ground state (H-maximizer = Paley tournament) embodies the QR structure directly.

### Coding Theory
Tournament identification is a "code" with persistent QR-controlled redundancy. The code rate approaches 1 as n → ∞, but the 2-adic structure v_2 = (n-1)/2 never disappears. This is an asymptotically rare type of redundancy — one that exists at the number-theoretic level rather than the information-theoretic level.

### Representation Theory
The class function f(σ) = 2^{orb(σ)} × [all odd] decomposes into S_n irreps. The trivial-representation coefficient is T(n), and its v_2 = (n-1)/2 constrains the Fourier coefficients of f. The n-cycle's dominance means T(n)'s 2-adic structure is controlled by a SINGLE conjugacy class.

### Dirichlet Characters
T(p) mod p decomposes into the Euler quotient f_p and the Wilson quotient w_p, both classical number-theoretic objects. The Legendre symbol (2/p) = χ_8(p) (the unique quadratic character of conductor 8) appears explicitly. This connects tournament counting to the L-function L(s, χ_8).

### Necklaces (Counter-example)
Binary necklaces use the SAME base (2) but a DIFFERENT group (Z/nZ instead of S_n). Their v_2 pattern is completely different: no clean (n-1)/2 formula. This shows the QR resonance depends on BOTH the base and the group action — it's not just about the number 2.

## The Fundamental Trichotomy

Every Burnside counting problem falls into one of three categories:

1. **Square-base problems** (b = k²): No QR information. b^{C(p,2)} ≡ 1 mod p² always. Examples: digraphs (b=4), functions (b=n²).

2. **Non-square base, unrestricted**: QR information exists at the labeled level (b^m ≡ (b/p) mod p²) but gets scrambled by even-cycle contributions when computing unlabeled counts. Example: graphs.

3. **Non-square base, parity-restricted**: QR information survives all the way to a clean formula for the unlabeled count. Example: tournaments.

## Open Questions

1. Are there OTHER parity-restricted Burnside problems (besides tournaments) that give clean v_2 formulas?

2. The q-adic valuations v_q(T(n)) for q > 2 have no clean pattern. Is there an explanation beyond "the q-cycle doesn't dominate for q > 2"?

3. The odd-part sequence T(n)/2^{(n-1)/2} = 1, 3, 57, 11971, 28242289, ... has terms that are sometimes prime. Does this sequence have any independent significance?

4. At Wieferich primes (1093, 3511): the Euler quotient f_p ≡ 0 mod p. What special structure does T(1093) have?
