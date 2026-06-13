# Number Theory Connections

## What is Number Theory?

Number theory is the study of the integers — 1, 2, 3, 4, ... — and the deep patterns hiding in their structure. It asks questions like:
- Which numbers are prime?
- How do primes distribute themselves?
- Which integers can be written as sums of squares?
- What happens when you do arithmetic on a "clock" (modulo some number)?

What makes number theory surprising is how often its patterns turn out to control other, seemingly unrelated problems. Tournament path counting is one of those "other problems" — number theory shows up throughout, in ways that range from the practical (building special tournaments, speeding up computation) to the philosophical (why these formulas have the shape they do).

---

## Clock Arithmetic: The Foundation

Before diving in, it helps to understand **modular arithmetic** — arithmetic on a clock.

A 12-hour clock "wraps around" at 12. 10 + 5 = 3 on a clock (not 15), because 15 mod 12 = 3. We write this as: 10 + 5 ≡ 3 (mod 12).

Similarly, you can do arithmetic mod any number n. Multiply, add, subtract — when the result exceeds n, wrap around. This is the arithmetic of integers modulo n, written Z/nZ or Z_n.

When n is prime, something special happens: every non-zero number has a multiplicative inverse (you can "divide" by any non-zero number), and the system behaves like a complete number system (a "field"). This is F_p — the field with p elements.

---

## Quadratic Residues: The First Connection

A **quadratic residue** mod p is a number that is a perfect square in the arithmetic of F_p.

In standard arithmetic, 9 is a perfect square (3² = 9). In mod 7 arithmetic: 1² = 1, 2² = 4, 3² = 2 (because 9 mod 7 = 2), 4² = 2, 5² = 4, 6² = 1. So the quadratic residues mod 7 are {1, 2, 4} — exactly half of the non-zero numbers {1, 2, 3, 4, 5, 6}.

This is always true for primes p ≡ 3 (mod 4): exactly half the nonzero elements are quadratic residues, and they distribute "evenly" in a precise sense.

### Paley Tournaments = Quadratic Residue Tournaments

The Paley tournament T_p uses this exactly: player A beats player B if A − B is a quadratic residue mod p. Since exactly half of differences are residues and half are non-residues, every player beats exactly (p−1)/2 other players. The tournament is perfectly balanced.

The mod-4 condition (p ≡ 3 mod 4) ensures this construction works: when p ≡ 3 mod 4, the Legendre symbol (−1/p) = −1, which means −1 is NOT a quadratic residue. This prevents the construction from being "anti-symmetric" in the wrong way. It's the same mod-4 condition that appears throughout Grinberg-Stanley's mod-4 refinement of Rédei's theorem.

---

## Gauss Sums: Counting Patterns with Complex Numbers

The **Gauss sum** mod p is: g = Σ (a/p) · ζ^a, summed over a = 0, 1, ..., p−1, where (a/p) is the Legendre symbol (±1 depending on whether a is a quadratic residue) and ζ = e^{2πi/p} is a complex number on the unit circle.

Gauss proved that g² = p (for p ≡ 1 mod 4) or g² = −p (for p ≡ 3 mod 4). In either case, |g| = √p.

This is why the eigenvalues of the Paley tournament's adjacency matrix include (−1 ± √p)/2 — the √p comes directly from the Gauss sum, which measures how "evenly spread" the quadratic residues are around the clock.

The same Gauss sum controls:
- The weight distribution of quadratic residue codes (Hamming, Golay)
- The eigenvalue spectrum of Paley graphs and tournaments
- The distribution of H(T) across Paley tournament families

**In plain language**: The "personality" of a prime p — how its quadratic residues are distributed — shows up directly in the spectrum of the Paley tournament and in the count H(T_p). Number theory determines tournament structure.

---

## The Mod-4 Condition: A Recurring Pattern

The number 4 keeps appearing:

1. **Rédei's theorem**: H(T) is always odd — equivalently, H(T) ≡ 1 (mod 2).

2. **Grinberg-Stanley's refinement (2023)**: H(T) mod 4 is controlled by the number of odd cycles. Specifically, H(T) ≡ 1 + 2·(# odd cycles) (mod 4).

3. **Paley condition**: T_p is self-complementary (and has maximum H) when p ≡ 3 (mod 4).

4. **Self-complementary tournament sizes**: A tournament on n vertices can be self-complementary only when n ≡ 3 (mod 4) or n ≡ 0 (mod 4).

5. **Forbidden values**: The permanently forbidden values H = 7 and H = 21 are both ≡ 3 (mod 4). Any value ≡ 1 (mod 4) seems achievable.

This isn't coincidence. The "mod 4" pattern is fundamental:
- The parity of H(T) (always odd) is a mod-2 statement
- Self-complementarity requires a square-root of −1 in the automorphism group, which exists only for sizes ≡ 3 (mod 4) — the same condition as for Paley tournaments
- The forbidden values being ≡ 3 (mod 4) reflects deeper congruence constraints on independence polynomial evaluations

Number theory's arithmetic mod 4 is the structural backbone of tournament path counting.

---

## Quadratic Reciprocity: A Deeper Connection

**Quadratic reciprocity** is one of the most celebrated theorems in all of mathematics (Gauss proved it at age 19 and called it the "golden theorem"). It says:

For two distinct odd primes p and q, (p/q) · (q/p) = (−1)^{(p−1)(q−1)/4}

In words: whether p is a perfect square mod q is related to whether q is a perfect square mod p, in a specific alternating way.

This enters tournament theory through the structure of Paley tournaments at different primes:
- Why T₇ and T₁₁ are the maximizers for 7 and 11 players, but the cyclic tournament takes over at p ≥ 13 — this crossover depends on how the quadratic residues of larger primes distribute, which is governed by quadratic reciprocity.
- The "crossover" between Paley and cyclic interval tournaments happens precisely because the quadratic residue structure becomes "less efficient" at packing vertex-disjoint odd cycles as p grows.

---

## The Chinese Remainder Theorem: A Practical Speedup

The **Chinese Remainder Theorem (CRT)** says: if you know a number x modulo several pairwise coprime moduli m₁, m₂, ..., m_k, you can reconstruct x exactly (as long as x is smaller than m₁ × m₂ × ... × m_k).

**In practice**: instead of computing H(T) directly as a huge integer (which can have hundreds of digits), you compute H(T) mod p₁, H(T) mod p₂, H(T) mod p₃, ... for several small primes. Each modular computation is fast. Then CRT combines the results to recover the exact value.

**Why this helps for tournaments**: H(T₂₃) = 15,760,206,976,379,349 — a 17-digit number. Computing this directly requires tracking huge intermediate values. But computing H(T₂₃) mod 10007, mod 10009, mod 10037, ... requires only small-integer arithmetic. CRT then stitches the small results back into the exact answer.

**The number theory behind it**: CRT works because the integers modulo m₁×m₂ factor completely as a product of rings. This is the "Chinese" structure: independent "clocks" at different prime sizes can be combined perfectly.

This speedup is why the project could compute H(T₁₉) = 1,172,695,746,915 and H(T₂₃) = 15,760,206,976,379,349 exactly — numbers too large to compute by brute force, but tractable via modular arithmetic.

---

## Sparse Matrix Representation: A Different Speedup

For larger Paley tournaments (like T₁₉ with 19 players), the matrix used to compute H(T) naively would require 42 gigabytes of memory — too large for any ordinary computer.

The key insight: the constraint matrix for Hamiltonian path counting is extremely **sparse** — most of its entries are zero. The matrix has 2^18 × 18! potential entries, but the vast majority are zero because most combinations of orderings violate the tournament structure.

By storing only the non-zero entries (in a format called CSC — Compressed Sparse Column), the representation shrinks from 42 GB to approximately 1.2 MB. The computation becomes feasible on a laptop.

This is the same idea behind all large-scale scientific computing: weather simulations, financial models, protein structure prediction. Most real-world matrices are sparse. Sparse representation is standard in computational science, but its specific application to tournament Hamiltonian path counting at Paley primes is new.

---

## The Krawtchouk Connection to Coding Theory

**Quadratic residue codes** are a family of error-correcting codes built from — exactly — the quadratic residues mod a prime p.

The two most famous:
- **Hamming [7,4,3] code**: transmits 4 bits using 7, can correct 1-bit errors. Built from QR mod 7.
- **Golay [23,12,7] code**: transmits 12 bits using 23, can correct 3-bit errors. Built from QR mod 23.

These are "perfect" codes — they use redundancy with optimal efficiency. No other codes of those lengths can correct as many errors.

The Paley tournaments T₇ and T₂₃ are built from the same quadratic residues. In the **Krawtchouk coordinate system** (a coordinate system for patterns on binary strings, standard in coding theory), the tournament matrices for T₇ and T₂₃ have the same "spectral fingerprint" as the parity-check matrices of the Hamming and Golay codes.

This means: **the same number-theoretic structure that makes a code perfect is what makes a Paley tournament maximize H(T)**. The optimality in each domain — error correction for codes, Hamiltonian path count for tournaments — springs from the same source: the perfectly balanced distribution of quadratic residues in F_p.

The Krawtchouk coordinate system was developed in the 1940s by the Ukrainian mathematician Mykhailo Krawtchouk specifically for problems about binary strings. It was later adopted by coding theorists (Delsarte 1972, MacWilliams-Sloane 1977). Its appearance in tournament theory is new.

---

## Primes Modulo 4: A Table of Connections

| Prime p | p mod 4 | T_p self-complementary? | T_p maximizes H? | Corresponding code |
|---|---|---|---|---|
| 3 | 3 | Yes | Yes (trivially) | — |
| 7 | 3 | Yes | Yes (verified) | Hamming [7,4,3] |
| 11 | 3 | Yes | Yes (verified) | — |
| 19 | 3 | Yes | Likely yes | — |
| 23 | 3 | Yes | Likely yes | Golay [23,12,7] |
| 5 | 1 | No | N/A | — |
| 13 | 1 | No | N/A | — |

The pattern is exact: Paley tournaments exist exactly for primes p ≡ 3 (mod 4), which are also the primes for which the Gauss sum is pure imaginary (g² = −p), which is the condition for the quadratic residue code to be self-dual.

---

## Summary: How Number Theory Drives Tournament Theory

| Number-theoretic object | What it controls in tournaments |
|---|---|
| Quadratic residues mod p | Which player beats whom in Paley T_p |
| Gauss sums (sum of Legendre symbols × roots of unity) | Eigenvalues (−1 ± √p)/2 of Paley tournament |
| p ≡ 3 (mod 4) condition | Self-complementarity of T_p; existence of quadratic residue codes |
| Quadratic reciprocity | Crossover between Paley and cyclic tournament optimality |
| Chinese Remainder Theorem | Speedup for computing exact H(T_p) at large primes |
| Krawtchouk polynomials | Spectral "fingerprint" connecting T_p to the Hamming/Golay codes |
| Mod-4 arithmetic structure | Forbidden values (7, 21 ≡ 3 mod 4); self-complementary tournament sizes |

Number theory is not a tangential tool borrowed for convenience. It is native to the problem: the tournaments with the most Hamiltonian paths are built by prime arithmetic, their spectra are governed by Gauss sums, their computational speedups rely on CRT, and their connection to perfect codes flows through the same quadratic residue structure. The integers are running the show.
