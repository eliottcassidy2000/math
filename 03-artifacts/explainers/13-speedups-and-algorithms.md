# Speedups and Algorithms: Making the Math Computable

## The Problem of Scale

Tournaments grow fast. With n players, the number of possible tournaments is 2^{n(n-1)/2}:

| Players | Possible tournaments | Distinct types |
|---|---|---|
| 5 | 1,024 | 12 |
| 7 | 2,097,152 | 456 |
| 10 | ~35 billion | ~9.5 million |
| 19 | a 41-digit number | incomprehensibly many |

Computing H(T) by brute force — literally listing all paths — requires checking every possible ordering of n players: that's n! orderings (over 3 million for just 10 players, 10^17 for 19 players). This approach quickly becomes impossible.

This project developed several techniques that make the computation tractable where brute force cannot.

---

## Speedup 1: The OCF Formula Reduces the Problem

**The old approach**: To compute H(T), check all n! orderings of players and count those that form valid paths. Scales as n! — grows like a factorial.

**The OCF approach**: Compute H(T) = I(Ω(T), 2) — evaluate the independence polynomial of the conflict graph. The conflict graph Ω(T) has one vertex per odd cycle, not per player. For highly structured tournaments (like Paley), there are far fewer odd cycle types than there are player orderings.

**The speedup**: Instead of checking 10^17 orderings for a 19-player tournament, you enumerate the odd cycles (a much smaller set for structured tournaments like Paley), build the conflict graph, and evaluate the independence polynomial.

For the Paley tournament T₁₉:
- Brute force: 19! ≈ 10^17 operations — infeasible
- OCF on Ω(T₁₉): the circulant structure of T₁₉ means many odd cycles are identical up to rotation, dramatically reducing the conflict graph size

The theoretical speedup depends on how dense Ω(T) is — for general tournaments it can still be hard, but for symmetric families like Paley it's dramatically faster.

---

## Speedup 2: Sparse Matrix Representation (42 GB → 1.2 MB)

**The problem**: Computing H(T₁₉) (Paley on 19 players) via matrix methods requires building a constraint matrix. In dense format, this matrix takes approximately **42 gigabytes** of memory — far beyond what a typical computer can hold.

**Why the matrix is mostly zeros**: The constraint matrix encodes which sequences of players can be extended to full Hamiltonian paths. Most sequences can't be extended (most partial orderings violate the tournament arc directions), so most entries are zero. The matrix is **sparse**.

**The CSC format**: CSC stands for Compressed Sparse Column. Instead of storing all n² entries (mostly zero), store only:
1. The values of the non-zero entries
2. Which row each non-zero entry is in
3. Where each column starts

For a matrix with 42 GB of zeros but only a few megabytes of non-zero entries, CSC reduces the storage to approximately **1.2 MB** — a factor of 35,000× compression.

**The result**: H(T₁₉) = 1,172,695,746,915 was computed exactly on ordinary hardware, by exploiting the structure of the tournament (most transitions are zero in the Hamiltonian path constraint matrix).

**Why this generalizes**: Sparse matrix methods are a cornerstone of scientific computing — used in weather simulation, structural engineering, machine learning, financial modeling. The specific application to tournament H(T) computation at Paley primes is new, but the technique is broadly applicable wherever matrices have special structure.

---

## Speedup 3: Chinese Remainder Theorem for Exact Large Integers

**The problem**: H(T₂₃) = 15,760,206,976,379,349. This is a 17-digit number. Computing it requires tracking very large integers throughout the computation. Large integers are slow: arithmetic on a 17-digit number takes many times longer than arithmetic on a 6-digit number.

**The CRT approach**:
1. Pick several prime moduli: p₁ = 10007, p₂ = 10009, p₃ = 10037, ...
2. Compute H(T₂₃) mod p₁, H(T₂₃) mod p₂, H(T₂₃) mod p₃ separately. Each computation uses only small integers (all results fit in 5 digits). Fast.
3. Use the Chinese Remainder Theorem to combine these small results into the exact large value.

**Why this works**: The Chinese Remainder Theorem guarantees that if you know a number modulo several coprime moduli, and the number is smaller than their product, you can reconstruct the number uniquely. H(T₂₃) has 17 digits, so using about 3-4 primes each of size ~10,000 gives a product ~10^16 > H(T₂₃), enough to reconstruct exactly.

**The speedup**: Each modular computation runs roughly 10,000× faster than the corresponding exact-integer computation (since all numbers stay small). With 4 primes, you do 4 fast computations instead of 1 slow one — still a massive net speedup.

**Connection to number theory**: CRT is a theorem about the structure of rings. The insight that H(T) can be recovered from its residues mod small primes comes from knowing the *size* of H(T) beforehand (via the Hadamard-type bound H(T) ≤ n!/2^{n-1}).

---

## Speedup 4: Circulant Structure and Block Diagonalization

**The structure**: Paley tournaments and cyclic interval tournaments are **circulant** — they look the same if you rotate all player labels by 1. This means they commute with the cyclic permutation matrix C (shift all labels by 1 mod p).

Any matrix that commutes with C can be **block-diagonalized** using the **Discrete Fourier Transform** (DFT). Instead of one large matrix computation, you get p independent small matrix computations (one per DFT frequency), each of size roughly 1/p of the original.

**For Paley T_p**: The path constraint matrix block-diagonalizes into p independent pieces. Each piece can be solved independently and in parallel. The solutions are then combined by inverse Fourier transform.

**The speedup**: A computation that would take p³ operations (for a p × p matrix) becomes p · (p/p)³ = p operations — a factor of p² improvement. At p = 23: roughly 500× speedup.

**The same technique works for path homology**: A February 2026 paper (arXiv:2602.04140) independently discovered that circulant digraphs' path homology computation block-diagonalizes the same way. This confirms the technique is fundamental to circulant tournament analysis, not a coincidence.

---

## Speedup 5: Real-Rootedness of I(Ω(T), x)

**The observation**: The conflict graph Ω(T) is always claw-free (no vertex has three mutually non-adjacent neighbors). Chudnovsky and Seymour proved in 2007 that independence polynomials of claw-free graphs have only real roots.

**Why this helps computationally**: A polynomial with all real roots p₁, p₂, ..., p_k can be written as:

I(G, x) = (1 + p₁x)(1 + p₂x)···(1 + p_k x)

This **factored form** is much more numerically stable to evaluate than the standard coefficient form. When you want I(Ω(T), 2), you:
1. Find all roots (which are real, so they're ordinary numbers — no complex arithmetic needed)
2. Evaluate (1 + 2p₁)(1 + 2p₂)···(1 + 2p_k)

Each factor is a simple product. No catastrophic cancellation of large intermediate values. The computation is stable and fast.

**A practical consequence**: You can compute H(T) approximately (by computing roots numerically) much faster than exactly, and the real-rootedness ensures the approximation is reliable. For very large tournaments where exact computation is impossible, this gives a reliable estimate.

---

## Speedup 6: Deletion-Contraction Recursion

**The idea**: H(T) satisfies a recursion. Pick any arc (A→B) in T. Consider two smaller tournaments:
- T₀: T with arc (A→B) as is — count paths through that arc
- T₁: T with arc (A→B) reversed — count paths through the reversed arc

H(T) = H(T₀) + H(T₁)

By applying this recursion, you can build a binary tree of smaller and smaller tournaments, eventually reaching base cases (trivially small tournaments) where H is known.

**The complexity**: This gives an O(2^n) algorithm — exponential, but with base 2 rather than the n! factorial of brute force. For n = 23: 2^23 ≈ 8 million operations vs. 23! ≈ 10^22 operations. An improvement of 10^16.

**Why it's useful**: The recursion can be combined with memoization (caching repeated subproblems) and symmetry exploitation (T₀ and T₁ may be isomorphic to already-computed cases). For Paley tournaments, the circulant symmetry dramatically increases the frequency of repeated subproblems.

---

## What These Speedups Enable

| Tournament | H(T) | Feasible without speedups? | Technique used |
|---|---|---|---|
| T₇ | 189 | Yes (n=7, brute force works) | — |
| T₁₁ | 95,095 | Marginal | OCF reduction |
| T₁₉ | 1,172,695,746,915 | No (42 GB dense matrix) | Sparse CSC |
| T₂₃ | 15,760,206,976,379,349 | No (17-digit exact integer) | CRT + modular |

The combination of these techniques moves the frontier of computability from n ≈ 11 (marginal) to n = 23 (exact), an increase of 12 in the tournament size.

---

## Broader Applicability

Each technique applies far beyond tournament math:

- **Sparse CSC matrices**: Standard in computational physics, machine learning, large-scale optimization
- **CRT for exact computation**: Used in computer algebra systems (Maple, SageMath), polynomial resultant computation, elliptic curve cryptography
- **Block diagonalization via DFT**: Signal processing, quantum chemistry, lattice QCD simulations
- **Deletion-contraction**: Graph coloring, reliability polynomials, flow polynomials in operations research
- **Real-rootedness for stable evaluation**: Numerical methods for characteristic polynomials, stable computation of permanents

The tournament setting — with its explicit algebraic structure from number theory — provides a clean laboratory where these techniques can be understood precisely. The insights transfer to any domain with similar structure.
