# The Barvinok–BIBD–QR Triangle

**opus-2026-04-05-S27**

## Three Frameworks, One Structure

Three apparently independent mathematical frameworks converge on the same object:

### 1. Barvinok (Analytic)
**Theorem (Barvinok 2014):** For matrices with ε ≤ a_ij ≤ 1:
  (1/2r)(ε/n)^r × per A ≤ ham A ≤ per A

The permanent (counting all cycle covers) approximates the Hamiltonian permanent (counting Hamiltonian cycles) up to a sub-exponential factor. The proof uses **cycle merging** — a 2-edge swap that combines two cycles into one.

For tournaments: a cycle cover requires all cycle lengths ≥ 3 and odd (no self-loops, no 2-cycles by antisymmetry). At n ∈ {3, 5, 7}, the ONLY partition of n into odd parts ≥ 3 is (n) itself, so **per = ham exactly** — every cycle cover is a Hamiltonian cycle.

Starting at n = 9 (where 3+3+3 = 9), multi-cycle covers become possible and per can exceed ham.

### 2. BIBD (Combinatorial)
For the Paley tournament P_p (p ≡ 3 mod 4 prime), the c₃ = p(p²-1)/24 directed 3-cycles form a **2-(p, 3, (p+1)/4) balanced incomplete block design**.

Every pair of vertices participates in exactly λ = (p+1)/4 directed 3-cycles. This uniformity follows from the Jacobi sum J(χ,χ) being constant over all pairs — a consequence of the QR structure.

Verified: P_7 has 14 3-cycles with λ = 2 uniform across all 21 pairs.

The BIBD property connects to H-maximization through the OCF: uniform 3-cycle distribution → maximum independence number in the odd-cycle conflict graph Ω(T) → maximum H(T) = I(Ω(T), 2).

### 3. QR Resonance (Number-Theoretic)
**THM-305:** v_2(T(n)) = (n-1)/2 = |QR_p| for odd n.
**THM-306:** 2^{C(p,2)} ≡ (2/p) mod p².

The Euler criterion 2^{(p-1)/2} ≡ (2/p) mod p creates a resonance between the binary structure of tournaments and the multiplicative group of F_p.

## The Triangle

```
           QR resonance
          /            \
    (2/p) via           v_2 = (n-1)/2
    Euler criterion     via Burnside
        /                    \
    BIBD ——————————————— Barvinok
    (Jacobi sum)        (per ≈ ham)
    (uniform λ)        (cycle merging)
```

**BIBD ↔ QR:** The Paley construction uses QR_p as the arc set. The BIBD property is the Jacobi sum J(χ,χ) being constant. Both are consequences of (Z/pZ)* acting on F_p.

**BIBD ↔ Barvinok:** BIBD = uniform local structure ≈ "δ-balanced" in Barvinok's sense. The Paley adjacency matrix scales to J/n (the uniform doubly stochastic matrix), giving the tightest possible van der Waerden bound.

**QR ↔ Barvinok:** The permanent per(A_T) = Σ_{σ: all odd, ≥3} Π a_{iσ(i)} and the Burnside sum = Σ_{σ: all odd} 2^{orb(σ)} / n! both sum over all-odd permutations. E_T[per(A_T)] = N_≥3(n) / 2^n connects them.

## Cycle Merging ↔ Wiggly Lines

Barvinok's proof technique and our wiggly framework are the same operation viewed from different angles:

| | Barvinok | Wiggly |
|---|---------|--------|
| **Space** | Cycle covers of T | Tournaments on n vertices |
| **Operation** | Merge 2 cycles (2-edge swap) | Flip 1 tile (1-arc swap) |
| **Effect** | Reduces cycle count by 1 | Changes tournament class |
| **Bound** | P(τ) ≥ ε² P(σ) | Each tile flip is weight-neutral |

Both are elementary moves in the edge space of K_n acting on binary structures.

## Key Computational Findings

1. **per = ham at n ≤ 7:** For n ∈ {3, 5, 7}, the only odd partition with parts ≥ 3 is (n). Every cycle cover is Hamiltonian. The per/ham distinction only matters starting at n = 9.

2. **P_7 is BIBD and H-max:** The Paley tournament P_7 has BIBD 3-cycles (λ=2) and achieves H_max = 189. The unique n=5 regular tournament does NOT have BIBD (5 ≡ 1 mod 4, no Paley).

3. **The BIBD defect predicts H:** Among regular tournaments, smaller BIBD defect δ(T) = max|λ(u,v) - (n+1)/4| correlates with larger H.

## Open Questions

1. At n=9 (first odd n where per can exceed ham): does the Paley tournament P_9... wait, 9 is not prime. What about n=11? Does P_11 minimize per/ham among regular 11-tournaments?

2. Is the BIBD→H-max implication provable using the OCF? The argument would be: uniform λ → maximum α(Ω) → maximum I(Ω, 2) → maximum H.

3. The partition constraint (per = ham for n ≤ 7) is a purely combinatorial fact about integer partitions. Is there a number-theoretic interpretation? The values n = 3, 5, 7 are the first three odd primes — coincidence?
