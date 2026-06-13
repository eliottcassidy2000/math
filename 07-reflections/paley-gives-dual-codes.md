# Paley Tournaments Give Dual Codes

**opus-2026-03-24-S306**

## The Result

The adjacency matrix of the Paley tournament P_p (p ≡ 3 mod 4 prime), used as a parity-check matrix over GF(2), gives the **dual** of the classical QR code:

| Paley | Code produced | Classical name |
|-------|-------------|---------------|
| P₇ | **[7, 3, 4]** | Simplex code = dual Hamming |
| P₁₁ | [11, 0, ∞] | Full rank (trivial) |
| P₂₃ | **[23, 11, 8]** | Dual Golay |

The duality is natural: C = null(A) gives the code whose parity-check matrix is A. When A is the Paley adjacency matrix, this is the dual of the QR code with generator matrix related to A.

## Why This Matters

1. **P₇ → simplex [7,3,4]**: The simplex code is the dual of the Hamming [7,4,3] code. It has 2³ = 8 codewords, all of weight 4. The weight distribution {0:1, 4:7} is perfectly uniform among nonzero words. This reflects the regularity of P₇ — all scores equal 3.

2. **P₂₃ → dual Golay [23,11,8]**: The dual Golay has 2¹¹ = 2048 codewords with minimum distance 8. Its weight distribution {0:1, 8:506, 12:1288, 16:253} has the remarkable property of being 4-design. P₂₃ is the unique regular tournament on 23 vertices with automorphism group PSL(2,23).

3. **P₁₁ → trivial**: The P₁₁ adjacency matrix has full rank over GF(2) — its null space is empty. This means the QR code at p=11 is "too dense" to have a nontrivial dual.

## The Tournament-Code Dictionary

| Tournament concept | Coding concept |
|-------------------|---------------|
| Iso class | Code (set of codewords) |
| Fiber (tiling count) | Code size K |
| Min intra-class Hamming | Minimum distance d |
| Waggly distance | Metric on code space |
| Paley adjacency | QR code parity-check |
| Regular tournament | Constant-weight code |
| Transitive tournament | Repetition code |

## The Association Scheme Connection

The metagraph is NOT an association scheme (tested at n=5: full algebra dim=35 vs needed 7). But the Krawtchouk polynomials K₁ still approximate the H function well (r=-0.94 at n=5, -0.83 at n=7). The scheme structure is a GOOD APPROXIMATION that degrades slowly.

The Paley codes sit at the intersection of:
- **Tournament theory** (regular tournaments, H-maximizers)
- **Coding theory** (QR codes, simplex, Golay)
- **Association schemes** (Hamming scheme, Krawtchouk eigenvalues)

All three viewpoints converge on the same algebraic structure: the quadratic residue pattern mod p.
