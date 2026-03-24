# Eigenvalues of the Merged Metagraph

**Session**: opus-2026-03-23-S268

---

## The Principal Finding: H IS the Second Eigenvector

The Hamiltonian path count H(T) is approximately the **second eigenvector** of the adjacency matrix of G_n/Z_2:

| n | |⟨v₁, H⟩|² | Correlation r |
|---|-----------|---------------|
| 3 | 100.0% | 1.000 |
| 4 | 99.5% | -0.998 |
| 5 | 72.3% | -0.850 |
| 6 | 78.6% | 0.888 |
| 7 | 73.3% | -0.856 |
| 8 | **87.8%** | 0.940 |

H captures **73-88% of the variance** in the second eigenvector of the metagraph's adjacency matrix. This is remarkable: the metagraph's geometry (encoded in its adjacency) "knows" about Hamiltonian path counts without any explicit computation.

**Why?** The second eigenvector of a graph captures the largest-scale variation after removing the degree effect (first eigenvector). In the metagraph, the main axis of variation IS the H-gradient: tournaments range from H=1 (transitive) to H=max (regular/Paley). The metagraph's edge structure reflects this gradient — nearby nodes in H tend to be connected.

The sign alternation (r flips between +/−) reflects the eigenvector orientation, which is arbitrary.

---

## The Markov Spectral Gap

The spectral gap of the random walk on G_n/Z_2:

| n | Gap = 1−μ₁ | 2/n | Ratio |
|---|-----------|-----|-------|
| 3 | 2.000 | 0.667 | 3.000 |
| 4 | 1.500 | 0.500 | 3.000 |
| 5 | 0.480 | 0.400 | 1.201 |
| 6 | 0.331 | 0.333 | **0.994** |
| 7 | 0.206 | 0.286 | 0.720 |
| 8 | 0.149 | 0.250 | 0.595 |

At n=6, the gap is almost exactly 2/n (ratio 0.994). But for n ≥ 7, the gap falls below 2/n and the ratio decreases. The metagraph becomes a **worse expander** at larger n.

**Conjecture (revised)**: The Markov gap ~ c/n for some c < 2, perhaps c → 1 as n → ∞.

**Mixing time**: τ_mix ~ 1/gap ~ n/c. Random walks on the metagraph mix in O(n) steps.

---

## Spectral Density: Approaching Semicircle

The eigenvalue distribution at n=8 (3528 eigenvalues) shows a clear **semicircle-like** shape, concentrated near 0 with tails extending to ±max. This is characteristic of:

1. **Random graphs** (Wigner semicircle law)
2. **Expander graphs** (eigenvalues near 0 imply good expansion)
3. **The GUE ensemble** (random Hermitian matrices)

The fact that the metagraph's spectrum looks semicircular suggests that at large n, the metagraph behaves like a **random graph on V_merged vertices with average degree d̄ ≈ C(n,2)**.

---

## Not Ramanujan (n ≥ 6)

A d-regular graph is Ramanujan if |λ₁| ≤ 2√(d−1). The metagraph is not regular, but we can check the analogous bound:

| n | |λ₁| | 2√(d̄−1) | Ramanujan? |
|---|------|----------|-----------|
| 5 | 1.85 | 3.58 | ✓ |
| 6 | 5.49 | 5.44 | ✗ (barely) |
| 7 | 12.75 | 7.64 | ✗ |
| 8 | 22.08 | 9.96 | ✗ |

The metagraph is NOT a good Ramanujan expander for n ≥ 6. The second eigenvalue grows faster than the Ramanujan bound. This is because H creates a strong gradient that prevents uniform expansion — the transitive tournament (H=1) is isolated from the regular tournaments (H=max).

---

## Triangle Count

| n | Triangles in G_n/Z_2 |
|---|---------------------|
| 3 | 0 |
| 4 | 1 |
| 5 | 12 |
| 6 | 139 |
| 7 | 1159 |
| 8 | 14184 |

The triangle count from tr(A³)/6 grows rapidly. Ratios: −, 12, 11.6, 8.3, 12.2. This measures the local clustering: how often two neighbors of a node are also neighbors of each other.

---

## Lie-Theoretic Interpretation

From the A_{n-1} Lie algebra identification (S265):
- Tournament arcs = positive roots
- S_n = Weyl group
- Metagraph = quotient of Coxeter complex

The eigenvalues of the metagraph should relate to the **representation theory** of S_n. Specifically:

1. The **trivial representation** (constant eigenvector) gives λ₀ (the spectral radius)
2. The **standard representation** of S_n on R^{n-1} (= Cartan subalgebra) should give eigenvalues related to H
3. Higher representations give the remaining eigenvalues

The fact that H ≈ v₁ suggests that the **standard representation dominates** the spectral structure. H is essentially the projection of the tournament onto the Cartan subalgebra (score sequence), and the 2nd eigenvector captures this.

---

## Key Eigenvalue Ratios

| n | λ₀ | λ₁ | λ₁/λ₀ | λ_{min} | λ_{min}/λ₀ |
|---|-----|-----|--------|---------|------------|
| 3 | 1.0 | −1.0 | −1.000 | −1.0 | −1.000 |
| 4 | 2.0 | −1.0 | −0.500 | −1.0 | −0.500 |
| 5 | 4.64 | 1.85 | 0.399 | −2.90 | −0.625 |
| 6 | 9.78 | 5.49 | 0.561 | −5.35 | −0.547 |
| 7 | 17.62 | 12.75 | 0.724 | −11.93 | −0.677 |
| 8 | 26.71 | 22.08 | 0.827 | −19.03 | −0.713 |

λ₁/λ₀ **increases** toward 1: the gap between the two largest eigenvalues narrows relative to their magnitude. This confirms the metagraph becomes a worse expander.

|λ_{min}|/λ₀ stays around 0.5−0.7: the most negative eigenvalue is about 50−70% of the spectral radius. This indicates a **bipartite-like tendency** (bipartite graphs have λ_{min} = −λ₀).

---

## Open Questions

1. Does H/v₁ alignment continue improving at n=9,10,...?
2. What is the exact asymptotic Markov gap? Is it c/n with c < 2?
3. Can the triangle count be computed from the Burnside tower?
4. Does the spectral density converge to a semicircle (Wigner) or some other universal distribution?
5. What representation of S_n generates the 2nd eigenvector?
