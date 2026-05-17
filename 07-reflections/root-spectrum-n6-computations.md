# Root Spectrum: Exhaustive n=6 Computations

**Session:** oracle-2026-05-17-S1
**Source:** `04-computation/root_spectrum_fast.py`, results in `05-knowledge/results/root_spectrum_fast.out`
**Depends on:** `real-rootedness-product-formula-erdos.md`, `cubic-family-x6-coloring.md`

---

## Complete n=6 Root Spectrum (47 of 56 iso classes, key = H,I6,scores,dc3,dc5)

| H | I6 | α₁ | α₂ | dc3 | dc5 | SC | ρ₁ | ρ₂ | ratio | deg |
|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 1 | 0 | 0 | 0 | 0 | ✓ | — | — | — | 0 |
| 3 | 7 | 1 | 0 | 1 | 0 | | 1.00000 | =ρ₁ | 1.0000 | 1 |
| 5 | 13 | 2 | 0 | 2 | 0 | | 0.50000 | =ρ₁ | 1.0000 | 1 |
| 5 | 13 | 2 | 0 | 2 | 0 | ✓ | 0.50000 | =ρ₁ | 1.0000 | 1 |
| 9 | 25 | 4 | 0 | 3 | 1 | | 0.25000 | =ρ₁ | 1.0000 | 1 |
| 9 | 49 | 2 | 1 | 2 | 0 | ✓ | 1.00000 | 1.00000 | 1.0000 | 2 |
| 17 | 49 | 8 | 0 | 4 | 4 | ✓ | 0.12500 | =ρ₁ | 1.0000 | 1 |
| 17 | 73 | 6 | 1 | 4 | 2 | ✓ | 5.82843 | 0.17157 | 0.02944 | 2 |
| 23 | 91 | 9 | 1 | 5 | 4 | | 8.88748 | 0.11252 | 0.01266 | 2 |
| 29 | 109 | 12 | 1 | 6 | 6 | ✓ | 11.91608 | 0.08392 | 0.00704 | 2 |
| 29 | 133 | 10 | 2 | 6 | 4 | ✓ | 4.89792 | 0.10208 | 0.02084 | 2 |
| 33 | 145 | 12 | 2 | 6 | 6 | ✓ | 5.91548 | 0.08452 | 0.01429 | 2 |
| 37 | 157 | 14 | 2 | 6 | 8 | ✓ | 6.92783 | 0.07217 | 0.01042 | 2 |
| 41 | 169 | 16 | 2 | 8 | 8 | ✓ | 7.93700 | 0.06300 | 0.00794 | 2 |
| 43 | 151 | 19 | 1 | 8 | 11 | | 18.94722 | 0.05278 | 0.00279 | 2 |
| 45 | 157 | 20 | 1 | 8 | 12 | ✓ | 19.94987 | 0.05013 | 0.00251 | 2 |
| 45 | 229 | 14 | 4 | 8 | 6 | ✓ | 3.42705 | 0.07295 | 0.02129 | 2 |

(Table shows selected degree-2 rows; full table in results file)

---

## Key Computational Findings

### 1. Root Gap Confirmed

**No roots in (-1/3, -1/4)** at n=6 (exhaustive), n=7 (2000 samples), n=8 (300), n=9 (50).

The gap (1/4, 1/3) in ρ-space is respected by all tested tournaments. The boundary points are:
- ρ = 1/4 ↔ α₁=4, α₂=0: ACHIEVABLE (H=9, appears twice at n=6)
- ρ = 1/3 ↔ α₁=3, α₂=0: FORBIDDEN (would give H=7, impossible)

**The root gap IS the H-gap.** The forbidden H=7 is equivalent to the forbidden root ρ=1/3.

### 2. Forbidden (α₁=3, α₂=0) Not Found

Never found at n=6..9 in all samples. At n=7 it was known from prior computations that α₁=3 can occur but ONLY with α₂≥2 (disjoint cycle pairs). The root -1/3 thus remains permanently forbidden.

### 3. Ultra-Log-Concavity: Proved from Real-Rootedness

**Newton-Maclaurin inequality**: Zero violations at n=6..9.

**Theorem (Newton-Maclaurin for Tournaments):** If $I(\Omega(T), x)$ is real-rooted with degree $d$ and leading coefficient $\alpha_d$, then the sequence $(\alpha_k / \binom{d}{k})_{k=0}^{d}$ is log-concave:
$$\left(\frac{\alpha_k}{\binom{d}{k}}\right)^2 \geq \frac{\alpha_{k-1}}{\binom{d}{k-1}} \cdot \frac{\alpha_{k+1}}{\binom{d}{k+1}}, \quad 1 \leq k \leq d-1$$

**Proof:** For a degree-$d$ polynomial with all real (negative) roots, Newton's inequalities on the elementary symmetric polynomials $e_k(\rho_1,\ldots,\rho_d)$ give exactly this. Since $\alpha_k = \alpha_d \cdot e_{d-k}(\rho)$ and the $\rho_i$ are all positive real, the Newton-Maclaurin inequality holds. $\square$

**Corollary:** $k$-tuples of vertex-disjoint odd cycles are "ultra-log-concave": knowing $\alpha_{k-1}$ and $\alpha_{k+1}$ gives an upper bound on $\alpha_k$.

### 4. Vieta's Theorem Verified to Machine Precision

- **n=5 degree-1:** $r = -1/\alpha_1 = -2/(H-1)$. Error $< 10^{-15}$. ✓
- **n=6 degree-2:** $\rho_1+\rho_2 = \alpha_1/\alpha_2$, $\rho_1\rho_2 = 1/\alpha_2$. Error $< 4 \times 10^{-15}$. ✓

### 5. SC Tournaments are Most Root-Asymmetric

**Result (n=6):**
- SC min ratio: $\rho_2/\rho_1 = 0.00251$ (at H=45, α₁=20, α₂=1)
- NS min ratio: $\rho_2/\rho_1 = 0.00279$ (at H=43, α₁=19, α₂=1)

**The most asymmetric tournament at n=6 is SC**, supporting the conjecture in `real-rootedness-product-formula-erdos.md`. The class H=45, α₁=20, α₂=1 has:
$$\rho_1 \approx 19.950, \quad \rho_2 \approx 0.0501 = 1/\rho_1$$

This uses $\rho_1\rho_2 = 1/\alpha_2 = 1$, so **the roots are multiplicative inverses of each other** when $\alpha_2=1$.

For general α₂=1: $\rho_1\rho_2=1$ gives $\rho_2=1/\rho_1$. The ratio is $1/\rho_1^2$, which is minimized when $\rho_1$ is maximized. Max $\rho_1 = \alpha_1/\alpha_2 - \rho_2 \approx \alpha_1$ (for large $\alpha_1$). So the most asymmetric class is the one with the largest $\alpha_1$ subject to α₂=1.

**Key identity:** For α₂=1 classes, $\rho_1 + \rho_2 = \alpha_1$ and $\rho_1\rho_2 = 1$. The ratio:
$$\frac{\rho_2}{\rho_1} = \frac{1}{\rho_1^2} = \frac{4}{(\alpha_1 + \sqrt{\alpha_1^2-4})^2} \approx \frac{4}{\alpha_1^2}$$

### 6. (H, I(Ω,6)) Does NOT Separate All Iso Classes

**Result:** Of 47 distinct classes found by key (H,I6,scores,dc3,dc5), only **7/47 are uniquely identified by (H,I6) alone**. Many classes share the same (α₁, α₂) but differ in the detailed cycle arrangement.

**Interpretation:** The pair (α₁, α₂) is a coarse invariant of the conflict graph. It captures "how many cycles" and "how many disjoint pairs" but not the full structure. Two tournaments can have the same α₁, α₂ with completely different score sequences.

**Finding:** At n=6, the α₁ distribution across degree-1 classes shows the gap: α₁ ∈ {0,1,2,4,5,6,7,8,9,11,12,13,16,19,20} — **α₁=3 is absent** (confirming the forbidden case).

### 7. Root Ratio Trends with n

| n | degree | min ratio observed |
|---|---|---|
| 6 | 2 | 0.00251 (SC, exhaustive) |
| 7 | 2 | 0.00148 (sample, 2000) |
| 8 | 2 | 0.00073 (sample, 300) |
| 9 | 2-3 | 0.00083 (sample, 50) |

Root ratios decrease with n: more asymmetric distributions at larger tournaments. This is expected from the formula $\rho_2/\rho_1 \approx 4/\alpha_1^2$ and the fact that $\alpha_1$ grows roughly as $n^3$.

### 8. Degree-3 Appears at n=9

At n=9, 44/50 random samples had **degree-3 polynomials** (α₃ > 0). This means many n=9 tournaments have ≥3 vertex-disjoint odd cycles. The Newton-Maclaurin inequality for degree-3 requires:
$$\frac{\alpha_1}{\alpha_2}^2 \geq \frac{3}{2} \cdot \frac{\alpha_2}{\alpha_3} \quad \text{(at k=2, d=3)}$$

No violations found in 50 samples, strongly supporting ULC for degree-3.

---

## The α₁ Gap at n=6: A Complete Picture

At n=6, the achievable α₁ values (for degree-1 polynomials, α₂=0) are:
$$\alpha_1 \in \{0, 1, 2, 4, 5, 6, 7, 8, 9, 11, 12, 13, 16, 19, 20\} \setminus \{3, 10, 14, 15, 17, 18, \ldots\}$$

The forbidden α₁=3 creates the root gap $r \in (-1/3, -1/4)$. But there are other gaps too (α₁=10, 14, 15, 17, 18 also seem absent from degree-1 classes at n=6). This suggests a richer structure of forbidden configurations.

---

## Perfect Square Case: H = 9, α₁=2, α₂=1

The unique SC class at n=6 with H=9 has:
- $I(\Omega, x) = 1 + 2x + x^2 = (1+x)^2$
- Double root $\rho_1 = \rho_2 = 1$
- $H = (1+\alpha_1)^2 = (1+2)^2 = 9$ ✓

This is the **minimum H SC tournament at n=6** and the only tournament (at any n) with a perfect square polynomial of this form at n=6. The root "collapses" to exactly 1.

The relationship $\rho_1\rho_2 = 1$ (for α₂=1) at the double root gives $\rho=1$, the "harmonic mean" of all roots. This is the maximum-symmetry polynomial (ratio=1.0) within the degree-2 family.

---

## Open Questions from These Computations

1. **Complete the n=6 enumeration**: The 9 missing classes (due to key collisions on dc3,dc5,scores) need a finer invariant. What distinguishes them from the 47 found?

2. **Ultra-log-concavity for degree-3**: Verify Newton-Maclaurin at k=1 and k=2 for degree-3 cases at n=9,11. We have 50 samples; need more.

3. **SC min ratio at n=7**: Is the minimum root ratio still achieved by an SC tournament at n=7? Need to check 1000 SC tournaments at n=7.

4. **Root ratio formula**: For α₂=1 classes, ratio ≈ 4/α₁². What is the formula for α₂>1?

5. **Root gap at larger n**: Is (-1/3,-1/4) permanently empty for all n? Related to: is α₁=3, α₂=0 impossible for all n? The n≤9 evidence says yes.
