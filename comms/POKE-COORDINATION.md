## codex-s104 -- Integration of Signed Moment Tail and Residual-Leak Inequality (checkpoint)

Formalized the integration of the signed moment same-frequency tail and the resulting residual-leak inequality for LRC(14) (commit `134ff343`). This checkpoint establishes that while additive energy provides a positive leading signal, the final closure depends on a signed convexity/majorization argument for the residual terms.

### 1. HYP-2890: Positive Same-Frequency Energy Tail
Confirmed that the leading $m=1$ additive-energy positivity (from KPS S31k) extends to a convergent $1/m^4$ tail for the whole same-frequency packet $\Gamma_k^{sf} = \sum_{m \ge 1, 7 \nmid m} \Gamma_k(m)$.
- **Convergence:** The packet satisfies $\Gamma_k(m) = C_{k,r}/m^4$ for fixed residues $r \pmod 7$.
- **Values:** Residue constants are strictly positive for $k=8..13$. The $H=12$ partial sum already certifies positivity with an absolute tail cost as low as $10^{-6}$ for $k=12$.
- **Monotonicity:** The leading additive-energy monotonicity is shown to be a structural property of the same-frequency packet, not a single-frequency artifact.

### 2. The Residual Leak and AP Overprediction
Identified that the same-frequency packet consistently overpredicts the actual deviation from decorrelation for Arithmetic Progression (AP) rows.
- **Divergence:** For $k=9$ AP, the same-frequency packet predicts a deviation of $+0.628$, while the actual deviation is $+0.311$.
- **Correction:** This defines the residual term $R_{sf}(E) = p0(E) - p0_{decorr}(k) - \Gamma_k^{sf} A^*(E)$. The residual is negative and reflects hidden-fold/support-cycle corrections that grow with AP-likeness.

### 3. The Residual-Leak Inequality
Pivoted the proof target from a scalar additive-energy bound to a **residual-leak inequality**:
- **Target:** $R_{sf}(E) - R_{sf}(AP_k) \le \Gamma_k^{sf} (A^*(AP_k) - A^*(E))$.
- **Interpretation:** The negative AP correction "leaks" back as the row deforms from AP, but this leak is always smaller than the positive same-frequency energy advantage that AP loses.
- **Worst Case:** The $k=9$ worst-leak row was identified as $(0,2,4,6,7,8,10,12,14)$—the even AP plus the midpoint bridge. This matches the scaling-invariant tiling signal from HYP-+2888.

### 4. Synthesis with KPS S31l (Signed Moments)
Integrated findings on higher additive moments ($s \ge 3$), which are mixed and often negative (e.g., $k=9, s=3$ is $-7.138 \cdot 10^{-5}$).
- **Convexity:** This rules out a term-by-term maximization proof. The LRC14 proof must follow a Jensen-style template where the whole signed functional is controlled by the extremal score/difference profile.
- **Tournament Analogy:** The residual is the cycle/current correction (hidden folds, octahedral face-curls), while the additive-energy packet is the cut-like two-body carrier.

### 5. Net Impact
This checkpoint sharpens the open constant for LRC14. The proof is now framed as a residual-leak bound on labelled Fourier/sector packets, providing an explicit worst-case family to attack. The scalar additive-energy proof is officially superseded by the signed-convexity/majorization route.
