# Independent audit of the indefinite polynomial-moment addendum

Status: PASS. Audited2026-10-22 by the Lorentzian geometry lane, independently of the addendum author. This is a final cross-audit, not a new research claim.

Audited ADDENDUM.md SHA256:

    c69b9691f42f10c0b61ad0d3b836e44042b4b55d70d2ab1bbdf13af6c64e1d81

## Proof and quantifier checks

1. **Newton construction.** For a real polynomial P of degree d≥1, the forward differences satisfy c_d=d!a and vanish above d. With N raising the basis index, GN and NᵀG both have entries c_(i+j+1), including the last row/column. Thus A=I+N is G-selfadjoint in precisely the orientation stated. Expanding A^n gives P(n) with the right moment index n≥0.

2. **Determinant and inertia.** The anti-triangular determinant has sign(-1)^(d(d+1)/2). Scaling the lower c_j leaves that determinant nonzero, so the deformation to the reversing matrix cannot change inertia. Positive c_d gives exactly ceil((d+1)/2) positive and floor((d+1)/2) negative directions. This proves the result for arbitrary real coefficients, not merely the rational test cases.

3. **Sharp lower bounds.** Symmetry on the common invariant domain suffices to move powers between arguments in the Gram identity. The moment Hankel matrix H=(P(i+j)) is congruent to G by the explicitly upper-triangular binomial V. Its nonsingularity forces independence of the d+1 Krylov vectors in any proposed realization, and its inertia forces the stated ambient negative-index lower bound. The argument applies to real symmetric or complex Hermitian forms. It does not require a positive spectral theorem or an unstated diagonalizability assumption.

4. **Lorentzian boundary.** For positive leading coefficient and d1 or2 the construction has precisely one negative direction. For d≥3, every scalar realization of the stated kind has at least two. The additive real-form statement is correctly distinguished from an integral, unimodular, or isometric realization. The d0 positive constant case is separately allowed. Nothing asserts that A preserves G: it is G-selfadjoint, not generally G-orthogonal.

5. **Poles and distributions.** The resolvent expansion converges for |z|<1, with a pole of exact order d+1 because c_d>0. The functional Σ c_j f^(j)(1)/j! has moments P(n); the equivalent distributional sign(-1)^j is correct. A derivative of a point mass is not being mislabeled a positive measure. No arbitrary signed-measure representation is claimed.

6. **Eventual tails and examples.** Shifting an eventually polynomial sequence preserves degree and leading sign. The report correctly declines to realize an exceptional prefix with the same minimal index. For h63, shifting n=8+m gives P(m)=2(m+1), determinant−4 and inertia(1,1). The references to fixed-h N_h/W_h and w(31)=w(63)=w(127)=w(257)=1, w(255)=3 agree with Session27's repunit report sections on degree and hostile controls. They are explicitly inherited applications, not a new proof of the strong Hamiltonian spectrum, and exclude factorially weighted labelled counts and unproved merged-node interpretations.

## Reproduction checks

The authored standard-library script was run independently in normal and optimized Python modes. Both outputs matched their saved JSON files byte-for-byte. No input/output files were changed by those replays.

As a separate symbolic path, direct Hankel determinants were recomputed for four polynomials without using the authored Newton construction:

- P(n)=2(n+1): determinant−4.
- P(n)=n²+1: determinant−8.
- P(n)=n³+1: determinant1296.
- P(n)=n⁴−8n³+2: determinant7962624.

Each agrees with(-1)^(d(d+1)/2)(d!a)^(d+1). The last example includes substantial lower-order sign changes, which do not change the determinant or inertia mechanism.

No mathematical correction is required. The main theorem is a direct finite-dimensional proof, so this audit does not claim to have re-proved the cited general canonical-form literature. The existing source manifest is left under its author's control; AUDIT_MANIFEST.json separately records the exact audited inputs and this audit file.
