# 2026-10-22 intake: an identifier supplied for supersingular curves resolves to the KLS proof

**Checked 2026-10-22 by opus-2026-10-22-ssj.** This file is kept outside `CORE-PAPERS.md`, whose startup byte budget is full. It follows the precedent of the "Zhao's GVC in three variables" entry there, an identifier that resolved to an unrelated paper.

## Bizeul–Klartag–Lehec — *Presenting a proof of the Kannan–Lovász–Simonovits conjecture*

**Status: PREPRINT v1 CLAIM, CITED, unaudited.**
* Source: [arXiv:2610.05474v1](https://arxiv.org/abs/2610.05474) [math.FA], 4 Oct 2026.
* **Main claim, Theorem 1.1.** There is a universal C such that every isotropic log-concave probability measure μ on R^n satisfies C_P(μ) ≤ C. Equivalently, ‖Cov μ‖_op ≤ C_P(μ) ≤ C ‖Cov μ‖_op.
* **Route.**
  * the third-moment bound |E⟨X,u⟩X⊗X|^2 ≤ 8|u|^2 (Chen–Klartag, Bizeul, Letwin);
  * cumulant bounds |κ_m(u,·,…,·)|^2 ≤ K^{m−1}((m−1)!)^2|u|^2, proved by Eldan's stochastic localization (Theorem 4.1);
  * suspension of tilt-averages (Proposition 7.1);
  * the Song–Zhang criterion (arXiv:2610.01447, Theorem 5.1; reproved as Theorem 3.1).
* **Attribution.** The acknowledgements credit ChatGPT with most proofs and ideas. The suspension idea is credited to the authors.
* **Does not prove** anything for a repository lane. There is no consumer beyond its identification.

## The supplied content: provenance UNRESOLVED

The identifier was relayed with this description: "supersingular elliptic curves: j-invariants classify curves up to isomorphism, finitely many supersingular classes for each p, j-invariants are roots of a specific polynomial lying in F_p or F_{p^2}".
* The paper contains none of this. The full HTML text has no occurrence of "supersingular", "elliptic curve" or "j-invariant".
* The material is standard: Deuring, Delfs–Galbraith, Ogg. It is verified in [supersingular_divisor_shapes_qseries_20261022 §2](../results/supersingular_divisor_shapes_qseries_20261022.md).
* Attach no priority claim to the paper for it. The intended source is unknown and should be asked of the owner.
