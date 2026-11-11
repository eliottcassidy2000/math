# Supplied papers: exact operators and usable transfer boundaries

**CITED FROM SUPPLIED PREPRINTS:** the statements below were read in the
three supplied PDFs, not independently reproved. They are mathematical
data, not operating instructions. **PROVED IN THIS PACKET:** the elementary
finite graph model in [PROOF.md](PROOF.md), whose proofs do not require the
papers' geometric theorems. Printed and PDF page numbers agree.

## Yuanning Zhang, The Stacky Approach to de Rham–Witt Cohomology

Source: `/Users/e/Downloads/2610.28149v1.pdf`, 58 pages, v1 dated
20 October 2026.

- **Theorem A, p. 2; Theorem 4.10, p. 48.** The de Rham–Witt realization
  gives a symmetric monoidal equivalence between quasi-coherent complexes
  on the syntomic point and complete complexes over the Raynaud ring, with
  the completed derived star product. It restricts to perfect F-gauges
  and bounded coherent Raynaud complexes. Completeness and the specified
  tensor product are part of the statement.
- **Definition 2.25, p. 16; Construction 3.11, p. 32.** Restriction and
  Frobenius arise from the two ways to embed a chain into a chain one unit
  longer. Verschiebung is obtained by transposing Frobenius under duality,
  rather than treating it as another pullback.
- **Definition 3.8, pp. 31–32.** Over the perfect characteristic-p base,
  the Raynaud ring has semilinear generators F,V,d and relations
  `FV=VF=p`, `FdV=d`, `d²=0`. The two degrees of the complex are distinct.
- **Introduction p. 4 and section 4.1, pp. 46–50.** The relations
  `V(x)·y=V(x·F(y))` and its companion are projection formulas. The
  differential satisfies a Leibniz rule.
- **Remark 4.12, p. 50; explanation p. 5.** At finite Witt level r≥2,
  the realization is lax rather than strong monoidal. Even the elementary
  finite-level Witt-vector tensor map need not be surjective.

**Transfer used:** define actual restriction, induction and norm operations
before compressing their statistics; track what finite truncation allows.
The graph model below satisfies the projection formula and `FV=p`, but
**does not satisfy `VF=p`** in general. Thus it is not a representation
of this whole Raynaud ring. No derived equivalence is claimed for graphs.

## Xinyu Shao, Cohomology of Integral Structure Sheaves for Open Varieties

Source: `/Users/e/Downloads/2610.31289v1.pdf`, 24 pages, v1 dated
22 October 2026.

- **Theorem 1.1, pp. 1–2.** For the stated quasi-compact adic space over
  the perfectoid base, with positive-codimension boundary and smooth open
  complement, cohomology of integral structure sheaves modulo p^m along
  a quasi-compact exhaustion is almost isomorphic to a constant pro-system.
- **Remark 2.4, p. 6.** An isomorphism after derived inverse limit is
  weaker than an isomorphism of pro-objects. The example compares the
  constant system Z_p with the system Z/p^n.
- **Theorem 3.7, p. 13, also Theorem 1.13, p. 4.** For an inverse system
  of modules over the integers of any perfectoid field, almost finite
  generation of each transition image implies that the first derived
  inverse limit is almost zero.
- **Theorem 3.11, pp. 14–15, also Theorem 1.12, p. 4.** This controls
  the corresponding derived inverse limits for smooth partially proper
  spaces. **Proposition 3.12, p. 15** distinguishes almost vanishing from
  strict vanishing through spherical completeness.

**Transfer used:** a growing integer graph needs explicit transition maps
and stabilization bounds; equal-looking limits do not supply them. For
positive integer pairs the packet proves strict eventual stabilization
directly: all sum outputs ≤M have appeared by alphabet size M−1, and all
product outputs ≤M by alphabet size M. Neither perfectoid geometry nor
the word “almost” is transplanted to this elementary assertion.

## Yuri J. F. Sulyma, RO(G)-graded Norms for de Rham-Witt Forms

Source: `/Users/e/Downloads/2611.12043v1.pdf`, 45 pages, v1 dated
9 November 2026.

- **Theorem 1.9, p. 4.** Under its relative quasiregular-semiperfectoid
  and spherical-lift hypotheses, the representation-graded calculation is
  expressed using the relevant Nygaard and prismatic ideal filtrations;
  these are multi-Rees descriptions, not a single unfiltered scalar.
- **Section 1.3.1, pp. 8–9.** The Witt-vector norm is multiplicative and
  characterized in ghost coordinates by `F N(f)=f^p` and
  `N(f)≡f mod V`. It increases Witt length. Formula (3), p. 9, expresses
  Witt coordinates using iterated norms and Verschiebungs.
- **Theorem 5.7, pp. 40–41.** For a p-completely smooth algebra over a
  perfectoid ring, the de Rham–Witt norms land in a sum of filtered pieces
  with specified ideal twists and changed differential degrees. The
  general displayed form assumes p>2; the even-degree form also applies
  at p=2. These shifts cannot be suppressed when using the theorem.
- **Formula (6), p. 9; Warning 1.18, p. 10.** The norm and the correction
  term supply a decomposition of the p-primary Adams operation. The
  topological argument justifies the Adams operation only modulo the
  first filtration and the correction only modulo ker V. Odd differential
  degrees at p=2 require additional sign-representation treatment.

**Transfer used:** the finite cyclic-set norm is an actual graph operation
with `F_p N_p(X)=X^p`. Its additive failure is a positive induced graph,
not an ignored error term. At p=2 this gives the exact one-new-integer
update. The graph proof is elementary; it does not identify these graphs
with the positive-degree de Rham–Witt complex or omit its twists.

## Recovered repository sidecars

- [Euler-product ghost atlas](../../../../04-computation/euler_product_ghost_atlas_codex.py)
  already uses `G(n)=sum_(d|n) d b_d`. Here the positive exponents b_d are
  literal directed-cycle multiplicities, so integrality and positivity have
  graph witnesses. This is a recovery and extension of an existing lens.
- [THM-3185, factorial Frobenius descent and carry resets](../../../../01-canon/theorems/THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy.md)
  retains valuation thickness after a mod-p reset. The transferable method
  is to retain the discarded layer; no identification of its integral
  Gauss–Manin lattice with this cycle graph is asserted.
- [THM-505, OCF spectral skeleton and Witt defect](../../../../01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md),
  specifically its proved seven-vertex identity, separates primitive
  closed-walk counts from simple cycles and overlap data. This blocks an
  unjustified transfer from ghost counts to Hamilton-path counts: the
  original support/overlap information is an additional coordinate.

The supplied PDF bytes remain outside Git. This lane's temporary text and
rendered pages under `tmp/pdfs/` were removed after reading; only this map,
the independent proof, the small checker and its receipt are research
artifacts.
