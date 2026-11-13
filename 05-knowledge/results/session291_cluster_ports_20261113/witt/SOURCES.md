# Source map and scope

## Attached source actually consulted

Yuri J. F. Sulyma, *RO(G)-graded norms for de Rham–Witt forms*,
arXiv:2611.12043v1, 45 pages. User-supplied attachment. Relevant passages:

- Introduction and Theorem1.9, pp3–5: representation grading retains filtration
  and prism twists absent from the integer-graded diagonal. Hypotheses involve
  delta-rings, spherical lifts and relatively quasiregular semiperfectoid rings.
- Warning1.7 p4: sum over conjugates is trace; norm is product over conjugates.
- Theorem1.15 pp7–8: the TR formula has a restricted representation range, not
  an unrestricted classification. The paper explicitly discusses extension
  difficulties outside that range.
- Sections1.3.1–1.3.2 pp8–9: Witt norm satisfies FN(f)=f^p and N=id modulo V;
  Witt coordinates use sums of V^k N^(n-k) terms. For forms, N changes degree
  and twists. The equation psi^p(omega)=N(omega)+V delta(omega) requires its
  proper mixed-degree carrier, not a naive scalar formula.
- Theorem5.7 pp40–42: multiplicativity N(omega tau)=N(omega)N(tau) coexists
  with N(omega+tau)=N(omega)+N(tau)+V s1(omega,tau), where
  (a+b)^p=a^p+b^p+p s1(a,b). Twisted target modules matter.

We use these as cited inspiration, without re-auditing the paper. None of its
homotopical results is asserted for tournaments. Our C_p is an independently
proved classical necklace carry, not the paper's N or delta.

## Repository lineage recovered

- `04-computation/euler_product_ghost_atlas_codex.py` and
  `07-reflections/euler-product-ghosts-and-zeta-side-cancellation-gates.md`:
  exponents, ghost coordinates and ordinary coefficients are different carriers;
  sparsity need not survive their transforms. Already an explicit old lesson.
- `04-computation/witt_necklace_s90at.py`: prior classical necklace/Mobius
  exploration. The necklace identity here is not claimed as new.
- `01-canon/theorems/THM-1926-tournament-zeta-euler-product-over-strong-core.md`:
  ordinary spectral ghosts enumerate primitive closed walks, factor over strong
  components, and ignore acyclic components.
- `01-canon/theorems/THM-505-the-ocf-non-spectral-defect-H-equals-spectral-skeleton-plus-witt-defect.md`:
  overlap defects obstruct replacing simple cycles by primitive closed walks;
  inherited OCF uses vertex-disjoint odd-cycle collections. Square-zero grading
  in this packet is an explicit classical way of enforcing that constraint.
- `THM-3185-iterated-factorial-frobenius-descent-and-witt-carry-reset-hierarchy.md`:
  retained state coordinates matter in the factorial recurrence. This is a
  different carrier, useful analogy only.
- `THM-3178-squarefree-resultant-tangent-cone-and-first-witt-norm.md`:
  its first conormal coefficient is a determinant norm in a collision algebra.
  It does not identify resultants, necklaces and de Rham–Witt forms. The shared
  discipline is specifying the algebra and readout before using the word norm.

No claim of comprehensive literature priority is made. The elementary identities
are independently derivable; the packet's contribution is the precise comparison,
physical-realization obstruction and explicit carrier-preserving test.
