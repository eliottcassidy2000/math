# Recovered mechanisms and map ledger

**Status:** PROVED toy bridges / CONDITIONAL Yang–Mills targets. The repo files below are read for their precise scoped mechanisms, not imported as physical results.

## Inheritance and concept board

- Anchor: low-energy spectral measure ↔ heat transform ↔ uniform visibility; obstruction is soft low-energy weight, not vacuum degeneracy.
- Niche: two-dimensional Wilson parity ↔ two-channel transfer; invariant is the full flux channel and composition law.
- Wildcard: tournament switching ↔ flat Z₂ lattice gauge fields; obstruction is a cohomology/toron sidecar lost on quotient.
- Cross-check: Pfaffian/Wick notation ↔ Gaussian characteristic function; obstruction is unrecorded higher cumulants and path ordering.
- Particle lane: energy marginal ↔ invariant mass measure; obstruction is loss of spatial momentum and mass-shell isolation.

Closest proved inherited mechanism: THM-2842's explicitly measured variance-jet decoder. Canonical hostile: S2's shrinking overlap makes fixed-word estimates uniform while normalized excitation gaps vanish. Corrected near miss: scalar decay and scalar commutant do not supply quantitative visibility. Least-used sidecar: torus flat holonomy, invisible to all local plaquettes.

## 1. Local-to-global observability

Read [THM-2842: ordered positive-cone Vandermonde multiplier observability](../../../../01-canon/theorems/THM-2842-ordered-positive-cone-vandermonde-multiplier-observability.md). It distinguishes external radial-variance jets from information already in an ordinary Gaussian scalar trace, and states a sharp finite-dimensional decoder count. That distinction, not a claimed transfer of the theorem to fields, drives S3 of PROOFS.md.

- Source: a finite quotient with scalar readouts and independently accessible multipliers.
- Target: positive-energy PVM with word-vector diagonal measures.
- Map: replace decoder coordinates by probe overlaps `||P_(0,u] x_j||²`.
- Preserved predicate: a lower visibility bound rules out an invisible nonzero target.
- Lost information: uniform conditioning, energy dependence, normalized vectors.
- Sidecar: explicit K in `||P|| <= K sum ||P x_j||²`.
- Cheapest decisive test: S2 gives K>=ε^−2. A proposed phase-uniform proof must prevent this blowup.

The positive bridge is an independently proved lemma; THM-2842 is provenance for the mechanism, not a dependency for its functional-analytic proof.

## 2. Transfer kernels compose; one scalar need not

Read [THM-1975: path-cover polynomial/refined compositional invariant](../../../../01-canon/theorems/THM-1975-the-path-cover-polynomial-is-the-refined-compositional-invariant.md). Its status is VERIFIED in the finite universe recorded there; its exact cyclic transfer kernel is explicitly left open. Do not promote its discussion of a general kernel to a new inherited theorem.

- Source: retain the path-component profile instead of only its Hamiltonian slice.
- Target: retain even/odd flux probabilities instead of only one selected transition entry.
- Map: both form a boundary-state vector; gluing is contraction by a fixed kernel.
- Preserved predicate: exact composition once channels are retained.
- Destroyed information: selecting a single entry forgets intermediate boundary states.
- Sidecar: both parity channels, their boundary condition, and torus closure.
- Cheapest test: the exact 2×2 matrix T in PROOFS.md; check powers against direct 2×2 lattice enumeration.

This target is completely elementary and proved here. No theorem about tournament complexity is needed to obtain it.

## 3. Switching, Wilson loops, and cohomology

Read [THM-479: level permutations via pair-orbit holonomy](../../../../01-canon/theorems/THM-479-level-holonomy-and-branch-split.md). Its affine tournament action and cut-space quotient retain orbit holonomy precisely because local orientation coordinates alone forget it.

- Source: F₂ edge cochains modulo vertex cuts; holonomy evaluates on cycles.
- Target: link signs modulo vertex gauge flips; Wilson signs evaluate on lattice cycles.
- Map: identify −1 with 1∈F₂; product becomes mod-two sum, gauge flips become coboundaries.
- Preserved predicate: cycle pairing is invariant under cuts.
- Lost information under plaquette-only projection: classes in H¹ of the surface.
- Sidecar: a cycle basis/torons, plus the product-of-all-plaquettes constraint on a closed torus.
- Cheapest test: same plaquette pattern with opposite fundamental-cycle holonomy; exact loop average zero when both sectors are sampled.

The tournament affine permutation cocycle is **not** itself a Yang–Mills connection or proof of confinement. The honest common structure is cochain/cut/cycle algebra. Center-vortex terminology would require specifying an embedding of Z₂ into the chosen gauge-group center, admissible lattice fluxes and their probabilities; no nonabelian claim is made here.

## 4. Pfaffian and Gaussian: keep the law, not just the determinant

Read [THM-174: Pfaffian determinant identity](../../../../01-canon/theorems/THM-174-pfaffian-determinant-identity.md). The even skew determinant square illustrates how a signed antisymmetric object can be hidden by its square.

- Source: skew kernel S and Pf(S), with det(S)=Pf(S)².
- Target: Gaussian covariance and characteristic functional, or fermionic Gaussian integration only when actually specified.
- Map: pair contractions are determined by a two-point kernel **under Gaussianity**.
- Preserved predicate: exact Gaussian expectations once the law and orientation/sign conventions are specified.
- Lost information: Pfaffian sign under squaring; higher cumulants under covariance truncation; ordering for noncommuting holonomies.
- Sidecar: field law, cumulants, representation, ordering and normalization.
- Cheapest test: Gaussian versus Rademacher of equal variance produces different Wilson characteristic functions.

No general Yang–Mills Gaussianity follows from a Pfaffian identity. A successful finite skew-algebra analogy must carry its measure hypothesis along.

## Scope and reproduction

Run `python3 05-knowledge/results/session272_yang_mills_bridges_20261111/spectral_loops/check.py` from the repo root. `check.out` records exact enumeration of all 4096 disk and 256 torus link configurations, all 16 face subsets in each case, 4088 fixed-word/phase controls, and parity transfer powers through 20. The continuous spectral values printed last are floating diagnostics; their formulas are proved in PROOFS.md. No reconstruction claims, source theorems, continuum theories or open Yang–Mills conclusions are computationally reverified here.
