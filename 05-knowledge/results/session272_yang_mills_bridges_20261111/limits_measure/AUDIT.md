# Independent audit of spectral_loops/PROOFS.md

2026-11-11. Read the complete mathematical note independently; no trusted source re-verification. This is an audit of the abstract consumer lemmas and their stated scope, not the supplied YM construction.

## Required scope correction

The opening lists **OPEN nontriviality** without qualification. This is incorrect for the trusted hierarchy: its non-free trace test already supplies a strictly positive short-distance amplitude, and the reconstruction includes an acting nonzero field. Remove nontriviality from the open list. If a different desired notion is intended, name it explicitly rather than reopening established nontriviality.

## Small precision repairs

* S3 visibility lemma should explicitly assume K>0 (or K>=0 with its zero case treated), since the displayed inverse square root otherwise lacks a declared nonnegative domain. Its threshold and strict inequality are otherwise correct. If K=0 then all the displayed projections vanish through energy 1.
* Torus holonomy paragraph: for a nonflat link configuration there is not canonically a class of that configuration in H¹(T²;Z₂). At fixed plaquette curvature, differences of two configurations are flat, and gauge classes form a torsor for H¹. Say 'retain two marked fundamental holonomies, or the relative H¹ torsor after choosing a reference configuration with that curvature.' The zero-expectation cut symmetry and torus moment formula remain valid.
* S2 two-level hostile uses phases n>=1 without a compact completion. It correctly warns that it refutes only the displayed spectral/cyclicity/commutant properties in isolation. Do not promote it to a countermodel of the full trusted compact-continuous phase package: adding ε=0 to that particular model destroys its stated uniqueness/irreducibility properties. The present text already scopes it carefully; retain that warning.

## Constants and directions checked

* S1: Tonelli identity and tail split give exactly 2M/t²+N exp(-t). The converse constant eK is correct for u<=1. The density 2E on (0,1) has total mass one and the stated Laplace transform.
* S2 hostile: products of σx and σz remain signed permutation matrices, so a fixed monomial with an off-diagonal entry contains an A and its magnitude is <=ε. Thus c_Pε and the optimized constant 4 exp(-2)c_P²/t² are valid. Uniformity is per polynomial, not per normalized vector.
* S3 fixed-word example: finite matrix units give finite support per word, a dense cyclic span and scalar commutant; spectral gaps nevertheless accumulate at zero. Common-rate dense-test lemma correctly uses bounded spectral projections, not a uniform prefactor.
* Visibility threshold min(1,(K sum M_j)^(-1/2)) is a valid lower spectral-edge bound for K>0, with zero-M case handled as in the note. It does not exclude spectrum exactly at the threshold. The note does not claim that it does.
* Finite disk Z₂ formula and torus conditioned formula are correct for the stated disk/torus cellulations. Relative-error control indeed depends on F-2a; pointwise fixed-a volume convergence is weaker than uniform macroscopic-loop control.
* Flux-parity transfer eigenvalues 1,q and its selected matrix element (1+q^a)/2 are correct. Covariance-only Gaussian/Rademacher hostile is valid.
* Particle section correctly distinguishes total-energy atoms from invariant-mass-shell components. The simple energy-spectrum examples do not purport to be relativistic local QFTs.

**Verdict:** abstract estimates and constants checked; no mathematical failure found beyond the two minor formulation repairs above. Correct the headline nontriviality status before synthesis/publication. This audit does not establish any additional YM gap, uniqueness, area law or particle classification.

## Resolution

All three requested edits were applied before the final checkpoint: the opening no longer lists nontriviality as open, K>0 is explicit, and fixed-curvature holonomy is described as a relative H¹ torsor. The narrow scope of the shrinking-gap countermodel remains unchanged.
