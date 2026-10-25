# Session31 synthesis: what the mirror, forest, and MRB structures actually share

**Final research synthesis.** Began 2026-10-22 at20:56 UTC. The statements below distinguish classical identities rederived here, proved session deductions, finite certificates, and open mechanisms. No literature-wide novelty claim is made.

## Answer to the narrow MRB question

We have not established a modular/CM mechanism generating the actual number A(8)−C_MRB. The exact honeycomb and singular-modulus identities survive, but they do not by themselves supply that missing global identity.

Two complementary results sharpen the question:

1. **An obstruction:** adding an explicit entire perturbation to the root kernel preserves every negative-axis branch jump A(m), any prescribed finite prefix of its integer values, and its positive-axis inverse-power asymptotics, yet changes the MRB alternating functional by a nonzero explicitly computable amount. Thus those data alone are insufficient in the stated Abel–Plana-admissible function class. The full analytic profile or stronger global growth conditions are essential.
2. **A positive full-object map:** the rooted-tree coefficient sequence determines a compact positive moment measure. Its Mellin transform reconstructs the complete analytic MRB kernel in the right half-plane. This is stronger than a special-value analogy, but substituting the known Mellin formula merely returns the original root series. It is a representation, not an evaluation of the tiny defect.

The simplest unweighted positive even-theta tail is rigorously excluded: it exceeds q⁴, whereas the certified defect is below q⁴. Arbitrary normalized or signed theta combinations are not excluded. A single CM specialization cannot establish an order of vanishing without an analytic family.

## 1. The coefficient sequence and the mirror equation

The supplied prefix is forced by the tree quotient

R(z)=−W₀(−z)/z=1+z+(3/2)z²+(8/3)z³+(125/24)z⁴+(54/5)z⁵+(16807/720)z⁶+⋯.

Its coefficient a_m=(m+1)^(m−1)/m! counts rooted labelled forests after multiplying by m!. For integer n>=3, the distinct positive mirror is R(log n/n). The n=2 case changes branches: this series returns2, not the distinct mirror4.

The rational mirror pairs are exactly2↔4. Gelfond–Schneider then proves every other integer-indexed mirror transcendental. In particular, 7 and21 are not exceptional in this arithmetic classification. No transfer of the tournament Hamiltonian gaps follows.

A precise natural interpretation replaces that loose analogy: for n>=3, r_n/n is the extinction probability of a Poisson(log n) Galton–Watson process. Conditioning on extinction changes its offspring mean to log r_n. This is the classical supercritical/subcritical duality in explicit mirror coordinates, not a new branching theorem.

Details and independent controls: [mirror packet](mirror/README.md).

## 2. A genuine tournament bridge, with exact losses recorded

For a uniformly random labelled tournament, let F_T(t)=det(I+tL_out) count directed rooted forests and κ=2^(n−1)H/n!. With a=1+nt/2,

E F_T(t)=a^(n−1),
E[κF_T(t)]=Σ_l binom(n−1−l,l)a^(n−1−2l)(t²/4)^l.

Conditioning on a fixed Hamiltonian path proves this exactly; it is not an illicit interchange of determinant and expectation. At t=2 the first expectation is (n+1)^(n−1)=n!a_n, directly recovering the user's coefficient family. The second formula gives strict positive covariance with every forest coefficient of degree2 through n−1, and a fixed-t asymptotic ratio through order n⁻³.

It also gives exact weighted inventory constraints for the merged metagraph. A class has h=H/|Aut| fixed-path presentations; converse pairs merge with weight2h and require averaging in/out forest observables. The resulting polynomial mass and signed-score-energy mass are fixed formulaically for every n.

This does not recover complete wiring. At order5, two explicit tournaments have identical complete in/out forest polynomials, but H=15 and13 and presentation multiplicities5 and13. The full polynomial is still a projection that loses path incidence.

Proofs, exact checks, and independent audits: [ensemble law](HAMILTONIAN_FOREST_SIZE_BIAS.md), [basic audit](mirror/HAMILTONIAN_FOREST_AUDIT.md), [addenda audit](bridge-audit/HAMILTONIAN_ADDENDA_AUDIT.md).

The [blue/black refinement](blue-forests/BLUE_FOREST_CONSTRAINTS.md), with its [independent audit](blue-audit/AUDIT.md), adds exact forest-weighted color masses and an orbit-compatible forest sum. Existing triangle-flow, blue parity and loop laws were recovered as prior canon, not promoted as new discoveries. A minimal even-order blue equal-energy pair is separated by its fourth forest coefficient, illustrating where the extra statistic helps without becoming complete.

## 2a. The tree-moment reconstruction and endpoint repair

The positive measure has fractional moments M(s)=(s+1)^(s−1)/Γ(s+1) for Re s>−1, a classical Nuttall–Bouwkamp/Stieltjes identity rederived with explicit domains. The full right-half-plane root kernel is w²/[Γ(1/w)M(1/w−1)]. This gives an exact route from the forest coefficient sequence to the full MRB analytic object.

Finite positive quadrature alone fails: it loses the divergent negative-moment endpoint and the reconstructed integer summands grow instead of tending to zero. An explicit continuous endpoint model repairs the leading singular terms; its difference from the true kernel is −(log k)²/(2k²)+O(log k/k²), giving an absolutely summable correction. No numerical speedup or modular evaluation is claimed for this repair.

For the mirror problem itself, exact rational Gauss bounds are useful: at the hardest integer input log3/3, the64-node bracket has width about3.26×10⁻⁸, versus about0.0906 for the equal-moment raw-series bound. The independently checked rational endpoints, not floating-point Lambert-W evaluation, supply the certificate.

Details: [tree moment bridge](tree-moments/TREE_MOMENT_BRIDGE.md), [independent moment audit](moment-audit/INDEPENDENT_MOMENT_AUDIT.md).

## 3. Certified MRB inverse clocks

Independent finite Euler transforms, Cauchy remainder bounds, and interval arithmetic certify the defect and the first four nearest-amplitude residue indices:

8,
105719004,
77401040908104410,
24845630402060641702669818431121505.

The first choice overshoots; the next three undershoot their positive residual magnitudes. Nearest amplitude is not generally the nearest integer to the continuous inverse, nor a fixed floor/ceiling rule. The signs must be retained.

For the stated dictionary, the signed greedy expansion either terminates exactly or converges absolutely, with a proved geometric residual bound and eventual quadratic scale. The next index is strictly larger whenever the process continues. These general convergence facts do not explain why the initial cancellation at8 is so small.

A separate entire deformation tests sensitivity. Its defect has a unique zero in a certified interval of width2×10⁻⁵⁷ near1.00000021113585. This proves neither a modular identity nor global uniqueness.

Details: [contour audit](mrb/MRB_CONTOUR_AUDIT.md), [inverse clocks](GREEDY_INVERSE_CLOCKS.md), [deformation](ENTIRE_MRB_DEFORMATION.md), [independent interval audit](analytic-audit/AUDIT.md).

## 4. Categories: exact scalar maps and definite obstructions

The exact scalar identity is

A(m)=sqrt(2m)m^(−1/m)S₀₀(SU(2)_(m−2))/π.

But this does not identify the loop model, its phase, a VOA character, or the MRB functional. At m=8 the polymer loop trace is0 whereas the SU(2)₆ fundamental dimension is sqrt(2+sqrt2); a pivotal map cannot preserve both. A(m) is transcendental and cannot be a fusion-category dimension. Even πA(8) is not an algebraic integer: its primitive minimal polynomial is2x⁸+8x⁴−1.

The elliptic theta nome, the usual character nome, and the quantum root of unity are different parameters. Dense and dilute n=0 theories also require a branch/phase choice. Ordinary characters and fixed finite-order twined characters cannot literally encode the unbounded-degree radical coefficient sequence under the stated finite-dimensional-grade hypotheses. The supertrace version additionally requires parity-preserving twining; the audit supplies a counterexample without it.

Details: [category note](categories/CATEGORICAL-HONEYCOMB-MRB-AUDIT.md), [independent audit](bridge-audit/CATEGORICAL_AUDIT.md).

## 5. Supplied-paper audits and useful repairs

The supplied Collatz manuscript does not prove convergence. In addition to a false halving transformation, explicit all-order families refute its cycle-contraction and post-collapse-infertility assertions. The repaired finite-word theorem independently prescribes an ascent length, a collapse valuation and the next ascent length. Arbitrarily long finite growing prefixes survive simple mod15/mod35 unit guards on the starting value. These are not native source-address guards and do not settle uniform153N coverage or construct an infinite divergent orbit.

The square-block MRB paper has a genuine printed sign error, confirmed visually. The corrected identity yields a general nonnegative curvature kernel, explicit rational tail bounds, and a sharp leading curvature asymptotic. A Galois conjugate counterexample prevents importing small real curvature as universal conjugate control. The Many Worlds paper's mathematics survives within its own idealizations and externally imposed alternating aggregation.

Details: [Collatz audit](collatz/AUDIT.md), [square-block audit](square-blocks/AUDIT.md), [independent structural audit](structural-audit/AUDIT.md).

## Remaining obligations

- Produce a target-preserving analytic transformation that evaluates or tightly explains the MRB defect, rather than re-expressing its ingredients.
- Preserve the singular Mellin endpoint in any numerical use of the tree-moment reconstruction; fixed finite quadrature alone destroys convergence of the reconstructed MRB series.
- Determine which additional local/path-incidence data make forest-weighted metagraph constraints sufficient for realizability. Current mass constraints are necessary only.
- Keep prior tournament loop remainder and uniform Collatz coverage problems open. This session supplies no missing global certificate.

All exact finite checks state their universe. Numerical diagnostics are not interval certificates unless explicitly identified. Primary-source trails, limitations, and reproducible scripts are provided in the linked notes.

Specific next proof obligations and their cheapest tests are recorded in [Next proof targets](NEXT_PROOF_TARGETS.md).
