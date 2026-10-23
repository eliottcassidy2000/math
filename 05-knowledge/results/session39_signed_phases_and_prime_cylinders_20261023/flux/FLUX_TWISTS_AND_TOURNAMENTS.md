# Flux attachment, quadratic refinements, and tournament constructions

23 October 2026. Status: sourced physical conventions; proved algebraic constructions and obstructions; small exact executable controls. No microscopic physical equivalence, phase-existence theorem, or Yang–Mills claim.

Repository inheritance is pinned to `e5d18e55`, not mutable HEAD. Relevant packet: `05-knowledge/results/session38_incidence_ports_electrical_20261023/`. The working checkout's HEAD was `e41a2c82b3` when inspected. No Git state was changed for this note.

## 1. What “seven flux quanta, two electrons, three attached to each” can mean

Write `m` for the number attached to each electron, to avoid confusing it with a Chern–Simons level. Use electron flux quantum Φ₀=h/e (equivalently hc/e in Gaussian units), not the superconducting h/(2e). Choose the statistical field to subtract from the applied field.

The uniform mean-field rule is

    B* = B − m Φ₀ ρ,       1/ν* = 1/ν − m.

Thus the schematic bulk ratio ν=2/7 and m=3 give ν*=2 and the bookkeeping 7−3·2=1. This is coherent as bulk flux accounting. It does not by itself specify a finite-system wavefunction or topological phase.

The exchange factor changes by exp(iπm): electrons have factor −1, so odd m produces composite bosons and even m preserves composite fermions. This follows directly from the exchange half-turn in the flux-attachment phase, or from the exchange parity of ∏(z_i−z_j)^m. Jain's original composite-fermion construction specifically uses even attachment [1]; the composite-boson field-theory route is different [2].

On a sphere, a homogeneous Jastrow factor ∏_{i<j}(u_i v_j−v_i u_j)^m adds degree m(N−1) to each particle coordinate. Therefore its exact flux-degree relation is

    N_Φ* = N_Φ − m(N−1).

For literal N=2, N_Φ=7, m=3, that is 4, not 1. The usual even-attachment spherical formula is explicitly stated in [3]; the odd-m version above is the same polynomial-degree derivation. A sphere lowest Landau level has N_Φ+1 orbitals. If “seven” instead counts orbitals, it corresponds to N_Φ=6. Geometry and counting convention matter. The thermodynamic filling relation can have a finite-size shift N_Φ=ν⁻¹N−𝒮; its orbital-spin data are additional to filling alone [4].

### A positive conditional completion of the user's m=3 idea

A bosonic integer Hall state can have

    K_B = [[0,1],[1,0]],     t = (1,1)^T,     ν_B = 2.

This is the interacting bosonic Hall state described by Senthil and Levin, not two filled free-fermion Landau levels [5]. Under the standard rank-one flux-attachment update, three units give

    K = K_B + 3 t t^T = [[3,4],[4,3]],
    K⁻¹ = (1/7)[[-3,4],[4,-3]],
    det K = −7,       t^T K⁻¹t = 2/7.

The inverse identity is immediate by multiplication. More generally, the rank-one inverse identity gives ν=ν_B/(1+mν_B) whenever the denominator is nonzero. The update changes the parity of each unit-charge local particle by m, as required. This provides a consistent effective topological theory conditional on the composite bosons occupying the specified gapped state.

The same matrix is

    K = −I₂ + 4 t t^T,

the two-channel reverse-field Jain form with four attached vortices. These are compatible descriptions of this K-matrix, not a reason to say three is even or that flux attachment numbers are uniquely intrinsic. Identifying full rotationally invariant phases would also require matching orbital-spin/shift and symmetry data. No finite N=2 phase claim follows.

For this K, the discriminant group is Z/7. With the representative ℓ=(1,0), charge in units of e is 1/7 and the representative exchange phase is exp(−3πi/7). Since K is odd, replacing ℓ by a local electron changes the exchange phase by a sign; the spin-theory qualification in §2 is essential.

A recent primary article, Yang et al., published 22 June 2026, includes precisely a schematic of two electrons sharing seven flux quanta in a ν=2/7 parton discussion [6]. It describes collective excitations in a many-electron Hall fluid, not a proof that a literal two-electron sphere is at the thermodynamic phase. More specifically, its Methods and Extended Data Fig. 1 use four vortices per electron for conventional composite fermions, with residual reversed flux 7−4·2=−1. In its alternative parton drawing, two charge-3e/7 fermionic partons share three flux quanta, while two charge-4e/7 bosonic partons share four. It does not assign three flux quanta to each electron. Our m=3 plus bosonic-integer-Hall construction is an alternative conditional EFT, not the paper's stated convention.

## 2. The algebra required for topological twists

For compact U(1) gauge fields a_I, an Abelian bulk Chern–Simons action has symmetric nonsingular integral K. In a fixed sign/orientation convention write the internal term as (1/4π)K_IJ a_I∧da_J, and use the associated algebraic conventions

    A_K = Z^r / K Z^r,
    b([ℓ],[m]) = ℓ^T K⁻¹m mod 1,
    Q_ℓ/e = t^T K⁻¹ℓ mod Z,
    ν = t^T K⁻¹t.

The sign of the physical Hall response also follows the chosen charge/current and spacetime orientation conventions. These formulas fix ours. The K-matrix Hall formalism is developed in [7]. Ordinary bosonic integral CS requires even diagonal K; odd diagonal K requires spin structure. The discriminant bilinear form, quadratic refinement, and appropriate framing/signature data enter the quantum classification [8]. A fusion group alone is insufficient.

Here are direct algebraic checks, independent of any physical realization. If ℓ is replaced by ℓ+Kn, the bilinear expression changes by the integer n^Tm. Thus b is well-defined and nondegenerate. For even K,

    q([ℓ]) = (1/2)ℓ^T K⁻¹ℓ mod 1,
    θ_[ℓ] = exp(2πi q([ℓ]))

are well-defined, because the change is n^Tℓ+(1/2)n^TKn, an integer. Polarization gives q(x+y)−q(x)−q(y)=b(x,y). The even-diagonal condition is exactly what makes this elementary descent work.

For odd K the same change can be a half-integer. Example K=[3]: representatives 0 and 3 are the same element of Z/3, but the naive phases are +1 and −1. One must retain the local-fermion/spin information rather than call this a bosonic twist function on the quotient. One algebraic repair uses a characteristic covector w∈Z^r with w_i≡K_ii mod 2:

    q_w([ℓ]) = (ℓ^T K⁻¹ℓ − ℓ^T K⁻¹w)/2 mod 1.

The change under ℓ→ℓ+Kn is integral since n^TKn−n^Tw is even. This is a normalized quadratic refinement, but it is not permission to identify it with the unqualified exchange phase above. Spin-theory modular/framing conventions carry additional information. For K=[[3,4],[4,3]], take w=t=(1,1), which matches its odd diagonal modulo two. Then q_w([e₁])=−2/7 mod 1, whereas the representative exchange phase remains exp(−3πi/7). They must not be silently identified.

## 3. A graph Laplacian gives a discriminant pairing, not automatically a bosonic spin assignment

Let G be a connected integer-weighted undirected multigraph. Drop loops in its electrical Laplacian, ground one vertex, and call the reduced positive-definite matrix L. Then

    Crit(G) = Z^(v−1)/L Z^(v−1),
    b_G([x],[y]) = x^T L⁻¹y mod 1

are a finite group and a nondegenerate symmetric linking pairing. The proof of descent is the same integer calculation above. Its order is det L, the spanning-tree count. With multiple components, ground one vertex in each and take the direct sum. Changing ground represents the same graph discriminant structure in different lattice coordinates.

Three distinct constructions are available:

1. Set K=L. This preserves the critical group and pairing. It is bosonic if every retained loop-free weighted vertex degree is even, and otherwise is a spin CS lattice. A physical charge vector is still an extra choice.
2. Set K=2L. This guarantees an even lattice but changes the group and pairing; |det(2L)|=2^(v−1)det L. It does not preserve the original electrical anyon count.
3. Set K_BF=[[0,L],[L,0]]. This is always even, integral, and nonsingular, of signature zero. Its group is Crit(G)⊕Crit(G), and its quadratic function is q([x],[y])=x^T L⁻¹y mod 1. This is an explicit nonchiral BF construction with the original graph pairing as electric–magnetic pairing. It doubles the sectors rather than asserting the original critical group already was a full anyon theory.

In Session38 the proven input is

    Crit(G₈)[2^∞] ≅ (Z/2)^3358 ⊕ (Z/4)^4 ⊕ Z/256.

Those Smith invariants specify a group, not a complete linking form or quadratic refinement. The divided energy forms used in its dyadic lifting tower are not automatically topological-spin functions. In particular, a Smith decomposition need not be orthogonal for b_G, so the abstract Z/256 summand is not already an independent orthogonal topological sector.

A very small exact counterexample to group-only identification is K=[256] versus K=[−256]. Both have discriminant group Z/256, but q_±(x)=±x²/512. Any group isomorphism is x↦ux for odd u. An isometry would require u²≡−1 mod 512, impossible already modulo 8 since odd squares are 1. These are distinct quadratic data on the same fusion group. Their generator twists have order 512, illustrating why “group exponent 256” is not a twist-order theorem.

## 4. A direct tournament skew-matrix bridge, with a sharp parity obstruction

For a tournament on n labeled vertices let S_ij=+1 for i→j, S_ji=−1, S_ii=0. S is skew-symmetric, not a valid symmetric CS K-matrix. Indeed the antisymmetric part of a constant coefficient matrix contributes only a boundary term to a_I∧da_J, since a_I∧da_J−a_J∧da_I=−d(a_I∧a_J). It does not define the desired bulk braiding matrix.

Nevertheless there is a canonical doubled construction when n is even:

    K_S = [[0,S],[S^T,0]].

**Theorem.** For every even-order tournament, S is nonsingular, det S is an odd square, and K_S is even, symmetric, nonsingular with signature zero. Its discriminant group has odd fourth-power order. Consequently this direct construction cannot contain the Session38 dyadic sector.

**Proof.** Over F₂ every tournament S reduces to I+J, where J is all ones. For even n, (I+J)²=I since J²=nJ=0 in F₂. Therefore det S is odd and nonzero over Z. The determinant of an even-sized skew-symmetric matrix is its Pfaffian squared. The block matrix K_S is symmetric with zero diagonal; its eigenvalues are plus and minus the singular values of S, so its signature is zero. Its cokernel is coker(S)⊕coker(S^T), and |det K_S|=(det S)², an odd fourth power. □

For odd n, det S=0 by skew symmetry, while any principal (n−1)-minor is nonsingular by the even case, so rank S=n−1. A choice of deleted vertex gives a nonsingular even-sized principal S, but that introduces a root choice; the undeleted K_S is singular and not a finite-anyon topological theory.

Relabeling vertices gives an integral congruence of K_S. Taking the converse S↦−S also gives an integral congruence, implemented by diag(−I,I). Thus this construction cannot distinguish a tournament from its converse. Seidel-type switching S↦DSD for diagonal sign D likewise induces a congruence. It is a genuine functorial algebraic construction, but deliberately forgets some tournament data.

The executable control enumerates all 64 labeled four-vertex tournaments and checks determinant, symmetry, parity, Smith factors, and converse congruence. The n=2 construction is unimodular, hence has no nontrivial discriminant sectors: a two-vertex tournament does not secretly supply seven anyons or three attached flux quanta.

## 5. Cycle holonomy is a different precise meaning of “twist”

Give the edges of an undirected support graph an arbitrary reference orientation, and let B be its vertex-edge incidence matrix (−1 at tail, +1 at head). A U(1) connection is an antisymmetric edge phase α; a gauge change α→α+B^Tf leaves every cycle holonomy unchanged, because for a cycle z∈ker B,

    z^T(α+B^Tf) = z^Tα.

Over F₃ the identical statement holds for cochains and sums rather than real phases. The number of independent graph cycle variables is e−v+c. A magnetic Laplacian with off-diagonal entries −w_uv exp(iα_uv) is Hermitian, and α→α+df acts by diagonal unitary conjugation. This constructs gauge-invariant spectral/holonomy observables, not automatically CS anyons.

Session38 establishes, with aligned vertices and reversed arcs,

    color_T(u→v) − color_Top(v→u) = c_v−c_u mod 3.

Treating this difference as a cochain on a fixed orientation, it is exactly dc. Its sum on every closed cycle is zero. Exponentiating by ω=exp(2πi/3) gives product one around every cycle, and any connections differing by it are gauge equivalent. This is a rigorous pure-gauge statement. It cannot itself produce a nontrivial Z/3 flux or anyonic topological twist. Nor does it imply that either original color assignment, separately, is a flat connection: only their aligned difference has been proved exact.

One can instead choose α_uv=φ S_uv for a tournament. On a directed triangle its holonomy is exp(±3iφ), while a transitive triangle traversed around its underlying 3-cycle has exp(±iφ). At φ=2π/3, directed triangles have trivial holonomy and transitive triangles nontrivial holonomy. This is an explicit orientation-sensitive connection, but φ is a freely chosen parameter, not a theorem that “three colors” fixes a statistical angle. Reversing all arrows complex-conjugates the holonomies.

Finally, plain incidence Laplacians B B^T of tournament edges forget every orientation: all tournaments of order n have underlying graph K_n. To retain tournament structure one must use S, an explicitly phased connection, orientation-dependent weights, or the authenticated Session38 metagraph. These are different inputs, not interchangeable names for the same object.

## 6. What would close the physical bridge

The clean next mathematical deliverable is the actual discriminant linking pairing of a specified authenticated grounded metagraph Laplacian, followed by a stated even-lattice or spin-lattice choice and a charge vector. Its consistency can be checked by exact descent, nondegeneracy, signatures and Gauss sums. That would produce a valid algebraic topological theory associated to the graph. A microscopic Hall interpretation would further require a Hamiltonian, conserved charge, dimensional/locality structure, and evidence of the claimed gapped phase. None of these are supplied by tournament enumeration, a Z/256 Smith factor, or the ratio two/seven alone.

## Sources

[1] J. K. Jain, “Composite-fermion approach for the fractional quantum Hall effect,” PRL 63, 199 (1989). https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.63.199

[2] S. C. Zhang, T. H. Hansson, S. Kivelson, “Effective-Field-Theory Model for the Fractional Quantum Hall Effect,” PRL 62, 82 (1989). https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.62.82

[3] A. C. Balram et al., “Phase Diagram of Fractional Quantum Hall Effect of Composite Fermions in Multi-Component Systems,” spherical-geometry discussion. https://arxiv.org/html/1410.7447

[4] X. G. Wen and A. Zee, “Shift and spin vector: New topological quantum numbers for the Hall fluids,” PRL 69, 953 (1992). https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.953 (The page at 69, 3000 is the erratum, not the original article.)

[5] T. Senthil and M. Levin, “Integer quantum Hall effect for bosons: A physical realization,” arXiv:1206.1604; PRL 110, 046801 (2013). https://arxiv.org/html/1206.1604

[6] Z. Yang et al., “Emergent partons in fractional quantum Hall systems,” Nature Physics 22, 1230–1237 (2026), published 22 June 2026. https://www.nature.com/articles/s41567-026-03338-9

[7] X. G. Wen and A. Zee, “Classification of Abelian quantum Hall states and matrix formulation of topological fluids,” PRB 46, 2290 (1992). https://journals.aps.org/prb/abstract/10.1103/PhysRevB.46.2290

[8] D. Belov and G. W. Moore, “Classification of abelian spin Chern-Simons theories,” arXiv:hep-th/0505235, especially §§1 and 6. https://arxiv.org/html/hep-th/0505235v1

## Reproduction

Run `python check_flux_twists.py` and `python -O check_flux_twists.py`. The saved receipts must agree exactly. Tests are algebraic controls, not numerical simulation of an electronic material or a proof of physical realization.
