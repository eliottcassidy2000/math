# Fourth-order code potential, representation sectors, and a disappearing full gap

Status: PROVED finite-dimensional spectral formulas and asymptotics, with exact symbolic certificates. This model is independently built from Session33's public Fano/Hamming computation, not from the attached Yang–Mills manuscripts. It is not a construction of continuum Yang–Mills theory or a physical mass-gap theorem. The word Hamiltonian below means a quantum energy operator; H(T) in the tournament notes counts directed Hamiltonian paths. They are different objects.

## 1. The model and its exact fourth-order potential

Index seven bits by nonzero a in F2³. Let X=F2^7 and s(x)=Σ_a x_a a. The code C=ker(s) has16 elements. On the128-dimensional space of complex functions on X, let X_a flip bit a, let P0 multiply by1_{s(x)=0}, and define for J>=0, h>0

    H_{J,h}=h Σ_a(I−X_a)+J(I−P0).                         (1)

Both terms are positive semidefinite. The diagonal potential penalizes precisely the112 nonzero-syndrome masks. For the corresponding Fano tournament family this potential equals J(189−H(T_x))/18, using the inherited171/189 identity; this is an explicitly chosen potential, not a derivation of quantum dynamics from tournaments.

Let Z_b(x)=(-1)^{Σ_a(b·a)x_a}. Character orthogonality gives

    J(I−P0)=7J/8−(J/8)Σ_{b≠0} Z_b.                       (2)

Every Z_b is a four-bit character. Thus the potential is exactly a sum of seven quartic spin interactions plus a constant. The quartic terms commute; the transverse flip terms do not commute with the potential. Fourth-order Fourier degree is not fourth spacetime dimension and is not the differential order of an operator square.

Translations by C commute with (1). This is a precisely specified finite code symmetry. It is NOT the vertex-switching gauge group of the simplicial edge model, and is not identified with the physical gauge redundancy of Yang–Mills theory.

## 2. Reynolds projection and the exact invariant-sector gap

The average R_C=(1/16)Σ_{c∈C}U_c is an orthogonal projector onto C-invariant functions, an8-dimensional space identified with syndromes. Bit flips give each nonzero syndrome displacement exactly once, so the kinetic operator becomes h(8I−J_8), where J_8 is the all-ones matrix. The potential is diag(0,J,...,J).

The subspace with value0 at syndrome0 and sum0 on the other seven points has dimension6 and eigenvalue J+8h. On the span of the origin and the normalized sum of the other points, the operator is

    [[7h, −sqrt(7)h], [−sqrt(7)h, J+h]].

Therefore the remaining two eigenvalues are

    λ±=(J+8h ± sqrt((J−6h)^2+28h²))/2.                    (3)

The invariant-sector ground is λ−; λ+ is the next level (equal to J+8h only at J=0). Its gap is

    Δ_inv=sqrt((J−6h)^2+28h²)>=2sqrt(7)h.                (4)

This is a genuine exact finite spectral gap. As J→infinity at fixed h it grows asymptotically like J. It is only a gap in the specified invariant sector.

## 3. The full space has16 sectors, not one

Let D=C^perp be the8-element simplex code, whose nonzero words have weight4. Fourier characters indexed by q∈F2^7 diagonalize the kinetic term with eigenvalue2h·wt(q). Restricting a character to C labels a sector; the sector is the coset q+D, with8 basis characters. In that Fourier basis its operator is exactly

    diag(J+2h·wt(v):v∈q+D) − (J/8)J_8.                 (5)

The16 sectors fall into four types, with weight enumerators and numbers of sectors:

- 1+7z^4: one sector, the invariant one;
- z+4z^3+3z^5: seven sectors;
- 3z²+4z^4+z^6: seven sectors;
- 7z³+z^7: one sector.

These enumerators are independently checked by explicit linear-code enumeration. Within a group of equal weights w and multiplicity m_w, the sum-zero subspace has eigenvalue J+2hw and dimension m_w−1. The remaining block has size2 or3. With t=h/J and z=E/J, its characteristic polynomial is

    ∏_w(z−1−2tw) · [1+(1/8)Σ_w m_w/(z−1−2tw)].         (6)

The rational expression simplifies to a polynomial; (6) concerns the small remaining block, not the extra repeated-weight eigenvalues already listed.

## 4. Rigorous large-penalty asymptotics of the full gap

At t=0 each sector has a simple eigenvalue0 and seven eigenvalues1. Its unique low branch is analytic near0 by the simple-root implicit function theorem. Expanding the four exact polynomials (6) gives

    E/J = 7t−7t²+c t³+O(t^4),

with c=−42,−6,+6,+42 respectively. The standard exact checker derives these coefficients symbolically, not from floating-point fitting. Since there are only four types and all remaining branches approach1, these expansions order the full low spectrum for all sufficiently small positive t.

Thus the full ground is the invariant-sector branch and the next level belongs to the seven sectors of type z+4z³+3z^5. The full excitation gap obeys

    Δ_full=36 h³/J² + O(h^4/J^3)                         (7)

as J/h→infinity at fixed h>0. In particular Δ_full→0 while Δ_inv→infinity. At J=0 the distinction is already visible: the cube has full gap2h, whereas its invariant K8 quotient has gap8h.

For every finite J and h>0 the full ground is unique: after adding a sufficiently large scalar multiple of I, minus H is an irreducible nonnegative matrix on the connected cube, so Perron–Frobenius applies. Its positive ground vector is invariant under all commuting code translations. Uniqueness at each finite parameter does not imply a uniform lower bound as the parameter changes.

### Why cubic tunneling emerges from a quartic potential

The closest distinct zero-potential configurations differ by a weight-three codeword. One or two flips cannot connect them. Each three-flip connection has3!=6 ordered paths through penalty-J states. The seven weight-three codewords generate a degree-seven Cayley graph on C with spectrum7,1,−1,−7 and multiplicities1,7,7,1. This explains the third-order coefficients −6 times those eigenvalues. Equation (7) is certified by the exact sector polynomials, so this perturbative interpretation is explanatory, not a substitute for proof.

This is an instructive mismatch of orders: quartic interaction degree, third-order tunneling, and a second-order inverse-penalty gap scale are three different measurements.

## 5. What this permits and what it blocks in Yang–Mills reasoning

A projection can be exact for selected invariant observables while omitting slower sectors. Therefore a spectral inequality on an observable quotient is not automatically a full-space spectral inequality. But the converse warning is equally essential: if C were independently justified as an actual gauge redundancy and the invariant subspace were the physical state space, the other sectors would not be physical excitations. This toy model does not decide that identification for any continuum theory.

The parameter limit J/h→infinity is an internal finite-model limit, not a continuum lattice-spacing limit or thermodynamic limit. We prove no reflection positivity of a spacetime field measure, Euclidean invariance, renormalization, continuum existence or Yang–Mills vacuum estimate. The official four-dimensional target remains [the Clay problem](https://www.claymath.org/millennium/yang-mills-the-maths-gap/).

The usable bridge is exact and modest: the same averaging operation behind representation-theoretic Burnside projects this quartic tournament-derived model, and the omitted representation sectors can be solved and compared. The result tells us precisely which sector a proposed gap estimate controls.

## 6. Reproduction

Run `python3 check_code_quantum_spectrum.py > code_quantum_spectrum.json`, also under `-O`. Exact SymPy polynomials plus integer character sums verify all16,384 Fourier potential matrix entries, all kinetic characters, all sector enumerators and characteristic factorizations. The output records the check count. No numerical eigenvalue threshold is used.
