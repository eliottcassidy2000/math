# Session 28: Hessian channels, positive moments, and the Weil boundary

Status: exact elementary results and source audit; no RH claim. Research snapshot: repository main 63012e65f0dbca29442e8f8544240eb5a7e44ff2, 2026-10-22.

## 1. An exact signed representation/Hessian bridge

Let E be the edges of K_n, oriented i→j for i<j. Put x_ij=2P(i→j)−1 and let F_n(x) be the expected number of directed Hamiltonian paths under independent edge orientations. This is a specified, unique multilinear polynomial, not an interpolation of finitely many adjacency spectra. Let B be the oriented incidence matrix with column e_i−e_j. For n≥2 set c_n=(n−2)!/2^(n−2). Then

**PROVED:** D²F_n(0)=c_n(2I−BᵀB).

Proof: the path formula is F_n=2^(1−n) Σ_(π∈S_n) ∏_(r=1)^(n−1)(1+ε_(π_r,π_(r+1))x_{unordered edge}). Multilinearity gives diagonal zero. For disjoint edges, reversing either two-vertex block pairs the contributions with opposite sign. Two edges meeting at v occur together in exactly 2(n−2)! permutations, as one consecutive three-vertex block. Their signs multiply to minus the incidence-column inner product. This gives the stated off-diagonal coefficient and hence the matrix.

Since BBᵀ=nI−11ᵀ, the edge space splits orthogonally into im(Bᵀ) and ker(B). The eigenvalues are

- c_n(2−n), multiplicity n−1, on the score/gradient space;
- 2c_n, multiplicity (n−1)(n−2)/2, on the cycle space.

For n≥3 it is a nondegenerate saddle, with trace zero. This is the differential Hessian at the unbiased interior point. It is not the discrete flip Hessian at a Boolean tournament and not the fixed-spine tiling Hessian at the transitive vertex. Their indices need not agree.

Relabelling acts on these orientation coordinates by a **signed** permutation. As an S_n-module this is Λ²R^n=V⊕Λ²V, where V is the standard representation. The unsigned permutation module on two-subsets instead has three Specht components, 1⊕V⊕S^(n−2,2), for n≥4 (for n=3 it is1⊕V, and for n=2 just1). These are two different representations of the same dimension; they must not be substituted for one another. This sharpens the repository's orientation-gauge correction recorded in MISTAKES.md on 2026-10-18 and the historical `07-reflections/representation-lens-fundamentals-2026-05-30.md`.

The exact quadratic Taylor/Walsh part is

F_(≤2)(x)=n!/2^(n−1)+c_n||x||²−(c_n/2)||Bx||².

On Boolean x, ||x||²=binom(n,2), while Bx=2(s−(n−1)1/2). Thus the quadratic component is determined by the squared score norm. This proves why a score/lattice lens captures the second-order part while higher-order information can distinguish tournaments at the same lattice point. It does not identify the Hessian with a Laplacian: BᵀB is positive semidefinite, whereas 2I−BᵀB is indefinite.

Verification: `check_structural.py` independently sums every Hamiltonian-path polynomial for n=2,...,7 (5,912 permutations in total) and checks every matrix entry exactly. The proof, not that finite census, establishes all n.

## 2. Polynomial enumeration cannot be a positive self-adjoint moment sequence

**PROVED:** Suppose a_n agrees eventually with a real polynomial P(n) of degree d>0 and positive leading coefficient. There is no self-adjoint operator A and vector v in the domains of all powers such that a_n=⟨v,A^n v⟩ for every n≥0.

Proof: b_m=a_(2m)=||A^m v||². Self-adjointness and Cauchy–Schwarz give

b_(m+1)²=|⟨A^m v,A^(m+2)v⟩|² ≤ b_m b_(m+2).

But if b_m=Cm^d+O(m^(d−1)), direct polynomial expansion gives

b_m b_(m+2)−b_(m+1)²=−d C²m^(2d−2)+O(m^(2d−3))<0

for large m. Contradiction. The statement also excludes any positive spectral measure with these moments, bounded or unbounded, when all displayed moments exist.

Application, conditional only on the inherited fixed-H counting theorem: session27's ordinary/automorphism-weighted fixed-H counts are eventually polynomials. Every case with positive padding degree is therefore excluded from this positive scalar spectral-moment interpretation. The h=63 sequence beginning at order8 is a_m=2(m+1); its even subsequence 2,6,10 already gives Hankel determinant 2·10−6²=−16. Ordinary consecutive determinant 2·6−4²=−4 already excludes ANY self-adjoint scalar moment model, by Cauchy–Schwarz on v and Av. The even-subsequence argument above is useful for eventual polynomial tails.

**Strongest survivor:** polynomial counts are naturally generalized-eigenvector/Jordan observables, not positive self-adjoint moments. For J=[[1,1],[0,1]], u=(1,0), v=(1,1), uᵀJ^m v=m+1 and Σ_(m≥0)(m+1)z^m=(1−z)^−2. This exact realization is nonnormal and uses two different vectors. It is not forbidden by the theorem.

A second survivor is a trace or determinant with different singularity rules. For finite self-adjoint A, the vector resolvent ⟨v,(I−zA)^−1v⟩ has simple poles with positive spectral weights in the partial fractions; an eigenvalue's multiplicity does not raise those pole orders. But det(I−zA)^−1 can have high-order poles even when A is self-adjoint. Repeated poles alone exclude the former interpretation, not the latter. Our Hankel proof states the precise prohibited map.

This is an application of classical moment positivity, not a new Cauchy–Schwarz theorem. Closest proved repository mechanisms are `01-canon/theorems/THM-2810-factorial-hankel-faithfulness-and-bounded-radial-carrier-no-go.md` (faithful factorial Hankel forms prohibit finite carriers) and `01-canon/theorems/THM-3056-product-gamma-reciprocal-hypergeometric-and-hankel-reversal.md` (coefficient reciprocals leave the positive moment cone). The new connection here is to eventual fixed-H metagraph padding polynomials and their repeated-pole generating functions. It is not a no-go theorem for all spectral encodings of tournaments.

## 2A. The actual decimal prime-density sequence fails a positive-moment test

**FINITE-EXACT, parent input table independently arithmetically checked:** Let δ_k=[π(10^k)−π(10^(k−1))]/(9·10^(k−1)), k≥1, with π(1)=0, and m_n=δ_(n+1). Thus δ begins4/9,21/90,143/900. For

q(x)=x(x−1)^2(x+1)(2x−1),

q(x)^2=x²−6x³+11x⁴−21x⁶+18x⁷+5x⁸−12x⁹+4x¹⁰.

Define the linear functional L(x^n)=m_n. Exact arithmetic gives

L(q²)=−2447633/22500000000<0.

Hence the **actual** density sequence cannot be the moments of a positive measure on the real line, or moments ⟨v,A^n v⟩ of one self-adjoint operator/vector. No limiting argument or RH assumption enters this finite obstruction. The input π(10^k), k=1,...,11, is4,25,168,1229,9592,78498,664579,5761455,50847534,455052511,4118054813; the parent lane owns its external source verification. Independent reproduction: `check_density_moment.py` and JSON. The source polynomial was found by the parent; this lane independently expands and evaluates it.

The control m_n=1/(n+1)=∫_0^1x^n dx gives L(q²)=1/396>0. Thus a logarithmic generating-function singularity is compatible with a continuous positive spectrum, but it does not ensure that the actual arithmetic coefficients have a positive spectral realization. This is a rejection of a particular observer/realization, not of RH or of all prime-spectrum models.

## 3. A concrete error in the proposed Conway-cusp audit

Primary source inspected: SRFP311T1 Collaboration, *A Structural Audit of Conway-Cusp Scattering and the Limits of Scalar Zero-Confinement Models*, October2026, https://srfp311t1.com/conway_cusp_scattering . Its “Cayley Transforms and Contractivity” section states that |(n−i)/(n+i)|≤1 is equivalent to Re(n)≥0. This is false:

|n+i|²−|n−i|²=4 Im(n), not 4 Re(n).

Exact hostile controls: n=−1+i has negative real part but Cayley norm²=1/5; n=1−i has positive real part but norm²=5. The right-half-plane transform (n−1)/(n+1) uses Re(n)≥0. The upper-half-plane transform (n−i)/(n+i) uses Im(n)≥0. Multiplying a boundary operator by i changes the convention and must be carried through Green's identity.

**Correct survivor:** Green's identity controls a sign-sensitive imaginary part of a Dirichlet-to-Neumann form at a complex energy. A real-part inequality requires a Rayleigh estimate if one chooses the real-half-plane Cayley normalization. The paper's claimed automatic mismatch is not proved with its displayed i-Cayley normalization.

**Correct survivor, independent of that error:** For any Re(a)>0, B_a(z)=(z−a)/(z+conj(a)) is contractive on Re(z)>0, has unit boundary modulus, and has an interior zero at a. This follows directly from |z+conj(a)|²−|z−a|²=4Re(z)Re(a). Contractivity cannot exclude interior scattering zeros. A scattering zero, a resolvent pole, and an L² eigenvalue require different definitions; they are not interchangeable.

The paper correctly labels its four-zeta coefficient as a candidate needing a double-coset/local-density calculation. The Leech identity Θ=E12−(65520/691)Δ gives the exact Epstein identity

Z_Leech(s)=2^(−s)(65520/691)[ζ(s)ζ(s−11)−L(s,Δ)]  (Re(s)>12),

with Z=Σ_(v≠0)||v||^(−2s). This is an identity for the Epstein Dirichlet series, not automatically for a particular scattering matrix entry. A channel map identifying representation-number coefficients with that entry remains necessary. Do not elevate the paper's schematic L(s−11,Δ) expression to an established scattering formula.

## 4. What a rigorous Guinand–Weil coercivity program must retain

Primary references checked:

- Connes–Consani, *Weil positivity and Trace formula, the archimedean place* (2021), https://arxiv.org/abs/2006.13771 . This establishes a scoped archimedean positivity result, not positivity of the full adelic Weil form.
- Burnol, *The Explicit Formula in simple terms* (1998), https://arxiv.org/abs/math/9810169 . It explains Weil positivity and local conductor operators; local self-adjointness does not prove global Weil positivity.
- Connes–Consani–Marcolli, *The Weil proof and the geometry of the adeles class space* (2007), Proposition6.2/Corollary6.3, https://math.jhu.edu/~kc/papers/Manin_Fest.pdf . Positivity of the full properly normalized trace pairing is an RH-equivalent criterion, not an already established positive form.
- Connes–Consani–Moscovici, *Zeta Spectral Triples* (2025), https://arxiv.org/abs/2511.22755 . The authors explicitly distinguish finite self-adjoint models/numerical agreement from the required convergence to the Xi function and all relevant zeros.
- Groskin, *A finite Guinand–Weil dictionary and archimedean tail order for the truncated Weil quadratic form* (2026), https://arxiv.org/html/2607.02828v1 . Its finite-cutoff sign-correction principle is relevant; the paper expressly makes no RH claim. This session did not independently reproduce its large interval computations.

The operative target is positivity on the full test space (with the proper involution, pole terms, normalization and domain), not necessarily a uniform strictly positive constant in an arbitrary norm. The positive operator diag(1,1/2,1/3,...) has no uniform coercivity constant; therefore positivity and coercivity must be distinguished.

**Elementary certification lemma (PROVED):** For finite Hermitian matrices, if Q∞=QT+R and 0≤R≤βI, then λ_j(QT)≤λ_j(Q∞)≤λ_j(QT)+β. Hence QT≥0 certifies the target nonnegative, and λ_min(QT)<−β certifies a target negative. A computed value in [−β,0) is inconclusive. Conversely, an arbitrary positive finite compression does not certify the omitted space: diag(1,−1) compressed to its first coordinate is the cheapest hostile.

**Complement-debt lemma (PROVED):** For a block Hermitian form Q=[[A,B],[B*,C]], assume A≥aI, C≥γI with a,γ>0, and ||B||≤β. Then

Q≥[(a+γ−sqrt((a−γ)²+4β²))/2]I.

In particular aγ>β² certifies strict positivity; equality gives nonnegativity. This follows by bounding the cross term by 2β||x||||y|| and diagonalizing the 2×2 comparison matrix. With C positive and boundedly invertible (automatic for positive-definite finite matrices), the exact condition is A−BC^−1B*≥0. This is a useful mathematical form of “pay the cross-term debt without recycling credit”: positivity of both diagonal sectors alone is insufficient.

For an unbounded form, the corresponding route needs a stated closed form, a dense form core, and controlled approximation. Positivity on every member of an increasing dense form-core family then passes by form convergence. Finitely many truncations do not. No such uniform completion is supplied here for the Weil form.

## 5. Prior-work recovery and correction ledger

Local physically present search for Guinand/Weil-coercivity/Hessian/self-adjoint themes returned447 files (paths in search_paths.txt). No local Guinand match and no connector Guinand result was found. Connector snapshots were saved for Weil(100 results, cap reached), self-adjoint(42), Hessian unbiased(3), cycle-space Hessian(61), and representation lens(100, cap reached). These are broad recovery evidence, not an exhaustive all-history claim.

Important inherited objects:

1. `07-reflections/n-2-negative-eigenvalues.md`: fixed-spine/transitive quadratic signatures and reflection sectors; separate carrier from §1.
2. `07-reflections/H-delta-field-gradient-hessian-and-the-7-21-holes-s699.md`: discrete differences and corrected even/not-always-divisible-by4 Hessian; does not supply PSD.
3. `07-reflections/polarized-delta-fields-band-gaps-and-numerical-semigroups-s699.md`: “Hessian=Laplacian” is an analogy needing the missing positivity/row-sum structure.
4. `07-reflections/representation-lens-fundamentals-2026-05-30.md`: unsigned two-subset module is valid as such; oriented sign coordinates use a different representation (§1).
5. `04-computation/tournament_L_function.py`: the defined finite product ∏(1−Q_k p^−s)^−1 has **poles**, not zeros, at s=(log Q_k+2πik)/log p for Q_k>0. The script's “zero” language is false. Paley Q=(p+1)/4 is not sqrt(p); its poles need not lie at realpart1/2. The finite rational function has no finite zeros. The survivor is a finite spectral determinant, not an arithmetic Euler product over primes or an RH criterion.
6. `04-computation/weil_rh_tournament.py`: a broad historical analogy involving Paley flatness, finite-field character sums and H maximization. A finite-field Gauss-sum estimate is not the zeta RH, and finite matrix trace identities are not a realized arithmetic trace formula.
7. `01-canon/theorems/THM-134-paley-local-max-parseval.md`: finite spectral interpolation and critical-point uniqueness need separate audit; another session28 lane owns this. Our §1 avoids both issues by specifying the independent-edge polynomial exactly.

No historical file was edited by this lane. Parent can attach precise warning repairs after review, retaining provenance rather than rewriting experiments.

## 6. Connection map and honest stopping point

Tournament orientation → signed edge module → score/cycle Hessian split preserves exact second derivatives, loses higher-order path interactions. Enumeration → positive spectral moments fails by an exact Hankel inequality; nonnormal Jordan realization survives. Leech theta → Epstein zeta preserves all shell coefficients; transfer to a scattering channel is not automatic. Weil truncation → full positivity requires a complement/tail budget and form-domain convergence. All are examples of identifying the operator and measure **before** using the word spectrum.

Open after this work: a genuinely global Weil positivity/coercivity completion; actual Conway scattering matrix coefficients; any arithmetic identification of tournament Hessian eigenvalues with zeta zeros. No numerical coincidence among dimensions supplies any of these missing maps.
