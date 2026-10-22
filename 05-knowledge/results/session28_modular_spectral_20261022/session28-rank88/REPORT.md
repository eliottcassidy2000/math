# Rank 88: Hecke closure, CM gluing, and a fixed Hessian sidecar

Session 28, 2026-10-22. Status labels apply section by section. No RH, Collatz, tournament-spectrum, or optimal-packing conjecture is claimed solved. New calculations below are reproducible; no priority claim is made for classical constructions or their elementary consequences.

## 1. Source update and the actual factorization of 88 (CITED)

[Sakamoto, arXiv:2609.15633v2, 27 September 2026](https://arxiv.org/html/2609.15633v2) constructs two extremal rank-88 even unimodular lattices using a Hermitian tensor product over Q(sqrt(-23)). The dimension mechanism is **88=2·4·11**: restriction of scalars of Hermitian ranks 4 and 11. The rank-4 trace lattice is E8; the rank-11 factor has a PSL2(23) action. The paper proves minimum 8 by controlling every tensor rank, including a separate rank-four obstruction, and establishes perfection and eutaxy, hence local packing extremeness. These are source claims, not a fresh independent certification of every computational sublemma. They do not assert globally optimal packing or universal energy optimality. The field is **-23, not -163**.

[Nebe–Sloane catalogue ideal115L](https://www.math.rwth-aachen.de/~Gabriele.Nebe/LATTICES/dim88b.html), updated October 17, supplies a Gram matrix. `check_catalogue.py` verifies all 7,744 entries form a symmetric even integral matrix, all 88 exact LDL pivots are positive, and determinant is exactly 1. Minimum 8 remains CITED; determinant and positivity alone do not certify it. Raw source HTML is a local input and should not be republished.

## 2. Exact scalar theta and full Hecke closure (PROVED on standard modular-form inputs)

Use q=exp(2πiτ) and theta=sum q^(|x|²/2). Weight is 44. The basis

    B_j=Δ^j E4^(11−3j),  j=0,1,2,3

has leading terms q^j, so prescribing coefficients q^0 through q^3 determines a form uniquely. The extremal theta is

    F=E4^11−2640 E4^8Δ+1813680 E4^5Δ²−244992000 E4²Δ³
     =1+168498000q⁴+2480127344640q⁵+6298997202432000q⁶+… .

`check_hecke.py` computes exact integer q-expansions through q^48. For prime p,

    a_n(T_p f)=a_(pn)(f)+p^43 a_(n/p)(f),

with the second term zero unless p divides n, and a_0(T_p f)=(1+p^43)a_0(f). See the [LMFDB reviewed Hecke definition](https://www.lmfdb.org/knowledge/show/cmf.hecke_operator).

T2F is not a multiple of F: its q² coefficient is 168498000 whereas F's is zero. More strongly, the coefficient-coordinate determinant of

    [F, T2F, T2²F, T2³F]

in this basis equals

    −41616498129939502499983528489711554439151616000000000,

so the smallest T2-stable space containing F is ALL of M44, of dimension four. Any exact linear Hecke representation carrying F and respecting this action therefore needs four dimensions; this does not prohibit a scalar-valued modular function with a nontrivial Hecke orbit.

The characteristic polynomial is

    (x−8796093022209)
    ·(x³+2209944x²−15663522502656x−19976984434430705664).

The cubic is irreducible over Q and has three real roots. Thus the rank-88 scalar theta is genuinely modular but not a Hecke eigenform. Milne's [Modular Functions and Modular Forms](https://www.jmilne.org/math/CourseNotes/MF.pdf), Proposition 5.5, gives the general obstruction: a nonzero simultaneous eigenform of positive weight has nonzero q coefficient. Our explicit q² witness detects failure already for T2.

The program also constructs T3, checks T2T3=T3T2, and verifies the represented action against extra q-coefficients through q^24 and q^16, respectively. Root independently reproduced the theta and cusp polynomial using a different E4/E6 basis and product expansion for Δ.

## 3. An honest self-adjoint realization, and what it is not (PROVED)

Let P be the cubic above, and C its multiplication-by-x companion matrix. The rational trace Gram matrix G_ij=Tr(C^(i+j)), 0≤i,j≤2, obeys GC=CᵀG. Its three leading principal minors are

    3,
    103748839982208,
    19104694120990804233835415658184468070400.

They are positive, so C is self-adjoint in this exact positive inner product. The normalized trace resolvent is

    m(z)=Tr((C−zI)^(-1))/3=−P'(z)/(3P(z)).

It has three real simple poles with residues −1/3 and positive imaginary part in the upper half-plane. This is an actual finite spectral realization of the weight-44 cusp Hecke spectrum, not a realization of q-plane theta singularities or zeta zeros. `check_spectral.py` supplies the complete integer certificate. Milne, section 5, supplies the classical Petersson self-adjointness context; our trace Gram is an independently constructed positive Gram, not asserted equal to the Petersson Gram.

For any self-adjoint A and finite-norm vector v,

    |<v,(A−zI)^(-1)v>| ≤ ||v||²/Im(z).

But F(iy)=y^(−44)F(i/y)~y^(−44) as y→0+. Therefore F(τ) itself cannot be such a scalar resolvent in spectral coordinate z=τ. A nonvanishing holomorphic prefactor or a nonsingular real boundary chart does not repair the exponent mismatch. This obstruction does NOT cover distributional vectors, singular reparametrizations, regularized transforms, determinants, logarithmic derivatives, or unrelated spectral constructions.

## 4. Recovering the real {2,11} hostile from prior work

The maintained repository note `05-knowledge/results/collatz_level22_oldform_bridge_20261005.md` gives U2 on the level-22 oldspace of the level-11 weight-two form:

    U=[−2  1; −2  0],   P_U=x²+2x+2,   U⁴=−4I.

Its eigenvalues −1±i are NOT real. This does not contradict good-prime Petersson self-adjointness: 2 divides the ambient level 22. An independent exact test in `check_cm_family.py` shows a symmetric Gram [[a,b],[b,d]] satisfying GU=UᵀG must have a=−2b−2d and determinant −(b+d)²−d², never positive definite. Thus a label 'Hecke operator' alone does not supply a positive self-adjoint model.

Separately, `THM-3183-factorial-hecke-lattice-square-and-oriented-wedge-continuant.md` concerns elementary modifications of integral lattices. It is not an identification with the modular operators above. These are distinct meanings linked only after specifying the object and action.

## 5. Explicit D=-163 CM gluing into E8 (PROVED on classical E8 uniqueness)

Let A have Gram [[2,1],[1,82]]. It is the binary CM norm lattice with

    Q(a,b)=a²+ab+41b²,  det(A)=163.

A generator of A*/A is e1*=(82,−1)/163, with discriminant quadratic value 41/163 modulo Z. Thus the discriminant group of A⁴ is F163⁴ with quadratic form 41Σx_i²/163. The plane generated by

    v1=(1,0,1,18),   v2=(0,1,−18,1)

is totally isotropic: each squared length is 326=2·163 and their pairing is zero. Adjoin representatives of this plane to A⁴. The resulting overlattice has index 163², is integral and even, and has determinant 163⁴/(163²)²=1. It is positive definite of rank eight, hence isometric to E8.

This is a genuine map with the missing gluing data supplied. Four copies of a determinant-163 binary lattice are NOT already E8 before gluing.

Multiplication by α=(1+sqrt(-163))/2 has matrix [[0,−41],[1,1]] on each original binary factor. On the discriminant group it acts as the scalar 82, so it preserves the gluing plane. It therefore descends to an integral endomorphism on the E8 overlattice, satisfying

    α²−α+41I=0,   α*+α=I,   α*α=41I.

Consequently J=2α−I is integral, skew-adjoint, and satisfies J²=−163I and J*J=163I. J/sqrt(163) is a real orthogonal complex structure, but it is generally not an integral lattice automorphism. The full explicit Gram and endomorphism are in `check_cm163_glue.json`; all identities were checked exactly.

### General family and its hostile to Heegner exceptionalism

The same construction works for every prime p≡3 mod4. Put c=(p+1)/4 and A_p=[[2,1],[1,2c]]. Choose a,b with a²+b²≡−1 modp, which is always possible: the two sets of squares and −1 minus squares each have (p+1)/2 elements and hence intersect. Use planes (1,0,a,b),(0,1,−b,a). This gives E8 with integral α²−α+cI=0 and J²=−pI.

In two copies the discriminant form has no isotropic line, since −1 is nonsquare; four copies are the first possible positive even unimodular rank (rank 8). The four-to-two isotropic-plane reduction is structural. Class number one is NOT needed. Prime 23 is an explicit non-Heegner control, and `check_cm_family.py` verifies all 24 primes p≡3 mod4 below 200. The count 24 in this finite test is only a cutoff count, not the Niemeier classification.

This does not produce an extremal rank-88 lattice over -163: the compatible rank-11 factor and tensor minimum proof are additional obligations. Taking three independent E8 copies gives the rootful Niemeier E8³, not Leech; root-removing gluing is again missing.

## 6. CM modular level is not erased by E8 gluing (CITED / PROVED consequence)

[Primitive Hecke-character theorem, Molin–Page Theorem 37](https://antsmath.org/ANTSXV/papers/ANTS-XV_molin-page.pdf), states that a type-(k,0) algebraic character of imaginary quadratic discriminant −D and conductor m induces a CM newform of weight k+1 and level D·Norm(m). All CM newforms arise this way. Hence a CM newform from -163 is not a level-one cusp eigenform. This concerns automorphic induction, not evaluation of a level-one modular function at a CM point.

The binary theta of A is naturally weight one with level 163 and quadratic character, not a level-one scalar modular form. Gluing A⁴ to E8 changes the lattice and sums discriminant cosets; it is not an operation preserving the original theta series or level. The original binary module, CM endomorphism, and gluing plane are the necessary sidecars.

Prior `ramanujan_heegner7_trivial_cycle_20261006.md` and its repair MISTAKE-568 already give the CM j-value −640320³, the distinction between nine Heegner fields and thirteen class-number-one orders, and the -7/-163 split-prime contrast. Those are inherited, not discoveries of this session. In particular a numerical CM value is not a modularity theorem for a tournament or Collatz generating function.

## 7. A fixed quartic sidecar controls ALL Gaussian Hessians (PROVED on harmonic-theta modularity)

Let L be any rank-88 extremal even unimodular lattice, and X its minimum shell, with |x|²=8 and N=168498000. Set F(t)=Theta_L(it). On the Euclidean space V=Sym_0(R^88), use inner product tr(HK). Parameterize determinant-one lattice deformations by exp(εH/2)L.

For harmonic homogeneous P of degree d, Theta_(L,P) has weight 44+d. It vanishes through q³. For d=2, weight 46 has no nonzero form vanishing to order at least four (valence bound); thus every shell is a 2-design, hence a 3-design by antipodality. For d=4, division by Δ⁴ gives a holomorphic weight-zero modular form, so

    Theta_(L,P)= (Σ_(x∈X)P(x)) Δ⁴.

This standard harmonic-theta mechanism is explained in [Nebe, Venkov's theory of lattices and spherical designs](https://arxiv.org/abs/1201.1834); the rank-88 specialization below follows that classical mechanism. The general fixed-shell Hessian method is also developed in [Heimendahl, thesis (2023), Lemma 6.2.1 and section 6.3](https://kups.ub.uni-koeln.de/65264/1/thesis.pdf); it is not claimed new here.

For H trace-free, decompose (xᵀHx)² into harmonic degrees 4,2,0. Its scalar spherical part is 2|x|⁴tr(H²)/(n(n+2)); its degree-two part sums to zero on every shell. Define the self-adjoint operator A on V by its quadratic form

    <H,AH> = Σ_(x∈X)(xᵀHx)² − 2723200 tr(H²),
    2723200 = 2·8²·168498000/(88·90).

Differentiating the convergent Gaussian sum twice gives exactly

    Hess_t = a(t)I + b(t)A,
    a(t)=2t² F''(t)/(88·90)+tF'(t)/88,
    b(t)=(πt)² Δ(it)^4 >0.

All these Hessians commute. They admit common eigenspaces, and their eigenvalue ordering is independent of temperature t. The smallest eigenvalue is a(t)+b(t)lambda_min(A). Thus a vast family of analytic positivity questions reduces to ONE finite angular eigenvalue and ONE scalar inequality.

This is not a claim that A is scalar or zero: a 3-design need not be a 4-design. The scalar theta determines a,b but loses A. Perfection gives Σ_X(xᵀHx)²>0 for H≠0, hence lambda_min(A)>−2723200, which by itself is too weak to guarantee positivity at every t.

Poisson duality gives a(t)/b(t)=a(1/t)/b(1/t). The initial coarse numerical grid found a(1)/b(1)≈11919.9492158637 and misleadingly suggested a minimum at t=1. **That suggestion is disproved: the self-dual point is a strict local maximum**, certified by exact symbolic reduction and rational intervals in the independent scalar follow-on described below. The original grid is retained as `probe_hessian_coarse_historical.json`, explicitly not a certificate. The updated probe samples the neighborhood of 1. The actual minimum-shell operator of either Sakamoto lattice has not been enumerated in this session. No universal energy optimality is inferred from local packing extremeness.

## 8. Connection contract and remaining work

- Scalar theta → four-dimensional Hecke module: preserve modular action, lose lattice isometry and angular shell geometry. Add harmonic theta and minimum-vector data.
- Binary CM lattice → E8: preserve the norm after overlattice inclusion and retain α, but add an explicit isotropic discriminant plane and change theta/level.
- Minimum shell → all Gaussian Hessians: retain its fourth-moment operator A; scalar radial theta alone is insufficient. Cheapest next test is calculate lambda_min(A), exploiting the finite automorphism group, then certify the one-variable lower bound a/b.
- Theta cusp → scalar resolvent: fails the growth bound in the stated coordinate. Replace the observable and prove its positivity/measure representation rather than merely renaming singularities eigenvalues.

No obstruction here excludes all modular transformations, all spectral realizations, or any open conjecture globally.

### 7a. Further rigorous consequences of the fixed operator

The dimension of V is 3915. For an orthonormal basis H_j of V,

    Σ_j (xᵀH_jx)² = |x|⁴(1−1/88).

Consequently tr(A)=0. If A is nonzero, it has both positive and negative eigenvalues; this angular anisotropy is invisible to scalar theta. A commutes with conjugation by every orthogonal automorphism of L, so finite-group decomposition of Sym_0 is a legitimate route to its eigenvalues. Irreducibility of the original 88-dimensional representation alone does not establish irreducibility of Sym_0 or make A scalar.

For t→∞ the explicit first shell gives

    a(t)/b(t)=2723200−15318000/(πt)+O(exp(−2πt)).

If L is perfect, lambda_min(A)>−2723200. Therefore every Hessian is positive definite for all sufficiently large t and, by duality, for all sufficiently small t. Possible Gaussian instability is confined to a compact temperature interval, despite local packing extremeness.

Moreover, each eigenmode has only finitely many zero crossings: a/b is real analytic on t>0; the limiting gap makes zeros impossible near either endpoint; and an analytic nonconstant function has finitely many zeros on a compact interval. Transition temperatures occur in reciprocal pairs. A zero at t=1 has even order in log(t), by reciprocal symmetry. These statements do not bound or compute the number of transitions and do not claim there are any.

### 7b. The dimension-residue classification behind this phenomenon

For an extremal even unimodular rank n=24m+r, r∈{0,8,16}, the minimum is 2(m+1). Division of a harmonic degree-four theta by Δ^(m+1) leaves weight r/2−8. For r=0 or 8 this weight is negative, so the harmonic quartic theta vanishes and all Gaussian Hessians are scalar. For r=16 it is zero, giving exactly the fixed-operator pattern above with Δ^(m+1). Thus rank 88 belongs to the **first non-scalar residue class**, not to an isolated numerical exception. This is a consequence of classical harmonic-theta modularity; `check_hessian_identities.py` audits the decomposition, trace, and rank-88 constants but does not compute A's spectrum.

### 5a. Composite discriminants, including 63 (PROVED)

The construction extends to every positive integer d≡3 mod4, not merely primes. Set c=(d+1)/4. For each odd prime power dividing d, solve a²+b²=−1 modulo the prime; at least one coordinate is nonzero, so its derivative 2a or 2b is a unit and the solution lifts to every prime power. Chinese remaindering produces a pair modulo d. The same two generators, with first two coordinates the identity matrix, generate a free isotropic subgroup (Z/dZ)² of order d². The discriminant form is cΣx_i²/d, and gcd(c,d)=1. The resulting overlattice is again positive even unimodular of rank eight, with α²−α+cI=0 and J²=−dI.

For d=63, one obtains an E8 embedding of the order of discriminant **−63**, conductor 3 in Q(sqrt(−7)), and α*α=16I. This is not a distinct field from Q(sqrt(−7)), not class-number-one evidence, and not a Collatz or tournament theorem. `check_cm_composite.py` supplies exact tests for nine composite d, including 63. Four binary copies are the smallest possible number for a positive even unimodular completion simply because its rank must be divisible by eight. The prime-only two-copy anisotropy proof above is not asserted unchanged for every composite discriminant.

### 7c. The sidecar is genuinely nonzero in rank 88 (CITED premise, PROVED consequence)

Sakamoto Proposition 6.2 proves that the minimum shell of **every** extremal rank-88 even unimodular lattice fails to be a 4-design. Therefore A≠0: otherwise polarization of the quadratic form on Sym_0, together with the 2-design identity for the trace component, would give the full isotropic fourth moment. Since tr(A)=0, A is indefinite. In particular every Gaussian Hessian at finite t>0 is genuinely nonscalar. This is compatible with positive definiteness of a(t)I+b(t)A, and does not decide it.

### 7d. Certified correction: the self-dual temperature is a local maximum

The independent companion `session28-cm-audit/scalar-followon/ADDENDUM.md` and `reduce_ratio.py` give a symbolic and rational-interval certificate. Write r=a/b and y=4π²E4(i). The Ramanujan identities give exactly

    r(1)=8192(304325−17489844/y),
    d²r(exp(u))/du² at u=0
      =1024(23708225y²−1361323584y−69959376)/y.

A rational Machin-formula bound for π, Taylor bounds for exp(−2π), and an E4 divisor-sum tail bound certify

    57.471215393911309459 < y < 57.471215393911309460.

Even the upper endpoint of the numerator polynomial is less than −21413. Hence the log-temperature curvature is strictly negative and t=1 is a strict local maximum. This reverses the original coarse-grid hypothesis; the exact certificate, not the grid, determines the conclusion.

Because r tends to 2723200 at both ends and r(1) is much smaller, a global minimum is attained away from 1. Reciprocal symmetry supplies at least two distinct off-center global minimizers. Their uniqueness, locations, and a complete global lower bound remain open here. Conditional on an angular eigenvalue satisfying −r(1)<λ<−min(r), that eigenmode would be stable at t=1 and both tails but unstable near reciprocal off-center minima. No actual lattice eigenvalue satisfying this condition has been established.
