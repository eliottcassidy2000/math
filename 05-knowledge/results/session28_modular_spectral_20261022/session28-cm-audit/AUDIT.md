# Independent audit: CM gluing and rank-88 Gaussian Hessians

Status: ACCEPTED for the exact claims below, with cited classical inputs and explicit boundaries. Session28, 2026-10-22. No repository source files were changed by this auditor.

## 1. CM construction at discriminant −163

The supplied certificate is checked as data, without executing the producer's script. With A = [[2,1],[1,82]], let G=A⊕A⊕A⊕A. Its discriminant group is cyclic of order163 in each factor; the dual-first-basis generator has quadratic value41/163 modulo Z. The two rows

(1,0,1,18), (0,1,−18,1)

span a size163² isotropic subgroup: both squared lengths are326, their mutual dot product is0, and the first two columns form the identity. Its inverse image in the dual lattice is therefore even and integral, of determinant163⁴/(163²)²=1. Positive definiteness is inherited from G. Rank-eight even positive-definite unimodular classification identifies the resulting lattice with E8. This is an actual gluing map, not a numerical analogy.

The supplied basis matrix B is independently checked to contain the original lattice and both glue vectors, to have determinant1/163², and to give the recorded integral even Gram matrix. Positive principal minors, determinant1, and all CM matrix identities hold. An independent exact LDL norm-shell enumeration also checks the expected240 norm-two vectors spanning rank8.

The endomorphism α satisfies

α²−α+41I=0,  α*+α=I,  α*α=41I.

Consequently J=2α−I is integral and skew-adjoint, with J²=−163I. The star means adjoint for the positive Gram form, not ordinary transpose in a nonorthonormal basis. Neither α nor J is a lattice isometry: their norm multipliers are41 and163. Dividing J by √163 gives an orthogonal real complex structure but not, in general, an integral lattice automorphism. This distinction must be retained in representation-theory interpretations.

## 2. The construction does not require the Heegner property

The producer's extension to primes p≡3 mod4 is valid. In fact, the elementary construction extends to every positive integer d≡3 mod4, using the quadratic order of discriminant −d.

Set m=(d+1)/4 and A_d=[[2,1],[1,2m]]. Choose integers a,b with a²+b²≡−1 mod d. Such a pair exists: over each odd prime field the two sets of squares and −1 minus squares, both of size(p+1)/2, intersect; the solution is nonsingular because at least one coordinate is nonzero. Lift to prime powers and use CRT. The rows

(1,0,a,b), (0,1,−b,a)

give an isotropic subgroup isomorphic to(Z/dZ)², because their leading2×2 block is I. The same determinant argument gives an even unimodular rank-eight overlattice. Multiplication by α, with α²−α+m=0, acts on each cyclic discriminant factor by the scalar(d+1)/2, so it preserves the glue subgroup. Therefore the CM action extends.

The independent test includes composite d=15,27,35,39,51,55,63,75,87,91,99. Class number one is not needed. The composite d=63 control is an arithmetic example, not an explanation of tournament repunit-count rarity. This theorem does not identify the rank-88 construction over−23 with a construction over−163, nor does it predict its extremal minimum.

## 3. Harmonic quartic normalization

Let H be a real symmetric trace-zero n×n matrix, r²=x·x, and P2(x)=xᵀH²x−tr(H²)r²/n. Then

P4(x)=(xᵀHx)²−4r²P2(x)/(n+4)−2r⁴tr(H²)/(n(n+2))

is harmonic. Directly, Δ(xᵀHx)²=8xᵀH²x, Δ(r²P2)=2(n+4)P2, and Δr⁴=4(n+2)r². All terms cancel. The equivalent expanded formula is

P4=(xᵀHx)²−4r²xᵀH²x/(n+4)+2r⁴tr(H²)/((n+4)(n+2)).

The audit checks these identities symbolically for deterministic non-diagonal trace-free H in dimensions2 through6. The proof itself is all-dimensional.

This is classical methodology, not a new harmonic decomposition: see Arne Heimendahl, *Geometric optimization problems in quantum computation and discrete mathematics: Stabilizer states and lattices* (2023), Lemma6.2.1, printedpp118–120, and the shell-Hessian discussion in§6.3, printedpp126–128. The rank88 specialization is a new application within this session, not a priority claim. Primary PDF: https://kups.ub.uni-koeln.de/65264/1/thesis.pdf .

## 4. Rank88 fixed-operator theorem

Assume L is an actual even unimodular rank88 lattice of minimum8. Harmonic theta series of degrees2 and4 have respective weights46 and48, and vanish through q³. Dividing by Δ⁴ gives a holomorphic modular form of weight−2 in the first case and weight0 in the second. Hence Θ_(L,P2)=0 and Θ_(L,P4)=C(H)Δ⁴. The modularity and valence principles are classical inputs, not consequences of a finite q expansion.

If X is the norm-eight shell and |X|=168498000, its two-design property yields

C(H)=Σ_(x∈X)(xᵀHx)²−2723200 tr(H²),

because168498000·2·8²/(88·90)=2723200. Polarize C to obtain a self-adjoint operator A on Sym₀(R88), with inner product tr(HK).

Write F(t)=Σ_(x∈L) exp(−πt||x||²). Under the determinant-one deformation exp(εH/2)L, the second derivative is

[2t²F''(t)/(88·90)+tF'(t)/88] tr(H²) + (πt)²Δ(it)⁴ C(H).

This follows by differentiating the Gaussian and applying the harmonic decomposition shell by shell. The termwise differentiations are justified by Gaussian decay for t>0. Thus the Gaussian Hessian is a(t)I+b(t)A with b(t)>0. Its eigenvectors do not depend on t, and its ordered eigenvalues are a(t)+b(t)λ_j(A). This conclusion is correct.

For a perfect minimal shell, Σ_X(xᵀHx)² is positive definite, so λ_min(A)>−2723200. That bound alone does not show a(t)+b(t)λ_min(A)>0. In particular, no claim of Gaussian universal optimality follows from perfection or packing extremality. Determining λ_min(A), or a stronger rigorous bound, remains a meaningful finite geometric obligation.

## 5. Duality and scope

Poisson duality for the actual unimodular lattice gives F(t)=t⁻⁴⁴F(1/t). Applying it to determinant-one deformations gives the corresponding Hessian identity, and Δ(i/t)=t¹²Δ(it) implies b(1/t)=t⁴⁴b(t). Therefore a(t)/b(t)=a(1/t)/b(1/t). Symmetry makes t=1 stationary, not necessarily globally minimal. A numerical plot or finite sample cannot establish the latter.

The underlying spherical-design/energy framework is classical; for example Coulangeon, *Spherical designs and zeta functions of lattices*, https://arxiv.org/abs/math/0611735 . Statements here concern lattice deformations, not arbitrary periodic sets. Related periodic-set literature has errata; no broader periodic-optimality theorem is imported.

## 6. Audited fixed-operator consequences

The operator A has trace zero. For an orthonormal basis H_j of Sym₀, Σ_j(xᵀH_jx)²=||xxᵀ−||x||²I/88||²=||x||⁴·87/88. Summing over the minimal shell agrees exactly with2723200·dim(Sym₀), where dim(Sym₀)=3915. Sakamoto's Proposition 6.2 (https://arxiv.org/html/2609.15633v2) supplies a separate congruence proof that this shell is not a4-design; consequently A≠0 and has both positive and negative eigenvalues. The positive scalar b(t) means the Gaussian Hessian is genuinely nonscalar at every temperature. This is not the assertion that the Hessian itself is indefinite.

For completeness, the nonvanishing follows directly from the zonal test. For a minimal vector a, put H_a=aaᵀ−(8/88)I. Then C(H_a)=S4(a)−261427200, where S4(a)=Σ_X(a·x)⁴. Proposition6.2 proves S4(a)≡0 modulo9 whereas261427200≡6 modulo9. Thus C(H_a) cannot vanish. The cited congruence is a source theorem, not an independent rank88 enumeration.

As t tends to infinity, the minimum-eight shell gives

F(t)=1+168498000 e^(−8πt)+O(e^(−10πt)),

and therefore

r(t)=a(t)/b(t)=2723200−15318000/(πt)+O(e^(−2πt)).

For a perfect shell, r(t)+λ_min(A) is eventually positive. By reciprocal duality it is also positive for sufficiently small t. Any temperature where the Hessian is not positive definite thus lies in a compact subinterval of(0,∞). Each eigenvalue has the form b(t)(r(t)+λ_j), a real-analytic function which is not identically zero because of its positive tails. Its zeros in that compact interval are finite. There are finitely many eigenvalues, so there are finitely many transition temperatures, paired under t↔1/t. At t=1 a zero has even order in log t because r(e^u) is even. These statements locate neither the transition values nor their existence.

The quartic operator also commutes with the conjugation representation of Aut(L) on Sym₀, directly from invariance of the minimal shell. Representation blocks may aid its computation. Commutation alone does not force scalar action on repeated irreducible representations; multiplicity spaces must be retained.

## 7. Reproduction and limits

Run `python3 check_independent.py`. It consumes the producer's `check_cm163_glue.json` in the sibling `session28-rank88` directory and writes `check_independent.json`. The input hash is recorded. Tests cover exact rank8 lattice/CM certificates, root enumeration, composite-discriminant controls, and harmonic polynomial identities. They do not enumerate rank88 short vectors, independently reprove the cited existence/minimum certificate, compute A, or prove Gaussian Hessian positivity.
