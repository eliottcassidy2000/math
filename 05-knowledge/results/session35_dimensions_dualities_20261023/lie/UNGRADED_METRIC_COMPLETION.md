# Removing the grading hypothesis: minimal metric completion and automatic contact grading

Status: PROVED elementary Lie-theoretic deductions from standard characteristic-zero complete reducibility, Levi decomposition and invariant-form facts. E7 existence is CITED through the standard E8 contact grading. The classical classification is used only to name the resulting simple248-dimensional algebra E8. No novelty claim or general existence classification. This separate note strengthens, but does not overwrite, the earlier frozen graded-completion notes.

## 1. General statement

Let s be a finite-dimensional complex simple Lie algebra, let V be a nontrivial finite-dimensional irreducible s-module carrying a nondegenerate invariant alternating form ω, and set

    g=s⋉H(V),   H(V)=V⊕Cz,   [u,v]=ω(u,v)z.

Write d=dim s and n=dim V. Let L be any finite-dimensional complex Lie algebra with a nondegenerate invariant symmetric bilinear form B, containing g as a Lie subalgebra. No grading, no semisimplicity, and no condition on B|g are assumed.

**Theorem.**

1. dim L≥d+2n+3.
2. If equality holds, L is simple and the centralizer T=C_L(s) is sl2.
3. At equality, a suitable h∈T gives a five-step grading extending the original degrees0,1,2 on s,V,z, respectively. The metric is homogeneous for this grading.

Equality need not exist for every(s,V). For s=e7 and n=56, E8 realizes equality. Thus the unrestricted minimum dimension is248, and every248-dimensional complex metric ambient Lie algebra containing E7⋉H57 is simple E8. No uniqueness of a fully marked embedding or metric normalization is asserted.

## 2. Initial metric lower bound without a grading

Every invariant symmetric bilinear form on g factors through g/H(V)=s and is a scalar multiple of the Killing form there. Here is the general proof, also underlying the initial E7 note.

The algebra g is perfect, since [s,s]=s, [s,V]=V and [V,V]=Cz. Thus its central element z is orthogonal to all of g under every invariant form. Further, V is not the adjoint s-module: the adjoint module has its unique invariant bilinear form symmetric, whereas V has its invariant form alternating. Hence the pairing V×s vanishes by Schur's lemma. Finally

    B(u,[v,a])=B([u,v],a)=ω(u,v)B(z,a)=0

for u,v∈V and a∈s, so B(V,V)=0. This proves the assertion about forms on g.

If B|g=0, then g is isotropic in L and

    dim L≥2dim g=2d+2n+2>d+2n+3,

because a complex simple s has d≥3. Therefore all candidate dimensions at or below d+2n+3 must have B|s nondegenerate and rad(B|g)=H(V).

Since H(V)⊂g⊥, nondegeneracy gives

    dim L≥dim g+dim H(V)=d+2n+2.

Only this paired dimension and the next dimension need analysis.

## 3. The forced s-module profile

Put N=d+2n+2, and assume dim L=N or N+1. Pairing with H gives a surjective s-equivariant map L→H*, whose kernel is H⊥. As B(g,H)=0, it induces a surjection

    L/g → H*.

Its kernel has dimension0 or1, respectively. A one-dimensional representation of simple s is trivial. Complete reducibility of finite-dimensional s-modules therefore gives

    L ≅ s ⊕ (V⊗M) ⊕ T,
    dim M=2,   dim T=2 or3.

Here V≅V* using ω, and T is precisely the s-invariant subspace, hence the centralizer C_L(s). It is a Lie subalgebra. The three isotypic blocks are B-orthogonal. Each restriction is nondegenerate. In particular B|T is an invariant nondegenerate symmetric form.

On W=V⊗M, Schur's lemma gives

    B|W=ω⊗ε

for a nondegenerate alternating form ε on M. The two alternating signs make their tensor product symmetric. T acts s-equivariantly on W, and invariance places its multiplicity action in

    ρ:T→sp(M,ε)≅sl2.

No claim of faithfulness has yet been made.

## 4. The solvable radical can only lie in the centralizer

Let R=rad(L). It is stable under s. There is only one adjoint s-module in L, namely the embedded s itself, so solvability prevents R from containing an adjoint summand. Hence B(R,s)=0.

Because [R,L]⊂R, invariance gives

    B(R,[L,s])=B([R,L],s)=0.

But [L,s]=s⊕W, and B is nondegenerate on this subspace. Therefore R⊂T. Since R is an ideal and is s-trivial, [R,W]=0: the bracket map R⊗W→L is s-equivariant, its image lies both in the V-isotypic part and in the trivial module R. Also [R,s]=0 by definition of T.

Thus the only possible noncentral action of R is inside T itself.

## 5. Small quadratic centralizers

For any quadratic Lie algebra T, the tensor

    τ(a,b,c)=B([a,b],c)

is alternating. In dimension2, Λ³T*=0, so τ=0 and nondegeneracy forces T abelian.

In dimension3, either τ=0 and T is abelian, or τ is a nonzero volume form. In the latter case the contraction map Λ²T→T* is an isomorphism, so [T,T]=T. A nonzero perfect3-dimensional complex Lie algebra is sl2: Levi decomposition has a nonzero semisimple part, whose smallest possible dimension is3, leaving no radical. Thus a3-dimensional quadratic T is either abelian or sl2.

If T is abelian, the preceding section makes R central in L. Levi decomposition then gives L=S⊕R, with S semisimple and [L,L]=S. Invariance makes S and R orthogonal, and the restrictions of B to each are nondegenerate. As g is perfect, the embedded g lies in S. Consequently any nonzero central R would yield a smaller quadratic ambient algebra containing the same g.

If T≅sl2, the solvable ideal R⊂T is an ideal of T, so R=0 directly.

## 6. Excluding the paired dimension N, without classifying simple algebras

At dim L=N we have dim T=2, so T is abelian and R is central. If R≠0, the smaller quadratic algebra S containing g contradicts the initial lower bound dim S≥N. Hence R=0 and L is semisimple.

Now the map ρ:T→sl2 is injective: an element of its kernel commutes with W, with s, and with the abelian T, hence is central in L; semisimplicity makes it zero. But an abelian subalgebra of sl2 has dimension at most1. Indeed the centralizer of a nonzero traceless2×2 matrix within sl2 is its own one-dimensional span. This contradicts dim T=2.

Thus dim L≥N+1=d+2n+3. No grading or classification of simple Lie algebras was needed for this lower bound.

## 7. Equality forces T=sl2 and L simple

At dim L=N+1, the centralizer has dimension3. If it were abelian, R would be central. Nonzero R would leave a quadratic ambient S of dimension≤N, now impossible by the proved lower bound. Thus R=0; the same faithfulness argument would embed an abelian3-dimensional T in sl2, impossible. Therefore T≅sl2 and R=0.

To see that semisimple L is actually simple, decompose it into simple ideals. Projection of the embedded s into any ideal is s-equivariant, and a nonzero projection supplies an adjoint s-submodule. Since L has only one adjoint s-copy, s is contained in one simple ideal I. All other ideals commute with s, so both nontrivial V-copies lie in I. The original z=[u,v]/ω(u,v) also lies in I. Thus I contains s⊕W and at least one dimension of T. Any remaining ideals together have dimension at most2, whereas every nonzero complex semisimple algebra has dimension at least3. There are none. Hence L is simple.

For the E7 case this gives simple dimension248, which is E8 by the classical classification; the classical dimension formulas have no248 solution, and E8 is the sole exceptional algebra of that dimension.

## 8. The missing grading is recovered automatically

At equality, ρ:T≅sl2→sp(M)≅sl2 is nonzero. If it were zero, T would commute with s,W and itself except for its own bracket, but invariance would give B([W,W],T)=B(W,[W,T])=0. That contradicts the nonzero original z∈[V,V]⊂T and nondegeneracy of B|T. Since T is simple, ρ is an isomorphism. Equivalently, after simplicity of L is known, the subalgebra generated by s and W is T-stable and hence an ideal of L; it equals L. Any element of T acting trivially on W centralizes these generators and therefore is central in L, hence zero. This gives a second direct faithfulness check.

The original V⊂g is an s-submodule V⊗ℓ for a line ℓ⊂M. Its central element z belongs to T, is nonzero, and satisfies [z,V⊗ℓ]=0. Thus the nonzero traceless2×2 matrix ρ(z) annihilates ℓ. It is nilpotent, with kernel and image ℓ.

Choose a basis of M with first vector spanning ℓ. Complete ρ(z) to a standard sl2-triple, and choose h∈T so that h has eigenvalues+1 on ℓ and−1 on its complement, with [h,z]=2z. Then ad h yields

    L−2=Cz−,  L−1=V−,  L0=s⊕Ch,  L1=V,  L2=Cz.

The original g is exactly s⊕L1⊕L2. Invariance gives (i+j)B(Li,Lj)=0, so the metric is homogeneous of degree0. The contact grading assumed in the earlier theorem is therefore a consequence at equality, not an extra hypothesis.

## 9. Hostiles and boundaries

- This is an embedding theorem for the original Lie bracket, not for a190-dimensional vector space or a quotient of g.
- The form on the ambient L is nondegenerate, but its restriction to g necessarily has57-dimensional radical in the minimal E7 case.
- The nontrivial irreducible symplectic hypothesis is essential to the two-copy multiplicity argument. Reducible modules, nonsimple Levi factors and positive characteristic are not covered by this proof.
- The low-dimensional centralizer argument applies at the two candidate minimum dimensions. It does not classify larger quadratic envelopes; the cotangent algebra of dimension380 remains a valid nonsemisimple completion.
- The theorem identifies the complex ambient algebra at E7 equality and recovers a compatible grading. It does not by itself classify every marked embedding, invariant-form normalization, or real form.
- For the specified split real original g, complexification gives the same lower bound. The real equality-algebra classification is not asserted here.

CITED E8 construction: [Landsberg–Manivel](https://people.tamu.edu/~jml/LMsexpub.pdf). Classification/standard E8 context: [Garibaldi](https://arxiv.org/pdf/1605.01721). The radical and multiplicity proofs above are explicit deductions, not imported from those papers.

Reproduce elementary matrix controls with `python3 session35-lie/ungraded_checks.py`. They verify the sl2 action, its invariant alternating form, and the dimension ledger; they do not computationally prove the general Lie-theoretic theorem.


## 10. Constructive control: every odd Heisenberg dimension occurs

Take V of dimension2m with its standard symplectic form, m≥1, and let s=sp(V). The larger symplectic vector space W0=Ce+⊕V⊕Ce− has Ω(e+,e−)=1 and Ω|V=ω. Inside sp(W0), define N_u by

    N_u(e−)=u,   N_u(v)=ω(u,v)e+,   N_u(e+)=0,

and define Z(e−)=2e+, with Z zero on V and e+. Direct calculation gives N_u,Z∈sp(W0), [N_u,N_v]=ω(u,v)Z, and [N_u,Z]=0. The embedded sp(V) acts on N_u in its standard representation and fixes Z. This explicitly embeds sp(V)⋉H(V) into sp(W0).

The standard symplectic Lie algebra has a nondegenerate Killing form. Its dimension is

    dim sp(2m+2)=(m+1)(2m+3)
                 =dim sp(2m)+2dim V+3,

so the general lower bound is attained for every m≥1. The diagonal h with standard-space weights(+1,0,−1) gives the same contact grading. This is a constructive equality family, not a dimensional guess.

For the same abstract H57 (m=28), the full symplectic Levi algebra has dimension1596, its semidirect product has dimension1653, and the completion sp58 has dimension1711. The exceptional choice E7 instead gives190 and248. Thus the arithmetic property of57 does not select the exceptional Levi action: that action is additional structure. No identification with another occurrence of1596 is asserted.
