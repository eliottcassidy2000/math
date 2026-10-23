# Equality rigidity: a248-dimensional homogeneous metric completion is simple E8

Status: PROVED from the audited completion theorem and elementary characteristic-zero Lie theory; classification of finite-dimensional complex simple Lie algebras is a CITED input only for naming the resulting simple algebra E8. No new classification or novelty claim. Initial audited files are unchanged.

## 1. Statement and scope

Keep exactly the completion hypotheses of `COMPLETION_ADDENDUM.md`: a finite-dimensional complex Z-graded Lie algebra L with an injective graded embedding of e7⋉H57, and a homogeneous degree-zero nondegenerate invariant symmetric bilinear form B. If dim L=248, then L is simple and hence isomorphic, as a complex Lie algebra, to E8.

More generally, for simple s and nontrivial irreducible symplectic V, any equality completion of dimension dim(s)+2dim(V)+3 is simple. Existence of such an equality completion is not asserted for arbitrary(s,V).

This does not yet assert uniqueness of a fully marked graded embedding, normalization of B, or a real form. The result is sufficient to rule out nonsemisimple248-dimensional contractions satisfying all the hypotheses.

## 2. Forced degree-zero coordinate

The paired lower count is247. Every extra dimension in a nonzero degree requires its opposite partner. At dimension248 the only possible profile is therefore

    L0=s⊕Ch,   L±1=V±,   L±2=Cz±.

Complete reducibility of finite-dimensional s-modules lets us choose h centralizing s. The subalgebra L0 is the direct sum of s and its one-dimensional center. Invariance makes B(s,h)=0, while homogeneity and nondegeneracy make B|L0 nondegenerate. Thus B|s is nondegenerate and B(h,h)≠0. In particular the new coordinate cannot be isotropic in this profile.

Schur's lemma gives [h,u]=a u on V+. The pairing V+×V− forces [h,v]=−a v on V−. Jacobi and the retained positive Heisenberg bracket give [h,z+]=2a z+, and the opposite pairing gives [h,z−]=−2a z−.

If a=0, h is central in all of L. Its nonisotropic orthogonal complement h⊥ is an ideal and a homogeneous metric Lie algebra of dimension247 containing the original g. That contradicts the audited minimum theorem. Consequently a≠0, and a−1 ad h is precisely the grading derivation.

## 3. No solvable radical remains

Let R be the solvable radical. The Z-grading makes R graded; the homogeneous metric supplies the perfect opposite-degree pairings. These are separate roles of the two assumptions.

Since R0 is s-stable and solvable, R0⊂Ch. As in the lower-bound proof, invariance forces [V+,V−] to have nonzero projection to s: otherwise B(u,[v,x])=0 for every x∈s, contradicting the nontrivial action and perfect pairing. Therefore neither V+ nor V− can be contained in R; irreducibility gives R±1=0.

If h belonged to R, then [h,V+]=aV+ would lie in R1, contradiction. Hence R0=0. The positive bracket and invariance force [V+,z−]≠0, so R−2=0.

Finally [z+,z−] belongs to L0 and centralizes s, hence is a scalar multiple of h. Its coefficient is nonzero because

    B([z+,z−],h)=B(z+,[z−,h])=2a B(z+,z−)≠0.

Thus z+ cannot lie in R either. All radical components vanish, so L is semisimple.

## 4. Semisimple actually means simple here

Decompose L into simple ideals I_i. Each is stable under the inner grading derivation ad h, hence graded. Projection of s into I_i is an s-equivariant map into (I_i)0. If nonzero, it supplies an adjoint s-submodule. But L0 contains exactly one adjoint s-module and one trivial module. Thus s belongs to precisely one simple ideal I.

The other simple ideals commute with s. Since V± are irreducible and nontrivial s-modules, they lie entirely in I. The original positive bracket then places z+ in I. Distinct simple ideals are B-orthogonal: invariance and [I_i,I_i]=I_i prove this directly. The one-dimensional graded space L−2 must belong entirely to one simple ideal, and nondegeneracy of B(z+,z−) forces that ideal to be I. Thus z− lies in I as well. The nonzero bracket [z+,z−] then places h in I. Hence I=L, proving simplicity.

By the classical classification, the only complex simple algebra of dimension248 is E8. The classical-family dimension formulas n(n+2), n(2n+1), n(2n−1) have no positive-integer rank solution248, and the exceptional dimensions are14,52,78,133,248. E8's contact grading gives the matching profile. Primary source for the E8 realization and contact-grading context: [Landsberg–Manivel](https://people.tamu.edu/~jml/LMsexpub.pdf). A further author-written E8 overview is [Garibaldi, E8, the most exceptional group](https://arxiv.org/pdf/1605.01721).

## 5. Why a contraction cannot evade the result

Suppose a limiting bracket obtained from E8 still satisfies Jacobi, retains the embedded positive Heisenberg algebra and E7 action, keeps this grading, and has a nondegenerate homogeneous invariant metric. Its dimension remains248, so the theorem makes the limiting algebra simple E8 again. A genuinely nonisomorphic contraction must lose at least one of those hypotheses; typically the ambient metric degenerates or the original positive bracket disappears. This statement concerns the algebra isomorphism class, not whether a chosen coordinate family has a nontrivial limit.

A cheap independent diagonal probe makes the mechanism visible. Scale the grade spaces by powers t^(w_j), keep s and h unscaled, and transport the E8 bracket. Keeping a fixed nondegenerate homogeneous metric forces w_−j=−w_j. Keeping [V+,V+]=Cz+ nonzero in the limit forces w2=2w1. Thus w_j=jw1 on every grade. These are ordinary grading automorphisms; every nonzero bracket exponent w_i+w_j−w_(i+j) is zero. Such a diagonal metric-preserving contraction is trivial.

This probe is deliberately narrower than the theorem: it tests homogeneous diagonal contractions with fixed B. It does not purport to enumerate all degenerations or all metric families.

## 6. What is genuinely repaired

The extra coordinate h cannot remain a harmless central appendage. Nondegeneracy forces it to be nonisotropic, and the247 obstruction forces it to act. Its action supplies a genuine inner grading element; the opposite central directions must then bracket back to h. Together these conditions eliminate the solvable radical and force a simple completion.

So, under the stated compatibility conditions,248 is not only the first possible integer dimension. It is the point where metric compatibility forces the exceptional simple ambient algebra. Marked-embedding uniqueness remains outside the proved claim.

Reproduce numerical controls with `python3 session35-lie/equality_checks.py`. Those checks verify arithmetic profiles and the diagonal scaling equations, not the abstract radical proof.
