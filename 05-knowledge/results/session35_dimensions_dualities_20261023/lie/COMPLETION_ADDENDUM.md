# Addendum: the least homogeneous metric completion has dimension248

Status: PROVED elementary lower bound under explicit hypotheses, with standard semisimple Lie theory as CITED inputs. E8 realizes equality. No uniqueness assertion and no unrestricted minimality assertion. This addendum does not modify the initial audited note or executable.

## 1. Precise completion problem

Work over C. Let g=s⋉(V⊕Cz), where s=e7, V is its56-dimensional symplectic fundamental representation, [u,v]=ω(u,v)z, and s acts trivially on z. Give g degrees0,1,2 respectively.

A homogeneous metric completion is a finite-dimensional Z-graded Lie algebra L=⊕L_j with:

1. An injective graded Lie map g→L, retaining the specified s-action and bracket.
2. A nondegenerate symmetric invariant bilinear form B on L.
3. Degree-zero homogeneity: B(L_i,L_j)=0 unless i+j=0.

No simplicity or semisimplicity of L is assumed. No assertion is made that every metric completion admits this compatible grading.

**Theorem. The minimum dimension of such a completion is248.**

## 2. Pairing lower bound and the forbidden247 case

Nondegeneracy and homogeneity give perfect pairings L_j×L_−j, hence dim L_j=dim L_−j. Therefore

    dim L ≥ 133+56+56+1+1 =247.

If equality holds, its entire grading is forced:

    L_0=s, L_1=V+, L_−1=V−, L_2=Cz+, L_−2=Cz−,

where V+=V and z+=z. The s-module V− is dual to V+, hence irreducible56-dimensional; z− is trivial. B on s is nondegenerate, and the opposite-grade pairings are nondegenerate.

Let R=rad(L), the solvable radical. It is characteristic and therefore stable under the grading derivation D, so is graded. Each R_j is also s-stable.

**Claim1: R_0=R_1=R_−1=0.** The component R0 is an ideal of s. Since s is simple and R solvable, R0=0. To exclude either V component, observe that [V+,V−] has a nonzero projection into s (indeed L0=s). Otherwise invariance would give

    0=B([u,v],a)=B(u,[v,a])

for every u∈V+, v∈V−, a∈s, contradicting nontrivial s-action and the nondegenerate opposite pairing. Irreducibility makes each R±1 either0 or allV±; the latter would force a nonzero [V+,V−] inside R0, contradiction.

**Claim2: R_−2=0.** Choose u,v∈V+ with ω(u,v)≠0. Since B(z+,z−)≠0,

    B(u,[v,z−])=B([u,v],z−)=ω(u,v)B(z+,z−)≠0.

Thus [V+,z−]≠0 in V−. If z− belonged to R, it would force R−1≠0, contradicting Claim1.

Consequently R is either0 or Cz+.

**Case R=0.** Then L is semisimple. Every derivation of a finite-dimensional complex semisimple Lie algebra is inner, so D=ad(h). Because [D,ad(h)]=ad(Dh)=0 and the center of L is zero, Dh=0. Thus h∈L0=s. Since D vanishes on s, h centralizes s, so h lies in the zero center of s. This forces D=0, contradicting D|L1=1.

**Case R=Cz+.** Idealness and the grading imply [R,L±1]=0 and [R,L−2]=0 because the resulting radical components vanish; s already acts trivially on z+. Thus R is central. Levi decomposition then gives L=S⊕R as a direct sum of a semisimple Levi algebra and a central abelian algebra, so [L,L]=S has zero intersection with R. But z+=[u,v]/ω(u,v) is a nonzero element of [L,L]∩R, contradiction.

Both cases fail. Hence dim L≥248.

## 3. Equality and two constructive repairs

The E8 contact grading has dimensions1,56,134,56,1. The degree-zero component is e7⊕Ch. Its Killing form is nondegenerate, invariant and homogeneous, and the original g embeds as e7⊕L1⊕L2. Therefore E8 is a completion of dimension248, attaining the lower bound. This uses the same standard contact-grading structure and root-space input as the initial note.

The extra degree-zero direction is not cosmetic. It permits the grading derivation to be inner, pairs the ±2 directions through [z+,z−], and participates in [V+,V−]. The theorem proves necessity of at least one additional dimension beyond the formal paired247 count, but does not classify all248-dimensional completions or their real forms.

An unconditional, larger repair for any finite-dimensional Lie algebra g is its cotangent algebra T*g=g⋉g*, using the coadjoint action. It has invariant nondegenerate pairing

    B((x,α),(y,β))=α(y)+β(x).

Here its dimension is380. Assign g*_−j dual to g_j to make this a homogeneous metric completion. It retains the distinguished E7 action but is much larger than E8. This is a positive control, not a minimality claim.

## 4. Why adding only the grading element does not suffice

The parabolic p=Ch⋉g has dimension191, with [h,V]=V, [h,z]=2z and [h,s]=0. Hence

    tr(ad h|p)=56+2=58≠0.

Every finite-dimensional Lie algebra with a nondegenerate invariant symmetric bilinear form is unimodular: invariance makes each ad x skew-adjoint, and over characteristic0 its trace is zero. Thus the191-dimensional repair is still nonmetric. This obstruction applies independently of any grading requirement on its proposed metric.

## 5. Hostile assumptions and exact scope

- Merely adding dual56 and dualcentral1 produces the tempting vector-space count247, but the theorem proves no bracket can satisfy all the required hypotheses. A dimension count alone misses Jacobi/Levi compatibility.
- Homogeneity is load-bearing: without it the proof that the radical splits into degree components and that L0=s applies only if a compatible grading is supplied. The unrestricted minimum is not claimed here.
- Retaining the original bracket is load-bearing. Replacing H57 by an abelian57-dimensional module changes the problem; it is not a completion of g.
- Changing the real form is not an arithmetic operation on the dimension. The result extends to the specified split real situation by complexification, provided the real grading and metric satisfy the same hypotheses.
- E8 provides existence, not a proof that its bracket is the unique equality solution.
- A nondegenerate ambient form need not restrict nondegenerately to the original g. In fact it cannot, as proved in the initial note.

## 6. Strengthened VOA obstruction, distinct from rationality

CITED input: for a **simple VOA of strong CFT type**, the invariant bilinear form is nondegenerate, pairs distinct weights orthogonally, and therefore restricts nondegenerately to V1. The strong-CFT condition includes V0=C1, no negative weights, and L(1)V1=0. See Dong–Mason, Section2, pp3–4: https://arxiv.org/pdf/math/0201318 .

Applying the initial invariant-form theorem therefore rules out g as the full weight-one Lie algebra already for a simple strong-CFT-type VOA, without needing rationality or C2-cofiniteness. This does not rule out a nonsimple intermediate vertex algebra, a different conformal structure failing strong CFT type, or an embedding of g into a larger weight-one algebra such as e8.

The constructive conclusion is sharper than a prohibition: preserving the contact grading, opposite weights and an extra degree-zero direction repair the metric obstruction exactly at E8 dimension248. Which vertex construction realizes the intended intermediate character data remains a separate question.

## 7. General form of the lower bound

The proof uses no E7-specific identity. Let s be any finite-dimensional complex simple Lie algebra, and V any nontrivial irreducible s-module with an invariant nondegenerate alternating form. For g=s⋉H(V), with degrees0,1,2 and its original Heisenberg bracket, every homogeneous metric completion obeys

    dim L ≥ dim s + 2 dim V + 3.

The identical radical argument excludes the paired count dim s+2dim V+2. Existence at this lower bound is not automatic for arbitrary (s,V); E7,V56 is the verified equality case here. This separates a general obstruction from the special exceptional realization.
