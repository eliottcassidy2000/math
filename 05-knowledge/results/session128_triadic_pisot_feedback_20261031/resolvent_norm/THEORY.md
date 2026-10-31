# Rational resolvent feedback is exactly the multinacci unit-norm gate

**Status:** PROVED elementary reduction uniformly at k=3,4,5; complete cubic
gain classification is CITED input + PROVED deduction; exact finite checks
corroborate the formulas. No complete quartic/quintic gain list is claimed.
This is a new Session128 connection of the previously retained sources, with
no historical-priority claim.

## 1. Source and exact feedback law

Fix k∈{3,4,5}, the field K=Q(α), and the order A=Z[α], where

    f_k(t)=t^k−t^(k−1)−...−t−1,
    F_k(p,q)=p^k−sum_(j=0)^(k−1) p^j q^(k−j).

Irreducibility for these k is certified in Session127. Let C be multiplication
by α in the ordered power basis, with subdiagonal ones and last column ones.
Its determinant is (−1)^(k+1), so C∈GL_k(Z).

For γ=a/b∈Q in lowest terms with b>0, define

    R_γ=(I−γC)^−1,
    C_γ=R_γ C,
    λ_γ=α/(1−γα)=bα/(b−aα).                                (1)

This is an exact implicit feedback update, for example

    x_(n+1)=C x_n+γ C x_(n+1)
       iff x_(n+1)=C_γ x_n.

It is not the explicit additive update C+γB, and it is not a claim that an
arbitrary nonlinear or delayed feedback law is conjugate to the original
operator. R_γ absorbs this specified algebraic loop. Equation (1) is meaningful
for every rational γ: a pole would make α or a conjugate equal to b/a, contrary
to irreducibility of degree k>1. Irrational real poles can still exist.

Since C represents multiplication by α, R_γ and C_γ represent multiplication
by (1−γα)^−1 and λ_γ. No smallness assumption or Neumann series is used.

## 2. Eight exact integrality conditions

For this specific rational-gain family, the following are equivalent:

1. λ_γ is an algebraic integer.
2. λ_γ∈A.
3. C_γ has integer entries in the named power basis.
4. R_γ has integer entries in that basis.
5. F_k(b,a)=±1.
6. b−aα is a unit of A.
7. bI−aC∈GL_k(Z).
8. U_k(γ) is an integer, where U_k(z)=z^(k−1)/(1−z−...−z^k).

**Necessity of the norm gate.** Directly from multiplicativity of field norm,

    N_K/Q(λ_γ)=(−1)^(k+1) b^k/F_k(b,a).                      (2)

Reduction of F modulo b gives F_k(b,a)≡−a^k (mod b), hence
`gcd(F_k(b,a),b)=1`. If λ_γ is an algebraic integer, its norm is an integer;
(2) therefore forces F_k(b,a)=±1. This uses the coprime numerator/denominator
structure; an arbitrary field element of integral norm need not be integral.

**Sufficiency.** Multiplication by b−aα has integral matrix bI−aC and determinant
F_k(b,a). If the determinant is ±1, its adjugate supplies an inverse inside A.
Then λ_γ=bα(b−aα)^−1∈A and R_γ=b(bI−aC)^−1 is an integral matrix.
Conversely an integral multiplication matrix has algebraic-integer eigenvalues.
Also C_γ=R_γC and R_γ=C_γC^−1, with C,C^−1 integral. This proves conditions1–7.

For condition8, primitive a,b give

    U_k(a/b)=b a^(k−1)/F_k(b,a),
    gcd(F_k(b,a),ab)=1.

Thus its reduced denominator is |F_k(b,a)|. The zero gain a=0 has b=1 and
satisfies the same conclusion directly.

This is an exact preservation of the **integrality predicate**, not an inverse
from a scalar to a source vector, nor an equivalence between all norm values
and all feedback systems.

### Fixed-S localization: the same Thue–Mahler gate

For any fixed finite set S of rational primes, the same proof gives

    λ_γ integral over Z[S^−1]
      iff λ_γ∈A[S^−1]
      iff every prime divisor of F_k(b,a) belongs to S.

Necessity follows from the localized integral norm and gcd(F,b)=1;
sufficiency follows from the adjugate after inverting S. The corresponding
resolvent matrices and U_k(γ) are S-integral under exactly the same condition.
This is the existing Thue–Mahler norm gate applied to the specified feedback
family. Under these equivalent S-integrality conditions, C_γ and R_γ are
invertible over Z[S^−1] **only if b is also an S-unit**, and conversely,
because their determinant magnitudes are b^k/|F_k(b,a)|. S-integrality alone
is not S-unimodularity.

For each fixed S and each of the fixed irreducible forms, classical
Thue–Mahler finiteness gives finitely many such rational gains. This is a
**CITED finiteness consequence**, not an explicit list or computed height
bound. The definition and finiteness statement are inspected in the
introduction, equation(2), of Gherga–Siksek,
[*Efficient Resolution of Thue–Mahler Equations*, arXiv:2207.14492v2](https://arxiv.org/html/2207.14492v2).
The ordinary cubic unit-gain list(5) is unchanged.

## 3. The generating-function value is a named coefficient of the resolvent

Let ε_0 extract the constant coefficient in the power basis. Then

    ε_0(λ_γ)=e_0^T C_γ e_0=U_k(γ).                          (3)

A short proof uses formal series, so it is valid without analytic convergence.
Let s_m=ε_0(α^m). Then s_0=1, s_1=...=s_(k−1)=0, and s_m obeys the k-nacci
recurrence. Consequently its formal generating function is

    S(z)=(1−z−...−z^(k−1))/(1−z−...−z^k).

The constant coefficient of α/(1−zα) is (S(z)−1)/z, which equals U_k(z).
This is a rational-function identity, hence holds at every nonpole, including
outside the Neumann/Taylor convergence disk.

Thus Session127's scalar −323008 is now identified precisely: it is the
constant power-basis coefficient of the exceptional closed-loop generator,
not a Tribonacci term, a volume, or a complete encoding of that generator.
The whole coefficient vector remains necessary for reconstruction.

## 4. Field degree and the exact new minimal polynomial

The inverse Möbius formula is

    α=bλ_γ/(b+aλ_γ)=λ_γ/(1+γλ_γ).

Its denominator cannot vanish: b+aλ_γ=b²/(b−aα)≠0. Therefore
Q(λ_γ)=Q(α), so the degree remains k. The monic minimal polynomial is

    P_γ(Y)=F_k(bY,b+aY)/F_k(b,a).                            (4)

Indeed the numerator vanishes at λ_γ, has degree k and leading coefficient
F_k(b,a)≠0; unchanged degree makes it the minimal polynomial after scaling.
Its constant term is −b^k/F_k(b,a), consistent with (2). It has integer
coefficients precisely at the unit-norm gains.

The new generator is in the **same field**, but the dynamics are not generally
conjugate to C: their eigenvalues are changed by λ=α/(1−γα). The explicit
polynomials below make that failure visible at the two nonzero cubic gains.

### Change of generator is a rational categorical equivalence

There is nevertheless a precise equivalence of **operator categories**.
The source has finite-dimensional Q-spaces with operators T satisfying
f_k(T)=0, with all Q-linear intertwining maps. The target uses operators D
satisfying P_γ(D)=0, with the same definition of morphism. The functors are

    T ↦ T(I−γT)^−1,       D ↦ D(I+γD)^−1,

and act as the identity on the underlying linear arrows. The necessary
inverses exist by coprimality of the defining polynomials with the linear
denominators. If H T_1=T_2 H, then H intertwines every polynomial and
available rational function of T_1,T_2. Applying the inverse rational map
proves the converse. Thus all Hom spaces agree under this change, and the
two compositions are identities.

This redefines which field generator represents one time step; it does
**not** provide a similarity between T and its image on the same vector
space. The differing spectra above are compatible with this categorical
equivalence. Unrestricted integral essential-surjectivity must not be
inferred: the proper order Z[λ_γ] generally admits lattices not stable under
the recovered generator α. Sections7–8 retain that missing arithmetic data.

## 5. Complete cubic rational gain classification

Put β=α_3. The inherited complete cubic Thue classification is rechecked
against the primary author-hosted paper: Knapp, Lemos and Neumann,
*Integral Values of Generating Functions of Recurrence Sequences*, Lemma2.4
and Example4.1 (printed pages7 and10–11). Their Delone–Nagell bound, applied
to discriminant−44, limits F_3(P,Q)=1 to four integer solutions; the four
witnesses are

    (P,Q)=(0,−1),(1,0),(2,1),(−103,−56).

Odd homogeneity gives the solutions for right side−1 by negating these pairs.
The paper also states the corresponding generating-function input list.
[Primary source](https://math.loyola.edu/~mknapp/papers/General%20Quadratic%20-%202-10-23.pdf).

Filtering this complete signed list by P=b>0, Q=a yields exactly

    γ∈{0, 1/2, 56/103}.                                     (5)

Equation(5) is a **global** classification of rational gains for which λ_γ
is integral; its completeness is not inferred from the finite checker.
It is a newly proved identification with the pre-existing classified norm/GF
slice, not a new proof or claimed discovery of the cubic Thue theorem.

Using β^−3=2−β and β^−16=56β−103 gives:

| γ | λ_γ | coefficients in (1,β,β²) | minimal polynomial |
|---|---|---|---|
| 0 | β | (0,1,0) | Y³−Y²−Y−1 |
| 1/2 | 2β⁴ | (2,4,4) | Y³−22Y²−20Y−8 |
| 56/103 | −103β¹⁷ | (−323008,−498623,−594104) | Y³+3249959Y²+2875039Y+1092727 |

The two nonzero formulas are immediate from

    2β/(2−β)=2β⁴,
    103β/(103−56β)=−103β¹⁷.

The exceptional gain is the reciprocal of the continued-fraction convergent
103/56 that precedes305. Thus the very same affine unit identity explains
the large CF digit, the generating-function denominator collapse, and this
exact integral feedback resolvent.

## 6. Nonperturbative and Pisot boundaries, with exact inequalities

The cubic has one positive real root β and a nonreal conjugate pair z,z̄.
Its discriminant is−44 and normβ=1, so |z|=β^−1/2. Exact sign checks give
9/5<β<2.

- γ=0 gives the positive Pisot number β.
- γ=1/2 gives the positive Pisot number2β⁴. Its other conjugates have modulus
  2β^−2<50/81<1.
- γ=56/103 gives a **negative** expanding real eigenvalue −103β¹⁷, not a
  positive Pisot number under the standard definition. Its sign reversal
  103β¹⁷ is Pisot: the other conjugates have modulus103β^−17/2<1, because

      9^17=16677181699666569
          >103²·5^17=8094024658203125.

The last gain lies on the far side of the real pole: β^−16=56β−103>0 implies
56/103>1/β. Hence spectral radius(γC)=γβ>1, and the full Neumann series
sum_(n≥0)(γC)^n diverges, while the exact rational resolvent still exists.
The large negative expanding mode and small internal modes are different
spectral coordinates; internal contraction does not make the full map stable.

Neither transformed Pisot number alone supplies an unaltered substitution
word system, positive incidence realization of the negative operator, or a
Euclidean Rauzy tiling. In particular the norm magnitudes below are nonunit.

## 7. Same field, different orders, different image-lattice index

Let A'=Z[λ_γ] at a unit-norm gain. Then A'⊂A, but generally A'≠A. For the
conjugates α_i, put λ_i=bα_i/(b−aα_i). Their differences satisfy

    λ_i−λ_j=b²(α_i−α_j)/[(b−aα_i)(b−aα_j)].

Multiplying squared pairwise differences gives

    disc(λ_γ)/disc(α)
      =b^(2k(k−1))/F_k(b,a)^(2(k−1)).

At F=±1, both generators have degree k, so the order discriminant/index
formula proves

    [A:A']=b^(k(k−1)).                                      (6)

The index can also be read directly as the absolute determinant of the
power-basis transition matrix with columns1,λ_γ,...,λ_γ^(k−1). Neither
calculation needs A to be the maximal order. In the cubic cases, the indices
are1,64,103^6=1194052296529.

Separately, multiplication by λ_γ has image-lattice index

    [A:λ_γ A]=|N(λ_γ)|=b^k,                                 (7)

which is **not** (6). Thus bI−aC is unimodular, but C_γ and R_γ have determinant
of magnitude b^k: they are integral matrices, not integral automorphisms for
b>1. The new generator is a nonunit. The original integral bilateral unit
recurrence cannot be inherited unchanged.

Because the order quotient A/A' has pure b-power order,

    A'[1/b]=A[1/b].                                         (8)

Only primes dividing b can support this order difference. In this localization
λ_γ is a unit as well. This is an arithmetic interface for finite-place data,
not by itself a tiling or attractor theorem. At the cubic nonzero gains the
relevant rational primes are respectively2 and103.

## 8. Exact matrices and finite controls

For γ=1/2, the power-basis matrices are

    Rγ=[[2,2,4],[2,4,6],[2,4,8]],
    Cγ=[[2,4,8],[4,6,12],[4,8,14]].

For γ=56/103 they are

    Rγ=[[-175615,-323008,-594104],
        [-271096,-498623,-917112],
        [-323008,-594104,-1092727]],
    Cγ=[[-323008,-594104,-1092727],
        [-498623,-917112,-1686831],
        [-594104,-1092727,-2009839]].

The standard-library checker verifies the exact resolvent equations, norm and
constant-coefficient identities, all integrality equivalences, Möbius inverse,
minimal-polynomial annihilation plus full degree, order-basis determinant,
closed powers and the exact Pisot inequalities. Its finite universe is all
primitive a/b with −6≤a≤6 and1≤b≤6 at k=3,4,5, plus56/103 at k=3. It includes
both negative/zero gains and nonintegral hostiles, so the test is not filtered
only to successful units.

Run normal and `-O` modes; the same2,436 checks pass. The two replay files are
identical. The global statements rely on the proofs and named cited input,
not on expansion of this finite box.
