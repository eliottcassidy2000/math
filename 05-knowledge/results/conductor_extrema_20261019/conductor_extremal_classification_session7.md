# Exact endpoint classification for uniform rank-three conductor defects

Date: 2026-10-19. Status: **PROVED** in the stated scope, with independent **FINITE-EXACT** controls.

This is an addendum to Sections 2–3 and the conductor question in Section 8 of `unprompted_bridges_session7.md`. It determines both sharp endpoints, gives rational configurations attaining them for every branch count, and characterizes equality at the upper endpoint geometrically. It does **not** determine which intermediate integers occur. The arguments use classical Hilbert-function, interpolation, and generalized-Vandermonde ingredients; no literature-wide novelty claim is made.

## 1. Main theorem

Let k be an algebraically closed field of characteristic zero. Let P={p₁,…,pᵣ} be r≥3 distinct projective points in P²(k), with no three collinear. Equivalently, their nonzero vector representatives realize the linear matroid U(3,r). Let X_P be the reduced union of their r affine lines through the origin in A³(k), let A_P be its coordinate ring, let B_P be its normalization, and let S_P be its seminormalization.

Write

- δ_sn(P)=length_k(S_P/A_P), the seminormal defect;
- N_q=binom(q+2,2), the dimension of ternary forms of degree q;
- s=s(r)=min{q≥0:N_q≥r};
- c=c(r)=ceil((r−1)/2)=floor(r/2).

Thus N_(s−1)<r≤N_s and s≥1; explicitly, s=ceil((√(8r+1)−3)/2). Then the **exact sharp bounds** are

δ_min(r) ≤ δ_sn(P) ≤ δ_max(r),                                      (1.1)

where

δ_min(r) = (s−1)r − binom(s+2,3) + 1,
δ_max(r) = floor((r−2)²/4).                                         (1.2)

Both endpoints are attained by configurations defined over Q:

- **Minimum:** Gᵣ={(1,t,tʳ):t=1,…,r}.
- **Maximum:** Cᵣ={(1,t,t²):t=1,…,r}.

The complete Hilbert-rank and branch-conductor formulas for these witnesses are

h_q(Gᵣ)=min(r,N_q),                  d_i(Gᵣ)=s for every i;
h_q(Cᵣ)=min(r,2q+1),                d_i(Cᵣ)=c for every i.           (1.3)

For r≥6, upper equality in (1.1) holds **if and only if all r points lie on a smooth conic**. For every such maximum configuration, all branch-conductor exponents equal c.

Lower equality holds if and only if the points impose independent conditions in every degree, meaning h_q=min(r,N_q) for all q. Equivalently, it suffices to check the two adjacent degrees

h_(s−1)=N_(s−1) and h_s=r.                                          (1.4)

This lower equality condition does not, in general, force all branch-conductor exponents to equal s. The displayed equality d_i=s is a property of the explicit minimum witness Gᵣ; Section 6 gives a minimum-defect counterexample to the stronger assertion.

The normalization defect is obtained by adding r−1 to each endpoint. It must not be conflated with δ_sn.

## 2. Graded evaluations and conductor conventions

Choose representatives v_i of p_i. For q≥0, let

E_q = image of [F ↦ (F(v₁),…,F(vᵣ))],
h_q = dim_k E_q,

where F ranges over homogeneous ternary forms of degree q. Changing projective coordinates or rescaling individual representatives preserves h_q and all separator degrees.

Choose a linear form ℓ nonzero on every point and rescale so that ℓ(v_i)=1. Such an ℓ exists because k is infinite. Then E_q⊆E_(q+1). Let

d_i=min{q:e_i∈E_q},

where e_i is the i-th coordinate vector in kʳ. Thus d_i is the smallest degree of a form that vanishes on every other point and not on p_i.

The predecessor's all-line-arrangement formula reads

δ_sn(P)=sum_(q≥1)(r−h_q),
length_k(B_P/A_P)=(r−1)+δ_sn(P),
(A_P:B_P)=direct sum_i t_i^(d_i) k[t_i].                            (2.1)

These identities can also be read directly from the grading. The normalization is the product of the branch rings k[t_i]; its degree-q part is kʳ. The original coordinate ring has degree-q part E_q. The seminormalization is the ring of branch functions with a common value at zero, so its degree-zero part is the diagonal k and its positive-degree parts are kʳ. A homogeneous function supported on branch i first occurs in degree d_i, and multiplication by ℓ supplies it in every subsequent degree. This proves the conductor formula. The same lengths and exponents hold after completion at the common origin.

This scope is essential: these are cones over projective point configurations, or equivalently unions of straight lines. For arbitrary curved branches, higher jets are additional data.

## 3. The universal Hilbert-rank sandwich

For every q≥0,

min(r,2q+1) ≤ h_q ≤ min(r,N_q).                                    (3.1)

The right inequality is the dimension bound for the evaluation map.

For the left inequality, set m=min(r,2q+1) and choose any m points. Fix one of them, p. Partition the other m−1 points into at most q pairs or singletons. For each pair use its joining line; that line avoids p because no three points are collinear. For a singleton choose a line through it avoiding p. Multiply the line equations, padding to degree q by powers of ℓ if necessary. The resulting form vanishes at the other m−1 selected points and is nonzero at p.

Doing this for each selected point shows that degree-q evaluation onto those m coordinates is surjective. Therefore h_q≥m. At q=0 the statement is simply h₀=1.

In particular, h_c=r. Thus every sum below is finite and every branch has d_i≤c.

### 3.1 Summing the lower bound for the defect

By the upper Hilbert bound,

δ_sn(P) ≥ sum_(q≥1) max(0,r−N_q)
          = sum_(q=1)^(s−1)(r−N_q).

The hockey-stick identity gives

sum_(q=1)^(s−1) N_q = binom(s+2,3)−1,

which proves the formula for δ_min(r).

More precisely,

δ_sn(P)−δ_min(r) = sum_(q≥1)[min(r,N_q)−h_q].                        (3.2)

Every summand is nonnegative. Equality therefore holds precisely when all the ranks are maximal.

The two-degree test (1.4) is equivalent. If evaluation is injective at degree s−1, it is injective at every smaller degree: a lower-degree form vanishing on all points could be multiplied by a power of ℓ to contradict injectivity at degree s−1. If evaluation is surjective at degree s, it stays surjective in higher degrees because E_s⊆E_(s+1)⊆⋯. The converse is immediate.

### 3.2 Summing the upper bound for the defect

By the lower Hilbert bound,

δ_sn(P) ≤ sum_(q≥1) max(0,r−(2q+1)).                               (3.3)

If r=2m, the nonzero terms have q=1,…,m−1 and sum to (m−1)². If r=2m+1, they again have q=1,…,m−1 and sum to m(m−1). These are exactly floor((r−2)²/4).

Again the equality information is termwise:

δ_max(r)−δ_sn(P) = sum_(q≥1)[h_q−min(r,2q+1)].                      (3.4)

## 4. Explicit rational attainment of the minimum

Consider Gᵣ={(1,t,tʳ):t=1,…,r}. These are rational points with distinct positive parameters.

### 4.1 A classical sparse-polynomial lemma

A nonzero real polynomial with n nonzero monomial terms has at most n−1 distinct positive roots. For a short proof, divide by its smallest monomial. The result has a nonzero constant term and n−1 other terms. If it had n distinct positive roots, Rolle's theorem would give its derivative at least n−1 positive roots. The derivative has n−1 nonzero terms, contradicting induction. The one-term case is immediate.

Consequently, for distinct nonnegative integers e₁,…,e_n and distinct positive real nodes u₁,…,u_n, the matrix (u_i^(e_j)) is nonsingular: a nontrivial column relation would give an n-term polynomial with n positive roots. Every square submatrix has the same property. This is the classical generalized-Vandermonde fact; a modern proof and attribution are given in [Khare, Section 5, Theorem 5.1 and Lemma 5.2](https://www.math.iisc.ac.in/~khare/teaching/Math341-notes.pdf).

### 4.2 Distinct exponents up to the interpolation threshold

The degree-q monomials x^(q−a−b)y^a z^b, with a,b≥0 and a+b≤q, restrict to t^(a+rb).

For q<r all these exponents are distinct. Indeed, if a+rb=a′+rb′, then r divides a−a′, while |a−a′|≤q<r. Hence a=a′ and then b=b′.

Thus for every q<r the degree-q evaluation matrix is a generalized Vandermonde matrix with N_q distinct columns, so

h_q(Gᵣ)=min(r,N_q).

The real argument proves that appropriate integer minors are nonzero over Q; they remain nonzero over every characteristic-zero extension. It does not require an ordering on k.

Since N_(r−1)≥r, one has s≤r−1. In particular the preceding calculation applies at q=s and gives h_s=r. Multiplication by x, which takes value 1 at every displayed point, gives h_q=r for every q≥s. This proves the entire minimum Hilbert function, including degrees q≥r where exponent collisions do occur.

The same generalized-Vandermonde lemma at degree one shows that every triple of rows is independent, since its exponents are 0,1,r. Thus Gᵣ realizes U(3,r). Equivalently, the affine graph z=tʳ is strictly convex for positive t, so no three of these points are collinear.

It follows from (2.1) that δ_sn(Gᵣ)=δ_min(r).

### 4.3 Exact branch-conductor exponents at the minimum witness

At q=s, evaluation is surjective, so every d_i≤s.

For q<s, N_q<r and hence N_q≤r−1. Deleting any selected row leaves r−1 positive distinct nodes. The remaining generalized-Vandermonde matrix still has full column rank N_q. Therefore a degree-q form vanishing at the other r−1 points must be zero, and it cannot separate the selected point.

Thus d_i(Gᵣ)=s for all i. In particular,

(A_G:B_G)=direct sum_(i=1)^r t_i^s k[t_i].                          (4.1)

## 5. Conic attainment and the complete upper equality classification

For Cᵣ={(1,t,t²):t=1,…,r}, degree-q monomials restrict to a spanning set for all polynomials of degree at most 2q. Ordinary Vandermonde interpolation therefore gives

h_q(Cᵣ)=min(r,2q+1).

The triple Vandermonde determinants are nonzero, so Cᵣ realizes U(3,r), and (3.3) is attained.

The same Hilbert function holds for any finite set of distinct points on a smooth conic. Over k the conic has the standard parametrization [u:v]↦[u²:uv:v²]. Degree-q plane forms restrict to every binary form of degree 2q. After choosing a parameter value outside the finite set as infinity, evaluation becomes ordinary one-variable interpolation.

A separator on a conic restricts to a degree-2q binary form with r−1 distinct zeros. It exists if and only if 2q≥r−1: necessity follows from the root bound, while sufficiency follows by taking the product of the r−1 linear factors for the other parameters and padding the degree. Consequently,

d_i=c=ceil((r−1)/2) for every branch of every conic configuration.    (5.1)

Now suppose r≥6 and δ_sn(P)=δ_max(r). Formula (3.4) forces h₂=5. The space of ternary quadrics has dimension six, so a nonzero quadratic vanishes at every point of P.

That conic cannot be reducible: a union of two projective lines, including a repeated line, contains at most four points of a configuration with no three collinear. Over an algebraically closed field a singular plane conic is reducible. Hence the conic is smooth. Conversely, every smooth-conic configuration has already been shown to attain the upper endpoint.

This proves the upper equality classification. In fact the conic is unique up to scalar, since h₂=5 makes the quadratic kernel one-dimensional. For a rational configuration its equation can be chosen over Q, by rational linear algebra.

For 3≤r≤5, the two bounds coincide: δ_sn=r−3 for every U(3,r) realization. At r=3 all d_i=1; at r=4,5 all d_i=2. The latter follows because the remaining points span k³, excluding linear separators, and degree-two separators exist by Section 3. At r=5 there is a unique smooth conic through the five points, since h₂=5 and a reducible conic can contain at most four. Thus r=6 is the first nontrivial upper-locus distinction.

## 6. What the endpoint theorem does and does not say about conductors

Let ρ(P)=max_i d_i, the first degree in which the full evaluation map is surjective. Then

s(r) ≤ ρ(P) ≤ c(r),                                                (6.1)

and both endpoints are rationally attained by Gᵣ and Cᵣ respectively. The lower bound follows because surjectivity requires N_q≥r; the upper bound follows from the pair-line separators. The equality ρ=min{q:h_q=r} uses E_q⊆E_(q+1).

Every minimum-defect configuration has ρ=s, by (1.4), but individual d_i can be smaller. Here is a concrete safeguard. Take

D₇={(1,t,t²):t=0,…,5} ∪ {(1,6,37)}.

No three points are collinear. Its Hilbert ranks are

1,3,6,7,7,…,

so δ_sn(D₇)=(7−3)+(7−6)=5=δ_min(7). Yet its branch-conductor exponents are

(3,3,3,3,3,3,2).

Indeed xz−y² separates the displaced point in degree two. A separator for any old conic point must vanish at the other five conic points, requiring degree at least three; the universal upper bound supplies degree three. There is no linear separator because the other points span k³.

In contrast, G₇ has the same minimal defect and the same Hilbert-rank sequence, but all seven conductor exponents equal three. Therefore even the complete Hilbert function does not determine the branch-conductor vector. Its ranks are less information than the actual subspaces E_q.

This observation refines, rather than weakens, the predecessor's message: the abstract tangent matroid loses the higher evaluation data; even the sequence of dimensions of those evaluation spaces loses the locations of individual coordinate separators.

## 7. Size of the variation and exact small values

Since s=√(2r)+O(1), the closed formula gives

δ_min(r)=(2√2/3)r^(3/2)+O(r),
δ_max(r)=r²/4+O(r).

Thus rational arrangements with the same abstract U(3,r) tangent matroid can have seminormal defects of different growth orders, r^(3/2) and r². The difference between the largest and smallest possible defects is exactly

floor((r−2)²/4) − [(s−1)r−binom(s+2,3)+1].                         (7.1)

This is stronger than the predecessor's conic-versus-one-displaced-point gap, which was only linear in r.

| r | s | c | δ_min | δ_max |
|---:|---:|---:|---:|---:|
| 3 | 1 | 1 | 0 | 0 |
| 4 | 2 | 2 | 1 | 1 |
| 5 | 2 | 2 | 2 | 2 |
| 6 | 2 | 3 | 3 | 4 |
| 7 | 3 | 3 | 5 | 6 |
| 8 | 3 | 4 | 7 | 9 |
| 9 | 3 | 4 | 9 | 12 |
| 10 | 3 | 5 | 11 | 16 |
| 11 | 4 | 5 | 14 | 20 |
| 12 | 4 | 6 | 17 | 25 |

The interval between these endpoint values is a containing interval, not an assertion that every integer in it occurs.

## 8. Independent finite-exact controls

The companion `verify_conductor_extrema.py` uses only the Python standard library. Its checks use explicit exceptions rather than assertions, so optimization cannot disable them. It uses neither floating-point arithmetic nor random choices.

The controls include:

1. Both closed endpoint formulas versus their defining sums for every r=3,…,1000.
2. Both explicit rational families for every r=3,…,18, all triple determinants, every Hilbert rank from q=0 through c+1, and all defect formulas. Every such rank is independently calculated by Fraction Gaussian elimination and primitive-integer fraction-free elimination.
3. Every branch's sharp conductor threshold in those families: all row-deleted ranks at degree d−1 and full rank at degree d, again by both algorithms.
4. Distinct power-curve exponents in every degree q<r, and the explicit collision boundary q=r. The proof uses exponent injectivity only below r; it does not silently assume injectivity in all degrees.
5. The seven-point minimum-defect example with unequal branch exponents, computed from first separator degrees.
6. An exhaustive finite universe: all subsets of size at least three from ten displayed conic points and two displayed off-conic rational points, retaining precisely the no-three-collinear subsets. For each retained configuration the program checks the rank sandwich, both endpoint bounds, the two-degree minimum equality test, the upper conic-rank equality test for r≥6, and the forced small-r defects.

Normal and Python-optimized replays returned PASS with byte-identical outputs: 998 branch-count formula checks, 32 explicit family instances, 224 double-checked full evaluation matrices, 336 double-checked row-deleted matrices, and all 4,017 retained subsets of the 12-point pool. Of those subsets, 2,510 have r≥6 and were checked against the quadratic-rank upper-equality criterion. The accompanying output records the detailed counts and data. These controls corroborate the all-r proofs; they do not establish the infinite claims by extrapolation.

## 9. Attribution, scope, and the remaining gate

The ordinary-singularity and graded-evaluation framework predates this note. The predecessor cites Leslie G. Roberts, *Ordinary Singularities with Decreasing Hilbert Function*, Canadian Journal of Mathematics 34 (1982), 169–180, [DOI 10.4153/CJM-1982-010-1](https://www.cambridge.org/core/journals/canadian-journal-of-mathematics/article/ordinary-singularities-with-decreasing-hilbert-function/E9D1AB1F50426FB4891D7CD245750555). That paper's references also explicitly include earlier work of Geramita and Orecchia on line arrangements, generic point configurations, and conductors. The sparse-positive-root and generalized-Vandermonde arguments are classical, as made explicit in Section 4. No ingredient here is presented as newly discovered general interpolation theory.

The scoped contribution is the assembled exact endpoint theorem for the repo's seminormal-defect question, its explicit rational all-r witnesses, the conic upper-equality classification, and the branchwise limitation at the lower endpoint.

**Closed gate:** exact minimum and maximum of δ_sn over rational realizations of U(3,r), with both values attained over Q; the same bounds hold over k.

**Still open in this artifact:** the set of intermediate attainable defects, geometric classification of the minimum-defect locus beyond condition (1.4), and the full attainable branch-conductor vectors. Nothing here classifies curved branches with these tangent matroids or supplies a consequence for the unrelated global conjectures named in the predecessor.
