# Simultaneous circle antiperiodicity and a measurable-plane obstruction

## Status and scope

This note independently proves a refinement of the circle argument described in §7 of the supplied `proof(1).md`. It does not assume the blueprint's global plane-coloring theorems, its algebraic transfer, its cyclotomic lemma, or any claimed formal verification. The plane argument below uses ordinary Lebesgue differentiation, elementary Fourier theory on the circle, continuity of integrals under translation, and the inverse function theorem.

The principal conclusions are:

1. A measurable sign function can change sign under every rotation in a family precisely when one integer Fourier frequency changes sign under all of them.
2. For the rotations π/3 and γ, this is possible precisely when γ/π is a reduced rational with odd numerator and odd denominator. Incompatible rational cases already have finite odd-cycle obstructions.
3. A measurable plane coloring cannot weakly avoid all equal-leg isosceles templates at a fixed leg length when those angles fail this common-frequency test.
4. Odd relations give quantitative circle-defect lower bounds. For rational rotations the exact optimization is a finite maximum-cut problem. For π/3 and an irrational angle, however, the infimum circle defect is zero despite impossibility of exact simultaneous alternation.

Compatibility of the circle equations is only the absence of this particular obstruction. It does not construct an avoiding plane coloring.

## 1. Definitions

Write T = R/(2πZ), with normalized Haar measure m. Let A = {α_j : j ∈ J} be a family of angles. A measurable function f:T→{−1,1} is **simultaneously antiperiodic** for A if

    f(t + α_j) = −f(t) for m-almost every t, for each j ∈ J.

The exceptional null set may depend on j. The theorem below even permits an arbitrary index set J. Finite families suffice for the applications.

Our Fourier convention is

    f̂(n) = ∫_T f(t)e^(−int) dm(t),       n ∈ Z.

For a given sign function define its antiperiodicity defect at α by

    D_α(f) = m{t : f(t + α) = f(t)}.

Thus D_α=0 means antiperiodicity almost everywhere, and D_−α=D_α.

## 2. Exact common-frequency criterion

**Theorem 1.** The following are equivalent:

(i) There is a measurable f:T→{−1,1} simultaneously antiperiodic for A.

(ii) There is a nonzero g∈L²(T) simultaneously satisfying g(t+α_j)=−g(t) almost everywhere.

(iii) There is an integer n such that

    nα_j ≡ π  (mod 2π)     for every j ∈ J.

If J is nonempty, such an n is necessarily nonzero. The same equivalence holds for the empty family, with no constraints on n.

**Proof.** (i) implies (ii), since a sign function has L² norm one. For (ii), Fourier completeness gives some n with ĝ(n)≠0. Translation by α_j multiplies the nth Fourier coefficient by e^(inα_j). The antiperiodicity identities therefore give

    (e^(inα_j)+1)ĝ(n)=0

for every j. Using the same nonzero coefficient for the whole family proves (iii); no intersection of uncountably many null sets is needed.

Conversely, given n satisfying (iii), take

    f(t)=sign(cos(nt)).

Assign either sign at its finitely many zeros if n≠0. Since cos(n(t+α_j))=−cos(nt), each identity holds away from a finite set. If the family is empty, any sign function works. ∎

More precisely, every simultaneous L² solution has Fourier support contained in

    Λ(A) = {n ∈ Z : e^(inα_j)=−1 for every j}.

Conversely, that support condition implies all the L² identities by completeness. The extra requirement of taking only the values ±1 is satisfied by the square-wave construction whenever Λ(A) is nonempty.

### Finite rational-angle classification

Suppose A is nonempty and finite, and write

    α_j/π = p_j/q_j,       gcd(|p_j|,q_j)=1, q_j>0.

**Corollary 2.** The circle equations are compatible if and only if every p_j is odd and all q_j have the same 2-adic valuation. In that case, setting Q=lcm_j q_j gives

    Λ(A)=Q(2Z+1).

**Proof.** The condition n p_j/q_j∈2Z+1 forces q_j|n. It then forces both p_j and n/q_j to be odd. Thus every p_j is odd and v₂(q_j)=v₂(n), independent of j. Conversely, if these conditions on the rational data hold, Q/q_j is odd for every j, so Qp_j/q_j is odd. Finally, every admissible n is a multiple of Q, and its multiplier must be odd. ∎

An irrational α_j/π immediately makes Λ(A) empty. For an infinite rational family, one also needs a finite common multiple of all denominators; the exact integer criterion in Theorem 1 remains the simplest statement.

## 3. The pair π/3 and γ

**Corollary 3.** A measurable sign function satisfies

    f(t+π/3)=−f(t),       f(t+γ)=−f(t)  almost everywhere

if and only if

    γ/π=p/q in lowest terms, with p and q both odd.

In the compatible case let N=lcm(3,q). Then Λ={N(2k+1):k∈Z}. Every simultaneous sign solution is equivalently a measurable sign function satisfying

    f(t+π/N)=−f(t) almost everywhere.

Thus all solutions, modulo null sets, can be obtained by choosing any measurable sign function h on [0,π/N), and setting

    f(s+kπ/N)=(-1)^k h(s),       0≤s<π/N, 0≤k<2N.

**Proof.** Antiperiodicity for π/3 says n=3 times an odd integer. If γ/π is irrational, nγ cannot be an odd multiple of π. If γ/π=p/q is rational, Corollary 2 applies because the denominator of 1/3 has 2-adic valuation zero. Both p and q must be odd, and the common Fourier support is N(2Z+1). This support is precisely the support allowed by antiperiodicity under π/N. The displayed interval construction gives all sign solutions to that latter equation. ∎

### Explicit rational odd closed walks

These rational exclusions do not need measurability or Fourier analysis at the circle level. Suppose γ/π=p/q is reduced.

- If p is even, q is odd and qγ≡0 modulo 2π. Following q successive γ-edges is an odd closed walk. Alternating signs cannot return consistently to the starting sign.
- If q is even, p is odd and qγ−3(π/3)=(p−1)π≡0 modulo 2π. Following q γ-edges and three negative π/3-edges is an odd closed walk of length q+3.

Negative shifts also reverse signs whenever the corresponding positive shift does. For almost-everywhere equations, the finitely many translated exceptional sets still have measure zero, so the same contradiction applies almost everywhere. An odd closed walk contains an odd cycle, although the displayed walk need not itself be a simple cycle.

Consequently, all rational angles failing the odd/odd test are excluded by finite parity, whereas an irrational angle requires a measurability or regularity obstruction. For example, if γ/π is irrational, the subgroup generated by π/3 and γ is abstractly C₆⊕Z and has a sign character sending both generators to −1. Choosing representatives of its cosets in T gives nonmeasurable exact solutions. Theorem 1 shows that none can be measurable.

## 4. Quantitative defects from odd relations

**Proposition 4.** Suppose integers k₁,…,k_s satisfy

    Σ_j k_j α_j ≡ 0  (mod 2π),       Σ_j k_j is odd.

Then every measurable sign function obeys

    Σ_j |k_j| D_(α_j)(f) ≥ 1.

**Proof.** Expand the relation as a closed walk whose steps are |k_j| copies of sign(k_j)α_j. Its length L=Σ|k_j| is odd because |k_j|≡k_j modulo 2. For any starting t, not all L edges of that closed walk can reverse the sign: the product of all successive sign ratios is 1, while L successful reversals would give −1. In fact, the number of failed reversals is odd and hence at least one. Sum the failure indicators along the walk, integrate in t, and use translation invariance and D_−α=D_α. ∎

The explicit rational walks in §3 yield

    q D_γ ≥ 1                         if p is even,
    q D_γ + 3 D_(π/3) ≥ 1            if q is even.

There may be shorter or stronger relations. These bounds are explicit certificates, not an assertion that the displayed walk always gives the optimal unrestricted bound.

For comparison, Parseval gives the exact identity

    D_α(f) = (1/4)||f(·+α)+f||²₂
           = Σ_n |f̂(n)|² cos²(nα/2).

Therefore, for weights w_j≥0,

    Σ_j w_j D_(α_j)(f) ≥ inf_(n∈Z) Σ_j w_j cos²(nα_j/2).

This is a valid spectral lower bound. The binary-sign restriction can make the true minimum larger.

## 5. Exact finite optimization for rational rotations

Suppose all α_j/(2π) are rational and let H be the subgroup they generate. It is cyclic of some order M, with

    H={2πk/M : k∈Z/MZ},       α_j=2πa_j/M modulo 2π.

Form the finite Cayley multigraph with vertices Z/MZ and one edge instance (k,j) from k to k+a_j. Retain all edge instances, including duplicates, opposite directions, and loops. For a vertex signing σ∈{−1,1}^M set

    d_j(σ) = (1/M) Σ_(k∈Z/MZ) 1{σ(k+a_j)=σ(k)}.

**Theorem 5.** The set of all defect vectors

    (D_(α₁)(f),…,D_(α_s)(f))

attained by measurable circle sign functions is exactly the convex hull of the finitely many vectors (d₁(σ),…,d_s(σ)). In particular, for w_j≥0,

    min_f Σ_j w_j D_(α_j)(f)
      = min_σ Σ_j w_j d_j(σ)
      = Σ_j w_j − MaxCut_w(G)/M,

where edge instance (k,j) has weight w_j. Here a cut counts every edge instance whose endpoints have different signs; edge directions are bookkeeping, not the one-direction convention sometimes used for directed cuts. The minimum is attained.

**Proof.** Write t=s+2πk/M uniquely with 0≤s<2π/M and k∈Z/MZ, apart from endpoint conventions. A measurable f determines a measurable family of finite signings

    σ_s(k)=f(s+2πk/M).

With s uniform on [0,2π/M), the vector of circle defects is the average of the finite vectors d(σ_s). This belongs to their convex hull. Conversely, partition [0,2π/M) into measurable intervals with the desired convex-combination proportions and use the selected finite signing on each interval. This realizes every convex combination. A linear objective minimizes at a finite signing, and its uncut edge weight is MΣ_jw_j minus its cut edge weight. ∎

This is an exact finite reduction, not a claim of an efficient closed-form algorithm for every collection of generators. Ordinary maximum-cut enumeration or integer optimization can evaluate individual finite instances.

When all weights are positive, a zero minimum occurs exactly when this graph is bipartite. For these rational Cayley graphs, that is also equivalent to Theorem 1's common-frequency criterion. If there is no common frequency, an odd cycle gives a positive defect certificate.

### A sharp closed form with π/3 enforced exactly

The following useful slice of the optimization has a simple formula.

**Proposition 6.** Let γ/π=p/q be reduced, q>0, and constrain f by D_(π/3)(f)=0. Then

    min D_γ(f) = 0                       if p and q are odd,
                 gcd(3,q)/q             otherwise.

**Proof.** Put g=gcd(3,q), L=q/g, and b=3p/g. Then

    Lγ=b(π/3).

If p and q are not both odd, their parities differ. Since g is odd, L+b is odd. Proposition 4 applied to this relation gives L D_γ+|b|D_(π/3)≥1 and hence D_γ≥1/L.

To attain the bound, work in H=〈π/3,γ〉. The subgroup A=〈π/3〉 has order six, and H/A has order L and is generated by the coset of γ. Every element of H has a unique representation

    kγ+l(π/3),       0≤k<L, 0≤l<6.

Assign to that vertex the sign (-1)^(k+l). Every π/3-edge reverses the sign. A γ-edge also reverses the sign except possibly the edges from k=L−1 back to k=0. Their sign ratio is (-1)^(b−L+1), which is +1 exactly when L+b is odd. There are six such edge instances among 6L, so D_γ=1/L. Repeat this signing on the fundamental interval as in Theorem 5.

If p and q are both odd, L+b is even, all the wrap edges also reverse, and the same construction has D_γ=0. ∎

## 6. Irrational incompatibility does not give a uniform defect gap

**Proposition 7.** If γ/π is irrational, there is no exact measurable simultaneous solution for π/3 and γ, but

    inf{D_γ(f) : D_(π/3)(f)=0}=0.

Thus no inequality D_(π/3)(f)+D_γ(f)≥c with c>0 can hold for all measurable sign functions in this irrational case.

**Proof.** Let n range over positive integers of the form 3(2k+1), and put f_n(t)=sign(cos(nt)). Then D_(π/3)(f_n)=0. A direct comparison of two half-circle arcs gives

    D_γ(f_n)=dist(nγ, π+2πZ)/π,

where the distance lies in [0,π]. Since 6γ/(2π)=3γ/π is irrational, the residues

    nγ = 3γ+6kγ  (mod 2π)

are dense in T. Select a sequence tending to π modulo 2π; its defects tend to zero. Exact compatibility is still ruled out by Theorem 1. ∎

This sequence uses increasingly oscillatory measurable circle functions. It does not produce avoiding plane colorings, does not quantify the density-boundary argument below, and does not contradict qualitative plane exclusion.

## 7. An independently proved measurable-plane consequence

Identify the plane with C. A measurable coloring c:C→{−1,1} **weakly avoids** a template (0,a,b) if the set of monochromatic placements

    {(p,t)∈C×T : c(p)=c(p+e^(it)a)=c(p+e^(it)b)}

has zero measure for planar Lebesgue measure times m. Ordinary avoidance implies weak avoidance.

**Theorem 8 (fixed-leg circle obstruction).** Fix r>0 and a nonempty finite family of angles α₁,…,α_s. If there is no integer n satisfying

    nα_j≡π  (mod 2π)       for all j,

then no measurable two-coloring of the plane weakly avoids every template

    (0,r,re^(iα_j)),       j=1,…,s.

In particular, for any measurable coloring at least one listed template has a positive-measure set of monochromatic placements. If only nondegenerate triangles are wanted, take α_j∉{0,π} modulo 2π.

We prove the measure-theoretic ingredients, rather than importing them from the blueprint.

### Lemma 8.1: neither color is null

If one color class is null, then for any fixed template, almost every placement has all three vertices in the other class. Indeed, for each t the three exceptional sets of p are translates of a null set; integrate in t and restrict p to any finite positive-area window. Thus weak avoidance of even one template implies that both color classes have positive area.

### Lemma 8.2: a nonempty null density boundary

For a measurable color class X define its upper area density at p using the ratios |X∩B(p,h)|/h². Let

    S={p : neither color class has density zero at p}.

The two density-zero sets are disjoint, since their ball measures sum to πh². Outside S, one color has density one and the other density zero. The Lebesgue density theorem gives |S|=0. The set S is Borel: for each fixed h the ball integral is continuous in p, and density zero can be tested along h=1/n. The latter sequence is equivalent to all h↓0 by monotonicity and the ratio of consecutive radii tending to one.

If both colors have positive area, then S is nonempty. Here is a direct argument. Choose density points of the two colors. For some sufficiently small common side length, the axis-parallel squares centered at those two points have X-proportions respectively above and below 1/2. The X-measure of a translated fixed square varies continuously, so a square Q₀ between them has exactly half of each color.

Given a balanced square Q, partition it into four half-side squares. Their proportions average 1/2, so one has proportion at least 1/2 and one at most 1/2. The positions of all half-side squares contained in Q form a connected square of translation parameters. Continuity therefore supplies a balanced half-side square contained in Q. Repeating gives nested balanced closed squares Q_n, with side lengths ℓ_n→0 and intersection {p}. Each lies inside B(p,2ℓ_n), and each color has area ℓ_n²/2 in Q_n. Along these radii each color's ball-area ratio is at least 1/8. Thus p∈S.

### Lemma 8.3: apex exclusion at positive upper density

Suppose c weakly avoids (0,a,b), and X has positive upper density at p. There is no t₀ such that both

    x=p+e^(it₀)a,       y=p+e^(it₀)b

are density-one points of X.

To prove this, assume otherwise and choose ε>0 and arbitrarily small h>0 with |X∩B(p,h)|≥εh². Set L=max(1,|a|,|b|) and let I_h be the arc |t−t₀|<h/L. For q∈B(p,h), t∈I_h, the other two vertices lie respectively in B(x,2h) and B(y,2h). For each such t, the measure of bad q for which either of those vertices is outside X is at most

    |Xᶜ∩B(x,2h)|+|Xᶜ∩B(y,2h)|=o(h²).

This bound is uniform in t. For sufficiently small h from the selected sequence, at least (ε/2)h² m(I_h)>0 placements have q∈X and both other vertices in X. This contradicts weak avoidance.

### Lemma 8.4: a null set cannot support itself on a fixed-radius circle

Let S⊂C be Borel and null, and r>0. If

    m{t:p+re^(it)∈S}>0       for every p∈S,

then S is empty.

Suppose p∈S. The set E={t:p+re^(it)∈S} has positive measure. For each t∈E, the point q=p+re^(it) lies in S and sees S on its radius-r circle in positive arc measure. Therefore Fubini gives positive two-dimensional angular measure to

    {(t,u): p+r(e^(it)+e^(iu))∈S}.

But the smooth map F(t,u)=p+r(e^(it)+e^(iu)) has Jacobian determinant r² sin(u−t). Its critical set is the union of two one-dimensional diagonal subtori and is null. Away from that set the inverse function theorem gives local diffeomorphisms. A countable collection of such neighborhoods covers the regular locus, and the inverse maps are locally Lipschitz after restricting to compact subneighborhoods. Each therefore pulls the planar null set S back to a null set. Adding the null critical set shows F⁻¹(S) is null, a contradiction.

### Proof of Theorem 8

Suppose the coloring weakly avoids every listed template. By Lemmas 8.1–8.2, its density boundary S is null and nonempty.

Fix p∈S. Assume, for contradiction, that its radius-r circle meets S in null angular measure. On almost every direction t, p+re^(it) is then a density-one point of a unique color. Let f_p(t) be this density color, with arbitrary signs assigned on the exceptional null set. It is measurable because the density sets are Borel.

At p, both colors have positive upper density. For each j, Lemma 8.3 rules out a same-color pair of density points at p+re^(it) and p+re^(i(t+α_j)). Away from the boundary directions and their α_j-translates, the two density colors must therefore differ. Hence

    f_p(t+α_j)=−f_p(t) almost everywhere, for every j.

Theorem 1 contradicts the assumed absence of a common integer n. Thus every p∈S sees S in positive angular measure on its radius-r circle. Lemma 8.4 says S is empty, contrary to its construction. ∎

### Explicit extension of the blueprint's irrational-angle core

Taking α₁=π/3 and α₂=γ gives:

**Corollary 9.** A measurable plane coloring cannot weakly avoid both an equilateral template of side r and the isosceles template (0,r,re^(iγ)) whenever either

- γ/π is irrational; or
- γ/π=p/q is reduced and at least one of p,q is even.

For instance, γ=π/2, 2π/3, or π/4 are rational-angle exclusions. The test does not exclude γ=π/5, π/7, or 3π/5, whose reduced numerator and denominator are both odd.

The slightly stronger one-template consequence is also worth isolating: if α/π is irrational, or is a reduced rational with even numerator, then no measurable coloring weakly avoids (0,r,re^(iα)). A rational angle with odd numerator and even denominator is not excluded by the one-template circle test; adding π/3 supplies the incompatible second rotation.

These are sufficient nonavoidance criteria, not a classification of all plane-coloring possibilities. In particular, this note makes no conclusion in the compatible odd/odd cases and no extension to nonmeasurable plane colorings.

## 8. Integration notes

- Theorem 1 and Corollary 3 can replace the irrational-only conclusion of the supplied `circle_alternation` statement, once its null-intersection hypotheses have produced the two antiperiodicity identities.
- With the blueprint's original asymmetric hypotheses, R∩(R−π/3), Rᶜ∩(Rᶜ−π/3), and R∩(R−γ) null, the first two give m(R)=1/2 and R−π/3=Rᶜ modulo null sets. Equal measures then upgrade the third condition to R−γ=Rᶜ modulo null sets. Thus the exact conclusion is that γ/π has reduced odd numerator and odd denominator, not merely that it is rational.
- The rational obstruction has a finite parity proof. The measurable Fourier theorem also handles irrational angles and gives a uniform formulation for many rotations.
- The quantitative statements here concern circle defects. The plane density proof is qualitative; no positive numerical lower bound for plane monochromatic-placement density is asserted.
- No claim here verifies the rest of the supplied blueprint or any Lean declarations.
