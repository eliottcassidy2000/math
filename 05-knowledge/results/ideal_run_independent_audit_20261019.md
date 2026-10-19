# Independent audit: ideal-run moments and the rooted-five neighborhood

Date: 2026-10-19.

## Verdict and scope

The ideal-run realizability inequality is valid. It eliminates the proposed one-sided saturated survivor without using seven-thinness, the cardinalities of the exceptional sets, or the nonchain condition. Applying it independently at the two exteriors also proves the stronger statement:

> If a finite poset has incomparable roots x,y whose closed incomparability neighborhood is the induced five-element poset with covers b<z<a, x<a, and b<y, then at least one of x:y, x:z, z:y is a 1/3–2/3 balanced pair in the entire poset.

There is no degree hypothesis in this statement. It addresses this specified rooted neighborhood, not arbitrary posets and not the full seven-thin conjecture.

This audit checks the new proof independently using the existing exact core-count formulas, derives the new inequalities explicitly, and verifies the formulas against direct linear-extension counting for actual small exterior posets. It does not infer a universal theorem from finite enumeration.

## 1. Primary-source verification

Richard P. Stanley, “Two combinatorial applications of the Aleksandrov–Fenchel inequalities,” Journal of Combinatorial Theory, Series A 31 (1981), 56–65. DOI: https://doi.org/10.1016/0097-3165(81)90053-4.

Author-hosted original paper: https://math.mit.edu/~rstan/pubs/pubfiles/48.pdf.

Theorem 3.1 on printed pages 62–63 (PDF pages 7–8) counts linear extensions with fixed ranks of elements on a chain and proves log-concavity in one of the ranks. Its single-element specialization is stated immediately after the theorem: if N_i counts extensions placing a fixed element at position i, then N_i²≥N_(i−1)N_(i+1). The original paper and this specialization were directly inspected. No unverified strengthening is needed here.

## 2. Exact rank-count bijection for an ideal

Let B be a finite poset and K an ideal. In a uniformly chosen internal linear extension of B, let R be the length of the initial run in K. Put

T_j = #{extensions whose first j elements are in K},
q_j = T_j/e(B) = Pr(R≥j),
q_0 = 1,
q_j = 0 for j>|K|.

Adjoin a new element t, with t below every point of B\K and incomparable with every point of K. This is a valid poset extension: K being an ideal rules out any relation from B\K down into K, so taking transitive closure does not produce t<K. The point t is minimal.

Deleting t from an enlarged-poset extension in which t has rank j+1 gives an extension of B whose first j elements are all in K. Conversely, insert t immediately after those first j elements. Every point required to follow t is indeed later, and t has no required predecessor. These operations are inverse. Therefore the rank-j+1 count for t is exactly T_j, with no factorial, binomial, or additional extension-count factor.

Stanley's theorem now gives q_j²≥q_(j−1)q_(j+1). The sequence is positive for 0≤j≤|K|, because an ideal can be placed first in a linear extension. Its consecutive ratios decrease. Thus for all i,j≥0,

q_(i+j)≤q_i q_j.

For i+j≤|K|, compare the j ratios starting at i+1 with those starting at 1. For i+j>|K| the inequality is immediate.

## 3. The moment bound

Write m=E R=Σ_(j≥1)q_j. Since all sums are finite,

m² = Σ_(i,j≥1)q_i q_j
    ≥ Σ_(i,j≥1)q_(i+j)
    = Σ_(k≥2)(k−1)q_k
    = E binom(R,2).

Consequently

E[(R+1)(R+2)/2]
=1+2m+E binom(R,2)
≤(m+1)².

Now let H⊆K be ideals and let S be the initial H-run in the same extension. Pointwise 0≤S≤R and

(S+1)(2R−S+2)/2
= (R+1)(R+2)/2 − (R−S)(R−S+1)/2.

Define d=1+E R and v=E[(S+1)(2R−S+2)/2]. The subtracted term is nonnegative, and the original integrand is at least R+1. Hence

d≤v≤d².                                                   (A)

Order duality gives the identical inequality c≤u≤c² for terminal runs in nested upper ideals of a lower exterior.

The auxiliary augmented poset need not obey the original degree cap; Stanley's theorem has no degree cap. This causes no issue.

## 4. One-sided consequence with an explicit witness rule

Suppose the rooted-five poset has no lower exterior, or its lower exterior is an ordinal prefix. In the upper exterior B, let K=I(a)∩B and H=I(z)∩B. These are ideals with H⊆K. The actual, uniform-internal-extension average gives

E = 5d+3v+3,
Nxy = 3d+2v+3,
Nzy = 4d+3,

up to the same positive factor e(B), which cancels from all probabilities.

If d≥3, then

2Nxy−E=d+v+3>0,
3Nxy−2E=3−d≤0,

so x:y is balanced.

If 1≤d<3, then (A) gives

3Nzy−E=7d+6−3v≥7d+6−3d²
=−(3d+2)(d−3)>0,

2E−3Nzy=6v−2d−3≥4d−3>0.

Thus z:y is balanced. In particular the saturated one-sided survivor is impossible. None of its specialized degree or nonchain properties is needed.

## 5. Both-sided extension and the exact polynomial audit

For a general rooted-five neighborhood, every lower exterior element is below x,y and every upper exterior element is above x,y, so all lower exterior elements are below all upper ones. Fix internal linear extensions of the two exteriors independently. In each fiber, the existing core formulas depend on lower parameters c,u and upper parameters d,v. The fiber counts are bilinear across the two sets of parameters; summing over all pairs of internal words therefore substitutes their independent means and introduces the common factor e(A)e(B).

This is not an assertion that the exterior words have independent uniform marginals under a uniform extension of the whole poset. It is an exact partition of extension counts into fibers.

The mean parameters satisfy c≤u≤c² and d≤v≤d² by Section 3. The normalized count polynomials are

E=2cd+c+2cv+2ud+u+d+1+v,
Nxy=cd+c+cv+ud+u+d+1+v,
Nxz=cd+c+2cv+d+1+v,
Nzy=cd+c+2ud+u+d+1.

Put C=c−1, D=d−1, U=u−1, V=v−1. In particular C,D≥0,

C≤U≤C²+2C,
D≤V≤D²+2D.

Direct expansion, independently checked, gives

2Nxy−E=c+u+d+1+v>0,                                  (B)
3Nxy−2E=2−C−D−CD−CV−DU,                              (C)
E−3Nxz=−CD−4CV−7C+2DU−D+3U−6V−10,                   (D)
14E−11(Nxz+Nzy)=6CD+6CV+6DU+4C+4D+9U+9V≥0.           (E)

If x:y is unbalanced, (B) implies Nxy/E>2/3. Equation (C) then forces 0≤C,D<2. In (D), the coefficient of U is 2D+3>0 and the coefficient of V is −4C−6<0. Substituting U≤C²+2C and V≥D therefore gives

E−3Nxz
≤(2D+3)C²−(D+1)C−7D−10
=(3C+5)(C−2)+D(2C²−C−7)<0.

The final sign is rigorous: the first summand is strictly negative for 0≤C<2; the quadratic 2C²−C−7 is at most −1 on [0,2], as follows from convexity and its endpoint values −7,−1. Thus p(x<z)>1/3. Duality gives p(z<y)>1/3.

Finally (E) yields p(x<z)+p(z<y)≤14/11<4/3, so these two probabilities cannot both exceed 2/3. At least one of x:z or z:y is balanced. This proves the stated arbitrary-degree rooted-neighborhood theorem.

## 6. Why the earlier arbitrary mixture was not realizable

The raw count-measure mixture 5μ_Q+2μ_T(0,0,1,5) would require the upper profile law to have R=0 with probability 5/7 and R=6 with probability 2/7, with the latter profile having S=5. Thus

d=19/7,
v=59/7,
v−d²=52/49>0.

It violates (A). Already its initial-run tails q_1=q_2=2/7 violate q_1²≥q_0q_2. The mixture remains a valid abstract positive mixture, but cannot be the actual exterior law needed for this rooted poset representation.

## 7. Reproduction

Run `python verify_ideal_run_audit_20261019.py` in this directory. It is standard-library-only and imports none of the earlier discovery or verification scripts. Its receipt is `ideal_run_independent_audit_20261019.json`.

The verifier enumerates all naturally labeled exterior posets B of size at most five and all ideal pairs H⊆K. It independently constructs the augmented marker poset, counts its rank sequence, compares it with the initial-run tail counts, checks log-concavity and submultiplicativity, and checks the moment inequality in exact rational arithmetic. It also constructs the actual rooted-five poset and directly counts its total and marked linear extensions to check the three one-sided formulas and the explicit witness rule.

A separate two-sided pass uses every such collar of size at most three on each side, dualizing the lower collar. It directly constructs and counts each full poset, checks the product-averaged count formulas, the four polynomial identities, the moment constraints, and the three-pair witness. Some tested posets have z of incomparability degree eight, so the regression is not limited to seven-thin instances.

All files in this audit are new. Earlier packets and repository branches were not modified.
