# Deficits outside the midpoint grammar force unbounded one-apex loss complexity

**Status: PROVED, with same-agent proof review and an independent finite audit of the signed-power lemma. This is a theorem about one-apex presentations, not arbitrary tournament feedback-arc number or universal spectrum completeness.**

Let C be the set of positive odd numbers of the form

    c=E(S)-1,
    E(S)=sum_(s<t in S) 2^max(t-s-1,1),

where S is a finite integer set of at least two elements symmetric about some midpoint. Translation does not change E. Packet e proves that9 is not in C and gives a finite test for membership of any fixed c.

## The theorem

Fix c>0 outside C and a positive integer r. Suppose a tournament with raw Hamiltonian-path count2^k-c admits a one-apex presentation in which the apex loses to at most r chain vertices. Then k is bounded in terms of r,c. In particular, for deficit9, one-apex realizations along any unbounded sequence of exponents require unboundedly many apex losses.

A coarse explicit bound is available. Set

    M=r+r(r-1)/2+1+popcount(c+1),
    B=M+floor(log_2(c+1))+4.

Then it is sufficient to take

    k <= (r+1)! 2^(r+1) B.                            (1)

The bound is intentionally loose. For fixed c it implies the necessary growth

    r = Omega(log k / log log k)

along unbounded one-apex realizations. It asserts neither that every exponent is realizable in that stratum nor that the displayed lower growth rate is optimal.

## 1. A bounded-span lemma for signed powers of two

**Lemma.** In a minimal vanishing sum of m signed powers of two, each coefficient being+1 or-1, the largest exponent minus the smallest is at most m-2. Repeated exponents are allowed; minimal means that no nonempty proper submultiset already sums to zero.

**Proof.** For m=2 the two exponents coincide. For m>2, the smallest exponent cannot occur with both signs, since those two terms would form a proper zero sum. All terms at that exponent have the same sign, and their multiplicity is even: divide by the smallest power and reduce modulo2. Combine two of them into one term at the next exponent. The resulting(m-1)-term sum remains minimal; a proper zero subsum would expand to a proper zero subsum of the original. Its largest exponent is unchanged, and its smallest exponent increases by at most1. Induction gives the bound m-2. The carry chain

    1+1+sum_(j=1)^(m-3)2^j-2^(m-2)=0

attains it for m>=3. ∎

`vanishing_power_clusters.py` audits all signed-multiset pairs in its declared universe (at most8 terms, exponents0..6) by comparing independently computed subset-sum sets. This detects minimality without using the carry proof. It also checks the sharp examples through16 terms.

## 2. The exact count equation is a finite union of rational polyhedra

Remove transitive source/sink padding from the one-apex presentation. This preserves H and cannot increase the number of losses. For a nontransitive core, the apex beats chain vertex1 and loses to its final vertex. Write the distinct losing positions as

    2<=t_1<...<t_s=K,          s<=r.

The known one-apex formula gives

    sum_i 2^(t_i-1)
      -sum_(i<j)2^max(t_j-t_i-1,1)
      -2^k+(c+1)=0.                                 (2)

Expand c+1 in binary. Equation(2) is a vanishing sum of at most M signed powers, with unit signed coefficients.

First select, for each pair i<j, the branch t_j-t_i<=2 or t_j-t_i>=3. Its exponent is respectively the constant1 or the affine form t_j-t_i-1. These branches are given by integer linear inequalities.

Partition the indexed terms of(2) into minimal vanishing submultisets. In each part choose a minimal exponent. By the lemma, every other exponent differs from it by an integer in[0,M-2]. There are finitely many partitions and difference assignments. Keep only assignments whose signed powers cancel after these fixed relative shifts.

For each kept assignment, impose those exponent-difference equations, the chosen pair-gap branches, t_1>=2, increasing positions, and k>=0. This produces a rational polyhedron in the s+1 variables(t_1,...,t_s,k). Every integer point of it satisfies(2), since each selected component cancels after factoring out its common power. Conversely every actual solution determines at least one of these finitely many polyhedra.

Every coefficient in these constraints has absolute value at most2. Every right-hand side has absolute value at most B: the relative shifts are at most M-2, exponent constants are -1,0,1 or binary positions at most floor(log_2(c+1)), and the ordering/gap bounds are at most3. Thus the deliberately larger stated B covers every constraint.

This step keeps the exponent positions, their integer clock and the adjacency exception. It does not replace the power equation with a numerical approximation.

## 3. An unbounded exponent would produce an excluded affine family

Suppose one of these polyhedra has an integer point and is unbounded in k. A rational polyhedron unbounded in a linear objective has a rational recession direction increasing that objective. Scale this direction to an integer vector and add its nonnegative integer multiples to the existing integer point. This gives an exact family

    t_i(a)=t_i(0)+a v_i,
    k(a)=k(0)+a v_k,       a=0,1,2,...,     v_k>0.

Nonnegativity and ordering constraints imply0<=v_1<=...<=v_s. Each equal-slope group is therefore a fixed affine apex cluster. The inequalities retain the branch choices and order for every a. Equation(2) holds identically along the family by the component factorization in Section2.

The affine-cluster rigidity theorem in [TOURNAMENT-PRICES.md](TOURNAMENT-PRICES.md) now applies. It says that a pure power of two minus a fixed positive deficit in such an affine family must be a symmetric midpoint fan. Its deficit must belong to C, contradicting the hypothesis. Consequently every integer-feasible polyhedron in Section2 has bounded k.

The standard polyhedral fact used here is about a *finite rational system*, not arbitrary closed sets: unboundedness of a linear objective supplies an improving recession ray. Rationality permits scaling to integer coordinates. Starting from an integer feasible point is essential; no integer ray is inferred from an empty integer fibre.

## 4. An explicit effective bound

All variables in these polyhedra are nonnegative, so they contain no line. A feasible finite linear program with bounded objective attains its optimum. Restrict to the optimal face for k and minimize the sum of the nonnegative coordinates there. This second minimum is attained, and its minimizing face is bounded because all coordinates are nonnegative and their sum is fixed. A vertex of this bounded face is a vertex of the original polyhedron. Thus k attains its maximum at a vertex.

At such a vertex, s+1 independent active constraints form an integer linear system. Its matrix entries have absolute value at most2, and its right-hand side entries have absolute value at most B. The nonzero determinant has absolute value at least1. By Cramer's rule and the permutation expansion, every coordinate is at most

    (s+1)! B 2^s <= (r+1)! B 2^(r+1).

This bounds k at every integer point in every integer-feasible polyhedron. The transitive case H=1 has2^k=c+1 and is also below(1). This proves the theorem.

## Interpretation and next test

This converts a failed affine construction into a limitation on **all bounded-loss one-apex witnesses** for an excluded deficit, even if their positions were originally chosen non-affinely. The finite signed-power decomposition forces any unbounded bounded-loss set to contain an affine ray, and that ray is forbidden by the already proved coefficient rigidity.

It does not bound the number of feedback vertices needed by general tournaments. In particular the two-feedback changing-word bank is untouched. The finite exact witnesses for2^k-9 at9<=k<=36 are preserved in `deficit_nine_queries.json`, with actual adjacency counts. Their varying words, rather than a fitted short fan, remain the next constructive object to study.

The next useful audit is to sharpen(1) by enumerating only feasible cancellation partitions for small r, retaining the branch inequalities and source order. The present enormous bound is an effective structural theorem, not a practical search cutoff or a proof of universal acceptance.
