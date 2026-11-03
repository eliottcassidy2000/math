# Fixed-rank cut moments and universal reciprocity to first order

Status: elementary all-order identities and small independent checks, not an induction or a characterization of tournament responses. The construction builds on the Session 152 integral carrier and marked-break identity. No literature-priority claim is made.

Let T be an n-vertex tournament, n>=2, with H directed Hamiltonian paths. Add x with x→i iff i∈S, and let F(S) be the new path count. Let c be any nonnegative integral carrier from `cut_induction`, so d(S)=(F(S)−H)/2=Σ_{i∈S,j∉S}c_ij. Put k_ij=c_ij+c_ji, W=Σ_{i<j}k_ij, deg_i=Σ_{j≠i}k_ij, and b_i=Σ_j(c_ij−c_ji). Thus Σ_i b_i=0. Define P_r(q)=Σ_{|S|=r}q^{d(S)}.

## Fixed-rank mean and first-order reciprocity

For 1≤r≤n−1,

  Σ_{|S|=r}d(S) = binom(n−2,r−1) W,
  mean_{|S|=r} F(S) = H+2W r(n−r)/(n(n−1)).

Indeed each directed edge crosses a rank-r cut in its forward direction for exactly binom(n−2,r−1) subsets. Consequently P_r(1)=P_{n−r}(1) and P'_r(1)=P'_{n−r}(1), so

  (q−1)² divides P_r(q)−P_{n−r}(q).

This holds for every tournament, even when it is not self-converse. Self-converse parents have the stronger exact equality P_r=P_{n−r}. For arbitrary parents rank distributions can differ: n=4 code 2 has H=3, W=9 and P_1=3q+q^6, P_3=1+3q^3. Thus the squared divisibility generally cannot be improved to equality or cubic divisibility. The complementary polynomials here use the half-increment d, not F itself; replacing q by z² and multiplying by z^H gives the corresponding ear enumerators.

Using W=((n−1)H+B)/2 from `cut_realizability`, where B counts permutations with exactly one backward adjacency, this provides the explicit average H+((n−1)H+B)r(n−r)/(n(n−1)). Since d is integral, at least one child has

  F(S) ≥ H+2 ceil(W floor(n²/4)/(n(n−1))).

This is an attainment bound, not interval coverage; it does not identify which child has the bound or its subsequent response state.

## Exact second-moment asymmetry

For n≥3 let D=Σ_i b_i deg_i. For 1≤r≤n−1,

  Σ_{|S|=r} d(S)² − Σ_{|S|=n−r} d(S)²
   = [(n−2r)/(n−2)] binom(n−2,r−1) D.                (1)

Proof. Complementation rewrites the left side as Σ_{|S|=r}(d(S)−d(V−S))(d(S)+d(V−S)). The first factor is Σ_{i∈S}b_i, and the second is Σ_{ij crossing} k_ij. For each edge ij, the sum of b_v over rank-r sets in which ij crosses equals (b_i+b_j)[binom(n−2,r−1)−2 binom(n−3,r−2)], with out-of-range binomials interpreted as zero. Sum over edges and use Σ_{ij}k_ij(b_i+b_j)=D. Simplifying the coefficient gives (1). Since the first moments agree, the left side also equals P''_r(1)−P''_{n−r}(1).

In particular for r≠n/2, cubic divisibility by q−1 is equivalent to D=0, and D is determined by the response despite carrier nonuniqueness. For self-converse T, D=0; the converse is not asserted. Code 2 has D=12 and rank 1 versus 3 squared-sum difference 12, showing sharpness. These restrictions apply to every integral carrier of a parent, not to arbitrary lists of child values.

## All moments from marked crossing arcs

For any nonnegative integer m and 0≤r≤n, there is a finite exact formula

  Σ_{|S|=r} d(S)^m
  = Σ_{(e_1,...,e_m)} (∏_{t=1}^m c_{e_t}) binom(n−a−b,r−a),

where the sum runs over ordered m-tuples of directed arcs (repetitions allowed). Let A be the set of their tails and B the set of their heads, a=|A|, b=|B|. Terms with A∩B≠∅ are zero, and otherwise an out-of-range binomial is zero. For m=0 take the unique empty tuple, yielding binom(n,r). Proof: expand d(S)^m and count rank-r subsets containing all tails and no heads. These moments determine each finite polynomial P_r, but computing them can be expensive and the carrier does not give a closed next-parent induction.

## Verification

`python3 check.py` independently constructs literal ear tournaments and counts paths by subset DP for every labeled tournament at n=3,4 and both session 152 order-5 codes. It derives k, b, W from the resulting F values, verifies the rank means and (1) at every rank, and prints the nonzero-D n=4 witness. It additionally checks the all-moment formula for the published code-8 carrier at every rank and moments 0 through 4, including repeated and conflicting arcs. The attainment bound follows from the proof; the test is not an exhaustive high-order theorem check. `out.txt` records the run.
