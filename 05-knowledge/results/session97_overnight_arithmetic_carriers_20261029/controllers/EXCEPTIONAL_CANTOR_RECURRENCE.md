# Exact exceptional-address dimension, cell recurrence and integer sparsity

**Status: PROVED.** Scope: the infinite language of expanding cheaper-ternary first returns. This is **not** the set of all Collatz-divergent integers. Its intersection with the positive integers remains unproved empty.

## 1. The exact inverse system and normalized metric

Let X=1+2Z_2 be the compact odd2-adics, with normalized metric

    d_X(x,y)=2|x−y|_2.

Thus X has diameter1, and a residue class modulo2^(L+1) is a level-L ball of diameter2^(−L). Normalized odd Haar measure gives that ball mass2^(−L).

For each of the five failure words w_i, retain its exact chart `(A_i x+C_i)/2^(D_i)`, with A_i=3^(J_i). Its inverse is

    g_i(y)=(2^(D_i)y−C_i)/A_i.

Because A_i is odd, g_i is an exact similarity of ratio2^(−D_i). Its image is the complete guarded source cylinder of the word. These five images are pairwise disjoint. Their division clocks are

    D=(1,6,11,11,11).

Let K be the unique compact attractor `K=union_i g_i(K)`. Equivalently K consists of the sources realizing an infinite sequence of the five failure words. Each finite prefix remains a genuine native valuation word with its exact affine carry, so no pole or illegal guard is inserted by the construction.

For a positive odd integer n, membership in K means that its successive cheaper-ternary first returns **all strictly increase**. Universal first-tree termination is already proved; it is the unending sequence of expanding retries that is at issue here. The all-period integer certificate in the preceding note removes every eventually periodic positive itinerary, but not an aperiodic one.

## 2. Exact Hausdorff dimension

There is a unique s>0 such that

    2^(−s)+2^(−6s)+3·2^(−11s)=1.                 (1)

Then

    dim_H(K)=s.

Here is a self-contained proof. At block length k, the five-map cylinders cover K, and the sum of their t-th powers of diameters equals

    [2^(−t)+2^(−6t)+3·2^(−11t)]^k.

For t>s this tends to0, giving the upper bound. For the lower bound assign Bernoulli address weights p_i=2^(−sD_i). Separation makes the induced probability measure on K unambiguous, and a block cylinder with cumulative division clock D has mass2^(−sD).

Take any level-L ball meeting K, and follow an address in it until its cumulative division clock first reaches or exceeds L. The previous cylinder has depth D_prev<L, while D_prev≥L−11. That previous cylinder contains the ball: two intersecting ultrametric balls are nested. Hence the ball's K-mass is at most

    2^(−sD_prev) ≤2^(11s)2^(−sL).

For L=0 the same bound is immediate. Any cover therefore has sum of s-th powers of diameters at least2^(−11s), by summing its masses. This proves the lower bound and positivity of s-dimensional Hausdorff measure. Together the bounds give the dimension assertion. No unverified appeal to ordinary real digit concatenation is involved.

Put lambda=2^s. Its exact algebraic equation is

    lambda^11−lambda^10−lambda^5−3=0,

with the unique root greater than1. Exact rational checks give

    1.357705 <lambda<1.357706,
    441/1000 <s<442/1000.

Numerically, s≈0.441170998809. The numerical display is only orientation; the rational brackets are checked by integer/Fraction arithmetic.

## 3. Occupied dyadic cells: an exact all-level recurrence

Let N_L be the number of level-L odd dyadic balls intersecting K. Because an odd-unit affine map permutes dyadic balls at every precision, for L≥11 each first-level component contributes exactly N_(L−D_i) cells, and separation prevents overlap. Thus

    N_L=N_(L−1)+N_(L−6)+3N_(L−11),   L≥11.       (2)

The initial values, including their necessary truncation collisions, are

    (N_0,...,N_10)=(1,2,3,4,6,9,13,18,24,31,40).

They can be obtained without enumerating addresses. The first image has depth1, contributing N_(L−1) for L≥1. For1≤L≤5, the other four images are single truncated cells, with distinct-cell counts respectively1,1,1,2,3; this gives N_1,...,N_5. At levels6 through10, the depth6 image contributes N_(L−6), and the three depth11 images are distinct from it and from each other, contributing3. Hence

    N_L=N_(L−1)+N_(L−6)+3,    6≤L≤10.

The corresponding rational generating function is

    sum_(L≥0) N_L z^L
      =[1+z+z²+z³+2z⁴+3(z⁵+z⁶+z⁷+z⁸+z⁹+z¹⁰)]
       /[1−z−z⁶−3z¹¹].                           (3)

This calculates all finite precision cell counts with O(L) recurrence steps, rather than5^k address enumeration.

**Why this is not naive binary concatenation.** Inverse transport multiplies by A_i^(−1), an odd2-adic unit, and shifts by the full carry. These transformations preserve dyadic balls and the similarity ratio but do not append a single fixed ordinary binary string independent of the tail. Equation(2) uses ball transport and separation; the low-level numerator in(3) records the actual overlaps after truncation.

## 4. Completed-address counts and clocks are a different recurrence

Let b_D count finite failure-word histories with total division clock D, including the empty history at D=0. Then

    b_0=1, b_D=0 for D<0,
    b_D=b_(D−1)+b_(D−6)+3b_(D−11),  D≥1,
    sum b_D z^D=1/[1−z−z⁶−3z¹¹].

These are not the N_L counts: b_D counts completed word boundaries, whereas N_L counts all occupied precision cells, including cuts inside a block. Retaining the odd-step clock gives

    sum b_(J,D) u^J z^D=1/[1−u z−u⁴ z⁶−3u⁷ z¹¹].

The coefficient3 summarizes three distinct carries. To retain them rather than merely their counts, append each chart separately with

    (J,D,C)→(J+J_i,D+D_i,3^(J_i)C+2^D C_i),

and reconstruct the exact source guard. The clock recurrence is not a claim that those three ports or carries may be identified in a later operation.

## 5. Certified O(X^s) bound for possible positive integer exceptions

The eleven initial values satisfy

    lambda^L ≤N_L≤(17/8)lambda^L,    0≤L≤10.

The checker verifies these inequalities against outward rational bounds on lambda. Induction using(1) and(2) proves them for every L≥0.

Let E=K intersect {positive integers}. For any real X≥1 choose L=floor(log_2 X). Each occupied residue class modulo2^(L+1) contains at most one integer in[1,X], because its modulus exceeds X. Therefore

    # {n∈E : n≤X} ≤N_L≤(17/8)lambda^L≤(17/8)X^s.  (4)

In particular hypothetical forever-expanding retry sources form a zero-density subset, with an explicit exponent less than0.442. Equation(4) neither proves that E is nonempty nor that E is empty. The lower bound for N_L concerns2-adic cells, **not** positive integer witnesses in K. Every such finite cell has positive integer representatives, but these representatives may leave the failure language at a later precision.

Nor does(4) count all divergent Collatz orbits: a divergent orbit could contain locally descending cheaper returns and never stay wholly within this failure alphabet. This is a theorem about a precisely specified exceptional-address set, not a global Collatz density statement.

## 6. Exact reproduction and stopping boundary

`check_exceptional_cantor.py` verifies:

- the rational lambda and s brackets by exact arithmetic;
- both initial growth inequalities;
- the recurrence through level1000;
- independent inverse-IFS cell sets through level24, including invariance at each precision;
- the rational generating function by coefficient convolution through level100.

Normal and optimized Python outputs agree in `cantor_checks.json`. No high-level cell count enumerates5^k addresses. This closes the requested dimension/recurrence/integer-sparsity extension; the remaining aperiodic arithmetic emptiness question is deliberately not pursued here.
