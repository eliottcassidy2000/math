# A Turán separator for fully covered cyclic reversal cubes

**PROVED:** an exact named two-port gap formula for a cyclic triple of blocks with identical ordered-path-cover profiles; unconditional strict separation for every such common profile; and an all-order sixteen-state reversal-cube theorem. Combined with incoming disjoint isolation, the theorem decodes every state of the two-level reversal cube. **FINITE-EXACT:** complete cubic/quartic response data for a pair-free nine-name covered example. **NOT CLAIMED:** log-concavity for every tournament profile, arbitrary prime quotient separation, or transfer of an integer gap through an arbitrary strong outer context.

## 1. Exact class and physical query

Let A_1,A_2,A_3 be nonempty tournaments on disjoint fixed named sets V_1,V_2,V_3. Assume their complete ordered path-cover profiles agree:

    q_k(A_1)=q_k(A_2)=q_k(A_3)=q_k,       k>=1.

Here q_k counts ordered spanning covers by k nonempty directed paths. Put q_0=0 and q_k=0 beyond the common block order m. This q_0 convention is for nonempty blocks, not the Hamilton count of an empty graph.

Let

    T=C3[A_1,A_2,A_3], with V_1->V_2->V_3->V_1,
    U=C3^op[A_1,A_2,A_3], with all child interiors unchanged.

Add named ports x->y. Take the actual attachment query

    a=V_1,       c=V_2.

Thus V_1 beats x and x beats V_2,V_3; V_2 beats y and y beats V_1,V_3. Every port sees each entire child uniformly.

**Theorem 1 (exact gap).**

    F_T(V_1,V_2)-F_U(V_1,V_2)
      = -sum_(k>=0) (k+1)(q_(k+2)+2q_(k+1))
                         (q_(k+1)^2-q_k q_(k+2)).          (1)

The sum is finite. It is an exact integer identity, not a congruence or an asymptotic statement.

## 2. Run-generating proof of the exact formula

Use commuting run markers X,Y,Z and put P=XYZ, R=1-P. A directed Hamilton path in a block substitution has maximal nonempty child-runs. For the cyclic quotient, their color word follows the cycle and may revisit a color. Its old-run path matrix is G=N/R, with

    N = [ X,  XY, P;
          P,  Y,  YZ;
          XZ, P,  Z ].

The entry G_ij counts nonempty old-run words starting in color i and ending in color j. Its coefficients are purely formal run counts. To turn a word with (r,s,t) runs into actual old vertices, apply

    L_q(X^r Y^s Z^t)=q_r q_s q_t.

Each ordered cover of a child fills its run slots in their order of appearance. This is a bijection, including multiple visits to a child. Since q_k vanishes for k>m, the functional uses only finitely many coefficients of every displayed series; no analytic convergence or substitution of actual vertices into a rational function is involved.

For a column selector a and its complement o=1-a, define the numerator factors

    L_a=R+1^t N a,
    R_o=R+o^t N 1,
    M_(o,c)=o^t N c.

For the query a=(1,0,0)^t and c=(0,1,0)^t, let p=1-c. The generating numerator for words containing each port exactly once is

    W(N)=L_a (R+M_(o,c)) R_p + L_c M_(p,a) R_o.             (2)

The denominator is R^3. The first term puts x before y, allowing their direct arc or a nonempty old middle segment. The second puts y before x, forcing a nonempty middle segment because y->x is forbidden. The prefix and suffix may be empty. This accounts for all actual port placements, not just the adjacent-port contribution.

For the converse quotient replace N by N^t. Let D=W(N)-W(N^t). A direct polynomial expansion gives

    sum_(pi in S3) D(pi(X,Y,Z)) = 2 R K,                    (3)

where

    K = X^2Y^2+X^2Z^2+Y^2Z^2
        +X^2Y+X^2Z+XY^2+XZ^2+Y^2Z+YZ^2
        -XYZ(X+Y+Z)-6XYZ.

Equation (3) is a finite integer polynomial identity; `verify.py` expands and checks every coefficient using only dictionaries of integer exponent triples. The same identity can be checked by multiplying the three explicit factors in (2).

The functional L_q is invariant under all six variable permutations because the profiles agree. Applying it to (3)/R^3 gives the gap as (1/3)L_q(K/R^2). Use

    1/R^2=sum_(k>=0)(k+1) X^kY^kZ^k.

For each k, one third of L_q(K P^k) equals

    q_k q_(k+2)^2 +2q_k q_(k+1)q_(k+2)
       -q_(k+1)^2 q_(k+2)-2q_(k+1)^3
    =-(q_(k+2)+2q_(k+1))(q_(k+1)^2-q_k q_(k+2)).

This proves (1). The division by three occurs in characteristic zero after an exact symmetric identity; no division modulo three or four is used.

### Independent cover interpretation

Let C_ij be the aggregate ordered two-cover count with the first path ending anywhere in V_i and the second starting anywhere in V_j. Its cyclic generating functions give

    C_12-C_21 = L_q( XY*((1+Z+YZ)(1+Z+XZ)
                        -(1+X+XZ)(1+Y+YZ))/R^2 ).          (4)

The same symmetric expansion shows that (4) is the negative of (1). This supplies a second count interpretation of the gap. The equality with the literal query gap uses equal-profile symmetry; it is not asserted for arbitrary unequal profiles.

## 3. An exact positive-square identity removes every profile hypothesis

For a finite sequence q_j>=0 with q_0=0 and q_1>0, extend q by zero at all other indices and write

    P(q)=sum_(j>=1) j(q_(j+1)+2q_j)(q_j^2-q_(j-1)q_(j+1)).

Individual Turán terms may be negative. Nevertheless the entire sum is always strictly positive, by the exact identity

    3P(q) = q_1^3
      +sum_(j>=1) j(q_(j-1)-q_(j+1))^2
                       (q_(j-1)+3q_j+2q_(j+1))
      +sum_(j>=1) (q_j-q_(j+1))^2
                       ((j+2)q_j+(2j+1)q_(j+1)).          (5)

All sums are finite. If m is the last nonzero index, the first square sum includes j=m+1. Omitting that boundary term would make the identity false. Every term on the right is nonnegative, and q_1^3 is strictly positive. In particular P(q)>=q_1^3/3>0 without log-concavity, unimodality, or real-rootedness.

Here is an explicit all-order verification. Put a=q_(j-1), b=q_j, c=q_(j+1). The two square summands expand respectively as

    j(a^3+3a^2b-6abc-3ac^2+3bc^2+2c^3),
    (j+2)b^3-3b^2c-3jbc^2+(2j+1)c^3.

The bc^2 terms cancel. After reindexing, the cube coefficient is 6j for each q_j^3, except that q_1^3 initially has coefficient five; the explicitly displayed q_1^3 repairs it. Also

    sum_(j>=1)(3j q_(j-1)^2 q_j-3q_j^2 q_(j+1))
        =sum_(j>=1)3j q_j^2 q_(j+1).

The remaining terms are -6jabc-3jac^2, exactly three times the expansion defining P. This proves (5) for every finite sequence. The checker independently expands both sides as cubic polynomials in algebraically independent q_1,...,q_m at cutoffs 1 through 20; those finite checks supplement, rather than replace, this reindexing proof.

**Theorem 2 (unconditional equal-profile noncancellation).** For any three nonempty children with the same ordered-cover profile, the actual physical gap in (1) is strictly negative. Thus the fully covered laminar operation consisting of the whole set and all three child modules cannot hide. Quantitatively,

    F_U(V_1,V_2)-F_T(V_1,V_2) >= H(A_1)^3/3 > 0.

All tournament Hamilton counts are positive. Apply (5) to their common finite nonnegative q profile and then (1). No property of that profile beyond q_1>0 and nonnegativity is used.

Conversing a child preserves its entire q profile and, because this query is uniform on that child, preserves the query value exactly. Hence the same two distinct root-query values work regardless of all three child-converse bits.

### Decode the complete two-level reversal cube

Suppose m>=2, so every child reversal is nontrivial. Parameterize a state by the cyclic quotient orientation theta and three child-converse bits epsilon_i. The uniform query first determines theta by Theorem 2. For that fixed quotient orientation, incoming [S1-S2 disjoint isolation](https://github.com/eliottcassidy2000/math/blob/b6135263ba67ee2908df7625c67a88a3df8fbe03/05-knowledge/results/simultaneous_module_reversal_20261104/PROOF.md) supplies one actual lifted query per child that determines epsilon_i, independently of the other two bits.

Thus all sixteen named full responses are distinct. Four adaptive queries suffice, or seven fixed queries: the root query and the three child decoders for each of the two possible quotient orientations. No optimality is claimed. The four original operations, the whole set and its three children, correspond bijectively to (theta,epsilon), since a whole-set reversal toggles the quotient and all three child orientations.

This adds a proved fully covered laminar class to the exposed-anchor theorem. It is not a theorem for arbitrary prime quotients or arbitrary unequal child profiles.

## 4. Transitive and cyclic examples

For a transitive m-name block, each nonempty vertex subset has exactly one Hamilton path. Therefore

    q_k = k! S(m,k),

where S(m,k) counts set partitions into k nonempty blocks. The ordered-partition recurrence

    q_k(m+1)=k(q_k(m)+q_(k-1)(m))

supplies convenient integer controls. For every m>=2, the full sixteen-state cube generated by the whole tournament and the three transitive children in C3[TT_m,TT_m,TT_m] has distinct named two-port responses. The root-query gaps for m=1,2,3,4 are respectively

    -2, -36, -2384, -407568.

Only the m=1 two-state root case applies when the children are singletons.

The theorem equally applies to three copies of any nontransitive tournament, or to nonisomorphic tournaments with the same complete q profile. Log-concavity of arbitrary tournament profiles is not established or needed here. The finite checks through order five are controls only. The ordinary q polynomial need not be real-rooted: for a directed triangle it is 3t+6t^2+6t^3, whose two nonzero roots are nonreal.

## 5. Pair-free nine-name fully covered control

Take all three children to be directed triangles, each with profile

    q=(0,3,6,6).

The formula gives the exact gap -2052. Both cyclic quotient orientations are strong and contain no homogeneous pair. They have H=3159, while the physical block-uniform query gives

    F_T(7,56)=70335,       F_U(7,56)=72387.

This is a fully covered laminar reversal, obtained by reversing the whole set and its three triangle children. It is not exposed by the homogeneous-pair sign obstruction. The complete difference has 54 nonzero quadratic, 432 nonzero cubic and 540 nonzero quartic coordinates; the constant and linear parts cancel.

`full_response.json` retains all 2,116 named coordinates of bidegree at most (2,2), including every zero, for T, U and the whole converse of T. It also retains the entire interpolation query battery for each. The whole-converse quartics cancel, whereas the internal quotient-converse difference retains them. No degree is dropped merely because the change is described using converses.

## 6. Ambient scope

Theorem 2 and its cube corollary are theorems about the actual cyclic block substitution T and U. Ordered earlier/later SCC padding preserves the separator by the exact S187 localization factor, a positive odd product of Hamilton counts.

No arbitrary strong outer quotient transfer is asserted. The integer gap can be divisible by four, as -36 already shows, and a scalar inequality for a port-extended core does not by itself prove coordinatewise domination of that core's entire ordered-cover profile. Applying incoming T1 outside its modulo-four scope would be invalid. The unrestricted covered-laminar residual remains open.

## 7. Reproduction

    g++ -std=c++17 -O3 -Wall -Wextra -Werror full.cpp -o full
    ./full > full_response.json
    python verify.py > receipt.json
    python -O verify.py > receipt.optimized.json

The Python verifier uses standard-library exact polynomial arithmetic for (3), independently computes all ordered-cover profiles through old order five, and checks literal blowups through twelve old names. It verifies both identities for artificial non-log-concave sequences too, including exhaustive nonnegative sequences of lengths one through five with entries zero through three and positive first entry. The positive-square identity is checked coefficientwise through cutoff 20. The transitive recurrence is checked through m=50 as a control, not an assumption of the theorem.

The C++ engine evaluates all 6,348 interpolation queries for the three nine-name responses by integer subset DP, with 96 additional dense controls. All counts are at most 11!, so uint64 arithmetic is exact. Python independently reinterpolates all retained coefficients and checks twelve additional literal dense counts by a different subset recurrence. Normal and optimized Python receipts agree byte for byte.
