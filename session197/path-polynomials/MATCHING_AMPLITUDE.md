# Maximum matchings determine the leading clone-response layer

Session 197 refinement of `CONSEQUENCES.md`, using the exact universal skeleton operator proved in `UNIVERSAL_CONTEXTS.md`. **All-order coefficient formula and amplitude lower bound; not an inverse-spectrum algorithm.** No literature-priority claim is made.

Let (Q,r) have N exterior vertices. As in S196, form the bipartite cross graph B_r on O={v:r->v} and I={v:v->r}, with an edge u-v exactly when u->v. Let nu be its maximum matching size and mu the number of its maximum matchings. When nu=0, the empty matching is counted once.

Write the exact context operator in shifts as

    A_Q=sum_(j=-nu)^nu p_j(x) S^j,       (Sf)(x)=f(x+1).

The companion consequence theorem shows that p_nu is nonzero and has nonnegative integer coefficients. Write x^(overline nu)=x(x+1)...(x+nu-1), with value one for nu=0.

**Theorem.**

1. The degree of p_nu(x) is exactly N-nu, and its leading coefficient is mu.
2. Coefficientwise,

       p_nu(x) >= mu x^(N-2nu) x^(overline nu).                  (1)

3. If N=2nu, then equality holds:

       p_nu(x)=mu x^(overline nu).                             (2)

4. The highest clone-run coefficient in S196 satisfies

       a_(nu+1)(Q,r)=p_nu(1) >= mu nu!.                        (3)

In particular, a physical target prescribing the highest exponential amplitude below mu nu! cannot be achieved by cloning this rooted template. Equality in (2) gives an exact amplitude in the balanced case. These are constraints on a specified template, not a converse realization theorem.

**Proof.** A skeleton contributing to shift nu has exactly nu marked clone-runs. Every outside path between consecutive clone-runs starts in O and ends in I, so has at least two vertices. If a crossing skeleton cycle has t runs and l exterior vertices, l>=2t. If a pure exterior odd cycle has length l, then l>=3.

For a skeleton with u unused exterior vertices, c_0 pure exterior cycles and c_1 crossing cycles, its contribution to p_nu is

    2^c_0 x^(u+c_0+c_1).

The highest forward coefficient of every crossing parity kernel is +1, which is why no other sign enters this expression. Since

    N-(u+c_0+c_1)
      = sum_(pure cycles)(l-1) + sum_(crossing cycles)(l-1)
      >= sum_(crossing cycles)t = nu,

the degree is at most N-nu. Equality requires no pure exterior cycle, and l-1=t for every crossing cycle. Together with l>=2t and t>=1 this forces t=1,l=2: every crossing cycle is one triangle consisting of a cross edge and one clone-run. The nu outside edges form a maximum matching, and all other outside vertices are unused. Conversely every maximum matching supplies exactly one such skeleton packing. This proves the degree and leading coefficient in part 1.

More generally fix any maximum matching M, regarding its nu directed O-to-I edges as nu distinct labelled pair objects. For any permutation of these objects, group the pair objects by the cycles of the permutation. Each permutation cycle gives a crossing cycle skeleton: traverse an O-to-I pair, then its marked clone-run, then the next pair, and so on. All boundaries are legal because I->r->O. Leave the N-2nu unmatched outside vertices unused.

The cutting inverse recovers both M and the permutation, so these skeletons are distinct as M and the permutation vary. A permutation with c cycles contributes x^(N-2nu+c) to p_nu. The classical elementary identity

    sum_(permutations pi of nu objects) x^(number of cycles of pi)
       = x^(overline nu)

follows by inserting the largest object: it either forms its own cycle, contributing x, or is inserted after one of the previous nu-1 objects. Summing over the mu matchings proves (1).

If N=2nu, every shift-nu skeleton uses all exterior vertices in its crossing paths. There is no room for a pure exterior cycle, unused exterior vertex, or crossing path of length greater than two. Thus every skeleton is exactly one of the matching/permutation constructions, proving (2).

Finally, evaluation at x=1 gives

    H(Q[r:Tr_m])=sum_j p_j(1)(1+j)^m.

The unique largest absolute base is nu+1, with coefficient p_nu(1). The inherited S196 Stirling clone-run expansion has the same largest base with coefficient a_(nu+1). Equality of the sequences, or their ratio to (nu+1)^m as m tends to infinity, identifies the two. Equation (3) follows from (1). The case nu=0 says a_1=H(Q-r)>=1 and agrees with the fact that r lies in no directed cycle. QED.

## Why this helps, and what it does not supply

This identifies a useful inverse-design parameter: the highest response is controlled first by a cross-matching problem and then by the lower-degree terms of p_nu. Balanced templates make the first layer completely explicit. For example a unique perfect cross matching at nu=2 gives p_2=x(x+1) and top clone amplitude two; at nu=3 it gives p_3=x(x+1)(x+2) and amplitude six. Multiple perfect matchings multiply these quantities by their exact number.

It remains necessary to control the lower exponential coefficients, including whether unwanted bases cancel and whether a desired offset is achieved. Maximum-matching data alone does not determine them. In particular a positive normal-form operator or a desired positive amplitude is not automatically an actual tournament context. The F/G physical port constructions in the companion packet escape particular smaller grammars, not this inverse obligation.

## Verification

Run `python3 verify_matching_amplitude.py` and the optimized version. The checker uses the already independently audited literal skeleton enumerator, a separate brute-force enumeration of bipartite matchings, and direct subset Hamilton DP for actual clone graphs. It verifies the coefficient lower bound, leading coefficient, balanced equality, and exact top amplitude. The deterministic finite universe and the exact rows for higher-rank controls are emitted in its JSON output. No finite table substitutes for the proof above.
