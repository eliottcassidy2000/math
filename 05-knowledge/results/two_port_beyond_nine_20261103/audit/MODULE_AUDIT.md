# Independent audit of the odd-module obstruction

Status: **PROVED — accepted**. Independent review of
[algebra/PROOF.md, A1-A3](../algebra/PROOF.md), 2026-11-03.
The all-order conclusion applies to uniform odd modules of size at least
three. No converse classification or universal two-port reconstruction
claim follows.

The proof's normalization is correct. q_k counts ordered nonempty path
slots. The slot-permutation action is free and preserves each monomial,
so every polynomial coefficient is divisible by k!. Homogeneous Euler
evaluated at a 0/1 adjacency gives

    sum_actual_arcs partial_uv q_k = (r-k)q_k.

The all-directed-arc derivative has a different normalization:

    sum_(u!=v) partial_uv q_k = k q_(k+1).

An unordered(k+1)-cover has k(k+1) ordered choices of paths to join and k!
orders of the resulting k slots. Comparing this with its(k+1)! ordered
covers gives factor k. The joining arc need not exist: differentiation
removes that variable. This confirms that absent directions must remain
in the all-arc sum.

For a uniform module M in Q[M], contracting maximal M-runs gives an
identity before evaluation,

    H(Q[M])=sum_k c_k q_k(M).

The exterior skeleton coefficients are nonnegative integers independent
of the internal arc variables. A single run contracts bijectively to
the replaced quotient vertex, so c_1=H(Q), which is odd. This step uses
the module hypothesis; it is invalid for an arbitrary odd subset.

Differentiating gives the claimed current

    J_M=sum_k c_k[2(r-k)q_k-kq_(k+1)].

For every k>=2, 2q_k is divisible by4 and k(k+1)! is divisible by4.
Only k=1 can survive modulo four. Complementary subsets in
P_M(2)=sum_A H(M[A])H(M[V\\A]) form2^(r-1) pairs, each contributing2
modulo4; for r>=2 this proves P_M(2)=0 modulo4. Since
q_2=P_M(2)-2H(M), the remaining bracket is2rH(M)-P_M(2), congruent to2
when r is odd. The odd coefficient c_1 preserves that residue. Thus
J_M=2 modulo4 in every ambient quotient. Singleton modules have current
zero, and the proof correctly excludes them.

The transpose step is also valid at the polynomial level:
q_k(X^t)=q_k(X). Differentiating and then evaluating at converse adjacency
transposes the internal gradient while leaving every exterior c_k fixed.
If the full ordered C matrices agreed, the internal current would be
zero, contradicting its residue. Other internal entries need not be
known and no sign-of-current adjacency rule is assumed.

For pairwise disjoint modules, the joint skeleton expansion has products
of their respective q_k polynomials. Reversing any other module preserves
each of its q_k values and hence preserves the first module's entire
internal gradient. Therefore one differing reversal bit is detected even
when other modules are reversed simultaneously. The resulting2^s choices
are distinct as labelled tournaments; the proof does not assert that all
choices are nonisomorphic.

Finally, the modulo-eight claim uses only an integer observable. Session27
proves J_uv=P_u-P_v-D_uv, where each P and D is an integer difference of
actual query values. For a module and its internal converse, fixed
reference orientations give current contrast2J_M=4 modulo8. No extraction
of C by division by two in a residue ring is used. Thus query values
modulo8 already detect the reversal. For regular modules the potential
sum cancels by equal indegree and outdegree, as claimed.

The normal and optimized algebra receipts agree and pin the current
source. The recorded controls include4,280 ambient pairs,97,144 internal
transpose entries,1,032 literal core-profile/derivative checks, and128
independent reversal variants. The odd nonmodule hostile with current-12
properly tests the essential hypothesis. These are controls of the proof;
the all-order statement rests on the polynomial argument above.
