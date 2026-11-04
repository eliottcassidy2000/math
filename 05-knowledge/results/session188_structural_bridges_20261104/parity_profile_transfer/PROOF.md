# Endpoint-parity transfer through arbitrary child tournaments

## 1. A universal one-query separator

Let Q be any fixed tournament on m quotient labels. Let A_i be arbitrary nonempty tournaments on disjoint physical child sets, and put

    T=Q[A_i],    U=Q^op[A_i].

More generally, the interiors on the U side may differ provided every corresponding complete ordered path-cover profile agrees. Physical names and the partition into child sets remain fixed. Let S_i(Q),E_i(Q) count quotient Hamilton paths starting and ending at i, and define

    p_i(Q)=S_i(Q)+E_i(Q) modulo2.

If p_i=1, the single actual query with all of A_i pointing to x, every other old name receiving from x, and y a universal source relative to old names satisfies

    F_T(1_(A_i),0)-F_U(1_(A_i),0)=2 modulo4.             (1)

Ports always have fixed edge x->y. The theorem allows singleton children, arbitrary quotient order and arbitrary child interiors. It is a sufficient separation criterion; p=0 is not a collision claim.

**Proof.** The independently audited uniform-block start identity (`uniform_block_starts/PROOF.md`, Section1) gives

    Delta(1_(A_i),0)-2 Delta(0,1_(A_i))
       =2[S_(A_i)(T)-S_(A_i)(U)].

Every augmented tournament has an odd Hamilton count by Rédei's theorem. Thus Delta(0,1_(A_i)) is even. In the maximal-child-run expansion of an old Hamilton path, q_k(A_j) is even for k>=2, because the k ordered cover slots can be permuted freely. The one-run count q_1(A_j)=H(A_j) is odd. Modulo two only words visiting every child exactly once survive. Consequently the two aggregate start counts reduce to S_i(Q) and S_i(Q^op)=E_i(Q). Reducing the displayed identity modulo four proves (1). No independent endpoint products are substituted for disjoint path covers.

## 2. The entire coarse contrast for nontrivial children

Assume now that every |A_i|>=2. Let a,c be arbitrary masks on quotient labels, expanded uniformly over the corresponding physical children. Then the **complete** block-uniform contrast obeys

    Delta(a,c)=2 sum_i p_i(Q)(a_i+c_i) modulo4.           (2)

Thus all mixed, quadratic, cubic and quartic coarse contrast coefficients vanish modulo four. This is a congruence for the full actual responses, not a truncation assumption and not a statement that the full physical named cubic/quartic terms vanish.

### Universal cover residues

For any tournament A of order n>=2,

    q_1(A) is odd,
    q_2(A)=2 modulo4,
    q_3(A)=2[n odd] modulo4,
    q_k(A)=0 modulo4 for k>=4.

Indeed q_k/k! counts unordered partitions into k nonempty vertex subsets with one Hamilton path on each. Every such subset has odd Hamilton count. Therefore q_k/k!=S(n,k) modulo two. Here S(n,2)=2^(n-1)-1 is odd, while S(n,3) has the parity of n for n>=2, by its recurrence and the base value S(2,3)=0. The factorial handles k>=4.

In a block-uniform run expansion modulo four, at most one child can be visited more than once. A word visiting a specified child three times has the form

    prefix, i, L1, i, L2, i, suffix,

where both return segments L1,L2 are nonempty and use disjoint remaining symbols. Exchange those segments. Every boundary incidence is unchanged, all port directions are retained within the segments, and the two nonempty disjoint symbol sets make the exchange fixed-point-free. Thus the number of these skeletons is even, so its q_3 contribution vanishes modulo four.

Let R_i,2(Q;a,c) count run words using every other child and both ports once and child i twice. The response therefore satisfies

    F_T(a,c)=(product_j H(A_j)) F_Q(a,c)
                         +2 sum_i R_i,2(Q;a,c) modulo4.

In the converse contrast the first factor is odd and F_Q-F_(Q^op) is even, so that odd factor can be removed modulo four. The remaining expression is exactly the contrast obtained when every child is a transitive pair. This proves a universal transfer to pair children, independent of all child orders and interiors.

### Pair-cycle calculation of the complete linear form

Use the classical odd-cycle Hamilton expansion of Grinberg and Stanley, [The Rédei–Berge symmetric function of a directed graph](https://www.cip.ifi.lmu.de/~grinberg/algebra/redeiberge-long.pdf), Theorem1.39. Modulo four it reads

    H(R)=1+2(number of nontrivial directed odd cycles of R) modulo4.

All unused vertices are singleton cycles. Consider a pair blowup with uniform port masks. If a used pair contributes one name, or its two names are not consecutive in a cycle, exchanging its physical names gives a factor two. Such cycle contributions vanish modulo four.

An old-only odd cycle and an odd cycle containing both ports each contain an odd number of old names, so at least one pair is free in this sense. Both kinds vanish. A surviving cycle containing exactly one port consists of consecutive whole pairs, each in its forced internal order. Collapsing them and deleting the port gives a simple directed quotient path P of at least two vertices. If P starts at s and ends at t, its first-port attachment factor is a_t(1-a_s). Its converse path contributes a_s(1-a_t). Their contrast is a_t-a_s, which is a_t+a_s modulo two. The same argument applies to c for the other port.

Hence the complete contrast is twice the linear endpoint incidence sum over all simple quotient paths. For every tournament Q,

    total simple paths starting at i = E_i(Q) modulo2,
    total simple paths ending at i = S_i(Q) modulo2,

including singleton paths. To prove the second identity induct on |Q|. The total ending-path count at i is 1 plus the sum of total ending-path counts at vertices w->i in Q-i. By induction replace them by S_w(Q-i). Their total over all w is H(Q-i), which is odd. The result is the sum over i->w of S_w(Q-i), precisely S_i(Q). The singleton base is immediate; converse gives the first identity. Singleton paths contribute twice to the endpoint incidence sum and disappear modulo two. This proves (2).

## 3. Ambient transfer and original-module consequences

A contrast equal to two modulo four survives every fixed ambient substitution when both ports use the original core module's exterior incidences. To see this directly, the two augmented cores have the same order. Their q_2 and q_3 residues above agree, and every higher q_k vanishes modulo four. In the occupied-run substitution formula, only their Hamilton difference remains. Its coefficient c_1 is the Hamilton count of the contracted ambient quotient and is odd. Thus the lifted physical query still has contrast two modulo four. No parity restriction on the augmented core order is needed.

For an original-module XOR family whose whole prime root is selected oddly, every proper selected module lies inside one root child. If the old Hamilton count increases, query (0,0) separates. Otherwise the independently established hereditary-equality theorem (`module_decomposition/PROOF.md`) preserves every child cover profile, so (1) applies whenever the root quotient has p_i=1. Therefore a fixed two-query battery, (0,0) and the indicated physical endpoint query, separates the original tournament from every such root-containing family. The same statement holds for families contained within a fixed ambient core, using its port lift.

More generally, in a scalar-equal original-module comparison, take the highest changed nodes in the persistent original strong-module tree. They form an antichain. If any highest prime quotient has a nonzero p_i, the physical query from (1) lifts through its unchanged ancestors. Other highest-node changes can be restored exactly because their profiles agree and the lifted ports see them uniformly. Thus a possible full-response collision must have p=0 at every highest changed prime quotient. This is an additional necessary condition, not a reconstruction theorem. It uses the antichain/restoration argument already proved in `highest_quotient_reduction/PROOF.md`.

## 4. A prime-five class with arbitrary child orders

The quotient with outgoing bit rows

    Q=(18,8,3,5,14)

is strong and prime. Its endpoint-parity mask is9: p_0=p_3=1 and every other p_i is zero. Thus either a=all of child0,c=0 or a=all of child3,c=0 is an all-order physical separator for arbitrary nonempty children, with arbitrary ambient lift.

The other two prime-five representatives from the finite classification,

    (20,1,2,7,14),    (24,17,3,6,12),

have p=0. Formula(2) proves coarse modulo-four blindness for their arbitrary nontrivial child blowups. It does not prove integer coarse equality or equality of full named responses. The existing exact prime-five triple certificate remains valid and nonzero on these classes.

Any quotient relabeling transports its child sets, endpoint vector and physical mask together; the two compared responses are never independently relabeled.

## 5. Full data and reproducibility

The verifier checks the universal cover residues and endpoint bridge on every labeled tournament through order five. It tests the three-return involution by a separate repeated-symbol automaton and compares the complete linear congruence with literal physical Hamilton counts on small quotients, arbitrary nontrivial children, and arbitrary ambient lifts. Singleton-child controls test the broader one-query theorem without claiming the full linear law there.

For the displayed prime-five quotient with five named pair children, it retains every one of the 3,136 Boolean coordinates of bidegree at most(2,2), including all zeros and all cubic/quartic coefficients of both full physical responses. Every one of its1,024 coarse queries is evaluated independently and checked against the complete named polynomial. The full named difference has nonzero cubic and quartic terms even though its coarse reduction modulo four is linear. A two-highest-prime-node control checks exact restoration of the unrelated changed node after physical port lifting.

Run `python3 verify.py > normal.json` and `python3 -O verify.py > optimized.json`, then compare the outputs. The standard-library verifier uses the sibling `uniform_block_starts/verify.py` primitives and an already installed C++17 compiler. The literal C++ engine counts exact integers, not residues; all tested augmented orders are at most17 and all counts are bounded by17!<2^64. No installation is performed. All checks remain active under optimization.

Frozen source allowlist: `PROOF.md`, `verify.py`, `literal.cpp`, `normal.json`, `optimized.json`. Exploratory `explore.py` and `exploration.json` are not required for this proof or reproduction.
