# Rooted merge envelopes close the single-triangle deletion sector

Session 193, 5 November 2026. **PROVED and independently reviewed.** The statements below have full proposed proofs; their finite checks are controls, not substitutes for proof. The unrestricted root-deletion minimum and the global Hamiltonian-path spectrum remain OPEN. No literature-priority claim is made.

## 0. Inheritance and the new bridge

Write H(T) for the number of directed Hamiltonian vertex words, with H(empty)=1. For a root a, let e_a(T) count Hamiltonian paths ending at a. Let R_a(T) count ordered spanning pairs (P,Q) of vertex-disjoint directed paths, with P possibly empty and Q nonempty starting at a. Set E_a(T)=R_a(T)+e_a(T); the two summands are distinct tagged types.

The inherited Busch minimum B(n), for strong tournaments, is

    B(1)=1,
    B(3q)=3*5^(q-1), B(3q+1)=5^q, B(3q+2)=9*5^(q-1).

There is no strong tournament of order two. The displayed formulas have their natural positive-order domains. The singleton is allowed separately throughout.

Sources:

* Arthur H. Busch, “A Note on the Number of Hamiltonian Paths in Strong Tournaments”, Electronic Journal of Combinatorics 13 (2006), N3, Lemmas 1–2 and Corollary 1, primary text inspected at https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/ . The source supplies the unrooted strong decomposition, path merging, and B. The root-avoiding refinement below is proved directly.
* THM-4744, `01-canon/theorems/THM-4744-interval-chains-attain-busch-and-are-the-extremal-strong-tournaments-through-order-twelve.md`: known Moon chain product and explicit endpoint response, also rederived in [Session 191 interval-chain growth](../../session191_tournament_growth_20261104/deletion_growth/INTERVAL_CHAINS.md).
* THM-4745, `01-canon/theorems/THM-4745-exact-gluings-of-interval-chains-and-the-reduction-of-busch-equality-to-a-one-way-lemma.md`: ordinary H merge inequality and its information loss.
* [Session 191 single-triangle reduction](../../session191_tournament_growth_20261104/deletion_growth/SINGLE_TRIANGLE.md): exact one-root-triangle side formula and weak-side floor. Its unresolved two-strong-side response optimization is the target closed conditionally on acceptance of the proofs below.
* [Session 192 weighted cuts](../../session192_surprising_approaches_20261105/weighted_cuts/PROOF.md) and [Session 192 cut sumsets](../../session192_surprising_approaches_20261105/cut_sumsets/PROOF.md): E_a(T) is exactly the exposed-edge capacity c_(a,b) in the ordered join T->b, with b a singleton. This is a concrete cut-boundary interpretation, not a claim that cut capacities are fixed under later growth.

Live concept board: marked path-pair response; tagged end-path response; strong gluing away from the mark; root-deletion transitivity; sharp chain support function. The inherited scalar R>=H loses precisely the e coordinate. Keeping the sum R+e repairs its support envelope.

## 1. Root-avoiding strong decomposition

**Lemma 1.** Let T be strong of order n>=5 and fix a vertex a. At least one of the following holds:

1. T is nearly transitive: it is a transitive tournament with its first-to-last arc reversed.
2. T-a is transitive.
3. There are induced strong subtournaments A,B covering T, each of order at least three, whose intersection is a single vertex v different from a. Thus a belongs to exactly one factor.

**Proof.** We prove (3) whenever (1) and (2) fail.

First suppose T remains strong after deletion of any vertex. Since T-a is nontransitive, it contains a directed triangle C disjoint from a. If T-C is strong, every vertex c of C has an in-neighbor and an out-neighbor outside C: deletion-strongness implies minimum indegree and outdegree at least two, and C supplies only one of each. Thus C and (T-C)+c are strong and meet at c, which is not a.

If T-C is weak, let W_1 and W_t be its source and sink strong components. At least two vertices of C have an out-neighbor in W_1. Otherwise deleting the sole such vertex would leave W_1 with no incoming arc from its complement, contradicting deletion-strongness of T. Likewise at least two vertices of C have an in-neighbor in W_t. These two subsets of the three-element C meet at some c. The tournament (T-C)+c is strong: c reaches the source component and is reached by the sink component. Again the two factors are C and (T-C)+c, meeting away from a. The complement T-C has at least two vertices because n>=5; these factors each have order at least three.

Now suppose T has a vertex whose deletion is weak. If T-a is weak, choose v=a. Otherwise choose any cut vertex v, necessarily different from a. Write the strong components of T-v as W_1->...->W_t, t>=2. Strongness of T gives an arc v->W_1 and an arc W_t->v.

If |W_1|>=3, choose w in W_1 with v->w. We can choose w different from a. This is automatic if v=a or a is outside W_1. In the only remaining case, T-a is strong; if a were the sole out-neighbor of v in W_1, then W_1-{a} would dominate every vertex of its nonempty complement in T-a, making T-a weak. This contradiction supplies w!=a. Take A=W_1 and B={v,w} union W_2 union ... union W_t. They are strong, cover T, meet at w, and have order at least three. For B, use v->w, w dominating every later component, and the sink component reaching v. If |W_t|>=3, the dual argument gives the same conclusion, choosing w!=a with w->v.

We may therefore assume W_1={w_1} and W_t={w_t}. If a middle component is nontrivial, choose a directed triangle C in that component and any c in C different from a. Take A=C and B=T-(C-{c}). The latter is strong because it retains the directed triangle v->w_1->w_t->v, and every retained middle vertex is reached from w_1 and reaches w_t. These factors cover T and meet at c!=a.

The remaining case is that T-v is transitive. If v=a this is alternative (2), contrary to the assumption. Thus v!=a. Let W be its middle vertices, and split W into W^-={x:x->v} and W^+={x:v->x}. If both are nonempty, the strong factors

    A={v,w_1} union W^-,  B={v,w_t} union W^+

meet exactly at v!=a: each contains its evident directed cycle through its two named vertices and all its other vertices. If W^+ is empty, v has sole out-neighbor w_1 and reversing v->w_1 makes T transitive. If W^- is empty, v has sole in-neighbor w_t and reversing w_t->v makes T transitive. Either is alternative (1), contrary to assumption. This exhausts the cases. QED.

This is not asserted for arbitrary digraphs. The ordered SCC structure of a tournament and the two-element neighbor sets in a triangle are used explicitly.

## 2. A tagged response merge inequality

**Lemma 2.** Let induced tournaments A,B meet at exactly one vertex v, with union T, and let a belong to A but not B. Cross-arcs are arbitrary. Then

    E_a(T) >= E_a(A) H(B).                              (1)

Neither factor nor T needs to be strong.

Recall the elementary path-merging fact: two disjoint directed paths in a tournament admit an order-preserving merge. Fix a deterministic choice of one merge for each pair, including empty paths. It follows by repeatedly inserting the next vertex of the second path into the appropriate suffix after its previously inserted vertex. Restricting the merged word to the two original vertex sets recovers the inputs.

Here is an explicit injection proving (1). Write a Hamiltonian path of B as C_- v C_+.

For an R-type input (P,Q), where Q starts at a:

* If v lies in P, write P=P_- v P_+. Output the R-type pair

      ( merge(P_-,C_-) v merge(P_+,C_+), Q ).

* If v lies in Q, write Q=Q_- v Q_+, with Q_- starting at a and nonempty because a!=v. Output

      ( merge(P,C_-), Q_- v merge(Q_+,C_+) ).

Every adjacency at v is valid because the merged prefixes end with an old predecessor of v and the merged suffixes start with an old successor of v.

For an e-type input, let L=L_- v L_+ be a Hamiltonian path of A ending at a. Merge it with the B path in the usual way:

    W=merge(L_-,C_-) v merge(L_+,C_+).

If W ends at a, output W as an e-type object. Otherwise split W immediately before a and output the resulting R-type pair. Since L ends at a, every vertex after a in W belongs to B, and this tail is nonempty.

**Injectivity and separation.** For every output, the original A object and B Hamiltonian path are recovered by restrictions; for the latter concatenate the B restrictions of the two bins when the output has R type. Distinct input pairs within any one case cannot coincide. An e-input giving e-output has a separate tag. An e-input giving R-output has v in its first bin and at least one B vertex in its second bin. An R-input with v in its first bin puts every B vertex in that first bin, while an R-input with v in its second bin keeps v there. Hence none of these images overlap. QED.

The condition a!=v is load-bearing in this proof. We do not silently use (1) for gluing at the marked root. Lemma 1 is exactly what supplies the missing hypothesis.

## 3. The universal rooted envelope

Always R_a(T)>=H(T): split a Hamiltonian path immediately before a. Also e_a(T)<=H(T)<=R_a(T). For strong T, e_a(T)>=1.

**Theorem 3.** For every strong rooted tournament (T,a) of order n>=4,

    R_a(T)+e_a(T) >= 2 B(n-1).                          (2)

For n=3 the exact minimum is four; for n=1 it is two.

**Proof.** Induct on n. The order-three case is the directed triangle. The unique strong order-four tournament is nearly transitive; alternatively its four-vertex path count and endpoints may be checked directly. For a nearly transitive T of any order n>=4,

    E_a(T)>=H(T)+1=2^(n-2)+2>=2B(n-1),

since a nearly transitive tournament of order n-1 has 2^(n-3)+1 paths and witnesses B(n-1)<=2^(n-3)+1.

If T-a is transitive, list it as x_1->...->x_(n-1). Strongness forces a->x_1 and x_(n-1)->a. There is exactly one Hamiltonian path ending at a, namely the chain followed by a. In an R-cover, Q={a} gives one choice. If the first chain vertex in Q is x_i, then a->x_i is required; all later chain vertices can independently be assigned to P or Q, while all earlier ones must belong to P. Therefore

    e_a=1,
    R_a=1+sum_(i:a->x_i) 2^(n-1-i),
    E_a>=2+2^(n-2)>=2B(n-1).                           (3)

In all other cases Lemma 1 supplies A,B meeting away from a; exchange their names so a belongs to A. Put p=|A|, q=|B|, so p,q>=3 and p+q-1=n. Lemma 2 and Busch give

    E_a(T)>=E_a(A)H(B)>=E_min(p) B(q).

If p>=4, induction bounds this by 2B(p-1)B(q)>=2B(n-1). The product inequality follows by concatenating optimal interval chains of orders p-1 and q, which gives a strong chain of order p+q-2=n-1. If p=3, it is at least 4B(q)>=2B(q+1)=2B(n-1), because B(q+1)<=2B(q) directly from the three residue formulas. QED.

Define c_1=1, c_3=1, and, for n>=4,

    c_n=2B(n-1)-B(n).                                  (4)

All these numbers are positive and at most B(n). The singleton order-two gap is never filled by a fictional strong tournament.

**Corollary 4 (sharp support function).** For every strong rooted n-tournament, n=1 or n>=3, and every real t in [0,1],

    R_a(T)+t e_a(T) >= B(n)+t c_n.                      (5)

The lower bound is attained at every t by one fixed rooted interval chain for that order.

The inequality is (1-t)R+t(R+e) and the two proved bounds. For attainability, mark the last spine endpoint of an optimal chain. Its R is B(n) and its e is the product of the factors preceding its last block. Choose the last block as follows:

* n=3: the single block 2, with (R,e)=(3,1);
* n=3q+1, q>=1: q blocks 3, ending in 3, giving e=5^(q-1);
* n=3q+2, q>=1: one block 4 and q-1 blocks 3, ending in 4, giving e=5^(q-1);
* n=3q, q>=2: one block 2 and q-1 blocks 3, ending in 3, giving e=3*5^(q-2).

These values are exactly c_n. For n=1 use the singleton. The interval [0,1] is essential; this does not classify the entire Pareto frontier for arbitrary weights.

## 4. Exact single-triangle deletion minimum

Suppose T is strong and its marked vertex v lies in exactly one directed triangle. Let A=N^-(v), B=N^+(v), and let b->a be the unique cross-arc from B to A. The inherited S191 identity is

    Delta(T,v)=H(T)-H(T-v)=R_a(A) R'_b(B)+e_a(A)s_b(B), (6)

where R'_b is the ending analogue of R. Strongness says a reaches every vertex of A and every B vertex reaches b.

**Corollary 5 (two strong sides).** For fixed side orders p,q, each 1 or at least 3, the minimum of (6) over strong A,B and all roots is exactly

    B(p)B(q)+c_p c_q.                                  (7)

Indeed, put R=R'_b(B), s=s_b(B). Then t=s/R belongs to [0,1]. Applying (5) on A gives (6)>=B(p)R+c_p s. Apply (5) to the converse of B with t=c_p/B(p), again in [0,1], to obtain (7). Attainability uses the last endpoint of the chain for A and the first endpoint of its converse counterpart for B. Adding v between them with only b->a reversed is exactly an interval chain with one marked triangle block between the two side chains.

**Theorem 6 (single-triangle sector, all orders).** For every n>=4, the minimum Delta(T,v) over strong n-tournaments whose root lies in exactly one triangle equals the interval-chain minimum g(n) proved in S191. Explicitly,

    n        4  5  6   7   8
    g(n)     2  4  6  10  16,

and, for every k>=0,

    g(9+3k)=26*5^k,
    g(10+3k)=46*5^k,
    g(11+3k)=78*5^k.                                   (8)

For two strong sides, Corollary 5 replaces the sides by endpoint chains of the same sizes without increasing Delta; their assembly is an interval chain of the same order. Thus Delta>=g(n).

If one side is weak, S191's proved exact floor is 2B(n-3) for n>=6. It equals g(n) at n=6,7 and exceeds g(n) for every n>=8, directly from the residue formulas (with n=8 checked as 18>16). At n=4 the only side sizes are 1+2, giving Delta>=2; at n=5 the possibilities 1+3 or 2+2 give Delta>=4, using the exact source-root weak minima r(2)=2,w(3)=4 and the C3 pair (3,1). These bounds are attained.

Conversely, every marked interior vertex in the S191 chain minimizers lies in exactly one directed triangle: a triangle through such an interior vertex must use the sole reversed edge of its own block, and its two endpoints then determine it. Hence every displayed g(n) is attained within this sector. QED.

For n>=4 a single-triangle root has strong deletion by the inherited normal form, so the sector lies within strong deletion pairs. We have not proved that a global minimum lies in this sector. The remaining unrestricted question is to exclude smaller increments at roots lying in at least two triangles. The inherited greedy arc-reversal counterexample still blocks that proposed reduction.

## 5. Why the methods now compose, and their limit

The scalar Busch theorem handles unmarked strong factors. The explicit tagged merge injection handles the missing rooted response, but only away from the root. The root-avoiding decomposition supplies precisely that condition; its two exceptions are handled by direct binary path counting. Convex combination gives the two-coordinate support function. The one-triangle structural normal form then converts that linear envelope into the exact bilinear deletion minimum. Finally Moon chain arithmetic evaluates the attained all-order result.

The capacity bridge is exact: apply S192's adjacent ordered-join formula to T->b, with b a singleton, to obtain c_(a,b)=R_a(T)+e_a(T). Thus Theorem 3 is also a sharp all-order floor for these particular boundary capacities. It is not a lower bound of the same form for every arbitrary internal capacity, and it is not a claim that a count alone determines all cuts.

The matching/real-rooted interval-cycle polynomial describes the attaining chains and their deletion cards; it is not needed to force arbitrary tournaments into that family. The new rooted envelope is the missing inequality that permits replacement only inside the proven one-triangle normal form. No spectral real-rootedness claim is exported to arbitrary tournaments.

## 6. Exact controls and reproduction

The standalone `verify.py` checks all 1,099 labelled tournaments through order five by literal Hamilton words and independent subset DP; verifies the tagged merge injection on 276,480 objects across all 61,440 relevant order-five marked gluing presentations, including weak factors; and independently evaluates the inherited complete 7,412-class bank through order eight. Every one of the 50,697 nonexceptional rooted strong cases at orders five through eight has a root-avoiding split; the 26 nearly-transitive and 52 transitive-root-deletion exceptions are counted separately. The observed sharp E minima are 2,4,6,10,18,30,50 at orders 1,3,4,5,6,7,8.

First run `g++ -O3 -std=c++17 generate_bank.cpp -o generate_bank` and `./generate_bank > inherited_tournaments_through8.txt` to regenerate the complete bank by vertex extension and degree-cell canonicalization. Then run `python3 verify.py > receipt.json`, `python3 -O verify.py > receipt.optimized.json`, and compare the files. Both runs passed and the outputs are byte-identical. All checks use explicit exceptions. The generator reconstructs the complete isomorphism-class bank; the independent audit also checked the structural induction without relying on this finite bank. These finite controls are not the proof of the all-order results.
