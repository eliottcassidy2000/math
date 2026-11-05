# The first symmetric 7p stratum is 4151

Session 196, 2026-11-05. **PROVED all-order classification, with separately recorded finite certificates.** No claim of literature priority is made. The full Hamiltonian-path spectrum remains open.

## 1. Statement and inherited boundary

Write H(T) for the number of directed Hamiltonian vertex words of a finite tournament T, and A(T)=|Aut(T)|. In this note, **automorphism-free** means A(T)=1; it does not mean the different property H(T)=A(T), sometimes called path-rigidity elsewhere in the repository.

**Theorem.** Among all tournaments for which H(T)=7p with p prime, the least Hamiltonian-path count admitting a nonidentity automorphism is

    4151 = 7 * 593.

It is attained at order 10. More precisely:

1. For every prime p<593, every tournament with H(T)=7p is automorphism-free. This includes empty strata, and is an all-order assertion, not an order-bounded census.
2. There are exactly two strongly connected, nontrivially symmetric isomorphism classes with H=4151. Both have order 10 and A=7. They are converses of one another and neither is self-converse.
3. For every n>=10, exactly 2(n-9) ordinary tournament classes of order n have H=4151 and A>1. Every such class has A=7. There are n-9 converse-merged classes of this kind and no self-converse ones. At n<10 there are none.
4. Order 10 is also the least order of any nontrivially symmetric tournament with H=7p, for any prime p.

The two strong cores are obtained by substituting the seven-vertex carousel C into either of two marked positions of the unique strong four-tournament Q. We use outgoing bitmask rows, vertex indices starting at zero:

    Q = [2,8,3,5],
    C: i -> i+1,i+2,i+3 modulo 7.

The marked positions are 0 and 2. One resulting ten-vertex witness is

    [142,156,184,240,225,195,135,512,255,383].

The known all-order omissions H=7 and H=21 are recorded in [THM-1370, h-spectrum-omits-7-21-all-n](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md). They remain part of the inherited background, and that source retracts an earlier novelty claim. For dependency independence, Appendix A supplies an elementary H=7 exclusion and proves that any hypothetical H=21 tournament would be automorphism-free. Thus the present symmetry theorem needs neither inherited omission as a black-box input, and does not assume that every other odd value occurs. All structural, symmetry, and counting steps used here are explicit.

The inherited finite-class multiplicity source is [Session 24, Actual Hamiltonian-path multiplicities](../../05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/PROOF.md). Its SCC rationality and eventual-polynomial architecture are not new here. The present new refinement specifies an entire automorphism stratum at the sharp target 4151, rather than asserting a new general SCC formula. In particular, this note does **not** determine the total number of classes with H=4151.

## 2. Three elementary path lemmas

**Lemma 1 (induced monotonicity).** If U is induced in T, then H(U)<=H(T).

Insert each exterior vertex z into a directed word just before its first vertex beaten by z, or append z if none is beaten. At an internal insertion the preceding vertex beats z, so the word stays directed. Fix the insertion order. Removing the added labels recovers the original word, proving injectivity.

**Lemma 2 (shuffle supermultiplicativity).** For disjoint vertex sets X,Y in a tournament,

    H(T[X union Y]) >= H(T[X]) H(T[Y]).                 (1)

This requires no uniformity of the cross arcs. Given directed words P on X and R on Y, repeatedly compare the two current first vertices and output the one that beats the other. Remove the chosen first vertex and repeat; when one word empties, append the other. The result is directed. If consecutive output vertices come from the same input word, their arc is one of that word's path arcs. If they come from different words, the later vertex was the unchanged competing first vertex when the earlier vertex was selected, so the earlier one beats it. Restriction of the output to X and Y recovers P and R. Thus different input pairs cannot collide. Iterating proves the corresponding product bound for any number of disjoint induced pieces.

This is a lower-bound injection, not a formula. Cross arcs can create many additional paths. It is also not an assertion that H(T[X]) divides H(T).

**Lemma 3 (strictness inside strong tournaments).** If U is a nontrivial strongly connected induced tournament properly contained in a strongly connected T, then

    H(U) < H(T).                                      (2)

Here are elementary details, including the only cycle fact needed. Every nontrivial strong tournament has a Hamiltonian cycle. For a longest directed cycle C, an exterior vertex mixed on C can be inserted at an in/out transition along the cycle. If there is no mixed exterior vertex, the exterior splits into A beaten by C and B beating C. Strongness forces both parts and an arc a->b with a in A,b in B. Inserting a,b consecutively in C makes a longer cycle. Hence a longest cycle spans.

A mixed exterior vertex z makes U+z strong. Let X be the vertices of U beaten by z, and Y those beating z. A Hamiltonian cycle of U has a Y->X transition. Cutting there gives a Hamiltonian word P starting in X and ending in Y. Both zP and Pz are valid extensions of P. Use one insertion for every other U-path. Deleting z separates the different old paths, so H(U+z)>=H(U)+1.

If there is no mixed exterior vertex on U, then U is a module: every exterior vertex is uniform on U. Strongness of T gives exterior a,b with

    U -> a -> b -> U.

The induced cyclic extension has at least 3H(U) paths, by the three cyclic placements of an intact U-path. In either case there is an induced enlargement with strictly more paths. Lemma 1 then proves (2). In particular, equality of path counts cannot hide additional vertices in a strong supertournament.

## 3. Automorphism orbits and the marked strong-four lemma

Every automorphism acts freely on Hamiltonian vertex words: fixing a word fixes every vertex in its ordered position. Consequently

    A(T) divides H(T).                                (3)

Let an automorphism g have odd prime order q. Its vertex orbits have sizes 1 or q. On a q-orbit O, label vertices by successive powers of g. The induced tournament is a circulant with some skew connection set D in Z/qZ: for each pair {d,-d}, exactly one is in D. Thus there are exactly 2^((q-1)/2) labelled connection sets to consider. The induced q-orbit is regular and strong. If there is exactly one nontrivial orbit, every other vertex is fixed by g and therefore has uniform adjacency to O. So O is a module.

The following small marked-subtournament fact controls every large quotient.

**Lemma 4.** Every strong tournament Q of order at least four has an induced strong four-tournament containing any specified vertex r.

First r lies on a directed triangle. Take a shortest directed cycle through r. If it had length at least four, the chord between r and the second following vertex either gives a shorter cycle through r or gives a triangle through r. Either way a triangle through r exists.

Let S be such a triangle. If an exterior vertex is mixed on S, adding it yields a strong four-set. Otherwise S is a module. Strongness of Q supplies a two-vertex ear S->a->b->S as in Lemma 3. Choose any two vertices x,y of S including r, ordered x->y. Then x->y->a->b->x is a directed cycle on four vertices, so the induced four-tournament is strong. This proves the lemma without an assumption of single-vertex growth from the originally chosen triangle.

If a strong module O is contracted in a strong tournament T, its quotient Q is strong. Indeed every directed T-path projects to a directed walk after repeated module labels are suppressed. Conversely, substituting strong blocks in a strong quotient yields a strong tournament. An induced marked four-set in Q therefore lifts to an induced strong subtournament consisting of O and three fixed vertices.

## 4. The two small orbit inventories and rooted substitution table

Define F2(U) to count **ordered** pairs of nonempty directed paths on disjoint sets partitioning U. Equivalently,

    F2(U) = sum_(empty != S proper-subset V(U)) H(U[S]) H(U-S).

For the cyclic two-vertex extension C(U)=U->a->b->U,

    H(C(U)) = 3H(U)+F2(U).                            (4)

If U appears in one run, there are three block orders. If it appears in two runs, the word is P,a,b,R and contributes one ordered two-path cover. There is no third run. This is the exact formula from the Session 195 strong-structure packet, reproved here for convenience.

The following finite values are certified from all skew connection sets, using the explicit subset recurrence in Section 8:

| prime orbit size | connection-set presentations | H(O) | F2(O) | H(C(O)) |
|---|---:|---:|---:|---:|
| 5 | 4 | 15 | 90 | 135 |
| 7, carousel type | 6 | 175 | 1638 | 2163 |
| 7, Paley type | 2 | 189 | 1638 | 2205 |

For q=7, the six carousel presentations are precisely the multiplicative translates of {1,2,3}, and the other two are {1,2,4} and its negative. Thus each row represents one ordinary isomorphism type. This finite classification is about prime-order orbits, not about all regular seven-tournaments.

There is exactly one ordinary strong four-tournament Q, represented by [2,8,3,5], and A(Q)=1. Substituting an arbitrary nonempty tournament U into a marked vertex of Q gives:

| marked vertex in this Q | H(Q[U]) |
|---|---|
| 0 or 2 | 5H(U)+2F2(U) |
| 1 or 3 | 5H(U)+3F2(U) |

For clarity, these are universal profile formulas, not interpolation from two test values. A Hamiltonian word can have at most two U-runs: any separator from one U-run to the next must start at a fixed out-neighbor of the marked vertex and finish at a distinct fixed in-neighbor, requiring at least two of the three exterior vertices. Two separators would require four exterior vertices. The one-run words correspond to the five Hamiltonian words of Q. For two runs, enumerate the orders of the three exterior labels and their placement before, between, and after the two U-runs, respecting Q's arcs. There are respectively two or three such patterns. Each pattern accepts every ordered two-path cover once. The verifier independently enumerates all 24 labelled strong four-tournaments with marked vertex zero and all eight seven-circulant presentations, and checks all 192 actual substitutions.

It follows that the exact possible counts for substituting a seven-orbit into a marked strong four-tournament are

    4151, 4221, 5789, 5859,                           (5)

with minimum 4151 attained only by the carousel type at marked positions 0 or 2.

## 5. Excluding the other prime divisor

Suppose H(T)=7p and p is prime. First p=2 cannot yield nontrivial symmetry, even without using the oddness of H: an involution would swap two vertices and reverse their arc, so no tournament automorphism has order two; an order-seven element already forces at least 175 paths. Cauchy's theorem and A(T)|14 exclude every nontrivial automorphism in this case.

Now take p odd. By the H=7 exclusion proved in Appendix A, the ordered strong-component factorization of H has exactly one nontrivial strong component S, and H(S)=7p. All other components are singletons. Automorphisms fix the ordered components, so A(T)=A(S). The case p=3 cannot have nontrivial symmetry by Appendix A, without needing the stronger known absence of H=21. We henceforth work with S and p>=5.

If p>=17, an automorphism of order p would produce a p-orbit O. For each d in its skew connection set D, the constant-step cycle

    0,d,2d,...,(p-1)d,0

is a Hamiltonian cycle. Different d give different directed cyclic successors, so these are (p-1)/2 distinct cycles, whose cuts are distinct Hamiltonian paths. Therefore

    H(S) >= H(O) >= p(p-1)/2 >= 8p > 7p,

a contradiction.

For p=11,13 the finite certificate file supplies **eight distinct Hamiltonian paths starting at zero for every skew connection set**, namely all 32+64 presentations. Their p translates have different starts and are distinct, again giving H(O)>=8p>7p. These are small explicit lower-bound witnesses, not a statistical count or an appeal to asymptotics. Thus an automorphism of order p is impossible for every prime p>=11. The independent session orbit-cubic inequality q(q^2+11)/12 supplies a shorter alternative proof, but is not needed as a dependency here.

At p=7, H(S)=49 is smaller than the 175 paths in any seven-orbit. At p=5, H(S)=35. An order-five automorphism cannot have two nontrivial orbits, by (1), since 15^2=225>35. With one orbit O it is a module. If S=O, H(S)=15; otherwise the strong quotient has at least three vertices and contains a triangle through the module vertex. Lifting it gives C(O), with H=135>35. Thus order-five symmetry is impossible too. Order-seven symmetry at H=35 is already excluded by its orbit minimum 175.

By Cauchy's theorem and (3), for every remaining possible symmetric S with H=7p, p>=11, there must be an automorphism of order seven.

## 6. Sharp minimum and equality

Assume H(S)=7p<=4151 and S has an order-seven automorphism g. Two nontrivial orbits would give

    H(S) >= 175^2 = 30625 > 4151

by the shuffle injection. Hence there is exactly one seven-orbit O, which is a module; the other vertices are fixed. Contract O to form the strong quotient Q.

- If |Q|=1, the counts are 175=7*25 or 189=7*27, neither seven times a prime.
- A strong quotient cannot have order two.
- If |Q|=3, it is a directed triangle. The counts are 2163=7*309 and 2205=7*315, again neither seven times a prime.
- If |Q|>=4, Lemma 4 supplies an induced marked strong four-tournament. Its lift is strong and has at least 4151 paths by (5). Therefore H(S)>=4151.

This proves the lower bound. At equality, the lift must have count 4151. Lemma 3 excludes any additional vertices of S. Thus S is exactly one of the two carousel substitutions at Q's marked positions 0 or 2.

The explicit rows in Section 1 have exactly 4151 Hamiltonian paths by both subset DP and independent directed-word DFS. The cyclic translation of its first seven vertices is an order-seven automorphism. The integer 593 is prime (trial division through 23 suffices). Equation (3) gives A(S)|7*593; since A(S) divides 10!, the factor 593 is impossible. Hence A(S)=7.

The two cores are not isomorphic. Each has a unique nontrivial orbit of its automorphism group: its seven-vertex module. Any isomorphism must preserve that module and induce a marked isomorphism of Q. But Q has no nonidentity automorphism, and its marked positions 0 and 2 are distinct. The map (0 2)(1 3) is an anti-automorphism of Q, and negation is an anti-automorphism of the carousel. Together they show that the two cores are converses. The marked-root argument also proves neither core is self-converse.

For the minimum-order assertion without a count restriction, a symmetric H=7p tournament cannot have an automorphism of order p by Section 5. A seven-symmetry needs at least seven vertices. On at most nine vertices there is at most one nontrivial orbit, and its strong quotient has order at most three. The same list 175,189,2163,2205 excludes H=7p. Thus ten is also the smallest possible order.

## 7. Exact contribution to A000568 and labelled counts

The strong-component quotient of every tournament is a uniquely ordered transitive word. With H=4151, the impossibility of H=7 means there is exactly one nontrivial strong component. For a symmetric tournament it is one of the two cores classified above. The remaining n-10 singleton components split into a run before the core and a run after it, in n-9 ways. All choices give different ordinary classes, by uniqueness of the ordered SCC decomposition. This proves

    N_n(H=4151, A>1) = 2(n-9)   for n>=10.            (6)

Every class counted has A=7, and none is self-converse. This is an exact all-order automorphism refinement of A000568, not a formula for A000568 itself or for the full H=4151 fibre.

Let N_n be the total number of ordinary classes with H=4151 and let L_n be the corresponding number of labelled tournaments. The other N_n-2(n-9) classes are automorphism-free. Consequently, for n>=10,

    L_n/n! = N_n - 12(n-9)/7.                         (7)

Likewise, the number of isomorphism classes of a tournament with one distinguished Hamiltonian path, restricted to this H fibre, is

    sum_[T]:H=4151 H(T)/A(T)
      = 4151 N_n - 7116(n-9).                         (8)

Thus the missing symmetry correction between ordinary-class, labelled, and path-pointed enumerators is known exactly at this target even though the total fibre size N_n is not computed here. Each symmetric core has path-orbit quotient H/A=593.

**Connection contract.** Source: ordinary tournament classes and their automorphism actions. Target: the exact H=7p arithmetic slice and its A000568 multiplicity correction. Map: take a prime-order automorphism, retain its circulant orbit(s), then contract the unique orbit module when the shuffle bound permits. Preserved data: raw H inequalities, exact prime divisibility, strongness, and marked quotient incidences. Lost data in contraction: the module's internal path counts and reentry paths. Restoring sidecar: the pair (H,F2), which is enough for the three- and four-vertex quotients used here. Hostile: treating H(U) as a divisor of H(T) would be false; the proof uses only a shuffle lower bound and the exact small-quotient formulas. The all-order threshold rests on those structural reductions, not finite-to-infinite extrapolation.

## 8. Reproduction and scope of the finite work

Run, from this directory:

    python3 verify.py
    python3 -O verify.py

The verifier uses exact integers and explicit exception gates, not assert. The subset recurrence is

    E({v},v)=1,
    E(S,v)=sum_(u in S-{v}, u->v) E(S-{v},u),
    H(S)=sum_(v in S) E(S,v).

It checks this against direct directed-word DFS on every labelled tournament through order five; tests every ordered bipartition's path-shuffle injection in that universe; checks marked strong-four containment and strict induced-strong monotonicity there; exhausts all 4+8 prime circulant presentations at orders five and seven; certifies eight anchored paths for every one of the 96 circulants at orders eleven and thirteen; checks all 192 marked strong-four/seven-orbit substitutions; and independently checks the ten-vertex witness, its automorphism order, the two equality-core classes and their converse pairing.

The proof is all-order because the group orbit, shuffle, module, marked-subtournament, and SCC arguments reduce every possible counterexample to those explicit finite configurations. The verifier alone is not claimed to prove all-order coverage. Neither (6) nor the sharp symmetry threshold proves that every prime target 7p is attained.


## Appendix A. Removing the spectrum-omission dependency

These reductions were independently supplied by the session's [prime-orbit auditor](../symmetry-audit/PROOF.md), and are included here so that the sharp symmetry classification is self-contained. They are not claimed as new global spectrum omissions.

**A1. Every strong tournament of order at least five contains an induced strong five-set.** Choose a strong four-set by Lemma 4. If a mixed exterior vertex exists, add it. Otherwise the four-set is a module and admits a two-vertex cyclic ear. Choose a directed triangle within the strong four-set and keep it together with that ear; the resulting five-set is strong.

The verifier exhausts all 1024 labelled five-tournaments. Exactly 544 are strong; their H histogram is

    9:240, 11:120, 13:120, 15:64.

Thus a strong tournament of order at least five has H>=9 by induced monotonicity. At the smaller strong orders, the only counts are 1,3,5. No strong component can have H=7, and SCC multiplicativity proves the all-order H=7 exclusion. This elementary bounded reduction does not use a general strong-order floor theorem.

**A2. Any hypothetical H=21 tournament is automorphism-free.** By free action its automorphism group has order dividing 21. An order-seven element is impossible because a seven-orbit already has at least 175 paths. Suppose instead there is an element of order three. Each of its nontrivial orbits induces a C3. There can be at most two, by the shuffle bound 3^3=27>21.

With exactly one nontrivial orbit it is a module. With two such orbits, label them by the synchronized powers of the chosen automorphism. After fixing the orientation of the first triangle, there are two choices for the second triangle's orientation and eight invariant cross-arc patterns. Directly enumerating all 16 cases gives H=9 in the four uniform-cross cases and H=45 in the twelve mixed-cross cases. Since induced counts cannot exceed 21, the cross arcs must be uniform. Each triangle is consequently a module relative both to the other triangle and to every fixed vertex.

Consider the strong component containing one of these triangle modules, and contract every nontrivial orbit it contains. If its strong quotient has at least four vertices, Lemma 4 gives a strong four through the selected module vertex. Keep the whole selected C3 and one vertex from each other selected quotient block. The universal rooted-four formula gives at least 5*3+2*6=27 paths, a contradiction. If the quotient has three vertices and contains one triangle module, its path count is 15 by (4), and 15 cannot divide the ambient count 21. If that quotient contains both triangle modules, the three blockwise cyclic paths already give at least 3*3*3=27 words. A strong quotient cannot have two vertices. Thus every nontrivial orbit would have to be its own strong component.

Two C3 components force 9 to divide 21, impossible. One C3 component forces the product of the remaining component counts to be 7; their ordered union would be a tournament with H=7, ruled out by A1. Therefore no order-three symmetry exists either, proving the claimed hypothetical rigidity at H=21. This last argument intentionally does not claim to prove the known absence of every H=21 tournament.
