# Several clone lengths: exact symmetry completion and orbit cycle indices

Session197. This separate extension uses the reviewed canonical transitive-module theorem and standard Burnside orbit counting. It gives an all-order construction/obstruction for making a fixed tournament self-converse by transitive cloning, and exact generating functions for all blow-ups of one fixed quotient. It does not assert preservation of Hamilton-path counts, unrestricted prime/7p coverage, or literature priority for the orbit-counting method.

## 1. Canonical blocks and the precise language

A transitive module is a nonempty vertex set that induces a transitive tournament and is seen uniformly by each outside vertex. The reviewed clone-length theorem proves that maximal transitive modules partition every nonempty tournament. Contract them to a tournament Q, and let their positive sizes be w_1,...,w_q. Thus

    T=Q[Tr_(w_1),...,Tr_(w_q)].

Either q=1, or Q has no transitive module of size at least two. We call Q the canonical quotient here; this is weaker than full modular indecomposability. In particular, Q may be weak.

For completeness, contracting a nontrivial transitive module of Q would unite its transitive blocks to a larger transitive module of T, contradicting maximality. Conversely, replacing the vertices of such a Q by arbitrary nonempty transitive chains leaves exactly those chains as the maximal transitive modules. Indeed, any transitive module meeting two chains can be united with every chain it meets, using the reviewed fact that overlapping transitive modules have a transitive modular union. Its image would then be a forbidden nontrivial transitive module in Q.

The permitted operation replaces a vertex by a nonempty transitive chain, with identical incidences to every other vertex. Repeated operations, even at newly created vertices, only increase the sizes of these original canonical blocks. Every vector v>=w is attainable: increase one chosen vertex of each block by v_i-w_i. Strongness is preserved in both directions: a blow-up is strong precisely when Q is strong, except that the one-vertex quotient gives transitive graphs, which are strong only at order one under the usual convention.

Let G=Aut(Q) and let A be the set of anti-automorphisms of Q. Every isomorphism between two blow-ups induces an element of G carrying their weight vectors, and has a unique lift: it preserves order within every chain. Every anti-isomorphism induces an element of A and has a unique lift reversing order within every chain. Therefore

    T(v) is self-converse iff there is a in A with v_i=v_(a(i)) for all i.  (1)

The automorphism and anti-automorphism maps of the actual graph are exactly these weight-preserving lifts, not just equal in number.

For q>=2, the reviewed converse-mass proof shows that G has odd order, every anti-automorphism has order twice an odd number, and every self-converse tournament has an anti-involution. The same argument applies inside the weight stabilizer: if a preserves v, its appropriate odd power is an anti-involution preserving v. Thus in (1) it suffices to check anti-involutions of Q. The singleton quotient is handled by the identity map, which also squares to the identity.

If specified original vertices are simultaneously replaced by chains of lengths m_j>=1, the resulting block weights are exactly

    v_i=w_i+sum_(j whose vertex lies in block i)(m_j-1).                (2)

Consequently the complete set of self-converse length vectors is a finite union of explicit systems of affine integer equalities v_i=v_(tau(i)), one system for each anti-involution tau. This is a quantifier-free Presburger description. No eventual-length approximation or unproved generic symmetry assumption is used.

## 2. Minimum cloning needed to obtain self-converseness

Fix a set M of blocks allowed to grow; blocks outside M must keep their original sizes. This includes any specified set of permitted clone roots, by taking the blocks containing those roots.

For an anti-involution tau of Q, its fixed block needs no growth. For each exchanged pair {i,j}, the smallest common final size is max(w_i,w_j). It is feasible exactly when every immutable block in that pair already has this maximum size. The added order for the pair is

    2 max(w_i,w_j)-w_i-w_j = |w_i-w_j|.                              (3)

Let D_tau be the sum of (3) over exchanged pairs if all pairs are feasible, and infinity otherwise. Then the exact minimum number of added vertices is

    D(T,M)=min_tau D_tau.                                             (4)

If there is no anti-involution, or none is feasible, no descendant in this language can be self-converse. Otherwise setting each pair to its maximum gives an explicit optimal descendant and an explicit reversing involution. Grow a single permitted root in each deficient block by the required amount.

Necessity is equally important: any self-converse descendant has a weight-preserving anti-involution of Q by Section1. Its two paired block sizes must agree and cannot be smaller than their original maximum. Hence it pays at least D_tau. This proves the lower bound in (4), including arbitrary repetitions and changing clone roots. There is no restriction to a bounded number of operations.

If every block is mutable, every self-converse quotient admits a completion, and (4) is simply the minimum of the sum of pairwise absolute size differences. A non-self-converse quotient is an all-order obstruction. For Q=C3 this specializes to

    D=min(|w_1-w_2|,|w_2-w_3|,|w_3-w_1|).

For example, C3[Tr_1,Tr_2,Tr_4] becomes self-converse after one added vertex, by changing the first size to two. If only the size-four block may grow, completion is impossible: the unequal sizes one and two cannot be matched, and neither can be increased to four. The construction changes H; no Hamilton-path conservation is claimed.

## 3. Exact ordinary and self-converse class generating functions

Now fix Q and allow every positive weight vector, grouped by total order n. Let U_Q(n) be the number of ordinary unmarked isomorphism classes of its blow-ups, and S_Q(n) the number of those classes that are self-converse. Write c for a cycle of a vertex permutation and |c| for its length.

Ordinary Burnside gives

    U_Q(z)=(1/|G|) sum_(g in G) product_(c cycle of g) z^|c|/(1-z^|c|). (5)

The exact self-converse counterpart is the coset trace

    S_Q(z)=(1/|G|) sum_(a in A) product_(c cycle of a) z^|c|/(1-z^|c|). (6)

If A is empty, (6) is zero. If it is nonempty, A is a coset of G. To prove (6) without treating its objects as ordinary Burnside orbits, consider one G-orbit O of weight vectors. If O is not self-converse, no anti-automorphism fixes any of its elements. If it is self-converse, for each v in O the set of anti-automorphisms fixing v is a coset of its automorphism stabilizer G_v, and has size |G_v|. Therefore

    sum_(v in O) #{a in A:a v=v} = |O| |G_v| = |G|.

After division by |G|, this orbit contributes exactly one. A weight vector fixed by a permutation is constant on its cycles, and its positive total-weight generating function is the product displayed in (5) and (6). This proves both formulas at every order, including all nontrivial stabilizers.

The ordinary automorphism mass is different and simpler:

    M_Q(z)=sum_[T] z^|T|/|Aut(T)| = z^q/[|G|(1-z)^q].                (7)

Indeed the number of weight vectors in an orbit is |G|/|G_v|. Thus the number of labelled blow-up tournaments of order n is exactly

    n! binom(n-1,q-1)/|G|,   n>=q.                                  (8)

These counts concern all unmarked outputs of this fixed quotient, not marked clone recipes. Distinct clone recipes inside one block need not give distinct graphs.

All three series are zero at n=0. For q=1 they are z/(1-z), so there is exactly one transitive class at each positive order and every output is self-converse. Quotients of order two cannot occur canonically.

## 4. Quantitative rarity of symmetry in a fixed nontrivial quotient

Suppose q>=3 and Q is self-converse. Choose an anti-involution tau and put C=|C_G(tau)|. All anti-involutions are conjugate in G by the reviewed converse-mass theorem, so their number is |G|/C.

An anti-automorphism has at most one fixed vertex: two fixed vertices would require an arc between them to equal its reverse. Its nontrivial cycles have lengths2 modulo4, because its order is twice an odd number and an odd cycle of length greater than one would make an odd anti power act identically on a nontrivial induced tournament. Thus an anti-involution has r transpositions and epsilon fixed points, where

    q=2r+epsilon,  epsilon in {0,1}.

It has the largest possible number r+epsilon of cycles. Every non-involutory anti-automorphism has strictly fewer cycles, since it contains a cycle of length at least six. Ordinary nonidentity automorphisms also have fewer cycles than the identity, as their order is odd.

The identity term in (5) yields

    U_Q(n) ~ n^(q-1)/[|G|(q-1)!].                                  (9)

For odd q=2r+1, the leading anti-involution contribution to (6) is

    (1/C) z^q/[(1-z)(1-z^2)^r].

The pole at z=1 has order r+1 and the pole at z=-1 only order r. All other anti terms have lower maximal pole order. Therefore

    S_Q(n) ~ n^r/[C 2^r r!],
    S_Q(n)/U_Q(n) ~ |G|(2r)!/[C 2^r r!] * n^(-r).                  (10)

For even q=2r, every anti cycle is even, so S_Q(n)=0 at odd n. On even n the anti-involution contribution is

    (1/C) [z^n] z^(2r)/(1-z^2)^r = (1/C) binom(n/2-1,r-1),

and the other anti terms have smaller order. Hence, along even n,

    S_Q(n) ~ n^(r-1)/[C 2^(r-1)(r-1)!],
    S_Q(n)/U_Q(n) ~ |G|(2r-1)!/[C 2^(r-1)(r-1)!] * n^(-r).        (11)

These are uniform ordinary-class sampling probabilities inside one fixed quotient family. They are neither Hamilton-weighted probabilities nor statements about the distribution on all tournaments. The singleton quotient and non-self-converse quotients have already been treated separately.

## 5. Explicit nontrivial group example

Take the Paley tournament on Z/7Z, with i->j exactly when j-i is1,2 or4. Its outgoing masks are

    [22,44,88,49,98,69,11].

It is a strong canonical quotient. Its21 automorphisms have cycle types1^7 once,7 six times, and1*3^2 fourteen times. Its21 anti-automorphisms have cycle types1*2^3 seven times and1*6 fourteen times. The standalone exhaustive permutation/backtracking control below checks these exact assertions; alternatively the affine maps i->a i+b with a in {1,2,4} give the21 automorphisms, and their remaining nonzero multipliers give the anti coset, with the same cycle types. Exhaustion verifies that there are no additional maps.

Thus

    U_Q(z)=(z^7/21)[1/(1-z)^7 +6/(1-z^7)
                         +14/((1-z)(1-z^3)^2)],
    S_Q(z)=(z^7/3)[1/((1-z)(1-z^2)^3)
                         +2/((1-z)(1-z^6))].                     (12)

Here C=3, not one, and the self-converse fraction is asymptotic to105/n^3. This explicitly checks the centralizer factor in (10). Every output is strong because this Q is strong. These formulas enumerate physical strong tournaments but impose no assertion that their H values cover primes or7p.

## 6. Reproduction, scope and inherited premises

The standalone C++17 verifier enumerates every pair-free quotient class through order five, plus a six-vertex ribbon quotient and the seven-vertex Paley quotient. A tournament has a nontrivial transitive module exactly when it has a two-vertex module somewhere inside one: adjacent vertices in the transitive module's order form such a pair. Thus this finite seed filter is exact.

For each seed it constructs every positive-weight blow-up through order eleven, compares the complete physical automorphism and anti-automorphism maps with their unique quotient lifts, checks actual maximal transitive modules through order nine, compares direct weight-orbit counts with both cycle-index formulas, and checks simultaneous affine clone changes. It also compares the minimum-completion formula, for every subset of mutable blocks, with exhaustive enumeration of all dominating weight vectors whose coordinates lie between one and three. The graph-level map checks are independent of the formal cycle-index calculation. The proof, not the finite range, establishes all-order completeness and optimality.

Run:

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug

Both outputs must be identical. The verifier has no file dependencies and uses explicit exceptions rather than disabled assertions for its checks.

Named reviewed mathematical dependencies:

- Canonical maximal transitive-module partition, union lemma and cloning invariance: session197/converse-clone-lengths/PROOF.md, published at5f51c4be6cf617bb98abe4f4173468d2a937f298, SHA256c3516df7522a100e692665f2d33e2427ff76a1dc83d0e873f072bbff6c66d890.
- Odd automorphism order, existence and conjugacy of anti-involutions: session197-orbit-wildcard/CONVERSE_MASS.md, SHA25607f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620.

Ordinary and coset Burnside counting is proved directly here and is standard group-action machinery. The previously reviewed binary-family enumeration is a particular different family, not a claimed proof of this general fixed-quotient statement. No finite census is extrapolated to an unrestricted arithmetic spectrum.
