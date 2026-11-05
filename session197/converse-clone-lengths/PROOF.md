# Finite or universal self-converse clone lengths, and eventual unit-clock rigidity

Session197. These are all-order statements about repeatedly replacing one fixed root by a nonempty transitive module. They do not classify arbitrary graph-changing constructions. The module arguments below are elementary and self-contained; the unit-clock consequences explicitly use previously reviewed algebraic classifications. No literature-priority claim is made.

## 1. The canonical transitive-block partition

A module M of a tournament T is a nonempty vertex set such that every exterior vertex either beats all of M or loses to all of M. A transitive module induces an acyclic tournament. Maximality below is by inclusion, not by cardinality.

**Lemma.** The maximal transitive modules partition V(T).

Every vertex lies in a singleton transitive module. If transitive modules A,B intersect, their union is a module: an exterior vertex has uniform incidence to each, and their intersection forces the two signs to agree. Their union is transitive as well. Any directed triangle in the union cannot be contained in either transitive module. If it uses two vertices from one module and one outside that module, the exterior vertex sees those two uniformly, so it is not a triangle. This covers all ways a triangle can meet A-B, B-A and their intersection. A tournament without a directed triangle is transitive. Thus intersecting transitive modules have a transitive modular union, and distinct maximal ones are disjoint. QED.

Let their quotient be Q, with positive weights w_i equal to block sizes. All arcs between two blocks are uniform, so T is exactly Q[Tr_(w_1),...,Tr_(w_q)]. The partition is canonical under both isomorphisms and converse isomorphisms. A transitive chain has only the identity order-preserving automorphism and only one order-reversing bijection to a chain of the same size. Therefore

    Aut(T) = Aut(Q,w),

and converse isomorphisms of T correspond exactly to arc-reversing permutations of Q that preserve the weights, each lifted by reversing the relevant chains. This is a bijection of actual maps, not just an equality of group orders.

If q=1, T is transitive. Otherwise q>=3: two transitive modules with a uniform arc between them have transitive union and could not be maximal.

## 2. Cloning changes just one weight

Fix a root r in a nonempty tournament T of order n, and let M be its maximal transitive module, of size l. Write T_m for replacement of r by Tr_m, m>=1. The replacement chain occurs at r's old position within M, so the expanded block M_m is transitive of size

    l+m-1.

**Lemma.** The maximal transitive-module partition of T_m is exactly the old partition, with M replaced by M_m.

All displayed blocks are transitive modules. Suppose a larger transitive module contains any part of M_m. Taking its union with M_m preserves transitivity and modularity by Section1, so a maximal such module contains all of M_m. Contract the inserted chain to r. The larger set would give a transitive module of T strictly larger than M, a contradiction. A putative larger module disjoint from M_m contracts unchanged and contradicts maximality of one of the other old blocks. Thus no new merge is possible. QED.

Hence Q is unchanged and only its marked weight varies:

    w_r(m)=l+m-1.

Define the finite set of possible weight-collision lengths

    D={ w_j-l+1 : j is another quotient vertex and w_j>=l }.

Outside D, the expanded block has a weight different from every other block. This is stronger than merely requiring it to be largest.

## 3. Exact finite/all self-converse dichotomy

Let S be obtained from T by contracting the whole block M to a single marked vertex s, leaving all other vertices intact.

**Theorem.** Exactly one of the following alternatives holds:

1. S has a converse isomorphism fixing s. Then T_m is self-converse for every m>=1.
2. S has no such converse isomorphism. Then every self-converse clone length belongs to D. For each m in D, self-converseness is decided exactly by whether (Q,w(m)) has a weight-preserving arc-reversing permutation.

In the second alternative the set can be empty. If it is nonempty, T is nontransitive, and

    every such m <= n-2l <= n-2.                         (1)

Thus for n>=2, a single self-converse clone at m>=n-1 implies self-converseness at every positive clone length. At n=1 the unique source is transitive and all lengths work. Equivalently, infinitely many self-converse clone lengths imply all positive lengths.

**Proof.** A converse isomorphism of S fixing s extends to T_m by reversing the replacement chain of size l+m-1. Conversely, any converse isomorphism of T_m that preserves M_m setwise contracts to one of S fixing s. If m is outside D, the unique size of M_m forces every converse isomorphism to preserve it, so a self-converse output outside D forces alternative1. Inside D, the exact weighted quotient test follows from Section1.

If alternative2 has a successful length, q>=3. For every other block j there is at least one additional nonempty block besides M and j. Hence w_j<=n-l-1, and m=w_j-l+1<=n-2l. This proves (1). QED.

This is a finite decision procedure for the entire infinite clone-length language. Compute the canonical blocks, test one root-fixed converse condition, and, if it fails, test at most the distinct lengths in D on the fixed finite quotient. No Hamilton-count cutoff, probabilistic test or bounded-order extrapolation is involved.

A root-fixed converse isomorphism can always be chosen to be an involution. Indeed a tournament automorphism group has odd order: an involution would exchange two vertices and reverse their joining edge, impossible; Cauchy's theorem then excludes even group order. If alpha reverses arcs and fixes s, alpha^2 has odd order d, so alpha^d is an arc-reversing involution fixing s.

### Sharp cutoff and a necessary warning

For every n>=4 take

    T=C3[K1,Tr_(n-2),K1],

with r in the first singleton. Its maximal transitive blocks are precisely the three displayed blocks. The converse permutations of the directed-triangle quotient are its three transpositions. A weighted directed triangle is self-converse exactly when at least two weights coincide. Thus T_m is self-converse exactly at m=1 and m=n-2. The bound in (1), and the universal threshold n-1, are sharp at every n>=4.

A particularly useful case is n=4, with outgoing rows[6,12,8,1] and root0. For all m>=1 it has

    H(T_m)=3*2^m-1,

but it is self-converse only at m=1,2. To see the count directly, in C3[Tr_m,Tr_L,K1] every directed cycle uses the singleton and nonempty increasing subsets of both chains, so a packing has at most one cycle. Its length is odd precisely when the two chosen subset sizes have even sum. The signed nonempty-subset sum is-1 on either chain. Therefore

    H=1+(2^m-1)(2^L-1)+1
      =(2^L-1)2^m+3-2^L.

At L=2 this is the claimed negative-unit binary clock. Its first self-converse counts are5 and11; the next count23 is not self-converse. A single self-converse source, or even two small self-converse clone lengths, cannot replace the infinite-length hypothesis.

## 4. Automorphism stabilization and unit-clock mass

Let G=Aut(T)_r be the stabilizer of the original marked vertex. Every automorphism fixing r preserves its canonical block M and fixes that transitive chain pointwise. Thus G is exactly the group of quotient automorphisms fixing M and preserving all other weights. It extends to every T_m, fixing M_m pointwise. Conversely the stabilizer of any one vertex in M_m is this same group.

For every m outside D, any automorphism of T_m must preserve its unique-weight block M_m. Hence

    Aut(T_m)=G   whenever m is outside D.                (2)

At collision lengths the full group may be larger. The equality in (2) identifies actual permutation actions after the natural exterior identification; in particular it identifies their orders. For nontransitive T, every m>=n-1 is outside D by (1)'s size bound, irrespective of self-converseness.

Now suppose an eventual physical clone identity is

    H(T_m)=A b^m+epsilon,
    A>0, b>1, epsilon in{-1,1}.

The reviewed general unit-clock theorem makes b,A integers, makes the equality exact for every m>=1, and proves (b-1)! divides A. Consecutive counts are coprime: their gcd divides b-1 by subtraction, whereas each count is epsilon modulo b-1 since b-1 divides A. Every element of G extends to both clone graphs and acts freely on their Hamilton words; a permutation fixing an entire Hamilton word fixes every vertex. Therefore |G| divides two coprime counts and G is trivial.

**Corollary.** Every eventual physical unit clone clock has trivial root stabilizer. Its physical outputs are fully rigid at every noncollision length, in particular at every m>=n-1. This statement has no symmetry hypothesis.

For each such rigid output of order N=n+m-1, the one represented unlabelled tournament class contributes N! labelled tournaments and H(T_m) ordinary Hamilton-word orbits. This is an exact automorphism-orbit consequence relevant to A000568's orbit accounting, not an evaluation of the whole enumeration sequence or a claim that different clone sources cannot coincide.

## 5. Consequences for self-converse unit clocks

Suppose additionally that infinitely many T_m are self-converse. Section3 gives a root-fixed self-converse contracted source S. Writing j=m+l-1, the contracted clone sequence has the eventual identity

    H(S_j)=(A/b^(l-1)) b^j+epsilon.                       (3)

This is an eventual identity for all sufficiently large j, not an unjustified assertion at deleted or nonpositive lengths. The reviewed finite-exponential theorem subsequently extends it to every positive j.

### Negative unit offset

Apply the reviewed marked negative-base rigidity theorem to S. Then

    b=3,
    A/3^(l-1)=2q^2, q odd,
    A/3^(l-1)>=242.                                    (4)

No attainment at242 is asserted. This extends the marked theorem to any fixed-root clone family that is self-converse infinitely often, without assuming the original root is fixed by an anti-involution. The four-vertex example above shows why the infinite-length hypothesis matters.

### Complete positive binary classification without a root-fixed assumption

Suppose instead H(T_m)=A*2^m+1 eventually and infinitely many T_m are self-converse. The reviewed marked binary theorem applied to S gives amplitude4^j and a single-backward-arc strong core of order2j+3, rooted at its middle, with equal transitive padding on both sides.

In that core the2j+1 internal vertices form its maximal transitive module containing the root. But s in S has a singleton maximal transitive module, by Section2 and contraction of all of M. Therefore j=0. Consequently

    A=2^(l-1),
    T=Tr_d -> C3[Tr_l,K1,K1] -> Tr_d,                   (5)

for some d>=0, with r at an arbitrary one of the l positions in the long transitive block. Conversely every source in (5) works for every m>=1, with H=2^(m+l-1)+1, and every output is self-converse. Thus all powers of two, not only powers of four, occur as amplitudes in this unmarked persistent-symmetry formulation.

For fixed l,d, there is one unlabelled unmarked source and exactly l rooted isomorphism classes. If l>1, the source is rigid: its unique large transitive core block fixes the cyclic quotient, and every chain and every outer singleton SCC is fixed pointwise. If l=1 the only automorphisms are the three rotations of the central triangle; its three possible roots lie in one orbit. In both cases every rooted stabilizer is trivial. At order n=l+2+2d there are exactly l*n! labelled rooted pairs, or l*(n-1)! with a prescribed root label. The unmarked labelled count is n!/3 for l=1 and n! otherwise.

An anti-involution fixes the specified root in (5) exactly when l is odd and r is the middle of its chain. This recovers the earlier powers-of-four marked classification as a strict subcase. Every output avoids divisibility by7; if an output is prime, m+l-1 must be a power of two. Neither assertion claims infinitely many Fermat primes.

## 6. Controls and inherited scope

The standalone C++17 checker enumerates all5405 rooted labelled sources through order5. It finds maximal transitive modules by testing every subset, checks the exact partition after a two-clone replacement, and directly backtracks all degree-compatible converse maps and automorphisms on the actual clone graphs. It separately tests the fixed weighted quotient. All32110 physical symmetry checks through m=n+1 agree with the finite/all criterion;64220 full/stabilizer automorphism checks agree with (2). There are1853 persistent sources and3552 finite cases, whose observed finite-length cardinalities are0:1664,1:1600,2:288. These small observed cardinalities are not asserted as an all-order bound of two.

The verifier also checks81 physical instances of the sharp family,105 actual positive binary clones with arbitrary root position, and eight actual negative-binary hostile clones with direct Hamilton counts. Optimized and assertion-enabled unoptimized builds have identical output. The structural proofs above, not this finite inventory, establish the all-order statements.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug

Named reviewed dependencies:

- General eventual-to-exact unit-clock, factorial and gcd reduction: published quaternary-obstruction PROOF.md, SHA256f0bf7a6bfb64a13b6b606bb346b7b28b64395b5d18ae56c4022b36adce59289c, at commitf633072db124b398e2563b9e3ba7b12e4fc2c7ec, path session197/quaternary-obstruction/PROOF.md.
- Marked binary classification: session197-orbit-wildcard/binary-rigidity/PROOF.md, SHA25665102468958b198f8538022288238c40d4e3ca33684cf0d419bc86b0993681a9.
- Marked negative-base rigidity: session197-prime-atoms/NEGATIVE_BASE_RIGIDITY.md, SHA256bea7513641935a9ad0b5af15551dac626ec981667bddc62dbae054bed3339b9d, which explicitly cites the reviewed square/amplitude242 theorem.

The finite/all clone-length decision, sharp threshold, stabilized automorphism action and the removal of the original root-fixed hypothesis are derived here. No assertion is made about arbitrary isolated self-converse sources beyond the exact finite candidate test.
