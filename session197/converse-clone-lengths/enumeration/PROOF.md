# Exact orbit enumeration and sampling bias in the persistent binary family

Session197. This is a separate corollary of the reviewed converse-clone-length classification, not an edit to that checkpoint. It distinguishes ordinary unlabelled classes, automorphism mass, rooted classes and Hamilton-normalized self-converse presentations. All formulas below concern this named family, not all tournaments or all self-converse tournaments. No literature-priority claim is made.

## 1. The classified physical family

The reviewed classification says that every fixed-root positive binary unit clock H_m=A*2^m+1 with infinitely many self-converse clone outputs has exactly the form

    T_(l,d)=Tr_d -> C3[Tr_l,K1,K1] -> Tr_d,
    l>=1, d>=0,

with its root at an arbitrary position in the l-chain. Its order is n=l+2+2d and amplitude A=2^(l-1). For l=1 any of the three triangle roots gives the same rooted class. For l>1 the l chain positions are pairwise nonisomorphic as marked roots. The graph has H=2^l+1; its full automorphism group has order3 if l=1 and order1 otherwise. Every rooted stabilizer is trivial.

Different pairs(l,d) give different unmarked classes: the unique nontrivial SCC has order l+2 and the two padding lengths are d. The quoted classification is complete, so the counts below omit no persistent self-converse binary unit-clock source. The empty and one-/two-vertex tournaments are not in this nonconstant family.

## 2. Ordinary and rooted class generating functions

Let u_n count unmarked classes, r_n rooted isomorphism classes, c_n rooted classes modulo converse, and f_n rooted classes fixed by converse. All four sequences vanish at n<3. Direct summation gives

    U(z)=z^3/[(1-z)(1-z^2)],
    R(z)=z^3/[(1-z)^2(1-z^2)],
    C(z)=z^3(1+z)/(1-z^2)^3,
    F(z)=z^3/(1-z^2)^2.                              (1)

Indeed each(l,d) contributes1 to u and l to r. Converse reverses the l root positions, giving ceil(l/2) rooted converse-orbits and one fixed rooted class exactly when l is odd. This gives C and F, and also R+F=2C coefficientwise.

For n>=3,

    u_n=floor((n-1)/2),
    r_n=floor((n-1)^2/4).

At n=2s+1 the rooted converse count is s(s+1)/2 and the fixed rooted count is s. At n=2s they are s(s-1)/2 and0.

The ordinary automorphism mass, summing1/|Aut(T)| over the unmarked family classes, is

    M_n=floor((n-1)/2)-(2/3)[n is odd],  n>=3,        (2)

and is0 at n<3. Its generating function is

    M(z)=U(z)-(2/3)z^3/(1-z^2).

Only the l=1 class at odd n needs the correction from weight1 to weight1/3. Consequently the number of labelled unmarked tournaments in this family is n! M_n. The number of labelled rooted pairs is n! r_n, or (n-1)! r_n if the root label is prescribed. These are different quantities; ordinary unlabelled enumeration such as A000568 counts classes with weight1, not weight1/|Aut| or by Hamilton paths.

## 3. Fixed Hamilton words and exact normalization fibres

Label T_(l,d) in its natural order0,...,n-1. All arcs go forward except the edge from the last vertex of the core to its first. The reflection tau(i)=n-1-i reverses every arc.

Every Hamilton word consists of the forced initial padding, a Hamilton word of the core, and the forced final padding. Within the core there is exactly one word avoiding the backward edge: its natural increasing order. It is fixed by the operation on Hamilton words that reverses the word and then applies tau.

Every other Hamilton word uses the backward edge. It assigns each of the l internal chain vertices either to an increasing prefix before that edge or to an increasing suffix after it. All2^l assignments are valid. If such a word is fixed, the unique backward edge must occupy its central pair, because tau reverses that edge setwise.

If l is odd, the central internal chain vertex is fixed by tau. It cannot be assigned to either side of the central edge in a fixed word. Thus no backward-edge word is fixed. If l=2t, each of the t mirrored internal pairs puts exactly one member in the first half; these2^t choices are all valid and distinct. Therefore the fixed-word count is exactly

    F_l=1                    if l is odd,
    F_l=2^(l/2)+1            if l is even.             (3)

Padding does not alter this count. The centralizer of tau in Aut(T_(l,d)) is trivial: for l>1 the full group is trivial, and at l=1 the three triangle rotations are inverted by tau, so only the identity commutes.

The reviewed converse-mass theorem identifies the number of distinct Hamilton-normalized self-converse presentations of an unmarked class as its fixed Hamilton-word count divided by this centralizer order. Hence the class T_(l,d) has precisely F_l such presentations. This is a per-class calculation; it is not inferred from positivity of a total mass.

Let S_n sum these normalization fibres over the family. Then

    S_(2s+1)=s,                 s>=1,
    S_(2s)=2^s+s-3,             s>=2,                 (4)

with S_n=0 for n<3. Equivalently

    S(z)=z^3/(1-z^2)^2
       +z^4(3-4z^2)/[(1-2z^2)(1-z^2)^2].            (5)

For example, at order14 there are six unmarked classes but132 normalized self-converse presentations. At order13 there are also six classes, each having exactly one such presentation.

The full normalized self-converse universe at order n has2^floor((n-1)^2/4) presentations by the named mass theorem. Thus S_n divided by that number is the exact probability that a uniformly sampled normalized presentation represents this family. This does not describe uniform sampling of unlabelled classes.

## 4. A concrete all-order discrepancy between two sampling measures

Condition on belonging to this family. Its strong class is exactly d=0.

At even order n=2s>=4, uniform sampling of unmarked classes gives strong probability

    1/(s-1),

whereas uniform sampling of Hamilton-normalized self-converse presentations gives

    (2^(s-1)+1)/(2^s+s-3).                           (6)

The first tends to0 and the second tends to1/2 as s grows. At odd order2s+1, both conditional probabilities are1/s, because every class's normalization fibre equals1.

This is an exact physical family in which Hamilton normalization produces a persistent parity-sensitive class bias. It is a counterexample to treating that finite normalization as an unbiased sampler of unmarked tournament classes. It makes no claim that either sampling measure is preferable for every purpose, and no claim of prime infinitude or unrestricted prime/7p coverage.

## 5. Reproduction and dependencies

The standalone C++17 verifier constructs every T_(l,d) of orders3 through14:42 physical graphs. It enumerates all degree-compatible automorphisms, checks the reflection and its exact centralizer, counts rooted/converse root orbits, and enumerates all10950 actual Hamilton words by adjacency DFS. It checks all282 fixed words, and independently forms the complete sets of normalized adjacency matrices from them. It verifies all displayed class counts, automorphism masses and normalization fibres in this range, including the nontrivial triangle automorphism case. Optimized and assertion-enabled unoptimized outputs are identical.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug

The all-order derivations above are independent of the finite control range.

Named reviewed dependencies:

- Complete physical classification and stabilizers: session197/converse-clone-lengths/PROOF.md, published at5f51c4be6cf617bb98abe4f4173468d2a937f298, SHA256c3516df7522a100e692665f2d33e2427ff76a1dc83d0e873f072bbff6c66d890.
- Exact normalized self-converse fibre and total-universe law: session197-orbit-wildcard/CONVERSE_MASS.md, SHA25607f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620.

The verifier has no input-file dependencies.
