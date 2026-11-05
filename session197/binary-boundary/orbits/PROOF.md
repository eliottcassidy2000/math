# Converse orbits of the quadratic boundary family

Session 197. **PROVED exact fixed-word count and rigid-family enumeration.** This is separate from the physical binary-boundary packet. It connects those prime/7p-capable boundary graphs to the exact weighted self-converse enumeration universe, rather than treating graph classes and path-pointed multiplicities as the same count.

## 1. A general fixed-word formula

Use the graph F_0 from `../binary-boundary/PROOF.md`: O is a marked tournament (S,a), I is its converse, all cross arcs point I->O except O_a->I_a. The designated anti-involution tau exchanges matching I/O labels. As there, set

    h=H(S), e=H_start_at_a(S),
    q=sum_(nonempty paths P ending at a) H(S-P).

Let f_tau be the number of Hamilton words fixed by reversal followed by tau. Then

    f_tau(F_0)=h+q-e.                                  (1)

This applies to every nonempty side S, including t=1. It refers to this designated anti-involution; it makes no general claim about the full automorphism group of F_0.

A fixed word is determined by its first half, which contains exactly one vertex from each tau-pair. The reversed marked arc O_a->I_a is itself fixed by reversal followed by tau. If it occurs in a fixed word, its position must therefore be central: any noncentral occurrence would require the same ordered pair of vertices at a second position, impossible in a Hamilton word.

If the central arc is O_a->I_a, the first half ends at O_a. Outside that arc no O-to-I transition exists, so the half-word consists of an optional I prefix followed by a nonempty O path ending at O_a. The I-prefix labels are exactly the complement of the O path's labels. Its directed orders correspond by converse to Hamilton paths on that complement. These are exactly q choices, and reflection determines a valid full word. Every cross join outside the centre is I->O and cannot be the missing I_a->O_a arc.

Otherwise the central arc must be I_b->O_b for b!=a. The first half ends in I and cannot contain any O vertex, since returning from O to I would require the already excluded marked edge outside the centre. It is therefore a full Hamilton path in I whose last label is not a. There are h-e such paths, by converse. Conversely each produces a valid fixed word. The two cases are disjoint and exhaustive, proving (1). At t=1 only the first case remains, giving one fixed word.

In particular f_tau is positive and odd: endpoint reciprocity gives q>=e and q=e modulo2, while h is positive and odd. These parity facts are also visible as consistency checks with the general anti-involution theorem; they are not used to infer the separate class weights below.

## 2. Rigidity of every quadratic-family graph

For k>=1 put t=k+2 and take S_k to be the transitive order 0,...,t-1 with only t-1->0 reversed, marking 0. The binary-boundary packet proves that its graph Q_k=F_0 has

    n=2t=2k+4, h=q=2^k+1, e=1,
    H(Q_k)=2h^2-1, C(Q_k)=1.

The outgoing degrees of its O vertices are

    deg(O_0)=t-1,
    deg(O_i)=t-1-i for 1<=i<=t-2,
    deg(O_(t-1))=1.

The paired I degrees are n-1 minus these. Thus every degree from 2 through n-3 occurs once, while degree1 and degree n-2 each occur twice. An automorphism must fix all unique-degree vertices. Within either two-element tied class, swapping its vertices would reverse their mutual directed edge, so each vertex is fixed there as well. Hence

    Aut(Q_k)=1                                       (2)

for every k>=1. The known paired anti-involution is consequently the unique anti-isomorphism of Q_k to its converse: composing any other with it would give a nonidentity automorphism. No numerical group census is needed for this all-order conclusion.

At every even order n>=6 this family contributes exactly one ordinary unlabelled tournament class and n! labelled tournaments. This is a statement about the specified family, not a classification of all rigid self-converse C=1 tournaments.

## 3. Exact path fibres and Burnside counts

Substituting h=q and e=1 into (1) gives

    f_tau(Q_k)=2h-1=2^(k+1)+1.                        (3)

The ordinary automorphism action on Hamilton words is trivial by (2), so its path-orbit mass is H(Q_k)=2h^2-1. The involution J_tau on those words has exactly 2h-1 fixed words, and therefore

    #orbits under J_tau = (H+f_tau)/2 = h^2+h-1,
    #nonfixed pairs = (H-f_tau)/2 = h^2-h.             (4)

For the 14-vertex 7p example k=5, H=2177=7*311, f_tau=65, and the converse-orbit count is1121. One explicit fixed Hamilton word in the O-then-I row labeling is

    1,2,3,4,5,6,0,7,13,12,11,10,9,8.

The independently reviewed self-converse normalization theorem assigns a class the weight f_tau/|C_Aut(tau)|. Here the centralizer is trivial, so Q_k contributes exactly 2^(k+1)+1 normalized presentations. The full normalized universe at order n=2t has size 2^(t(t-1)); this single class occupies the exact fraction

    (2^(t-1)+1) / 2^(t(t-1)).                         (5)

The denominator counts normalized Hamilton-word-pointed self-converse objects, not uniformly sampled ordinary A000568 classes. For the 14-vertex example the fraction is65/2^42, despite there being only one ordinary class and2177 unrestricted Hamilton-word presentations. These are different, explicitly identified enumerations.

This formula does not imply that the family covers arbitrary prime or7p targets. It gives an exact all-order relationship between a constructive arithmetic family, its ordinary class count, its automorphism mass and its converse-fixed fibres.

## 4. Dependencies and reproducibility

The physical graph/counting dependency is `../binary-boundary/PROOF.md`, SHA256 `a6c358fc78204c266d05ad04d458877abe2aa980e25403740511d1af2ab0abf7`. The positive odd fibre and normalized total are inherited from `../CONVERSE_MASS.md`, SHA256 `07f4b20dde9f9d23b7711f65bd484d7468bdbef3cc6da0cc2888aa91d402a620`, whose normalization and centralizer proof was independently reviewed. The general fixed-word decomposition and degree-based rigidity are proved here. No literature-priority claim is made.

Compile and run the standalone C++17 `verify.cpp`. It enumerates all pair-transversal first halves and all their orders, independently checking the complete fixed-word formula for every rooted labelled side through order five plus twenty deterministic order-six roots. It separately enumerates all degree-compatible automorphisms and anti-isomorphisms of six family graphs, checks their fixed words and Burnside values, and checks the displayed14-vertex fixed word. Optimized and unoptimized assertion-enabled receipts are compared. All-order claims follow from Sections1–3, not from these finite checks.
