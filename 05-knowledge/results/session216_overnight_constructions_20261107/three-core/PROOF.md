# Three fixed cores and a binary inverse problem

**Status:** PROVED construction formula; FINITE-EXACT retained coverage certificate; PROVED elementary all-length inverse bounds and two permanent template omissions; OPEN all-length inverse conjecture. Independently audited: PASS; see audit/REVIEW.md and audit/INPUT_FREEZE.md. No unrestricted coverage-record or literature-priority claim is made.

## 1. A small physical construction language

Use one of the following strong tournament cores, in outgoing-bitmask notation:

    C0 = (8,1,3,6),          H(C0)=5;
    C1 = (2,4,9,19,7),      H(C1)=9;
    C2 = (2,8,3,21,7),      H(C2)=9.

Choose nonnegative integers a,b, a cut on the core, an a-bit word p and a b-bit word r, allowing leading zeroes. Form

    L = TT_a ⇒ C_i ⇒ TT_b.

The chosen cut S consists of the selected prefix vertices, selected core vertices, and selected suffix vertices, reading each word from its most significant bit to its least significant bit along the corresponding transitive order. Add vertices v,x with

    L→v→x,       x→u for u in S,       u→x for u in L minus S.

Call this tournament T. Require the first prefix bit to be1 when a>0; if a=0 require a nonempty core cut. This condition makes T strong: x reaches the initial source component of L, hence all of L, and every L vertex reaches v→x. This condition is sufficient, not an extra promise inferred from the count. Vertex v has exactly one outneighbor, x.

There are only80 core-cut states:16 from C0 and32 from each of C1,C2. Arbitrary padding lengths and bit words give the unbounded family. The retained finite certificate uses a+b<=16.

## 2. The exact cut coordinates

For a tournament K with a cut S, let h=H(K), H(empty)=1, and

- q=sum_(U subset V(K)) H(K[U])H(K[V(K) minus U]);
- A count ordered spanning path covers(P,Q) with Q nonempty starting in S and P optional;
- G count ordered spanning path covers(P,Q) with P nonempty ending in S and Q optional;
- T count such covers with both paths nonempty, P ending in S and Q starting in S;
- S0 count Hamilton paths of K starting in S;
- E count Hamilton paths of K ending in S.

All paths are directed. Here the symbol T for a coordinate should not be confused with the final tournament; in the displayed state table it occupies the fifth numerical column. The ordered core-cut state is(h,q,A,G,T,S0,E). `CORE_STATES.txt` gives all80 states, first core cuts0..15, then second and third core cuts0..31. `verify.py` independently recomputes every entry from ordinary vertex-subset start/end tables.

For two nonempty cut-labeled tournaments K1⇒K2, splitting each cover path into its block portions gives

    h' = h1 h2,             q' = q1 q2,
    A' = A1 q2+h1 A2,       G' = G2 q1+h2 G1,
    T' = A1 G2+T1 h2+h1 T2+E1 S02,
    S0' = S01 h2,           E' = h1 E2.              (1)

For example, the four terms in T' distinguish whether the first path uses K2 and whether the second path uses K1. When both occur, the term is A1G2; when only the second occurs it is T1h2; when only the first occurs it is h1T2; and when neither occurs the paths are entire block paths, giving E1S02. Empty padding is handled by omitting the corresponding ordinal factor.

For a transitive block of length n with cut word z,

    h=1, q=2^n, A=z, G=rev_n(z), T=f(z),
    S0=first_n(z), E=z mod2,

where rev_n reverses exactly n binary digits and first_n is the first digit. For the empty word both rev and first are0. The function f is equivalently

    f(z)=sum_(i<j, bits i,j selected in z)
             (2^(j−i−1)+indicator(j=i+1)),

or recursively

    f(0)=0,
    f(2z)=f(z),
    f(2z+1)=f(z)+z+(z mod2).                       (2)

The pair sum prescribes the first vertex of the second path and last vertex of the first path. The extra adjacent term accounts for the partition whose two paths occur without overlap of their index intervals. Leading zeroes do not change f.

## 3. Closed counting formula

For the star closure L→v→x above, the two-run decomposition gives

    H(T)=2A_L+h_L−T_L.                            (3)

A direct proof is to remove x from a Hamilton path. The remaining one or two paths traverse L and the singleton v in order. Equivalently this is the R=empty case of the exact star factorization proved in the separate structural packet; the decomposition gives coefficients c2=A_L and c1=2A_L+h_L−T_L.

Applying(1) to the prefix, fixed core, and suffix gives the following standalone binary formula. For the selected core-cut constants(h,q,A,G,T,S0,E), put

    F(p,r,b) = (qp+A)(2^(b+1)−rev_b(r))
               +h(2r+1−f(r))
               −pG−hf(p)−T−(p mod2)S0−E first_b(r).       (4)

Then H(T)=F(p,r,b). In deriving this, the prefix-core coordinates are

    A_pre=qp+A,
    T_pre=pG+hf(p)+T+(p mod2)S0,

and suffixing contributes A_pre rev_b(r)+hf(r)+E first_b(r) to T_pre. The prefix length a is retained to build the graph and test strongness, but does not otherwise enter(4), since leading zeroes in p have no effect on its numerical terms.

There is also a returning bit compiler, used by the producer rather than formula(4). Appending a sink with cut membership ε to a current body maps

    (h,q,A,G,T,S0,E)
      ↦ (h,2q,2A+εh,G+εq,T+ε(A+E),S0,εh).

Prepending a source of membership ε maps it to

    (h,2q,A+εq,2G+εh,T+ε(G+S0),εh,E).

These follow immediately from(1). They are physical operations on the cut-labeled source, not scalar recurrences on unmarked output counts.

## 4. Exact finite coverage, with retained recipes

**Finite theorem.** With a+b<=16, the three-core construction above realizes every odd N in[23,1,000,000] except47 and63. There are499,987 represented targets. Every retained recipe has at most22 vertices, although the enumerated construction universe allows up to23 vertices.

`construct.cpp` exhausts all80 cut states, every split a+b<=16, every allowed prefix word, and every suffix word:88,735,747 parameter choices. It evaluates the returning compiler and retains the first recipe for each target. `RECIPES.bin` contains explicit little-endian unsigned64-bit entries indexed by floor(N/2). A zero entry means no retained recipe. Otherwise subtract1 and read:

    bits0..6:   core-cut index c;
    bits7..11:  prefix length a;
    bits12..16: suffix length b;
    bits17..32: prefix word p;
    bits33..48: suffix word r.

The file has500,001 entries, hence4,000,008 bytes. Its meaning is fully specified without native-endian assumptions.

The verifier independently recomputes the80 finite core states, decodes every retained recipe, checks its cut and size constraints, and evaluates the closed formula(4), using the direct pair sum for f rather than the producer's bit update. It additionally counts120 deterministically sampled small retained graphs by literal vertex-subset DP and counts every allowed template with at most3 padding vertices. The displayed strongness condition proves strongness for every retained graph without assuming a generic count is strong.

`EXCEPTIONS.json` supplies separate literal strong tournaments with counts47 and63, checked by subset DP and displayed Hamilton cycles. Consequently these recipes plus those two graphs also certify every odd N in[23,1,000,000] as a strong tournament count. This is a scoped recertification result; no assertion that this exceeds every prior unrestricted coverage range is made.

## 5. Elementary all-length inverse bounds and permanent omissions

The covers counted by T_L are a subset of those counted by A_L: the former impose the extra condition that the first path is nonempty and ends in the cut. Therefore

    0<=T_L<=A_L,    H(T)=2A_L+h−T_L>=A_L+h.       (5)

For a strong recipe, put A0=qP+A. Its strongness condition guarantees A0>0. Suffix padding gives exactly

    A_L=2^b A0+hr.

Consequently any target N satisfies

    qP+A<=N−h,
    2^b(qP+A)<=N−h,
    0<=r<=min(2^b−1, floor((N−h−2^b(qP+A))/h)).  (6)

For P>0 the prefix length is exactly its bit length, because the first prefix bit is1. For P=0 the legal strong prefix is empty and the core cut must be nonempty. There are thus finitely many legal P,b,r for each fixed target and cut state. Equations(4),(6) give an exact **all-length membership algorithm** for this three-core family, without assuming a padding cap. They do not prove membership for every N.

For N in{47,63}, the80 literal core states yield a particularly small exhaustive universe. The H5 core has q=32 and minimum positive A=5; hence P<=1. If P=1 then A0>=32 and(6) forces b=0. If P=0 then b<=3. The H9 cores have q=84 and minimum positive A=9; hence P=0 and b<=2. Therefore every possible strong recipe for either target has total padding at most3.

`verify.py` independently counts every legal template with a+b<=3 by literal vertex-subset DP, checks formula(4), and finds neither47 nor63. These two counts are therefore omitted at **all padding lengths**, by an elementary argument with no external tournament-order bound. This does not exclude them from the unrestricted spectrum; their separate strong witnesses are included.

## 6. The remaining inverse problem

**OPEN three-core inverse conjecture.** Every odd N>=65 is represented by formula(4) for one of the80 specified core-cut states and legal binary words p,r of arbitrary lengths.

Together with the already retained low targets23..63 and the two exceptional graphs, this conjecture would prove the requested all-odd Chan–Pak conclusion. This is a specific sufficient inverse-numeration statement about three literal cores, not an equivalent characterization of all tournaments and not a proved induction. The finite million-range certificate does not bound the next possible gap.

The formulation retains precisely what the arithmetic needs: two binary words, their lengths for literal reconstruction, and one of80 fixed states. It replaces a search over arbitrary marked graphs by a fully explicit binary expression. A denominator-removal/rebase operation, or an all-target inverse proof for(4), remains missing.

## 7. Inheritance and reproduction

The path-cover composition method and star closure law are inherited from the earlier Session216 structural packet, which in turn uses the established path-cover kernel. The new content here is this fixed three-core specialization, its explicit two-word formula, the retained finite certificate, the permanent restricted omissions, and the exact inverse conjecture. Candidate graphs were discovered in the existing small-tournament bank, but no bank completeness assumption is needed: all three cores, all constants, all recipes, and the separate exception graphs are local.

Reproduce from this directory:

    g++ -std=c++17 -O3 construct.cpp -o /tmp/three_core_construct
    /tmp/three_core_construct CORE_STATES.txt /tmp/three_core_recipes.bin
    cmp /tmp/three_core_recipes.bin RECIPES.bin
    python3 verify.py > /tmp/three_core_verify.json
    cmp /tmp/three_core_verify.json VERIFY.json

Normal and optimized Python receipts are retained separately and must agree. An independent audit should use its own table/count code and should not infer completeness from the producer's hit file alone.
