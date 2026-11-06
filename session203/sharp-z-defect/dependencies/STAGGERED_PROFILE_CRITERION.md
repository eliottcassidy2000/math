# A staggered endpoint-profile criterion for the odd spectrum

**PROVED conditional implication; finite certificates and guarded compiler independently audited.** The all-order profile lemma below is **OPEN**. No unconditional eventual threshold is proved. This packet isolates a concrete physical sufficiency condition and supplies a compiler that refuses to assume a missing certificate.

## 1. Actual profiles and inherited readout

For a nonempty rooted tournament (U,v), let h=H(U), s=start_v(U), and let L count ordered spanning two-path covers (P,Q), where P is nonempty and ends at v, and Q may be empty. Define

    B(U,v)=h+L-s.

Adding a vertex a with v->a and a->U-v makes an actual tournament E(U), and the inherited rooted-amalgamation formula gives H(E(U))=B(U,v). Thus B is always a positive odd integer. Its value is coupled to the actual rooted graph; it is not an independently chosen coordinate.

The separately audited positive-difference construction, frozen in dependencies/, gives another actual graph C(U). Its rooted core has outgoing rows(2,8,3,5), with root0 identified with v; every other core vertex dominates U-v, and all edges of U are retained. Its exact count is

    H(C(U))=3h+3L-s=3B(U,v)+2s.                (1)

This graph compiler, not formal cancellation, is the mechanism used below. No strength condition is imposed on U or required of the target witnesses.

## 2. The uniform finite-width condition remains open

For integers J>=0 and odd B0>2J, define the following all-order statement F(J,B0):

For every odd B>=B0 and each r in {1,2}, there exist an integer j with0<=j<=J and an actual nonempty rooted tournament (U,v) such that

    B(U,v)=B-2j,
    start_v(U)=r+3j.                            (2)

The two residues may use different j and different graphs. The fixed J is uniform over every B. Neither a finite list of successful bases nor positivity of a profile vector proves this condition.

The concrete target in this packet is

    F(8,137).                                  (OPEN)

In set language, if Gamma(b) denotes the start counts realized by actual rooted U with B(U)=b, condition(2) asks for at least one of

    1+3j in Gamma(B-2j),  0<=j<=8,

and at least one of the corresponding memberships with2+3j. This permits a small shift of the base instead of demanding prescribed endpoint counts1 and2 at exactly B.

Substituting (2) into the literal graph law(1) gives

    H(C(U))=3(B-2j)+2(r+3j)=3B+2r.             (3)

The correction terms cancel in the count of one explicitly constructed tournament. No subtraction of graph counts or arbitrary-base addition is invoked.

## 3. Conditional implication with a finite seed interval

**General implication.** Let N be an odd integer at least3. Suppose F(J,N) holds and every odd integer in[N,3N-2] is the count of an actual tournament. Then every odd integer n>=N is a tournament Hamilton-path count.

**Proof by strong induction on odd n.** The stated finite interval supplies the initial values. At the next and every later odd n, n>=3N.

If3 divides n, then n/3 is an odd integer in[N,n), so induction provides a tournament T with H(T)=n/3. Form the ordered join of a directed triangle before T. Every Hamilton path visits the first factor before the second, so the counts multiply and the join has3H(T)=n.

If n is2 modulo3, put r=1 and B=(n-2)/3. If n is1 modulo3, put r=2 and B=(n-4)/3. Since n is odd, B is odd. The first values in these two cases after3N are3N+2 and3N+4, respectively, so B>=N. Apply F(J,N) to B and r, and use(3). This supplies an actual n-count graph. The three residue cases exhaust every odd n, completing induction.

**Concrete consequence, still conditional.** The supplied seed file contains literal witnesses for every odd value23<=n<=409. In particular it contains[N,3N-2]=[137,409] for N=137. Therefore F(8,137), if proved, would imply every odd n>=137 occurs. The remaining supplied seeds23 through135 would then give every odd n>=23.

This is only positive existence above23. The exclusions of7 and21 and any other classification of smaller values are not asserted or used here. The OPEN hypothesis is the substantive infinite step; the verified finite seed interval does not establish it.

## 4. Physical certificate compiler and its failure boundary

SEEDS.tsv gives one literal adjacency-list witness for each of the194 odd integers23 through409. A seed line is

    target, order, outgoing rows.

PROFILE_CERTIFICATES.tsv contains862 actual rooted graphs, two certificates for every odd B from137 through997. A profile line is

    B, r, j, root, order, outgoing rows.

The uniform bounds are r in{1,2} and0<=j<=8. The verifier counts h,s,L directly and checks both equations(2). It then compiles the four-core graph and counts its Hamilton paths literally by vertex-subset DP, checking(3). These certificates prove F(8,137) only over the finite base interval137<=B<=997, not at larger bases.

compiler.py accepts an odd requested target n>=23. It proceeds as follows:

1. If n has a listed seed, validate its tournament and exact path count, then return its rows.
2. If n is divisible by3, recursively compile n/3 and prepend an ordered triangle factor.
3. Otherwise compute B and r as in the proof. Read an actual matching profile certificate, check its tournament, actual h,s,L values and equations(2), compile C(U), and verify its exact path count.
4. If that certificate is absent, return MISSING_OPEN_PROFILE_CERTIFICATE and stop. The code never treats the unproved all-order statement as an oracle answer.

The joins can be built at arbitrary graph order using Python integer bitmasks. Exact subset-DP counting is deliberately limited to small certificate graphs; large compiled joins are certified by the proved product bijection and their already checked factor witnesses. No literal enumeration of a large joined graph is claimed.

Examples:

    python3 compiler.py 413
    python3 compiler.py 1239
    python3 compiler.py 3001

The first uses an actual profile certificate; the second additionally uses a triangle join. The third deliberately stops because its required base999 certificate is not in this finite packet. That failure is not a nonexistence claim for3001.

## 5. Verification, provenance and scope

The source graphs were extracted from an existing tournament bank through order8; no new tournament enumeration was performed. Completeness or isomorphism uniqueness of that bank is not a dependency. Every published graph stands on its own adjacency rows and independent exact counting. The seed values themselves are inherited finite existence facts, not new spectrum claims.

The standalone verifier.cpp validates every orientation, counts all194 seed graphs, checks all862 actual source profiles and all862 compiled children, and checks the complete declared finite ranges. There are1,918 literal graph-count calculations. The largest compiled child count is2,995. Normal and debug receipts must agree.

    g++ -O3 -std=c++17 verify.cpp -o verify
    ./verify > NORMAL.json
    g++ -O0 -D_GLIBCXX_ASSERTIONS -std=c++17 verify.cpp -o verify-debug
    ./verify-debug > DEBUG.json
    cmp NORMAL.json DEBUG.json

The compiled physical map and the conditional induction are exact. The arbitrary-base normal form, the uniform finite-width hypothesis F(8,137), and the unconditional odd-spectrum threshold remain OPEN. This packet does not convert the finite interval of profile certificates into an infinite assertion.
