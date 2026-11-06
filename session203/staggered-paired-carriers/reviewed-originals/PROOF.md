# Two exact advances toward the staggered-profile criterion

**PROOF CANDIDATE, pending independent audit.** The two results below are all-order statements about actual rooted tournaments. Neither proves nor disproves F(8,137), and no eventual odd-spectrum threshold is claimed.

## 1. Definitions and inheritance

For a nonempty rooted tournament (U,v), let h=H(U), s=start_v(U), e=end_v(U), d=H(U-v), with H(empty)=1. Let L and R count ordered spanning two-path covers whose designated first path is nonempty and respectively ends or starts at v; the other path may be empty. Write b for the number of two labeled spanning paths on U-v, both allowed empty, and C for that count on U. Finally T counts three labeled spanning paths on U whose first path starts at v and is nonempty. Put

    B=h+L-s,  V=R+b.

The physical O/Z/J moves and their complete nine-coordinate transfers are inherited from `session203/rooted-cover-closure/PROOF.md`. O appends a global sink and retains the root. Z appends a vertex beating only the root and retains the root. J makes the same addition and marks the new vertex. For nontrivial positive-start inputs, all preserve s.

F(8,137), as defined in `session203/staggered-profile-criterion/PROOF.md`, asks for every odd X>=137 and each r in {1,2} for an actual source satisfying

    B=X-2j,  s=r+3j,  0<=j<=8.

The inherited physical four-vertex context then has count 3B+2s=3X+2r. That context has root-core rows (2,8,3,5), identifies root 0 with v, and makes every other core vertex dominate U-v.

Inheritance board: the finite staggered certificates; the complete nine-coordinate transfer; the bounded-order theorem for each requested (B,s); the start-two multi-entry source mechanism; and coupled rather than separate residue occupancy. The hostile is silently upgrading a finite seed orbit into an all-source classification. The least-used sidecar is simultaneous membership of the two staggered residue tracks.

## 2. Four physical paired carriers

The following eight actual rooted sources are specified by outgoing rows. The label A/B within each pair is not a claimed isomorphism class.

    677 A: root 2; rows (2,12,17,101,75,23,167,63)
    677 B: root 2; rows (2,12,17,165,43,71,31,119)
    759 A: root 1; rows (2,44,49,133,11,89,31,119)
    759 B: root 1; rows (2,52,25,35,105,197,15,95)
    1547 A: root 0; rows (2,148,97,7,45,75,27,125)
    1547 B: root 0; rows (2,148,41,67,45,75,23,125)
    821 A: root 2; rows (4,9,50,85,163,75,151,47)
    821 B: root 2; rows (4,9,50,85,99,139,167,31)

Their literal profiles in order (h,s,e,L,R,d,b) are

    677 A: (165,10,27,516,237,49,824)
    677 B: (175,11,27,507,249,49,812)
    759 A: (253,13,29,511,287,29,640)
    759 B: (253,11,29,511,287,29,640)
    1547 A: (259,13,101,1293,267,109,1424)
    1547 B: (259,11,101,1293,267,109,1424)
    821 A: (295,16,33,532,359,33,720)
    821 B: (293,17,33,535,359,33,720).

Use respectively j_A,j_B=(3,3),(4,3),(4,3),(5,5). In every pair A has s=1+3j_A and B has s=2+3j_B. The adjusted base X=B+2j and the other padding coordinates agree:

    (X,e,V,d) = (677,27,1061,49),
                 (759,29,927,29),
                 (1547,101,1691,109),
                 (821,33,1079,33).

**Theorem 1.** Apply any common finite word in O,Z to both members of any displayed pair. The resulting two actual tournaments supply both F(8,137) certificates at the same odd base X given by the exact row dynamics

    O: (X,e,V,d) -> (X,0,2V,d),
    Z: (X,e,V,d) -> (X+V-e,d,2V-d,d).            (1)

Every base in these four physical word orbits is therefore certified for both residues. Applying the inherited four-vertex context gives actual Hamilton-path counts 3X+2 and 3X+4.

**Proof.** The original rows are certified by literal subset path counting. The O/Z transfers give B_O=B, B_Z=B+R+b-e, V_O=2V, V_Z=2V-d. Since s and hence j stay fixed, the same identities hold with B replaced by X=B+2j. Both members therefore retain identical (X,e,V,d), and their fixed j values remain in [0,8]. These are physical graph operations, so no positivity cone is treated as a realization theorem. The inherited count context gives the final assertion. O leaves X unchanged, and Z increases it; thus all bases remain at least137. QED.

For example, repeated Z gives these four explicitly certified infinite base sequences, for k>=1:

    1012*2^k-313,
     898*2^k-139,
    1582*2^k-27,
    1046*2^k-225.                              (2)

The second and fourth formulas also hold for k=0. In general (1) gives V_k=2^k(V-d)+d, e_k=d for k>=1, and

    X_k=X_0+(2^k-1)(V_0-d)+d-e_0,  k>=1.

These sequences are a genuinely coupled residue result with bounded starts. They do not cover every odd base. They are not a proof that this four-orbit language has a tail.

The restriction to O,Z matters. For the 759 pair, C is respectively 4384 and4552 and T is2713 and2797. A common J changes the adjusted bases to4379 and4547, so the paired invariant is lost. Merely retaining (X,e,V,d) is insufficient for J.

## 3. No finite seed bank is complete for all bounded-start nine-states

Define the three nonnegative quantities

    e,  p=d-e,  f=R-h+e-d.                    (3)

Their nonnegativity has a direct counting interpretation. Removing the root injects root-ending Hamilton paths into paths of U-v, so d>=e. To interpret f, use covers (Q,P) with Q starting at v. Cutting each Hamilton path just before v gives h distinct covers. Independently, taking Q=(v) gives d covers. Their intersection has e covers. Thus f is precisely the number of remaining covers: Q has at least two vertices, P is nonempty, and the last vertex of P is beaten by v.

Every O output has e=0, every Z output has d-e=0, and every J output has f=0. The last identity follows directly from the all-input transfer:

    h'=h+R-e, e'=h-e, R'=h+R, d'=h,
    R'-h'+e'-d'=0.                            (4)

Now let U_n be the rooted-edge graph after n-2 repeated J moves, with each new root relabeled0 and the preceding root1. For n>=3, U_n is strong. Indeed U_3 is a directed triangle, and adding a vertex beating only the current root preserves strong connectivity. Both vertices0 and1 have outdegree1 and start count1: the new root starts through the old root, while the old root retains its unique starting path under the underlying Z operation. The ordinary one-vertex statements apply from the nontrivial rooted edge onward.

For every n>=4, form W_n by adding a new root v which beats precisely vertices0 and1 of U_n and loses to every other vertex. Then

    start_v(W_n)=start_0(U_n)+start_1(U_n)=2.   (5)

W_n is strong. Its new root reaches all of U_n through0; every old vertex reaches a vertex other than0,1, which beats v. Therefore e(W_n)>0. Also

    d(W_n)-e(W_n)=end_0(U_n)+end_1(U_n)>0.     (6)

Finally U_n-0 is U_(n-1), with root1, and is strong for n>=4. Choose a Hamilton path P of U_n-0 ending at1. Then Q=(v,0), together with P, is one cover counted by f: Q has two vertices, P is nonempty, and its last vertex1 is beaten by v. Hence

    f(W_n)>0.                                (7)

**Theorem 2.** Every nine-coordinate state of W_n, n>=4, lies strictly inside all three faces in (3). It is not the image under O,Z, or J of any physical rooted tournament. Consequently no fixed finite bank of physical seeds, extended solely by O/Z/J, is complete for the nine-coordinate states of all rooted tournaments with start count at most2.

**Proof.** Equations(5)-(7) and (4) exclude every last letter. Thus any representation of the full state in such a seed grammar would have to be an empty word, so the state must already occur among the seeds. There are infinitely many distinct states: d(W_n)=H(U_n)>=n, because U_n is strong and every vertex starts a Hamilton path. A finite bank cannot contain them all. QED.

This is stronger than a graph-isomorphism obstruction: even a different generated graph cannot have one of these full states, because (3)-(4) are state equalities. But the conclusion must NOT be strengthened to the pair (B,s). Another generated graph could have the same pair while differing in the retained coordinates. The theorem does not disprove profile-completeness for (B,s), and in particular it does not refute F(8,137).

The surviving route is a grammar that includes multi-entry constructions, or a proof that every needed (B,s) has an alternate representative in a chosen seed grammar. The order bound in `session203/rooted-profile-order-bound/PROOF.md` remains applicable to each requested profile, but is not an all-target occupancy theorem.

## 4. Reproduction and scope

`verify.cpp` directly validates every orientation, counts the seven displayed coordinates, checks all common O/Z words of length at most5 from all eight carriers, and checks the four pure-Z formulas. It also builds U_n,W_n for4<=n<=14 and checks strong connectivity, both unit-start inputs, the s=2 output, and all three strict inequalities. Normal and unoptimized debug outputs must agree. These finite controls validate the formulas and boundaries; the all-order arguments are Sections2-3.

The discovery bank was the inherited `session194-structural/inherited_tournaments_through8.txt`; its completeness is not used. `PAIR_PROBE.txt` and `pair_probe.cpp` record discovery only. The independent literal controls use just the eight displayed graphs and the stated physical constructions. No literature-priority claim is made. F(8,137) and the eventual threshold23 remain OPEN.
