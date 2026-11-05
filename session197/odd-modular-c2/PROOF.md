# Every odd-modulus residue occurs in both strong C1 and C2 classes

Session197. This is an all-order modular representability theorem, not a cofinite-coverage theorem. The physical constructions are ordinary start-one gadgets followed by one reviewed two-source C2 graft. No negative graph operation, formal source tuple, or finite-to-infinite extrapolation is used.

**Theorem.** For every positive odd integer m and every integer r, each of the classes C=1 and C=2 contains a finite strong tournament T with

    H(T)=r modulo m.

The proof is constructive but gives no useful small-order or bit-complexity bound. It does not say that a prescribed positive integer r itself occurs. It does not settle residue coverage for even moduli, except for the necessary oddness of every actual output.

## 1. Physical matrices and a fixed final readout

Use the reviewed rooted start-one source grammar. Its triple v=(D,E,U)^T counts Hamilton paths, paths ending at the root, and reserved-start ordered covers. Beginning with a triangle gives v=(3,1,3)^T. Every ordinary gadget preserves the one-start condition; applied to a nontrivial strong source it produces another strong source.

The gadget matrix is

    M_L=[[g(L),-(L mod2),L],
         [1,-[L=1],0],
         [Q(L),0,B(L)]].                            (1)

The named ROOTED dependency proves these physical transfers and the binary recurrences for g,Q,B. In particular

    A=M_1=[[1,-1,1],[1,-1,0],[1,0,1]], det A=1;
    M_5=[[5,-1,5],[1,0,0],[7,0,8]], det M_5=8;
    M_9=[[9,-1,9],[1,0,0],[15,0,16]],
    M_11=[[9,-1,11],[1,0,0],[11,0,16]],
    M_13=[[11,-1,13],[1,0,0],[13,0,16]].

The last three determinants are16. Thus A,M_9,M_11,M_13 are invertible modulo every odd modulus. Their inverses in a modular calculation are realized by positive powers: an invertible element of a finite group has finite order, and G^{-1}=G^{ord(G)-1}. Consequently **every finite group word in these matrices has a positive physical gadget word with the same matrix modulo m**. No actual graph deletion or inverse gadget is performed. Enumerating successive modular powers until the identity is one terminating way to obtain the needed exponents.

After constructing a physical source with triple(D,E,U), append one final M_5. Its four responses are

    (h,e,R,d)=(5D-E+5U,D,7D+8U,D).

Graft it with the nontrivial transitive two-vertex source using the reviewed N construction. Since the latter source has quadruple(1,0,2,1), the resulting strong C2 graph has

    H=3h-2e+6R-2d=53D-3E+63U=:F(D,E,U).           (2)

This is the physical final readout used throughout. The C2 graft theorem applies even though the left two-vertex source is weak; both graft sources are nontrivial.

## 2. Exact matrix identities, including a rejected shortcut

Let E_ij denote the matrix unit and T_ij(t)=I+tE_ij. All following identities are exact over the rationals, with powers-of-two denominators only. Write [X,Y]=XYX^{-1}Y^{-1}. Direct multiplication gives

    P=M_9 M_11^{-1}, Q=M_11 M_13^{-1},
    [P,Q]=T_12(3/4).                               (3)

For completeness, both P and Q have shape

    [[1,*,b],[0,1,0],[0,c,1]].

Their (b,c) pairs are(-1/8,4) and(-1/8,-2). The central commutator coefficient is bc'-b'c=3/4, proving(3) without a numerical matrix search.

Here are identities recovering every elementary direction from T_12(1) and conjugation by A. First,

    [T_12(1), A^2 T_12(1) A^{-2}] T_12(1)^2
       = T_32(-4).                                 (4)

Once T_32(1) is available, put

    S=A^{-1}T_32(1)A,
    T_31(1)=[T_32(1),S]T_32(1),
    R=T_31(1)T_32(-1),
    T_12(1)=A R A^{-1}T_32(-1),
    T_21(1)=T_12(-1)S T_12(1),
    T_13(1)=T_31(-1)[A T_32(1) A^{-1}]T_31(1),
    T_23(1)=[T_21(1),T_13(1)].                     (5)

These identities are independently testable exact matrices. For example S=I+(e1+e2)(e1-e2)^T, making the first two commutators transparent.

A tempting shorter identity is **false**:

    M_5 M_7^{-1} != I+2E_32.

The actual matrix is

    [[1,5/4,-1/4],[0,1,0],[0,2,1]].

Its first-row change cannot be dropped. The proof uses(3), not that incorrect shortcut. The verifier retains it as a negative control.

## 3. Relative elementary moves away from3

Factor the desired odd modulus as

    m=3^a q, a>=0, gcd(q,6)=1.

If q=1 there is no relative-adjustment step. Otherwise choose an integer k>=0 such that

    k is divisible by3^{max(a-1,0)},
    3k=4 modulo q.

Such k exists because q is coprime to3. Raising(3) to k gives a group word W_12 with

    W_12=I modulo3^a,
    W_12=T_12(1) modulo q.                          (6)

For a=0 the first condition is vacuous. All fractions have denominators invertible modulo m.

Apply(4) with W_12 and raise the result to an exponent congruent to(-4)^{-1} modulo q. It gives a relative W_32 which is identity modulo3^a and T_32(1) modulo q. Then(5) gives every relative W_ij. Conjugating identity modulo3^a remains identity, so every step preserves the first congruence. Powers W_ij^c realize arbitrary elementary parameters c modulo q, still acting identically modulo3^a. All these group words are replaceable by positive physical gadget words by Section1.

Thus any prescribed finite sequence of elementary integer row operations can be performed on a source triple modulo q, while leaving its entire triple unchanged modulo3^a. The proof below gives explicit such operations; no transitivity theorem over a composite ring is needed.

## 4. Set the3-adic residue using two actual primitive rays

For j>=1, apply gadget L=2^j to the singleton. Its actual strong source has

    v0=(z,1,z), z=2^j+1.

Its final readout(2) is

    F(v0)=116z-3=113+116*2^j.                       (7)

Alternatively first apply M_5 once more, giving

    v1=(10z-1,z,15z),
    F(v1)=1472z-53=1419+1472*2^j.                  (8)

Both vectors are primitive as integer triples: v0 has middle entry1, and the first two entries of v1 have gcd1 since(10z-1)-10z=-1.

For a>=1,2 generates every unit modulo3^a. One elementary proof is that4=1+3 has order3^{a-1}: cubing a number1+3^t u with3 not dividing u raises the exact3-adic valuation of its difference from1 by one. Together with2=-1 modulo3, this gives order2*3^{a-1}, the full number of units. Exponents may always be chosen positive, so j>=1 imposes no restriction on those residues.

The coefficient116 is a3-adic unit and113=2 modulo3. Hence(7) attains exactly those residues not equal to2 modulo3. Similarly1472 is a unit and1419=0 modulo3, so(8) attains exactly those residues not equal to0 modulo3. Their union is every residue modulo3^a.

Choose v0 or v1 and positive j so that F(v)=r modulo3^a. If a=0 simply choose v0 with j=1. This supplies a physical primitive source before the relative adjustment.

## 5. Adjust the other odd factor and finish physically

Choose the explicit primitive integer target vector

    w_r=(2r+3,35r+53,0)^T.

It satisfies F(w_r)=r and -35(2r+3)+2(35r+53)=1. It is used only as a modular target, never asserted to be a physical source.

Here are explicit elementary row operations mapping the chosen physical seed vector to e1, and then to w_r.

For v0=(z,1,z), successively apply

    T_12(-z), T_32(-z), T_12(1), T_21(-1).

The resulting vectors are(0,1,z),(0,1,0),(1,1,0),e1. For v1=(10z-1,z,15z), successively apply

    T_12(-10), T_21(z), T_31(15z),
    diag(-1,-1,1).

The first three produce(-1,0,0). The last matrix is elementary: it is the square of T_12(1)T_21(-1)T_12(1).

Finally the determinant-one matrix mapping e1 to w_r is

    G_r=[[2r+3,-2,0],[35r+53,-35,0],[0,0,1]]
       =[T_21(17)T_12(1)T_21(2)T_12(-1)]T_21(-r).

Products act in the ordinary right-to-left order. This gives a completely explicit finite elementary sequence; its parameters may be reduced modulo q.

Replace each elementary operation by its relative physical group word from Section3. On the q factor it maps the current source triple to w_r. On the3^a factor it leaves the original physical seed triple unchanged. The resulting positive-gadget source v' therefore satisfies

    F(v')=r modulo q and modulo3^a.

The Chinese remainder theorem gives F(v')=r modulo m. Append the final M_5 and perform the N graft as in Section1. This is a literal finite strong C2 tournament, with Hamilton-path count exactly F(v') and hence the desired congruence. The cases q=1 and a=0 require only their nontrivial factor; m=1 is immediate. This proves the C2 part of the theorem.

Only matrix residues, not signed numerical source vectors, were physically transformed. Every intermediate actual source comes from positive ordinary gadgets, preserving the source condition and strongness. Modular inverses were compiled into positive repetitions; they never mean deleting vertices or running a graph operation backwards.

## 6. The same construction gives C1 residue surjectivity

No final graft is needed for the C1 claim. Every strong source with exactly one starting Hamilton path has exactly one Hamilton cycle: each cycle gives a distinct rooted starting path after cutting at its root, so there is at most one, while strong tournaments have a Hamilton cycle (the standard Camion input already used by the inherited source grammar).

For the3-adic selection use the same primitive seeds, but read D rather than F:

    D(v0)=1+2^j, D(v1)=9+10*2^j.

These attain respectively all residues not equal to1 modulo3 and all residues not equal to0 modulo3, modulo every3^a. Thus choose a primitive physical seed with D=r modulo3^a. For a=0 any seed suffices.

The relative elementary operations from Section3 preserve this3-adic triple. Reduce the seed to e1 modulo q using the explicit operations of Section5. Then successively apply relative T_21(1) and T_12(r-1), sending e1 to(r,1,0) modulo q. The resulting actual strong start-one graph therefore has D=r modulo m and C=1. Again q=1 needs no relative adjustment. This proves the C1 part of the theorem without a C2 graft or an assumption that a modular target vector is physically realizable over the integers.

## 7. Checks, dependencies and limitations

The verifier checks(3)-(5) by exact rational arithmetic, retains the false M5/M7 identity as a hostile control, and checks the relative CRT transvection for every odd modulus through101. It separately uses positive physical generator words to obtain all residues in both scalar readouts for every odd modulus through41, records exact uncapped counts, and checks the3-adic ray residue sets for both classes through3^6. Literal graph counters check the final C2 readout for seven small source words, including the three larger generators. No finite result is used as a premise of the all-order argument.

Run `python3 verify.py` and `python3 -O verify.py` in a copy. It writes VALIDATION.json and WITNESSES.json. The unchanged core.py supplies the previously reviewed actual gadget constructors and path/cycle counters. The only mathematical inputs are the pinned ROOTED.md physical transfers and BILINEAR_C2_INVERSE.md physical C2 graft theorem. The new group identities, relative congruence argument, explicit elementary maps, and3-adic seed selection are given in full above.

This theorem excludes an obstruction to either individual class depending solely on residue classes modulo an odd integer. It does not exclude sparse nonmodular gaps, prove any prescribed count, imply cofinite C1/C2 coverage, or provide a bounded-width construction. No priority claim is made.
