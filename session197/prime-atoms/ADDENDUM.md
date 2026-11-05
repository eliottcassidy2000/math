# Exact component transfer and the first prime-amplitude gaps

Session 197. This addendum extends the independently checked bounded-slack normal form in `PROOF.md`. **All-order exact classifications at deletion counts 9, 11 and 13, and of odd prime unit amplitudes through 19. The full Hamilton-path spectrum remains open.** It does not replace or change the frozen first packet.

## 1. The full boundary sum factors through strong components

Let S=C1 -> ... -> Cs, and adjoin a root r with arbitrary out-neighbor set O. Put O_i=O intersect C_i. For each component define seven nonnegative quantities

    h_i=H(C_i), q_i=P_Ci(2),
    u_i=sum_(x in O_i) R_x(C_i),
    v_i=sum_(y in O_i) L_y(C_i),
    a_i=sum_(y in O_i) ends_y(C_i),
    b_i=sum_(x in O_i) starts_x(C_i),
    f_i=sum_(distinct x,y in O_i) B_yx(C_i).

R and L are the rooted starting/ending two-path-cover responses. B_yx has its first path ending at y and second starting at x. P_C(2) allows either path slot to be empty. These definitions do not assume that the root extension of C_i is strong.

For a prefix of the ordered components retain

    d=H(prefix), U=sum_(x in O intersect prefix) R_x(prefix),
    f=sum_(distinct x,y in O intersect prefix) B_yx(prefix),
    a=sum_(y in O intersect prefix) ends_y(prefix).

Start at (d,U,f,a)=(1,0,0,0). Appending a component with data (h,q,u,v,A,b,F) gives the exact nonnegative linear transfer

    d' = h d,
    U' = q U + u d,
    f' = h f + F d + v U + b a,
    a' = A d.                                       (1)

After all components, the actual root-extension profile is

    H(T)=d+U-f-a, e_r(T)=d-a, R_r(T)=d+U,
    H(T-r)=d.                                       (2)

**Proof.** Deletion count multiplies. The U formula is the SCC transport from the first packet: old starting covers see the new component in two independent path slots, and a start in the new component puts all old vertices in the first path. Only the last strong component contributes endpoint mass, giving a'=A d.

For f, pairs of out-neighbors both in the old prefix force the new component entirely into the second path, contributing h f. Pairs both in the new component force all old vertices into the first path, contributing F d. If the second path starts at an old out-neighbor and the first ends at a new one, their old parts form a rooted starting cover and their new parts a rooted ending cover, contributing U v. In the opposite order, the first path lies entirely in the old prefix and ends at an out-neighbor, while the second lies entirely in the new component and starts at an out-neighbor. This contributes a b. These four disjoint cases are exhaustive. Equation (2) is the first packet's U and slack identities. QED.

Equivalently, for x in C_i and y in C_j with i<j,

    B_yx(S)=(product_(k<i) h_k) R_x(C_i)
             (product_(i<k<j) q_k) L_y(C_j)
             (product_(k>j) h_k).

The reverse boundary B_xy vanishes unless j=i+1, in which case it is the corresponding prefix/suffix product times ends_x(C_i) starts_y(C_j). This also verifies (1) directly.

For a singleton beaten by r the local tuple is (1,2,1,1,1,1,0); for a singleton beating r it is (1,2,0,0,0,0,0). Thus the binary padding transitions are precisely

    root loses: (d,U,f,a) -> (d,2U,f,0),
    root wins:  (d,U,f,a) -> (d,2U+d,f+U+a,d).

These recover the inherited d=1 slack recursion and the stationary K/D operations, while (1) handles arbitrary nontrivial components exactly. It avoids subset-DP on a long padded deletion graph; it does not replace the full source geometry by a residue calculation.

A useful limitation is visible: beta=d-f updates as

    beta'=h beta-F d-v U-b a.

Once beta is negative, appending a component cannot make it positive. At beta=-1, a unit return requires h=1 and F d+v U+b a=0. This statement concerns this fixed-root ordered-component operation, not root-switching or general reconstruction.

## 2. Why the small core inventory is complete at every order

Use the [S195 strong-core theorem](https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/strong-structure/PROOF.md). A directed triangle has H=3 and F2=6, so its anchored-growth threshold is 3*3+6=15. Every strong tournament with H<15 therefore grows from any chosen triangle by a chain of induced strong single-vertex extensions. All intermediate Hamilton counts are at most the final count by induced monotonicity.

An exhaustive labelled order-six computation gives minimum H=15 among strong tournaments. The lower bound is sharp: rows

    [36,1,2,7,15,30]

give a strong order-six H15 graph. **The minimum is 15, not 17.** The verifier explicitly tests all 32,768 labelled order-six tournaments and records this minimum.

Consequently a strong tournament with H<=13 has at most five vertices. Otherwise its strong singleton chain would contain an order-six strong induced subtournament of count at most13, contradicting the certified floor15. This structural reduction, not a guessed graph-order bound, makes the following finite inventory complete:

    H=3: one strong isomorphism type, order3;
    H=5: one type, order4;
    H=7: none;
    H=9: two types, order5;
    H=11: one type, order5;
    H=13: one type, order5.

The exact canonical adjacency rows are in the receipt. The inventory is independently counted by the original subset-DP verifier and compared under all vertex permutations.

Since strong-component Hamilton counts multiply, a deletion tournament with d=9 consists, apart from singleton components, of either one of the two H9 cores or two ordered directed triangles. For d=11 or13 it contains just the corresponding one strong core. There is no deletion count7: a prime total seven would require a strong component of count seven, and none exists.

## 3. Exact unit-source table through deletion count 13

For beta=d-f in {+1,-1}, the following are all possible positive U in a rooted source. As in the first packet, each displayed set is multiplied by {2^j:j>=0}.

| d | beta=+1 | beta=-1 |
|---|---|---|
| 1 | {1} | {3,5} |
| 3 | {6} | {9,15} |
| 5 | {11,12,16} | {15,39} |
| 7 | empty | empty |
| 9 | {20,24,32,44} | {26,27,34,35,41,45,99} |
| 11 | {29,37,47} | empty |
| 13 | empty | empty |

The first three rows are in the frozen packet. Here is the exact finite reduction for the new rows. Use every deletion-core bank identified in Section 2. Insert at most k(f)+ell(f) singleton components in all possible positions between and around its nontrivial components. Choose every root incidence subset with at most k(f) out-neighbors. Retain precisely the reduced cores from Theorem 1 of the first packet, and use (1) to calculate their profiles. Filter for slack f=d-beta. The resulting base sets are exactly those displayed. The theorem guarantees that every possible physical source normalizes to one of these cores followed by a K-tail. Conversely every retained core is a literal tournament, and the K-tail supplies all dyadic multiples.

The enumeration retains representative adjacency rows and verifies them independently by subset DP. The new exact search sizes are:

    d=9, beta=+1: 18,092 reduced cores, 5 matching cores;
    d=9, beta=-1: 77,271 reduced cores, 15 matching cores;
    d=11, either sign: 8,493 reduced cores;
      beta=+1 has 3 matching cores; beta=-1 has none;
    d=13, either sign: 8,493 reduced cores, no matches.

The numbers of matching cores can exceed the numbers of U values. No isomorphism or value filtering precedes the completeness test. In particular the empty rows are all-order exclusions under the specified deletion data, not failures of a finite gadget bank.

The new positive entries give actual families U*2^m+beta for every m>=1. For example the prime amplitudes 29,37,47 all have positive unit offset with d=11, and amplitude41 has negative unit offset with d=9. The receipt contains graph witnesses and checks two successive actual K extensions for every new U entry. These scalar rays are not asserted to be new individually; the exact inverse classification is the contribution.

## 4. Complete odd prime unit-amplitude classification through 19

The S196 geometric completeness and prime-amplitude reduction apply to **general two-piece binary-star synthesis**, not merely the singleton-partner subclass. If its amplitude alpha=p is an odd prime and its offset beta is ±1, the factorization alpha=U_X L_Y forces L_Y=1, H(Y)=H(Y-c)=1, U_X=p, and beta=d_X-f_X. Thus the complete singleton-partner table applies to every such prime-amplitude amplifier.

The endpoint reinsertion bound from the first packet supplies

    3d_X <= 2p+2beta-1.                              (3)

For p<=19 the largest bound is d_X<=13. Because p is odd, no positive K-tail can occur in U_X=p; only the undoubled base entries of Section 3 need be checked.

**Theorem 4.** For odd prime p<=19 and beta in {+1,-1}, an actual positive binary transitive-cloning amplifier with

    H(T[r:m])=p*2^m+beta, for every m>=1

exists if and only if

    (p,beta) is (3,-1), (5,-1), or (11,+1).           (4)

All three positive cases have explicit physical witnesses in the first packet. The negative cases follow from the exhaustive inverse budget (3) and Section 3. In particular the primes7,13,17,19 admit neither unit offset, at every graph order. The same argument also excludes p=23 with beta=-1, since then d_X<=14 and only odd d_X<=13 need be tested. It does not decide p=23 with beta=+1.

For example alpha=19,beta=-1 would give H37 followed indefinitely by H->2H+1; the theorem excludes that cloning presentation. It does not exclude H37 itself or a beta=+1 amplifier whose composite amplitude is18. Likewise the prime-amplitude statement is not a forbidden-Hamilton-count theorem: the numbers 13,15,25,27,33,35,37,39 are not declared absent from the full spectrum.

For composite alpha, the second factor L_Y need not equal one. The deletion/slack table alone then does not settle the general two-piece determinant. This restriction is essential; no composite-amplitude exclusion is inferred merely from its absence in the U table.

## 5. Reproduction and a failed arithmetic shortcut

Run

    python3 verify_addendum.py > receipt-addendum.json
    python3 -O verify_addendum.py > receipt-addendum.optimized.json
    cmp receipt-addendum.json receipt-addendum.optimized.json

The checker imports only the frozen local `verify.py`. It verifies the SCC transfer independently against full subset-DP profiles on every labelled rooted tournament through order5 and larger random controls, checks the complete order-six strong floor, computes the finite-core table, and checks actual graph extensions. All checks remain active in optimized mode.

A tempting power-of-two shortcut was explicitly rejected: it is false that U a power of two implies f<=d. The actual rooted tournament

    rows=[52,1,2,7,14,30], root5,
    (h,e,R,d)=(19,1,73,9), U=64, f=46

is a counterexample. Thus factorization of a composite amplitude cannot be dismissed by this claimed inequality. The surviving all-order task is to solve the physical determinant with variable partners or use a richer operation such as the separately developed two-port grammar. Neither the tables nor modular compatibility prove full prime/7p coverage.
