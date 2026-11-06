# A closed physical nine-coordinate state for O/Z/J

**PROOF CANDIDATE, pending independent audit.** The transfer theorem is all-input. The bounded inverse is complete for the stated physical seed grammar. It does not prove F(8,137) or an eventual Hamilton-path threshold.

## 1. Coordinates and actual operations

Let (U,v) be a nonempty rooted tournament. Set h=H(U), s=start_v(U), e=end_v(U), d=H(U-v), with H(empty)=1. Let L and R be the ordered spanning two-path-cover counts with the designated first path nonempty and respectively ending/starting at v; the other labeled path may be empty. Let C=B_U(2) and b=B_(U-v)(2) count two labeled spanning paths, either allowed empty. Finally let T count three labeled spanning paths whose first path is nonempty and starts at v, with the other two allowed empty.

The coordinate order is

    (h,s,e,L,R,d,C,b,T).

O appends a global sink and keeps v. Z appends z with z->v and every old nonroot vertex pointing to z, keeping v. J is the same graph as Z but re-marks z as the root. These are literal vertex additions, not formal profile arithmetic. Define B=h+L-s and V=R+b, as in the staggered-profile criterion.

For every nontrivial U, the exact transfers are

    O: (h,s,0,L,2R,d,2C,2b,3T),

    Z: (h+R-e,s,d,L+b,2R-d,d,
        2C+2T-2L,2b,3T-b),

    J: (h+R-e,s,h-e,C-L,h+R,h,
        2C+2T-2L,C,T+C).                    (1)

At a singleton input all formulas remain valid except the s coordinate of Z: its correct value is zero rather than one. O and J have no singleton exception. The coordinates stay actual nonnegative integers because they count the displayed physical constructions.

## 2. Derivation, including the new closed coordinate

The inherited ordinary h,e,R,d laws follow by insertion at z. For clarity the relevant graph law is H(Z(U))=h+R-e, and re-marking does not change h. The old rooted cover laws are L_O=L, L_Z=L+b and L_J=C-L. The last one counts a cover ending at the new root: its final old vertex may be any vertex other than v, so it is precisely an arbitrary ordered two-cover of U except those whose first path ends at v. Similarly R_J=h+R, splitting its designated path into the singleton new root or the new root followed by v.

A global sink doubles the ordinary two-bin count C: append it to either labeled path, including an empty path. It similarly doubles b. Root deletion under Z still appends a global sink, whereas root deletion under J leaves exactly U, proving b_Z=2b and b_J=C.

For C_Z, split according to the bin containing z, giving a factor two. If the old vertices in that bin form A, then the Hamilton count on A+z is H(U[A]) when v is absent from A, and H(U[A])+R_v(U[A])-e_v(U[A]) when v is present. Thus

    C_Z=2 sum_(A subset U) H(U[A]) H(U-A)
        +2 sum_(A containing v) (R_v(U[A])-e_v(U[A]))H(U-A)
       =2C+2T-2L.                              (2)

The second convolution of R is T: its existing ordinary path and the new complementary path are two independently labeled ordinary bins. Empty ordinary bins are retained. Equation(2) also gives C_J because J only changes the mark.

To compute T_O, append the new sink to any one of its three labeled paths, giving 3T. For T_Z, the new vertex must again end its path: if it used its only outgoing edge to v, it would precede v in the designated root-starting path, which is impossible. Appending it to any of the three old bins gives 3T candidates. The only illegal candidates append it immediately after the singleton designated path (v), because z->v. These are counted exactly by b, hence T_Z=3T-b.

For T_J, the designated path either is the singleton z, leaving C choices for the other two bins on U, or begins z,v and continues as an old designated root-starting path, leaving T choices. Thus T_J=T+C. This closes the state at nine coordinates. No unbounded path-cover hierarchy is required for this operation alphabet.

## 3. Positive B growth

For every actual rooted U,

    C>=2L,     R>=d>=e,     h>=s,     L<=B.       (3)

For the first inequality, the C covers in which the root ends the first path and those in which it ends the second path are disjoint classes of L each. The root-singleton designated covers give R>=d, and deleting the root injects root-ending Hamilton paths into paths of U-v, giving d>=e. The remaining inequalities are immediate.

If U is nontrivial and s>0, the root has an out-neighbor. A designated path (v,w) and any Hamilton path on the remaining vertices show R>d. Since U-v is nonempty, its two-bin covers include the 2d covers with one empty bin, so b>=2d>=2.

From(1),

    B_O=B,
    B_Z-B=R+b-e>=2,
    B_J-B=R-e+C-2L>0.                          (4)

Every B is odd because it is the Hamilton-path count of the actual edge extension E(U): add a with v->a and a->U-v. Paths starting at a correspond to the h-s old paths not starting at v; every other path has v immediately before a and corresponds to one L-cover. Thus H(E(U))=h-s+L=B, and inherited tournament Hamilton-path parity applies. Consequently the final strict inequality improves to B_J>=B+2. All three moves preserve s for these inputs.

This is a source-height monotonicity fact, not a recurrence on arbitrary odd integers: all nine coordinates remain coupled to the actual graph.

## 4. A finite target-sensitive inverse for any finite seed bank

Fix a finite set of actual nontrivial rooted seeds with positive s, and target bound M. Consider exactly their O/Z/J extension grammar. Every represented pair (B,s) with B<=M has a representative word with no terminal O letters, because deleting that terminal run preserves B and s.

By(4), such a word has at most floor((M-B_seed)/2) letters in {Z,J}. Write it as successive blocks O^k Z or O^k J. For k>=1 the preceding O run transforms e to zero and R,C,b to 2^k times their old values, while B,L stay fixed. Therefore

    B_(O^k Z)=B+2^k(R+b),
    B_(O^k J)=B+2^k(R+C)-2L.                    (5)

If the resulting B is at most M, the first formula bounds 2^k(R+b)<=M-B, and the second gives

    2^k(R+C)<=M-B+2L<=M+B<=2M.                 (6)

Thus each O run has k<=floor(log2(2M)). Zero-length runs are permitted and use the original e, not the k>=1 formula. These bounds already give a finite exhaustive word universe.

There is a sharper order bound. On an n-vertex source, C>=2^n and b>=2^(n-1), because every ordered vertex bipartition has at least one path on each bin. A Z child has B_Z>=b>=2^(n-1). A J child has L_J=C-L>=C/2>=2^(n-1), so B_J>=2^(n-1) too. After deleting the terminal O run, every nonempty remaining word ends in Z or J. Therefore any B<=M grammar witness has a representative of order at most

    max(seed order, floor(log2 M)+2).             (7)

The seed-order alternative is essential for an empty or all-O word. This is a complete target-dependent cap, not a guessed order cutoff. State merging by all nine coordinates is valid because(1) determines every continuation. Positive s and nontrivial seeds are essential to the monotonicity used above; no arbitrary source-membership test is supplied by a relaxed cone.

For endpoint one, the separately inherited complete physical O/Z/J source theorem identifies all nontrivial rooted s=1 sources with the orbit of the rooted edge. Hence this procedure is also a complete all-order decision method for membership 1 in Gamma(B), with B=1 additionally admitting the singleton. This consequence inherits that source-classification theorem; it is not claimed for arbitrary s or all sources behind F(8,137).

No finite successful run of this inverse supplies an all-target induction. A target-dependent normal form or descent proving every required B exists remains open.

## 5. A nonnegative cover block

The added closure has a useful positive representation. Put

    A=C-2L>=0,   W=T-b>=0,   q=b-L>=0.

The first inequality was proved above. The root-singleton designated three-path covers give T>=b. Removing the terminal root from the designated path of an L-cover injects it into an arbitrary ordered two-cover of U-v, giving L<=b. Empty first bins after removal are retained. In coordinates (L,A,W,q), the exact cover maps are

    O: (L, 2L+2A, 3W+L+q, L+2q),
    Z: (2L+q, 2A+2W, 3W, q),
    J: (L+A, 2W+2L+2q, W+L+q, L).

All coefficients are nonnegative. This is a physical state change, not permission to treat every nonnegative vector as realizable. The inherited h/e/R/d cone can be retained separately, with s fixed as before. The new block alone already avoids signed cover coordinates in exact continuation and inverse work.

## 6. Guarded target compiler

`compiler.py` takes a positive odd requested B. It enumerates the physical rooted-edge grammar by the nine-coordinate transfer, pruning only states with B above the target (justified by monotonicity), and uses the order cap in(7). Breadth-first merging retains the shortest recipe for each full nine-coordinate state, so it cannot remove a continuation available within the order budget. It returns a literal source graph if one exists. Exhaustive failure is labeled `ABSENT_UNIT_START_PROFILE`, not a Hamilton-count exclusion. If the configured state budget is exhausted, it returns `INCOMPLETE` and permits no absence conclusion.

Examples:

    python3 compiler.py 137
    python3 compiler.py 139
    python3 compiler.py 997 --inventory

The first excludes endpoint one at B=137 at every order, using the inherited complete unit-start source theorem. It does not refute the staggered criterion, which permits shifted bases and other starts. The second supplies an actual unit-start source at B=139. The third records 330 attained unit-start bases through 997, with the finite exhaustive search bounded by the theorem rather than a guessed order. These are scoped finite-target conclusions, not occupancy of an infinite tail.

`check_compiler.py` directly counts every physical recipe through length 5, checks the resource-limit behavior, validates the positive target, and compares the complete B<=137 support with a separately generated canonical strong-unit grammar through the proved order cap 9. Weak unit-start graphs add only singleton component tails, which do not change B; the singleton value1 is included separately.

## 7. Exact controls and scope

    python3 verify.py > NORMAL.json
    python3 -O verify.py > OPTIMIZED.json
    cmp NORMAL.json OPTIMIZED.json

The verifier directly counts paths, two-bin covers and the designated three-bin cover on every labeled root-zero tournament through order5 and all three literal children. It checks the singleton Z exception and the B-growth statements for all positive-start nontrivial inputs. A separate three-vertex bin-allocation/permutation calculation checks the T coordinate. No removed assertions are used as validation gates.

The frozen `dependencies/PHYSICAL_SOURCE_GRAMMAR.md`, Theorem1, supplies the inherited complete unit-start O/Z/J grammar and ordinary h,e,R,d laws. Its public source route is `session198/bilinear-descent/PROOF.md` in eliottcassidy2000/math; the frozen snapshot was recovered from the audited Session200 signed-block-transport dependency packet. Other inherited ingredients are the root-amalgamation and B edge-extension conventions in Session197/203; and tournament Hamilton-path parity. The contribution here is the explicit C,b,T closure and its B-targeted finite-word bound. These scopes must be independently audited before promotion.
