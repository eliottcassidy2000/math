# Independent audit of the common retained-profile source pair

6 November 2026. **PASS in the stated constructive scope.** The two physical sources share one padding state, retain distinct start counts, and remain a common-state pair under arbitrary common rooted right amalgamation. This does not prove F(8,137), arbitrary base normalization, or an eventual odd-spectrum threshold.

## Reviewed inputs

ROW_NORMALIZATION.md has SHA-256 6424bba04da878d83705191df27434f3110a0be36fb4ea0f9fa8c4b68663f8fc. The producer verifier has SHA-256 8fc057ac001a8153ba6ec9b14564c8bf6fec4f659105d87b7e2d2d2c56eae111. BASE999_WITNESSES.json has SHA-256 42113a43d18592d444ab5a3ebd80d2d4d71be58b443f247c07d4c44eaf90e5a2. All entries of AUDIT_INPUTS.json passed. Exact reviewed snapshots and frozen inherited dependencies are included.

The only correction requested during this review was finite-universe wording: the producer checks all 75 labeled right tournaments through order 4 with root fixed at vertex 0. The current proof says this explicitly. This universe is sufficient up to relabeling; the independent checker additionally tests all 285 actual marked-root choices. No mathematical formula changed.

## Physical state and padding

For a nonempty rooted tournament U, L and R count two labeled paths covering U, with the designated nonempty path ending or starting at the root. The other path may be empty. The quantity beta is the ordered two-bin path-cover count of U-root, allowing either bin to be empty. With H(empty)=1, define

    B=h+L-s,  V=R+beta,  row=(B,e,V,d).

The hypothesis |U|>=2 is necessary for the claimed start-preserving Z move. For the rooted singleton, Z changes the start count from 1 to 0 and the displayed row update is not valid. This excluded boundary is tested explicitly.

O adds a universal sink. It leaves h,L,s,d unchanged, sets e to zero, and doubles R and beta. Z adds a vertex beating only the root. A root-starting Hamilton path must place it last, so s is preserved for |U|>=2. Z gives h'=h+R-e, e'=d and R'=2R-d. The subtraction of e removes old Hamilton paths ending at the root, since they cannot append the new vertex. Root-ending covers place the new vertex either last in the ordinary path or immediately before the root in the designated path, giving L'=L+beta. Root deletion makes the new vertex a sink, giving beta'=2beta and d'=d. These are disjoint physical path decompositions, proving the exact updates for every valid input.

Thus the row evolves by

    O: (B,e,V,d) -> (B,0,2V,d),
    Z: (B,e,V,d) -> (B+V-e,d,2V-d,d),

while s remains unchanged. Equality of rows is preserved by every common finite O/Z word, without any inference from a finite suffix search.

## The actual collision and positive ray

The source rows are exactly

    U6=(10,20,33,38,13,19), root 0,
    U7=(10,12,17,52,35,7), root 0.

Independent literal enumeration gives the full coordinates (h,s,e,L,R,d,beta):

    U6: (45,6,9,72,57,13,104),
    U7: (43,7,9,75,57,13,104).

Both are strong, and both rows are (111,9,161,13). Their common padding histories therefore give actual four-core child counts 3B_w+12 and 3B_w+14. At the empty word these are 345 and 347.

For k>=1 consecutive Z moves, V_k=148*2^k+13 and e_k=13. The first step gives B_1=263; summing the subsequent increments yields B_k=148*2^k-33. Hence the child counts are 444*2^k-87 and 444*2^k-85. These formulas do not apply at k=0. Every Z-only source remains strong, and the fixed strong four-core readout of a strong source is strong. An O move may instead introduce a sink component. No strongness claim for every general word is made.

## Arbitrary common right input

The common row is not a new independent coordinate system. If E(U) is the ordered amalgamation of the rooted edge root->a with U, retaining the old root, its standard profile (H,e,R,d) is exactly (B,e,V,d). This follows directly from the frozen rooted-amalgamation formulas and is also checked on literal edge extensions. Deleting the distinguished edge vertex a recovers U and its different start count. Equality of the edge-extension profiles does not preserve this deletion sidecar automatically.

Physical associativity gives E(U*W)=E(U)*W. For every rooted right input W of order at least 2, the inherited retained-profile law therefore gives precisely the four displayed row coordinates. The original left source is nontrivial, so its start endpoint scales as

    s(U*W)=s(U) H(W-root).

For the singleton right input, all coordinates are unchanged. In particular e'=e, rather than the nontrivial formula d*e_W. The proof handles this identity boundary correctly. The other formulas do survive it.

For the specified pair, the common output base is

    B'=13h_W+148L_W-4s_W-46d_W,

and the two start counts are 6d_W and 7d_W. Although this displayed linear form has negative coefficients, it is the count of an actual edge extension. Equivalently it is

    111d_W +148(L_W-d_W)+4(d_W-s_W)+13(h_W-d_W),

a positive retained-state expression. This is not a subtraction of separately constructed graphs. The coordinates of W must still belong to a single actual rooted tournament.

The endpoint gap becomes d_W, and is one exactly when d_W=1. There is no independent prescription of the resulting base. One additional scope observation is that 6d_W is always 0 modulo 3, so this pair and its common right lifts alone never provide both required endpoint residues 1 and 2 for F. This observation does not rule out different source pairs or additional operations.

## First missing-base certificates

Both archived constructions at the conditional compiler's first absent base B=999 are literal and valid:

- Starting at root 6 of (2,68,9,19,39,15,61), append O,Z. The source has order 9 and (h,s,L)=(485,16,520), hence B(U)=989=999-2*5. Its order-12 readout has H=2999.
- Starting at root 3 of (2,20,9,3,13), append O,O,Z,O,Z. The source has order 10 and (h,s,L)=(341,2,660), hence B(U)=999. Its order-13 readout has H=3001.

The archived original, final source and compiled rows agree exactly with independent physical construction. The two counts are recomputed by literal path enumeration. These are two new finite certificates, not a new universal claim. The published conditional compiler and its frozen missing-certificate controls remain unchanged.

## Independent verification

check.cpp enumerates actual directed prefixes, updating exact path counts on every visited vertex set. It does not use the producer's subset recurrence. Full child counts are independently obtained by nonmemoized Hamilton-path backtracking. Physical amalgamation explicitly identifies arbitrary root labels; no producer graph-building code is imported.

The independent normal and debug builds agree byte-for-byte and check:

- All 31 O/Z words of length at most 4 on both sources, including empty and all-Z boundaries.
- All 75 labeled right tournaments of orders 1 through 4 and all 285 choices of marked root, producing 570 right-amalgamation/child pairs.
- Literal graph equality in E(U*W)=E(U)*W for every such pair, not only equality of counts.
- Another 568 generic O/Z state checks on all nontrivial marked inputs through order 4.
- The excluded singleton Z input, direct retained edge profiles, and the two base-999 constructions.

There are 1,523 exact profile computations, 634 literal full Hamilton counts, and 54,389,277 directed-prefix visits. These numbers describe finite controls, not the source of the arbitrary-order theorem. The proof is the physical decomposition and retained-state law. All graphs counted have at most 13 vertices, so the unsigned 64-bit counters are well within the 13! bound.

The standalone Python witness-row checker passes normally and with optimization. The producer C++ verifier is also replayed with optimized and debug builds, agreeing with both frozen author receipts.

## Reproduction

In this directory:

    g++ -std=c++17 -O3 -Wall -Wextra check.cpp -o check
    ./check > NORMAL.json
    g++ -std=c++17 -O0 -g -D_GLIBCXX_ASSERTIONS -Wall -Wextra check.cpp -o check_debug
    ./check_debug > DEBUG.json
    cmp NORMAL.json DEBUG.json
    python3 check_witness_rows.py > WITNESSES_NORMAL.json
    python3 -O check_witness_rows.py > WITNESSES_OPTIMIZED.json
    cmp WITNESSES_NORMAL.json WITNESSES_OPTIMIZED.json
    g++ -std=c++17 -O3 PRODUCER_VERIFIER.cpp -o producer
    ./producer > PRODUCER_NORMAL.json
    g++ -std=c++17 -O0 -g -D_GLIBCXX_ASSERTIONS PRODUCER_VERIFIER.cpp -o producer_debug
    ./producer_debug > PRODUCER_DEBUG.json
    cmp PRODUCER_NORMAL.json PRODUCER_DEBUG.json

The reviewed packet includes no later arbitrary-module source-pair theorem. Such an extension needs its own proof and audit. Publication and status promotion remain separate from this review.
