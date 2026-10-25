# Two overlapping magic hourglasses: exact seam, unique completion, and AP-grid chords

Status: PROVED elementary structural equivalences; FINITE-EXACT bounded controls. No integer seven-distinct-square hourglass or nine-distinct-square magic square is constructed or excluded. No priority claim. All entries below are VALUES, not their chosen square roots.

## 1. Inheritance and primary-source boundaries

Cain, *Gaussian Integers, Rings, Finite Fields, and the Magic Square of Squares*, arXiv:1908.03236v2, §2, defines the horizontal seven-cell hourglass by retaining the top and bottom rows and center. The five constraints are those two rows, the middle column, and the two diagonals. Theorem 2.1 and Corollary 2.1 give total three times the center and opposite pairs summing to twice the center. The rotated copy supplies the other orientation. Primary text checked at https://arxiv.org/html/1908.03236v2 ; publication identity/version checked at https://arxiv.org/abs/1908.03236v2 . We use those elementary identities, reproved below, not the paper's later Gaussian parametrizations or unproved assertions.

The precise earlier repository `3+1+3=7` reference is `05-knowledge/results/FRAGMENT-DECODE-mu3-depth2-tree-klein-S691.md`. Its top-priority correction and the MISTAKE-319 heading “root audit of klein-S691” in `01-canon/MISTAKES.md` identify three successive inverse fibers of the sporadic Keller map, with cardinalities 3,1,3. The underlying theorem is THM-2473, `01-canon/theorems/THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy.md`. The missing two sheets escape to infinity; they do not collide at finite points. The visible fixing symmetry is an involution, not μ3. The former Arithmetic-Kakeya identification was withdrawn. This establishes a genuine prior decomposition but NOT an arithmetic bridge to the present hourglass. A cardinality-preserving arbitrary bijection does not preserve magic sums, squares, or inverse-map fibers.

Other historical repository uses of “hourglass” concern adjacent-layer cycle residues, not Cain's seven square entries. Shared terminology is not a theorem dependency. The closest operative mechanisms are Session57's AP-K3,3 labelled-grid equivalence and its fixed-gap doubled-elliptic-x carrier; exact connections appear below.

## 2. The exact gluing theorem, without division

Write a nine-cell array as

    A B C
    D E F
    G H I.

Let the horizontal hourglass have cell set {A,B,C,E,G,H,I} and lines ABC, GHI, AEI, BEH, CEG. Let the vertical hourglass have cell set {A,C,D,E,F,G,I} and lines ADG, CFI, AEI, DEF, CEG.

**Theorem.** Over any additive abelian group, compatible horizontal and vertical magic hourglasses glue to a magic square iff their values agree at all five shared labelled cells A,C,E,G,I. Their totals then agree automatically. Conversely restriction of a magic square gives exactly such a compatible pair. Positivity and square membership are checked on the two copies. For these compatible hourglasses, injectivity on both copies implies injectivity on their union, as proved below. These predicates cannot be checked on unaligned copies.

Proof: the shared cells are the four corners and center, so 7+7−5=9. The shared lines are AEI and CEG, so 5+5−2=8, precisely all magic-square lines. For either hourglass, summing its three center lines gives 3T while summing its outer rows and three copies of the center gives 2T+3E. Cancellation gives T=3E with no division. Therefore agreement at E already makes the two totals equal. Agreement at the other four labelled shared cells makes the union an actual array, and its eight line equations all hold. Restriction is the inverse operation.

The five shared cell values are not five independent scalar coordinates: opposite-pair identities force A+I=C+G=2E. Thus the core has just three free coordinates E,A,C. Agreement at these three named cells suffices, because I and G are then forced. Agreement at an unordered multiset or only at E and T does NOT suffice. Root signs are immaterial to equality of square VALUES.

Over Z or Q with square entries, the existence of two separate, unaligned hourglasses is not a construction. In fact seven distinct values in EACH compatible copy do guarantee nine distinct values globally. The only collision not already tested within a copy would equate one of E±(x+z) with one of E±(x−z). Such an equality forces 2x=0 or 2z=0, which already identifies two shared corners E±x or E±z. This contradiction works in any additive abelian group. Square membership remains independent. A pre-audit draft incorrectly denied this implication; the draft witness E=10,x=2,z=1 already repeats values inside its vertical hourglass, so it does not refute the corrected theorem.

## 3. Every scalar hourglass has a unique scalar completion

Put x=A−E and z=C−E. The five hourglass equations force

    E+x     E−x−z     E+z
      ·       E        ·
    E−z     E+x+z     E−x.

Indeed opposite-pair sums give I,G,H, and ABC=3E gives B; the bottom row then follows. The unique possible missing values are

    D=E+z−x,   F=E+x−z.

Substitution proves that these values always complete all eight sums. Thus a scalar hourglass is not a two-parameter family of extensions: its extension fiber is a singleton. Arithmetic completion can fail only because the forced cells fail square membership, positivity, or global distinctness.

The five equations in seven cell variables and T have rank5. Embedded in all nine cells plus T they still have rank5; the full eight equations have rank7. The other hourglass adds only TWO independent linear constraints, although it contributes three previously absent lines. This rank count is compatible with the three-parameter normal form and holds over any field: the preceding eliminations use only coefficients ±1. At fixed E the family has two scalar parameters.

In characteristic2, opposite values coincide because 2E=0 and −v=v. Hence a distinct seven-cell hourglass, and therefore a distinct nine-cell square, is impossible. In characteristic3 the total is0, but the parametrization and gluing proof remain valid; dividing by3 would incorrectly discard this case.

## 4. Exact square-gap completion gate and collision guards

For E in Q define D_E={t in Q : E+t and E−t are rational squares}. Require E itself to be a square. Then

- a square-valued horizontal hourglass is exactly x,z,x+z in D_E;
- its full square-valued completion exists exactly when x−z also belongs to D_E.

For integer-square entries use integer E,x,z and integer squares in this definition. Membership in D_E is symmetric under t→−t, but arbitrary independent changes of x and z can exchange which hourglass is being tested. Do not drop the signed sum relation.

Over Q or R the horizontal seven entries are distinct iff

    x z (x+z)(x−z)(2x+z)(x+2z) ≠ 0.

All nine entries are distinct iff

    x z (x²−z²)(x²−4z²)(4x²−z²) ≠ 0.

Proof: the respective offset sets are {0,±x,±z,±(x+z)} and {0,±x,±z,±(x+z),±(x−z)}. Comparing their signed pairs gives the displayed factors. For positive entries the full array requires E>|x|+|z|; the hourglass alone requires E>max(|x|,|z|,|x+z|). No order/positivity statement is imported to finite fields.

These are reformulations of the simultaneous square conditions, not proofs that the required gap configuration exists.

## 5. Exact AP-K3,3 connection: the overlap is a spanning tree

Use two indexed potential shores A_i=i x and B_j=E+j z for i,j∈{−1,0,1}; their edge values are L_ij=E+i x+j z. Potentials are formal indexed vertices; coincident numeric potentials do not identify graph vertices. Session57's scaling/gauge can supply positive disjoint integer vertex values when an eligible positive integer full-square witness is available.

The horizontal hourglass retains the seven edges except (1,−1),(−1,1). The vertical hourglass retains the seven edges except (1,1),(−1,−1). Their intersection is the five edges with i=0 or j=0: a spanning tree of K3,3. The central row and column AP guards already live on that tree. Each hourglass adds two chords; its additive edge consistency closes two fundamental cycles. Their union closes all four cycles of K3,3.

The map preserves every squared-entry value, its affine-grid address, center, gap relation, and all square-membership predicates. The original magic placement is recovered from Section3. Dropping the edge labels loses the arithmetic: the bare K3,3 graphs exist regardless of whether the labels are squares. Dropping which midpoint is the center loses the AP guard. Thus the exact missing object is TWO forced square-labelled chords, not an arbitrary second K3,3-minus-two-edges subgraph.

This is different from Session57's pair of Hamilton paths, whose five-edge supports share one edge and cover K3,3. Here each hourglass has seven edges and they share five. These are different decompositions of the same labelled carrier, not interchangeable counts.

## 6. Elliptic connection and its limit

Session57's fixed-gap set is C_d={c>0 rational : c−d,c,c+d are squares}. For positive gap d, c in C_d is equivalent to a rational right triangle of area d/4 with hypotenuse squared c, and to c=4x(2P) on the eligible non-torsion positive branch of Y²=X³−(d/4)²X.

The hourglass condition says E lies in C_|x|, C_|z|, and C_|x+z|, excluding zero gaps for distinctness. These are generally THREE DIFFERENT elliptic curves/areas. Completion requires E also in C_|x−z|, again with all collision guards. There is no asserted elliptic addition law between those points on different curves.

By contrast, the completed grid gives E−z,E,E+z all in ONE fixed-gap set C_|x|. This recovers Session57's exact fixed-curve three-center AP condition, including the doubled-x image restriction. The hourglass alone lacks one square at each of the two off-center fibers. A point count, independent curve ranks, or a group-theoretic progression of elliptic points does not fill these two missing square conditions.

## 7. The seven-cell incidence is not the Fano plane

The five hourglass lines form a linear triple system: no pair of vertices occurs in two lines. But its two outer rows are disjoint, whereas any two lines of the Fano plane meet. Hence this system cannot be embedded as five Fano lines on the same seven vertices.

Stronger: every pair touching the center is already used by a spoke. The used outer pairs consist of the two three-cliques (rows) and the three opposite-pair matching edges. Unused pairs are therefore K3,3 minus a perfect matching, a six-cycle. A new line would require a triangle of unused pairs, but this graph is triangle-free. Thus the hourglass is a MAXIMAL partial Steiner triple system with five lines, not an incomplete Fano plane awaiting two additional lines. This is a genuine incidence obstruction, but it supplies no integer-square obstruction.

The two appearances of K3,3 must not be conflated: here its six vertices are outer CELLS and its edges are available cell pairs; in Section5 its six vertices are POTENTIALS and its edges carry cell values. Only the latter carries the Session57 square-sum arithmetic directly.

## 8. Hypotheses, hostile controls, and bounded evidence

H1: compatible labelled horizontal+vertical hourglasses are equivalent to full magic arrays. PROVED in Section2.
H2: one square hourglass automatically completes to a square magic array. REFUTED over F13 by the seven distinct square values

    1 3 9
    · 0 ·
    4 10 12.

All five sums are0 mod13, but the forced missing entries are8 and5, both nonsquares. All nine completed scalar values are distinct. This does not refute the integer implication on seven-distinct-square hourglasses: no such integer premise is supplied here.

An integer repeated-entry hostile is E=25,x=24,z=−24: the hourglass entries are49,25,1;25;49,25,1, all squares, while D=−23,F=73 are not squares. Its repetitions are explicitly disallowed in the original distinct problem.

H3: three square AP spokes at one center suffice. REFUTED by Cain's §2 displayed center1105² example: opposite-pair sums agree but the top-row defect is nonzero. Our checker reconstructs the values and checks this without treating it as an actual hourglass.
H4: seven-cell injectivity in BOTH compatible copies implies full injectivity. PROVED by the shared-corner collision argument in Section2. Injectivity in only ONE copy is insufficient: E=10,x=2,z=1 has seven distinct horizontal entries but repeats in its vertical mate.
H5: the seven-cell layout extends to a Fano plane. REFUTED structurally in Section7.
H6: the earlier Keller 3+1+3 and the hourglass are the same mechanism. UNSUPPORTED, with explicit historical correction against finite-point-identification reasoning. No transfer theorem is claimed.

Run `python3 check.py`, and independently `python3 -O check.py`. The standard-library script uses explicit exceptions, not removable assertions. It verifies exact rational Gaussian elimination, a three-vector kernel basis, all1331 parameter triples E,x,z∈[−5,5], collision guards, the unused-pair obstruction, and the AP-grid tree map. It enumerates all prime fields p≤59 by TWO methods: (i) actual square values A,B,C,E and opposite sums, and (ii) offset square-gap sets D_E. Their complete labelled seven-distinct-hourglass sets agree. Counts include F13:12 hourglasses/0 full completions; F29:84/56; F47:276/0; F59:1044/232. Entries may include zero; squares are residues, never integer lifts. No isomorphism quotient or orbit census is claimed. Normal and optimized output are byte-identical. certificate.json retains the full counts and first noncompletion witnesses.

Stopping boundary: this resolves the carrier, gluing, rank, square-gap, and prior-work interfaces. Existence or nonexistence of distinct integer-square hourglasses and full squares remains outside these structural deductions.

Second-wave supplement: MITRE_AFFINE_BRIDGE.md identifies the classical mitre configuration from Brouwer1977, proves the minimal9-point Steiner host, reconstructs54same-center mitre pairs and432magic layouts, and gives the exact missing-line defects±3x,±3z. It also certifies a nonzero-center F37 hostile.
