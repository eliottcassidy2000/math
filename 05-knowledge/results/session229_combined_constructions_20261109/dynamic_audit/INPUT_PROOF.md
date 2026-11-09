# Context-stable gap-two switches and a universal three-value repair core

Status: PROVED exact all-order construction identities; FINITE-EXACT small-core interfaces and exhaustive search. No unbounded interval or Chan–Pak completeness claim. Independent audit: [PASS report](../dynamic_audit/AUDIT.md).

## 1. Inheritance, target, and preserved interface

The source is Session216 `three-core/PROOF.md`, its cut-state composition law and arbitrary binary prefix/suffix star closure. Session216 `numeration/PROOF.md` warns that a balanced inverse candidate does not establish occupancy. Session217 `run-entropy/PROOF.md` retains the reversed coordinate G and proves a logarithmic correction, not coverage. Session229 `session229-overlap/PROOF.md` supplies six H19 equal-alpha fibers and an H9 compensated carry. We combine these constructions rather than extrapolate their finite coverage. Canonical hostiles remain the old grammar's permanent omissions 47 and 63. Concept board: preserved interfaces; pure scalar defects; arbitrary binary contexts; cross-topology fibers; iteration divisibility. Anchor: context-stable adjacent odd values. Niche: minimal repair interface. Wildcard: vertex-transitive cores and larger scalar fibers.

For a nonempty tournament B and cut M, use the inherited state (h,q,A,G,T,S,E): h counts Hamilton paths; q counts ordered spanning two-path covers allowing empty paths; A counts covers whose second path is nonempty and starts in M; G counts covers whose first path is nonempty and ends in M; T imposes both conditions; S and E count full paths starting and ending in M.

The fixed interface for a repair is **(h,q,A,G,S,E)**, omitting only T. This retains more information than equal-alpha modules. For the old star closure, form L=TT_a⇒B⇒TT_b with prefix cut word P and suffix cut word R, add v,x with L→v→x, and x beating precisely the cut vertices of L. When P>0 it starts with1 and a is its bit length; when P=0 the prefix is empty and the core cut is nonempty. These rules ensure strongness. Every L vertex reaches x via v, and x reaches the initial block through its selected vertex; the cores used below are strong, while an initial singleton selected by prefix1 also suffices.

Let rev_b(R) reverse b binary digits and let

    f(R)=sum_(i<j, both selected) (2^(j-i-1)+1_(j=i+1)).

The inherited exact count is

    N=(qP+A)(2^(b+1)-rev_b(R))+h(2R+1-f(R))
       -PG-hf(P)-T-(P mod2)S-E first_b(R).                 (1)

For empty suffix use first=rev=f=0.

**Theorem 1 (pure-T repair).** Two physical cores with the same interface (h,q,A,G,S,E) and T-values differing by2 yield outputs differing by2 in every common legal binary context. A fiber with T-values t,t+2,...,t+2k yields k+1 consecutive odd output values in every such context.

Proof: subtract (1); every term except -T cancels. Each resulting graph is a tournament, hence its Hamilton-path count is odd. Strongness must still be supplied physically as above. No scalar inverse or density claim follows from this cancellation.

For fixed h,q, the retained interface is also necessary for a difference independent of **all** binary contexts, algebraically including P=0. Indeed P=R=0 and varying b first force ΔA=0. With b=0, comparing P=0,2 forces ΔG=0, then P=1 forces ΔS=0. Choosing suffix first bit1 forces ΔE=0. The remaining constant difference is -ΔT. This does not assert necessity in restricted context languages.

## 2. A five-vertex universal gap-two switch

Use outgoing masks

    B5=[24,17,3,6,12].

Its masks3 and5 have states respectively

    (15,120,42,42,10,6,6),
    (15,120,42,42, 8,6,6).                              (2)

Thus switching3→5 increases every output by exactly2, with arbitrary binary words on both sides. Without padding the outputs are89 and91. With empty prefix and all-zero suffix of length b, they are

    84*2^b+5, 84*2^b+7.                                (3)

This is an infinite family of explicit strong twin-count constructions, not an interval covering all odds. It is a context-stable pure-T difference: the interface does not change at all. The scalar gap2 is a unit difference in the half-count coordinate, not arithmetic unsaturation: gcd(15,120)=15, so this module is arithmetically saturated.

The graph is the vertex-transitive cyclic tournament on five vertices, in a different labeling. Vertex transitivity makes all single-vertex start/end and cover-start/cover-end weights equal. Therefore every two-element mask has A=G=42,S=E=6; the two pair orbits have T8 and10. Likewise three-element masks have the common interface A=G=63,S=E=9 and T26 or28. The identity is not an accidental equality of output counts.

Exhausting every labeled tournament and cut through four vertices finds no fixed-interface pair with T difference2. The displayed five-vertex switch is therefore minimal by deletion-core order **for this exact fixed-interface notion**, not for all possible context-dependent constructions.

## 3. Six vertices support universal triples, across two core topologies

Use

    B6a=[48,49,3,7,12,28],
    B6b=[40,37,17,6,11,28].

The following three choices share interface (45,432,129,129,15,15):

    B6a, mask9: T=20;
    B6a, mask3: T=22;
    B6b, mask3: T=24.                                  (4)

Each core is strong. With any common legal arbitrary prefix/suffix, these give three consecutive odd counts. Without padding they give283,281,279. With empty prefix and zero suffix of length b they give

    258*2^b+21, 258*2^b+23, 258*2^b+25.                 (5)

The highest and lowest values differ by4. The switch need not preserve deletion topology, but it preserves all coordinates that the old arbitrary-context formula sees except T.

`search.cpp` exhausts all labeled tournaments and all cuts on orders1 through6. It groups states across all these topologies by the full interface, retains each attainable T, and checks maximal consecutive step-two runs. The maximal run lengths are1 through order4,2 through order5, and3 through order6. `SEARCH.stats` reports cumulative interface groups across orders:2,6,18,64,342,3202. Per-order counts are2,4,12,46,278,2860, over respectively2,8,64,1024,32768,2097152 labeled rooted objects. Thus no four-value universal step-two fiber exists in this enumerated small-core universe. This is FINITE-EXACT and not an obstruction to larger cores. Reproduction instructions below specify the entire universe; no preexisting core bank or isomorphism-class completeness assumption is used.

## 4. Exact composition explains why these windows do not iterate for free

For ordinal composition of two nonempty cut-labeled blocks, the inherited law is

    h'=h1 h2, q'=q1 q2,
    A'=A1 q2+h1 A2, G'=G2 q1+h2 G1,
    T'=A1 G2+T1 h2+h1 T2+E1 S2,
    S'=S1 h2, E'=h1 E2.                                (6)

A pure-T replacement in a middle block therefore remains pure-T under arbitrary ordinal left/right contexts, with its difference multiplied by h_left*h_right. Empty context has multiplier1. Singleton contexts have h=1, explaining Theorem1 independently of the binary expression.

For t copies of B5, choosing either mask3 or5 independently, all cross terms in (6) are fixed. The attainable T changes are exactly

    -2*15^(t-1)*k,  k=0,...,t.

Hence the output count changes are +2*15^(t-1)*k. These are t+1 equally spaced values but are not consecutive odd values once t>=2. The same mechanism gives step2*45^(t-1) when composing the triple modules. Every such conclusion remains valid with arbitrary singleton padding inserted between cores. This explicitly blocks the unjustified iteration from one universal gap-two switch to unbounded consecutive intervals.

## 5. The H19 and H9 bridges identify the missing interface

For the H19 base B=[40,1,3,6,15,30] and its six equal-A fibers, the full differences in (G,T,S,E) from the first mask to the second are:

    17→40: (12,6,0,2)
    24→38: (78,26,-2,10)
    19→42: (12,0,0,2)
    21→44: (12,10,0,2)
    25→39: (78,44,-2,10)
    23→46: (12,4,0,2).

All have h19,q256 and equal A. Consequently, writing e=first_b(R),

    ΔN=-P ΔG-ΔT-(P mod2)ΔS-eΔE.                       (7)

Precisely the 19→42 pair permits |ΔN|=2: it requires P=0,e=1. Every other displayed pair has absolute difference at least4 (or6), and P>=1 makes each difference larger than2. This gives a context-dependent gap-two switch for arbitrary suffix tails; its smallest suffix is R=1,b=1, giving343 and345 on nine vertices. Among this six-fiber H19 universe, empty prefix plus a one-bit suffix is the minimal padding needed for gap2. This is not a minimality claim over all tournament constructions.

For the H9 block C=[16,1,3,7,14], compare the two-block cut words (2,28) and (4,0). The states are

    (81,7056,1344,3165,558,9,36),
    (81,7056,1344, 144,  0,18, 0).

They preserve A through a genuine carry, but the star-closure output difference is

    ΔN=3021P+558-9(P mod2)+36e,                       (8)

always at least558. Thus this carry does not supply a gap-two core repair. The H19 prefix coefficient and the H9 coefficient3021 show exactly why equal-alpha freedom alone was insufficient; pure-T fibers remove both defects.

A parallel H29 block [26,44,57,16,34,9] has equal-A81 masks9,10,17,18. Relative to9 their output offsets are

    0, 12P+6+2e, 36P+10+6e, 48P+8+8e.

At P=0,e=0 the last three give offsets6,10,8, so there are consecutive triples175,177,179 without padding, preserved by arbitrary suffixes starting0. Arbitrary prefix padding breaks this constant-width property. This is another distinction between equal-A fibers and the stronger pure-T interface.

## 6. Structural next step and hostile controls

Vertex-transitive cores are a principled source of fixed interfaces: for a fixed mask cardinality, A,G,S,E are constant, while T varies as a weighted pair statistic. But symmetry does not guarantee a consecutive T interval. Cyclic cores on7 vertices have h175,q1988, and two-element masks yield T68,82,84. On9 vertices, h3267,q49608, the corresponding values are1080,1196,1252,1258. The checker also records three- and four-element fibers. These explicit gaps reject a naive general interval claim.

The new operation creates universal width-two and width-four scalar windows in the old arbitrary-binary language. A proof of coverage would still need an occupancy theorem showing these moving windows meet every target, or a growing-interface construction whose T fiber contains long consecutive intervals without the multiplicative spacing of Section4. Neither is proved here. Small-core exhaustive maximality and the cyclic hostile examples narrow that next operation without ruling it out.

## 7. Reproduction and independent routes

Run from this directory:

    python3 verify.py > CHECK.json
    python3 -O verify.py > CHECK.optimized.json
    cmp CHECK.json CHECK.optimized.json
    g++ -O3 -std=c++17 search.cpp -o /tmp/dynamic-core-search
    /tmp/dynamic-core-search > SEARCH.txt 2> SEARCH.stats

`probe.py` supplies ordinary subset start/end DP and the inherited composition law. `verify.py` uses a separate literal permutation-and-split count to verify core states and whole-graph subset Hamilton DP to verify padded closures. No Python assertions or external repository files are required. The C++ exhaustive search has its own independently implemented start/end and pair-weight enumeration. Its graph encoding reads unordered pairs(i,j), i<j, lexicographically, with bit1 meaning i→j. The retained search records include this integer encoding and cut masks. The all-order quantifiers follow from (1),(6), not from any finite testing range.
