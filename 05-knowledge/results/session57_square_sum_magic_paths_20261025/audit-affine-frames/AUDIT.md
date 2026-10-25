# Independent affine-frame audit

Verdict: the AP theorem 6(q−9), its ordered factors and 432-layout correspondence pass. One path-scope sentence needs correction; the suggested all-frame strengthening also passes independent algebra and finite enumeration.

## Base packet verification

- Original verifier passes normally and under python -O, with byte-identical results: unordered AP counts 432, 4320, 13068 for q=81,729,2187. All recorded SHA256SUMS and source/origin pins match.
- Independent trial division by every monic polynomial of degree at most floor(n/2) confirms all three field moduli irreducible.
- Independent F3-vector coordinates (no imports from source verifier) confirm q=9,27,81 counts, including translated planes with C outside U. All six direction slices and all nine C representative changes were checked; the source verifier checks representative changes for only one direction pair, though its proof covers all.
- A complete 9! permutation search independently finds exactly 432 magic layouts, all with magic constant zero. This additionally verifies completeness, not merely correctness of the explicit 432 constructions.
- The simultaneous gauge is injective because L∩M={0}; individual basepoint freedom does not create a factor of 3 or 9. Changing C by l+m sends g to g+l and preserves the frame. Arbitrary ambient C is valid.
- AP-only restriction is necessary: the nonline transversal example is valid. Independent checks count all q−9 disjoint gauges for that example.
- q=81 and q=729 witnesses contain zero; q=2187 does not. The packet only promises zero-allowed square witnesses generally, so its statements remain valid, but mentioning both zero-containing examples would be clearer.

## Required wording correction: paths

The statement that repeated realizations can occur across frames is false for ordinary vertex-labeled simple paths. A six-vertex path determines its unordered alternating bipartition, hence its unique frame. Extra edges do not change this. Therefore 216(q−9) AP-frame/path incidences also count distinct six-vertex ambient paths in the AP-frame subclass. They need not exhaust all ambient paths, nor are these Hamilton paths of the entire q-vertex graph. Independent explicit deduplication confirms 15,552 distinct paths for q=81.

Suggested replacement: “There are 216(q−9) frame/path incidences. Each is a distinct six-vertex simple path in the ambient square-sum graph, because its alternating bipartition recovers its frame. This is not a census of all ambient paths.”

## Independent proposed all-frame strengthening

Assume disjoint three-element parts B,D and exactly nine distinct sums B+D=A=C+U. Fixing one vertex of either part shows each part lies in a coset of U. Their nonzero difference sets must be disjoint, or two cross sums collide. A nonlinear triple in U has exactly six nonzero differences; the complement in U\{0} is a punctured line. Thus at least one part is an affine line. If one is nonlinear, the other must have that complementary direction.

For each direction L, normalize the other part's transversal T of U/L by T∩L={0}. There are 3·3=9 transversals, of which three are lines (the three other directions) and six are nonlinear. Each mixed frame therefore has a unique normalized representation

B=g+L, D=C−g+T,

with four choices of L, six nonlinear choices of T, and q choices of g. If two representations agree, the unique AP part forces equal L and gauges differing in L; T∩L={0} then forces the difference zero. Disjointness excludes exactly nine gauges because T−L=U and 2g∈C+U characterizes collision. Consequently mixed frames number 24(q−9), total frames 30(q−9), and distinct six-vertex paths arising from these frames number 1080(q−9).

The independent verifier enumerates every ambient unordered triple B, prunes only the necessary condition B−b0⊆U, and enumerates every D⊆A−b0. It checks cross-sum equality and disjointness directly, without AP/type assumptions. Results (normal and -O identical):

- q9: 0 frames and paths.
- q27: 540 total, 108 AP/AP, 432 mixed, no nonlinear/nonlinear; 19,440 distinct paths.
- q81: 2160 total, 432 AP/AP, 1728 mixed, no nonlinear/nonlinear; 77,760 distinct paths.

This is a census only for a fixed plane A with nine distinct cross sums and three distinct vertices in each disjoint part. It is not a count of all square-sum bicliques without these restrictions, nor an induced-biclique count.

## Reproduction

Run in this directory:

    python verify.py
    python -O verify.py
    python independent_audit.py
    python -O independent_audit.py
    python check_pins_fields.py
    python all_frames_audit.py
    python -O all_frames_audit.py
    sha256sum -c SHA256SUMS

All work is confined to this audit directory. The original packet and remote repositories were not modified. The all-frame candidate's producer artifact was not yet available at the time of this report; the independent result above audits the mathematical claim, not a final producer implementation.

## Final producer-extension audit

The producer artifact `/tmp/session57-affine-bare-frame-counts` was subsequently inspected and copied to `bare-extension/` here for reproduction. Its proof correctly fixes an orientation convention before defining the gauge action, proves free unordered orbits, and obtains exactly nine bad gauges from distinct cross differences. Its implementation passes normal and -O execution with identical results for all three witnessed fields, including exhaustive F81 first-triple search, 30 explicit gauge orbits, and path collision checks. Disjoint totals 2160, 21600, 65340 and path totals 77760, 777600, 2352240 match independent reasoning.

The additional claim of 36 unordered two-Hamilton-path covers of each K3,3 also passes an independent six-vertex permutation and edge-set enumeration (`check_two_path_covers.py`, normal and -O). Each cover has exactly one common edge. The union graph recovers its unique bipartition, establishing injection across frames.

Only a minor wording qualification was suggested: a proportion “one fifth of disjoint frames” presupposes q>9, since q=9 has none. The count formulas themselves hold at q=9. No remaining substantive mathematical or verifier defect found in the supplement.

Final wording fix confirmed: the q=9 proportion qualification is now present. Final supplement SHA256SUMS all pass; verifier and field implementation remain byte-identical to the audited copies. Final theorem/README/pin/hash metadata refreshed in the audit copy. Verdict: PASS for final supplement.
