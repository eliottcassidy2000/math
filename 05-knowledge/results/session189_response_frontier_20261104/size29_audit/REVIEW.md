# Independent audit of the complete 29-name boundary

**PASS for the finite root and cyclic-two-exterior certificates, the source-relative exclusion through 29 old names, and the resulting lower bound of 30 old names.** The 24-query root battery is approved for corresponding profile-equal children. The separate paired battery of at most 48 queries is approved for arbitrary root-containing original-module families using the reviewed monotonicity theorem. These are different statements.

The approved producer `PROOF.md` SHA-256 is `dc50c092907d846e7baa195a43cb063f9e2d43e849cb21e270f3bcf1000080d7`. Frozen inputs, producer identity, and inherited reviewed dependencies are recorded separately. This approval supersedes the producer's pending-review status, not the scope restrictions in its statement.

## 1. Exact structural alternatives

The imported reviewed results say that every highest changed node in a full-response collision from an ORIGINAL-source-module XOR family has a prime quotient of order at least nine, and all its children are nontransitive. The previously reviewed size-28 result already excludes sources with at most 28 old names.

For a hypothetical 29-name collision, every highest changed node has size at least 27. Since such nodes are disjoint, there is exactly one, M, and every changed old edge lies inside M. This last fact is essential: an exterior u-v edge cannot change independently. It would require another highest changed node or a change above M, both excluded by the reviewed restrictions and the selection of M.

Writing k for the number of children gives 29>=|M|>=3k>=27, hence k=9. All possibilities are exactly:

- |M|=27, with nine C3 children and two unchanged exterior vertices;
- |M|=28, with eight C3 children and one nontransitive four-child, plus one unchanged exterior vertex;
- |M|=29, the root, with either one five-child and eight C3 children, or two four-children and seven C3 children.

There is no additional child-size partition, second highest changed node, or omitted exterior edge change. Strong-module terminology is structural and does not mean strongly connected.

For one exterior vertex, both ports imitate M's incidence to it, making the exterior vertex universally first or universally last in any augmented Hamilton path. Deletion is an exact bijection. For two exterior vertices, contraction of M produces one of the eight labeled three-name tournaments: six transitive possibilities or two cyclic orientations. Ports imitating M retain the augmented core as a module. The transitive case has forced block order and exact ambient count H(S). The cyclic case is the genuinely new boundary audited below.

## 2. Cyclic exterior formula and the obstruction to scalar lifting

Let the cyclic quotient be S->u->v->S, where S is the augmented 27-name core including its two ports. A Hamilton path has either one S-run, with block order S,u,v or u,v,S or v,S,u, or two S-runs, necessarily S,u,v,S. Therefore the exact integer identity is

    G(S)=3H(S)+q2(S).

Here q2 counts ORDERED covers by two nonempty directed paths. No factorial correction or division occurs: the two run positions distinguish the two paths.

The modulo-eight identity for q2 follows correctly from the separately reviewed odd-cycle expansion H(A)=1+2O(A) modulo four. Write q2 as twice the sum H(A)H(B) over unordered nonempty bipartitions. The constant is 2^n-2, equal to -2 modulo eight for n>=3. The remaining contribution is four times the number of odd cycles wholly within one bipartition side. A fixed odd cycle on r<n vertices contributes 2^(n-r)-1 times, an odd number; a spanning cycle contributes zero. Thus for odd n,

    q2(S)=2H(S)+4-4C(S) (mod 8),
    G(S)=5H(S)+4-4C(S) (mod 8),

where C counts directed Hamilton cycles once per cyclic ordering. For even n the C term is absent. Here n=29 is odd. In particular G=H modulo four, so the endpoint-nonzero core branch lifts directly without invoking a blanket arbitrary-ambient result.

The equation retains the exact obstruction: if Delta H=4 modulo eight, it survives this cyclic exterior precisely when cycle parities agree. If core residues agree and cycle parities differ, the ambient residue can instead separate. The retained bank has 7,624 quotient/query pairs with odd cycle-parity difference. Ignoring this coordinate would be an invalid scalar-lift argument.

For nine physical C3 children, anchor a Hamilton cycle at the singleton port x. Its maximal child runs have filling factors (3,6,6). Every repeated child has even factor, so modulo two only one run per child survives. Consequently the physical 29-vertex augmented core's Hamilton-cycle parity equals that of the eleven-name augmented quotient. Anchoring at x eliminates a possible cyclic-rotation multiplicity. Triangle family changes preserve the entire profile and do not alter this argument.

The proof uses this specific cyclic two-singleton exterior and the separately treated transitive/singleton exteriors. It proves no arbitrary-core scalar lifting or arbitrary-exterior theorem.

## 3. Two independent complete ambient calculations

### Cycle covers over F2

`cycle_covers.cpp` does not use the producer's anchored path/endpoint recurrence. For every subset A of an eleven-name augmented quotient it computes det(Adj[A]) over F2. In characteristic two this determinant counts directed cycle covers modulo two. Every cycle cover has a unique cycle containing vertex zero. Consequently

    cover(S)=sum over A containing0 of cycle(A) cover(S-A).

The term A=S is isolated, and smaller subsets recursively determine the Hamilton-cycle parity. The empty cover is one; loops are absent. This is a separately implemented linear-algebra/subset-convolution method, with cycles counted once as functional cyclic permutations, not once per starting vertex.

It independently reproduces all 28,176 inherited quotient/query source-target parity pairs and the extra query at node 56186, totaling 28,177 pairs. Every parity agrees. The extra (a,c)=(2,4) has cycle parities (1,1). The count of 7,624 odd contrasts is reproduced independently.

### Direct expanded ambient path words

`ambient_words.cpp` bypasses q2, the odd-cycle expansion, and the Hamilton-cycle identity entirely. It counts Hamilton paths through nine C3 child symbols plus four individually named singleton symbols x,y,u,v, with the actual cyclic exterior arcs and x->y. Raw directed words are partitioned by their exact set of zero, one, or two repeated children. A predecessor recurrence counts unweighted words, and only after completion are their triangle cover weights applied. All terms with at least three repeated triangles vanish modulo eight, by the exact C3 profile.

The resulting normalized count is multiplied by 3^9=3 modulo eight. All 2,348 chosen cyclic-ambient witnesses are independently recounted on both sides. Every source and target residue agrees, and every gap is four. For the unique inherited-battery miss, node 56186, the added query (2,4) gives physical ambient residues (7,3).

Thus every ambient selected separator has a direct recount independent of the proposed q2/C lifting formula, while every cycle-bank entry has an independent determinant-based calculation. The formula and the direct calculations agree entry by entry on the selected domain. This does not claim an independent recount of every redundant inherited core-H bank value; the inherited full bank remains a pinned reviewed input.

## 4. Complete five-child audit and the profile-changing obstruction

`local_checks.py` enumerates all 1,024 labeled five-name tournaments, independently canonicalizes them over all 120 permutations, and recovers exactly twelve isomorphism classes, of which eleven are nontransitive, with 904 labeled nontransitive sources. No classification by degree sequence alone is used.

The complete ordered-cover profile is computed by a different formula from the producer's direct cut enumeration or its subset-convolution verifier. For a permutation with b backward adjacent arcs, every such arc must be a cut; among the remaining n-1-b forward adjacencies, choose h-1-b cuts. Its contribution to q_h is therefore binomial(n-1-b,h-1-b). Summing over all permutations gives every exact profile.

The audit finds all original modules and reconstructs every representative-family target directly by edge XOR. All 1,122 representative family outputs are checked, corresponding to 130,288 labeled source/family outputs, of which 7,408 have nontransitive sources. Profiles are coordinatewise nondecreasing, but equality fails exactly for nontransitive representative codes 2,5,8,11.

The explicit hostile example is retained and confirmed: source code 8, rows [16,1,3,7,14], has original overlapping pairs {1,2} and {2,3}. Their simultaneous XOR gives rows [16,5,9,3,14] and changes

    (9,66,186,240,120) to (15,78,198,240,120).

Thus unconditional five-child family profile preservation is false. The main collision argument correctly invokes the inherited hereditary equality theorem before substituting target child profiles. It does not infer that equality merely from the local classification.

Every five-child first coordinate is odd. Modulo eight, its normalized profile is exactly one of

    A=(1,6,6,0,0), B=(1,2,2,0,0).

The second agrees with the triangle profile for all nonzero terms. Four-child profiles normalize to (1,6,4,0), as separately reviewed. For a single five-child the physical factor is 3^8 q1=q1 modulo eight. For two four-children it is 3^7 q1 q1', also odd. Thus every normalized gap four remains physical gap four for every actual type, not only the canonical retained interiors.

Modulo sixteen there are four five-child normalized profiles, with nonzero fifth-visit weight eight. Four-child fourth visits can likewise have weight eight. The new recurrence is only modulo eight and makes no unsupported modulo-sixteen truncation claim.

The independent local audit also recomputes every named response coefficient by enumerating seven-name Hamilton orders and expanding their actual Boolean port-edge factors. All 24 vectors of 1,024 coefficients and 1,024 values agree, including zeros, 2,184 nonzero cubic coefficients, and 2,160 nonzero quartic coefficients. These are complete SMALL controls, not a full 29-name response polynomial bank.

## 5. Independent root recurrence and every placement certificate

`split_words.cpp` is a separately written unweighted predecessor recurrence using exact repeated-set decomposition. It is not the producer's weighted forward state machine or its completion-ratio calculation. For each R of size zero, one, or two, members of R have capacities three and all other child/port symbols have capacity one. Completion requires members of R to occur twice or three times and all others once. Hence the cases partition admissible words without overlap.

Raw endpoint counts are stored as bitplanes and are computed at precisions eight, four, and two, respectively, because final cover weights are divisible by 2^|R|. After counting complete words, the audit directly multiplies the chosen normalized profile weights: triangle (1,2,2), five-A (1,6,6), and four (1,6,4). No weighted transition, even-residue division, distinguishable-repeat factorial, or producer state indexer is used. Fourth/fifth visits have zero weight modulo eight; three repeated children have a factor eight.

For all 8,769 distinct quotient/query pairs needed by the selected certificates, the independent counter computes all 46 configurations on both sides. All 403,374 configuration/query records, or 806,748 residues, agree with the retained bank. This includes every one of the 105,660 new placements:

- 21,132 five-A placements;
- 84,528 pairs of four-child positions.

The 21,132 five-B placements use the independently recounted baseline column, with their actual physical positions and factors separately checked. All 126,792 selected canonical witnesses have a nonzero physical gap four. The eleven initial misses and the two added queries (0,2), (1,17) are confirmed exactly, including first-witness consistency.

The root battery has sixteen residual queries; its union with the nine endpoint anchors has 24 members. The finite quotient domain and endpoint partition are inherited from the separately certified complete prime-nine census. This audit does not claim to recertify that entire census or to independently recount all 3,024,224 redundant initial bank residues.

## 6. Physical graphs, third counting method, and retained controls

Every one of the 126,792 root physical witness pairs is checked against a separately written construction, as are all 2,348 ambient pairs. This checks 258,280 complete augmented 31-row graphs, not only hashes or coarse masks.

Root old names are 0,...,28; ports are 29->30. A five-child receives both names27,28; two four-children receive one each, in the stated ordered pair of positions. The retained canonical interiors are the five-A, five-B, and source-plus-cycle four templates. Other actual types are covered by their proven normalized profiles and odd physical factors, not falsely presented as individually retained full graphs.

For the cyclic exterior the old core is 0,...,26, exterior names are u=27,v=28, and ports29,30. Physical port masks include v and exclude u. The audit checks every old edge, port incidence, unchanged child interior, quotient-converse relation, and fixed port edge through exact row equality with the independent constructor. A common quotient isomorphism/converse transport moves both sides, all whole child sets, and both masks together; no target-only response relabeling is used.

`physical_subsets.cpp` supplies a third method, derived from the ordinary physical-vertex endpoint recurrence. Only actual directed-triangle MODULES are compressed by their rotation automorphisms. All other physical vertices, including interiors of prime five-children and strong four-children, remain individually represented. No cover normalization or repeated-child truncation occurs. The standard size-one/size-two/size-three triangle-subset equations retain the terminal endpoint of the canonical two-vertex edge, so incoming internal arcs are counted exactly.

The focused recount checks 226 actual 31-vertex graphs:200 root graphs and 26 ambient graphs. The root controls cover all 54 positional configuration types, all eleven initial misses, all eleven nontransitive five-child isomorphism types, and all nine ordered four-type combinations. The ambient controls cover the selected query types and the exceptional node56186. The largest tested state space has seven compressed triangles and ten singletons, or 2^24 physical subset orbits. Every residue agrees.

Another 11,472 actual small graphs, of augmented order at most13, compare this orbit recurrence with uncompressed physical-subset counts over the integers. They include all quotient orientations and coarse queries through three quotient children, all twelve five-child classes including the transitive hostile control, and all nine ordered four-type combinations. Every seventeenth case receives a common deterministic physical relabeling before the literal count. All comparisons agree modulo eight.

There are also1,256 independent small q2/H/cycle identity checks, including even orders, and 672 exact integer tests covering all eight labeled three-block exterior quotients over small cores. These support the derivation but are not used to extrapolate it.

## 7. Profile-equal24 versus arbitrary-family paired 48

The 24-query root theorem assumes corresponding complete child profiles equal, as ensured under a hypothetical full-response collision by the reviewed hereditary equality theorem. For a prime root, every proper original module lies in a child. Selecting the root reverses the quotient; root-induced interior converse can be absorbed by the whole-child original module. Hereditary equality therefore justifies the exact run substitution needed for the contradiction. Together with all exterior alternatives, this proves the source-relative lower bound30.

The optional family-uniform paired theorem needs a separate argument. Let D(a,c) compare the source with the quotient-only converse retaining the source child profiles. At a child-uniform query, exact full converse-complement duality and invariance of complete profiles under child converse give

    D(511 xor c,511 xor a)=-D(a,c).

A nonzero certified contrast therefore has a strictly negative member in its paired set, over the integers, even though the finite bank records only residues. The reviewed original-module monotonicity makes every actual target child profile coordinatewise at least the corresponding source profile. All occupied-word coefficients are nonnegative integers, so the actual family target count is at least the quotient-only target at that same query. A strictly negative source-minus-baseline contrast cannot be canceled by this increase.

Hence the 24-query battery and its converse-complement partners, at most48 queries, separate every root-containing original family. This does not claim that the original unpaired24 queries separate profile-changing families. Child-only families remain outside this coarse family battery and are handled by inherited highest-node restrictions in the full-collision argument.

The distinction is challenged by explicit controls for all four profile-changing five classes. For each class, all 64 coarse queries of a three-child cyclic substitution are evaluated on the source, quotient-only target, and actual profile-increasing target. These 768 exact 13-vertex physical counts verify the sign-reversing duality, target monotonicity, and preservation of strict negative contrasts. Each source has 44 nonzero baseline query contrasts. This is a control for the general proof, not its replacement.

## 8. Reproduction, inherited inputs, and stopping boundary

Run `python3 verify.py` and `python3 -O verify.py`. The complete frozen-domain checks, independent-output comparisons, physical-row reconstruction, local classification, named controls, and identity controls remain active under optimization. The final receipts agree byte for byte. `python3 verify.py --recount` recompiles and repeats all independent C++ calculations. No network or software installation is needed.

The root and ambient three-query benchmarks also pass address/undefined-behavior sanitizer replays. Leak scanning is disabled and is not claimed as evidence. All source files, exact inputs/outputs, controls, frozen dependencies, and receipts are retained in the release allowlist; binaries and routine progress logs are excluded.

Approved conclusion: a nonidentity simultaneous XOR target of ORIGINAL source modules with the same full named response requires at least30 old names. This is a lower bound, not a collision construction at 30. General named-response injectivity, arbitrary prime-nine children, arbitrary ambient lifting, arbitrary pairs of family outputs, and later-created-module operations remain outside scope.
