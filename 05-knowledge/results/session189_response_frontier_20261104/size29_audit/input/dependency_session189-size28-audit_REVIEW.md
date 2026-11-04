# Independent audit of the 28-name original-module boundary

**PASS for the finite eight-C3/one-nontransitive-four theorem, the root-containing 22-query family battery, and the source-relative lower bound of 29 old vertices.** All independently recounted selected separators, focused physical-subset counts, complete small controls, and normal/optimized verification pass. The reviewed producer `PROOF.md` SHA-256 is `360d5bb08cf33c1ceba3fc8fcef90c86916b28e35b1d08a7dd0f00e800084ba5`. Its pending-review status is superseded by this separate approval. The frozen mathematical inputs and reviewed dependencies are identified by `INPUT_PINS.json` and `DEPENDENCIES.json`.

## 1. Approved scope and inherited domain

The finite theorem concerns a prime nine-name quotient with eight directed-triangle children and one nontransitive four-vertex child, on 28 fixed old physical names. Source and target compare the same physical child sets at the same two-port query, with the quotient conversed and ports fixed as x->y. Child interiors may then be replaced by targets of simultaneous XOR families of ORIGINAL modules because their complete ordered-cover profiles agree.

This audit imports the separately certified prime-nine census, its endpoint partition, the reviewed highest-prime-node restriction, the reviewed nontransitive-child restriction, and the reviewed exclusion of at most 27 old names. It does not claim to recertify that entire census here. The previously independent degree-bucket/orbit-mass census contains 69,445 prime classes modulo isomorphism and converse: 67,097 have an odd endpoint coordinate and 2,348 have every endpoint coordinate even. The retained residual row arrays are frozen unchanged.

The all-order endpoint lemma uses the actual coarse query (a,c)=(e_i,0), not its transpose, and gives a nonzero gap modulo four for the first branch. Its arbitrary-nonempty-child scope is a separately reviewed theorem, not an extrapolation of a finite prime-eight calculation. Only the residual 2,348 by nine placements require the new modulo-eight calculation.

The resulting bound is source-relative: a nonidentity simultaneous XOR target of ORIGINAL modules with equal full named response needs at least 29 old names. It neither supplies a collision at 29 nor proves unrestricted named-response injectivity, arbitrary-pair injectivity inside a module orbit, an arbitrary-child prime-nine theorem, or validity for modules created after earlier operations.

## 2. Four-child classification, family operations, and profile normalization

Independent enumeration decodes all 64 labeled four-name tournaments, finds every original module by its uniform outside-incidence definition, and constructs every selected-family XOR target directly edge by edge. Instead of the producer's permutation-and-cut cover counter, the audit first counts directed paths on every physical subset and then convolves disjoint nonempty subsets to count ordered spanning path covers.

The 40 nontransitive sources have exactly these profiles:

- 16 labeled sources have (3,18,36,24), split into eight cycle-plus-source and eight cycle-plus-sink realizations;
- 24 labeled sources have (5,22,36,24), the unique strongly connected four-name type.

Explicit permutations of the three canonical templates give disjoint labeled orbits of sizes 8, 8, and 24, whose union is exactly the 40 nontransitive sources. Equality of profiles alone is not used as the classification test. Every nontransitive source has two nonsingleton original modules, one proper and the whole child. All 160 source/family pairs retain the source's complete profile. The proper module is respectively the cycle or the transitive pair. The local whole-child operation also absorbs the internal converse contributed by selecting the global root, so the family enumeration covers every target interior needed for the root-family theorem.

Multiplying the profiles by their odd first-coordinate inverses in Z/8Z gives (1,6,4,0) in both cases. The triangle gives (1,2,2). A four-visit contribution vanishes modulo eight. Every repeated child contributes at least one factor two, so contributions with three repeated children vanish. Neither fact discards a potentially nonzero named cubic/quartic contribution: these are statements about complete actual-query evaluations, using all allowed run words.

For eight triangle children and one four-child the normalization factor is 3^8 q1, which equals q1 modulo eight. Thus factors three and five respectively restore the physical residues. Each is odd, so a normalized gap four remains physical gap four.

The canonical four-child templates are independently confirmed as [2,4,1,7], [10,12,9,0], and [6,4,8,3] on local labels [3j,3j+1,3j+2,27]. Their labels and orientations agree with the retained physical graphs.

## 3. Independent full residual recount by exact repeated-set expansion

`split_words.cpp` is separately written and does not use the producer's weighted forward recurrence, normalized transition weights, state indexer, or completed-state positional ratio. It partitions run words by the EXACT set R of repeated old children, with |R|=0,1,2.

For each fixed R, all remaining old children and both ports have capacity one. Each member of R has capacity three. An unweighted predecessor recurrence counts all directed words on this capacity box:

    D(s,i) = sum over j->i of D(s-e_i,j),

with the singleton-word base equal to one. The loopless quotient excludes consecutive equal child symbols. At completion, nonrepeat symbols have occurred exactly once and every member of R has occurred exactly twice or three times. Hence the boxes are disjoint by their declared repeated set; no word is omitted or counted in two boxes. No distinguishable copies or factorial division is introduced.

Only after counting complete words does the audit multiply by the actual normalized profile product, using (1,2,2) at triangles and (1,6,4) at the chosen four-child position. All nine placements are computed from these raw completed-word counts. The raw counts need only precision 8, 4, and 2 for |R|=0,1,2, because their final weight is divisible by 2^|R|. Endpoint coordinates are stored as bitplanes and predecessor sums use exact population counts. This precision reduction is justified by divisibility of the complete weight, not an assumption that individual endpoint counts vanish.

The independent run recounted 6,283 distinct quotient/query pairs needed by the selected certificates. Each gives all nine special-child positions on both sides, producing 56,547 placement/query rows and 113,094 independently computed residues. Every such entry agrees with the retained producer bank. In particular, all 21,132 selected placement certificates independently have source-minus-target residue four modulo eight.

This is an independent recurrence implementation and repeated-set decomposition of the exact occupied-run bijection. It is not misrepresented as an independent physical-subset recount of all 21,132 placements; the latter, materially different physical method is used for the focused checks below. Nor does the audit claim to independently recount every redundant entry in the producer's full twelve-query bank.

All 13 cases missed by the initial twelve queries are retained and independently separated by the added queries (65,9) and (3,32). The first-witness histogram for the fixed fourteen-query order is:

    10706, 5320, 2591, 1297, 580, 0, 303,
    164, 69, 52, 24, 13, 8, 5.

Adding nine endpoint anchors gives 22 distinct queries, since (1,0) is shared. The battery is in the certified representative coordinates and must be transported with the quotient and entire physical children; it is not claimed to be the same raw fourteen or twenty-two masks under every arbitrary labeling.

## 4. Independent physical-subset method and hostile controls

`physical_subsets.cpp` starts with the ordinary physical-vertex endpoint recurrence. It quotients only by rotations of actual directed triangles, not by internal path-cover weights. Every physical subset orbit remains, irrespective of how many child visits a Hamilton path makes.

For a triangle with occupied size d, store T_i, the total endpoint count in that triangle, and, when d=2, B_i, the count ending at the terminal vertex of its canonical directed edge. Let E be the external-predecessor sum on the smaller state. The physical recurrence is:

- d=1: T_i=E;
- d=2: T_i=2E+T_i(previous), and B_i=E+T_i(previous);
- d=3: T_i=3(E+B_i(previous)).

A singleton uses the ordinary predecessor sum. The singleton-subset base count is one. These equations follow by removing the last physical vertex and using the triangle's rotation symmetry; they do not use normalized run counts or the two-repeated-child truncation.

For cycle-plus-source/sink four-children, the cycle remains one genuine rotational block and the extra vertex is a singleton. With nine cycles and three singletons including ports, the recurrence has 2,097,152 physical subset orbits. For the strong four-child, all four physical vertices remain distinct singleton groups; only the other eight cycles are compressed. With eight cycles and six singletons including ports, it has 4,194,304 subset orbits. In both cases all actual intergroup edges are checked uniform before compression, and the four-child's internal arcs remain fully represented.

The focused large recount covers 31 quotient/position placements, including every initial-battery miss, every one of the nine possible four-child positions, and every query index used by a selected witness. All three four-child types and both quotient orientations are counted: 186 actual 30-vertex graphs. Every physical residue agrees with the independently recounted run-word residue after multiplication by three or five.

A further 1,644 small actual graphs compare the orbit method against uncompressed physical-subset counting over the integers. These include all quotient orientations and all coarse queries for one, two, or three children, and all three four-child types, with a specified placement for each case. The maximum augmented graph order is twelve. Every seventh case is also relabeled by a fixed common physical permutation before the uncompressed count, testing that compression does not depend on a hidden physical naming convention. All comparisons agree modulo eight.

Finally, literal permutation enumeration verifies all 256 named two-port values for each of the three canonical four-children and their converses, and independently reconstructs all 1,536 Boolean coefficients, including zeros. The controls contain 244 nonzero cubic and 180 nonzero quartic coefficients. They are complete SMALL named-response controls, not a full 28-name coefficient bank. Another 192 exact small tests verify the uniform-singleton exterior deletion identity used in the structural argument.

An address/undefined-behavior sanitizer replay of three complete independent word-count pairs also agrees. Leak scanning is disabled because of the execution environment; it is not claimed as evidence.

## 5. Every retained physical query and orbit transport

The audit checks all 63,396 retained physical witness pairs, comprising 126,792 complete 30-row augmented tournaments, against a separately written direct construction. It checks source and target orientations, all old names, the common child sets, both expanded port masks, the fixed physical edge 28->29, each local four-child template, the odd-factor residue restoration, and the nonzero physical gap.

For a selected child j, a coarse mask m selects the three names 3i,3i+1,3i+2 for every selected bit i, and also name 27 exactly when j is selected. The target reverses precisely the between-child relations. No target-only renaming is used.

An arbitrary named quotient is isomorphic to a certified representative or its converse. One common isomorphism transports the source, target, entire child sets, and both physical masks. In the converse case the compared sides exchange, preserving a nonzero difference. An arbitrary four-child in its classified type may likewise be transported by one common internal isomorphism. More generally, profile-equal family targets agree at the same child-uniform query by the exact occupied-run substitution identity, whether or not one chooses a matching internal isomorphism. This is equality of evaluated counts on the same names and masks, not equality of two independently relabeled full response polynomials.

## 6. Highest-node classification and the one-exterior-name boundary

Suppose a nonidentity original-module target U satisfies full response equality with source T and |T|<=28. The reviewed highest-node and nontransitive-child theorems give a highest changed original strong node M with prime quotient of order k>=9 and all children nontransitive. Therefore

    28 >= |T| >= |M| >= 3k >= 27.

Thus k=9 and |M| is 27 or 28. If M has 27 vertices, all its children are C3. With no exterior, the reviewed 27-name result applies directly. With one exterior vertex v, M remains a source module with an UNCHANGED boundary: changing that boundary would change a higher quotient, contrary to the chosen highest-node condition and the inherited exclusion of changed linear nodes. There can be no disjoint second nontrivial highest changed node inside the singleton complement.

If v dominates M, set its physical port bits a_v=c_v=1. Then v dominates M and both ports on both compared sides, so every augmented Hamilton path starts at v. If v loses to M, set a_v=c_v=0, making v lose to M and both ports, so every path ends at v. Deleting v is an exact bijection with augmented core paths. This transfers the reviewed core separator as an integer equality of each side's count, not by assuming a scalar modulo-eight gap survives arbitrary substitution.

The only remaining case has M=T of order 28, with eight size-three children and one size-four child. All are nontransitive, so precisely the finite domain above applies. There is no exterior and no other highest changed node.

At a prime root every proper original module lies within a child. A changed root quotient is its converse. The actual target interiors arise from original child-module XORs, with the global root's internal converse absorbed by the whole-child module. The exhaustive four-child and triangle profile invariance therefore shows that every root-containing family has exactly the quotient-only target's value at every child-uniform actual query. The same 22-query battery separates every such family.

Families changing only child interiors are NOT claimed to be separated by this coarse battery. They are excluded as full-response collisions by the inherited highest-node restriction. Combining these branches proves the lower bound 29 for the stated source-relative original-module collision problem.

## 7. Reproduction and release

Run `python3 verify.py` and `python3 -O verify.py` to validate frozen inputs, all independent output identities, every physical record, the complete local controls, and the exact census/placement domain. Their receipts agree. Run `python3 verify.py --recount` to compile the separately written C++ counters and reproduce all independent word and physical-subset counts. This needs only the already installed compiler and Python, no network or new software. `prepare_physical.py` reconstructs the focused large and small physical inputs from frozen inputs.

`PUBLICATION_FILES.json` is the exact publication allowlist. `SHA256SUMS` pins those release files; binaries and routine progress logs are excluded. The independent census and earlier structural results remain inherited pinned dependencies. The approved new conclusion stops at the finite size-28 boundary and source-relative lower bound 29.
