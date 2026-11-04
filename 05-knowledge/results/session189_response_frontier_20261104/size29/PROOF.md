# The complete 29-name original-module boundary

**Producer finite-exact theorem and source-relative integration, pending independent review.** The proposed consequence is that any nonidentity original-module XOR target with the same full named two-port response requires at least 30 old names. This is not unrestricted response injectivity, arbitrary prime-nine child separation, or arbitrary-core ambient lifting.

## 1. Exact structural alternatives

Write F_T(a,c) for the number of Hamilton paths after adjoining actual named ports x->y, with a_v=[v->x] and c_v=[v->y]. A target is obtained by simultaneous XOR of modules of the original source, never by modules created during intermediate operations.

Assume a collision F_T=F_U with T!=U and |T|<=29. The reviewed size 28 result already excludes |T|<=28. The reviewed highest-node restrictions give a nonempty antichain of highest changed original strong nodes; each has a prime quotient of order at least nine and every child is nontransitive. Each such node consequently has size at least 27. Two disjoint such nodes cannot fit inside 29 names, so there is exactly one, M, and every changed old edge lies inside M. In particular, exterior incidences and an edge between exterior names remain unchanged. If M has k children, then

    29 >= |M| >= 3k >= 27.

Thus k=9. There are exactly these alternatives:

1. |M|=27: its children are nine directed triangles and there are two exterior old names u,v. Contracting M gives a three-name tournament, either transitive or cyclic. M need not be the root.
2. |M|=28: eight triangle children and one nontransitive four-child, with one exterior singleton.
3. |M|=29: M is the root. Its child sizes are either 5+3+3+3+3+3+3+3+3 or 4+4+3+3+3+3+3+3+3.

No case of two changed exterior vertices or an extra changed highest node is omitted. A source strong module is a module which overlaps no other source module; it is not a synonym for a strongly connected subtournament.

## 2. The exterior cases

For |M|=28, give both ports M's incidence to the one exterior singleton. Every augmented Hamilton path puts that singleton first or last, so deletion is an exact bijection with augmented core Hamilton paths. The reviewed size 28 battery therefore separates.

For |M|=27 and a transitive exterior quotient, give both ports M's incidence to each exterior name. The augmented core S=M+x+y remains a module; all three blocks have a forced order, and the ambient count is exactly H(S).

For the cyclic quotient, name the exterior vertices so that M->u->v->M, and again have both ports imitate M. Decomposing ambient Hamilton paths into S-runs gives the exact integer formula

    G(S)=3H(S)+q2(S).

Here q2 counts complete ordered two-path covers. Since |S|=29 is odd, the odd-cycle expansion proves

    q2(S) = 2H(S)+4-4C(S) (mod 8),
    G(S) = 5H(S)+4-4C(S) (mod 8),

where C is the directed Hamilton-cycle count, counted once per cyclic ordering. In particular G(S)=H(S) mod 4, so the inherited endpoint contrast 2 mod 4 lifts directly. For endpoint-zero quotients the complete inherited H bank plus independently computed cycle parity gives a finite certificate: the inherited 12 queries separate 2347 of 2348 quotients; adding (a,c)=(2,4) separates the remaining node 56186. All 2348 chosen physical 29-name witness pairs are retained in ambient/physical_witnesses.jsonl.gz.

The reason cycle parity is tractable is exact: anchor a cycle at singleton x, decompose into child runs, and multiply by the complete triangle profiles (3,6,6). Every repeated triangle contributes an even factor. Modulo 2, only one run per child survives, so C(S) parity equals the cycle parity of the 11-name augmented quotient. Triangle original-module changes preserve their complete profiles.

See ambient/PROOF.md for the derivation, exact pinned input, complete cycle bank, original physical masks and graphs, and 1296 literal small checks. Scalar gaps alone do not universally lift:7624 inherited quotient/query pairs have odd cycle-parity contrast, hence nonzero q2 contrast mod 8. This cancellation is explicitly retained rather than assumed away.

## 3. Local four- and five-child profiles

Four-child complete profiles, already reviewed, are (3,18,36,24) for a source/sink extension of a triangle and (5,22,36,24) for the strongly connected four-type. Dividing by their odd first coordinate in Z/8Z gives (1,6,4,0). Every original-module family inside a nontransitive four-child retains its source profile.

The new complete five-child classification in five_local/PROOF.md has 12 isomorphism classes,11 nontransitive, accounting for 904 labeled nontransitive sources and 7408 original-family outputs with labeled multiplicities. Every nontransitive five-child has one of two normalized mod 8 profiles:

    A=(1,6,6,0,0),  B=(1,2,2,0,0).

Profile B agrees with the normalized triangle profile for all visits relevant modulo 8. Profile A is obtained by multiplying the triangle's second- and third-visit factors by 3.

Unlike size 4, nontransitive size 5 does not have universal family profile preservation. The explicit source C3[TT3,1,1], rows [16,1,3,7,14], has original overlapping pairs {1,2} and {2,3}. Reversing both gives rows [16,5,9,3,14], isomorphic to C3[C3,1,1]. The complete profile changes from (9,66,186,240,120) to (15,78,198,240,120). Four of the 11 nontransitive classes admit profile-changing families.

For the collision argument, this is handled by the reviewed hereditary equality theorem: F_T=F_U implies H(T)=H(U) at a=c=0, which forces equal complete profiles on every corresponding old-name subset. Thus the target's actual five-child profile equals its source's under the collision assumption. There is no reliance on unconditional profile invariance.

## 4. Exact normalized root recurrence

For a child-uniform actual query, every physical Hamilton path has a unique maximal child-run word, interspersed with the two singleton ports. A child visited h times has precisely q_h physical fillings, because its h run positions distinguish the ordered paths of its spanning cover. These factors account for complete physical counts and all named cubic and quartic terms; no endpoint-product approximation or coefficient truncation is made.

Normalize each profile by odd q1 modulo 8. The baseline triangle profile is (1,2,2). Every repeated child contributes a factor 2, so words with three or more repeated children vanish. Fourth visits to four-children and fourth/fifth visits to five-children vanish modulo 8. The retained state space has 94208 states: each of the nine child occupancies is 0,1,2,3; the two ports each have Boolean occupancy; at most two children have occupancy>=2.

root_bank.cpp counts baseline normalized words by multiplying by 2 at a child's second visit and by 1 on its other visits. At each completed state it records 46 outputs:

- column 0: baseline triangle count;
- columns 1..9: replace child j by normalized profile A, using completion multiplier t(1)=1, t(2)=3, t(3)=3;
- columns 10..45: replace the lexicographically ordered pair j<k by two four-children, using the product of t(1)=1, t(2)=3, t(3)=2 at those two positions.

All multipliers are integers; no division by an even residue occurs. Both ports are counted on their fixed names with x->y. The target reverses only between-child arcs. The scalar physical factor is the product of all odd child q1 values, so every normalized gap 4 mod 8 remains physical gap 4. Other normalized five profile B uses column 0. The complete local census, not a degree-sequence guess, establishes the two five profiles.

Modulo 16 this truncation would fail: five-child fifth visits and four-child fourth visits have nonzero weight 8. No modulo 16 result is claimed.

## 5. Complete finite root domain and battery

The reviewed prime-nine census has 69445 representatives modulo common isomorphism/converse. Of these 67097 have an odd endpoint parity coordinate and are separated by an inherited child-uniform endpoint query for arbitrary corresponding profile-equal nonempty children. The remaining 2348 endpoint-zero representatives form the exact new finite domain.

The inherited 14-query battery is

    (0,1), (1,0), (1,2), (1,3), (3,1), (3,3), (3,5),
    (0,7), (0,18), (40,256), (4,3), (130,128), (65,9), (3,32).

root_residues.tsv retains every source and target normalized residue for all 46 columns at every one of these 14 queries on every residual quotient:32872 rows and 3024224 individual residues. It covers all but 11 of the 105660 new configurations (21132 exceptional-five placements and 84528 two-four placements). The exact 11 misses are retained in root_initial_misses.json.

Adding two queries, (0,2) and (1,17), covers all 11. root_selected_residues.tsv retains both queries for all 46 columns on all 11 affected quotient nodes, including redundant entries. No result is inferred merely from an earlier C3 or size 28 witness. The normalized-B five placements inherit column 0 of the same complete bank, producing 21132 further retained witnesses. Thus root_witnesses.tsv has 126792 chosen configurations, each with a physical gap 4.

The fixed 16-query root battery together with all nine endpoint anchors has 24 distinct queries, since (1,0) is already present. This battery separates every root quotient-converse with corresponding complete source/target child profiles equal, for either 29-name root partition.

## 6. Physical labels and families

Every root physical record has old names 0,...,28 and ports x=29,y=30. Child i initially contains 3i,3i+1,3i+2. A five-child receives both 27,28. A pair j<k of four-children receives 27 in j and 28 in k. Coarse masks expand to exactly these fixed physical child sets on both sides.

The canonical retained interiors are: exceptional five rows [4,1,2,7,15], profile(3,42,162,240,120); inherited five rows [16,1,3,7,14], profile(9,66,186,240,120); and four rows [2,4,1,7], profile(3,18,36,24). Their normalized odd physical factors are respectively 3,1,3. root_physical_witnesses.jsonl.gz retains both complete 31-row augmented graphs, masks, normalized residues and physical residues for each of 126792 cases, totaling 253584 graphs. These are canonical realizations; other child isomorphism types are covered by their proved normalized profiles and odd physical factors, not falsely listed as individually retained graphs.

An arbitrary named quotient is transported by the same isomorphism and child permutation on both source and target and on both masks. Converse exchanges source and target. The positions and every literal interior remain paired; no independent response relabeling is permitted.

A proper original module meeting a prime root lies within a single child. Reversing the root oddly therefore converses its quotient; descendant operations change only child interiors. Under the assumed full-response collision their complete profiles match, so exact run substitution identifies their actual query count with the quotient-only target. The 24-query battery contradicts equality. Combined with the exterior alternatives, this excludes every nontrivial collision on 29 old names, subject to independent review of this packet and its pinned imports.

A separate family-uniform root statement is available without assuming profile equality: pair the 24 root queries with their core converse-complements J(a,c)=(1-c,1-a). The reviewed original-module monotonicity and exact augmented converse-complement duality make at least one baseline contrast strictly negative; increasing the target child profiles cannot cancel it. Hence at most 48 paired queries detect every root-containing original family relative to its source. This is distinct from the 24-query profile-equal collision argument and does not detect child-only families by itself.

## 7. Verification and scope

Run python3 verify.py and python3 -O verify.py for complete producer replay and physical validation. Five-local and ambient certificates have their own normal/optimized verifiers and receipts. root_crosscheck records a distinct full-occupancy recurrence on 33 physical pairs: it retains fourth/fifth visits and arbitrarily many repeated children, counts unweighted run words first, then multiplies complete unnormalized profiles. All 66 counts agree. This bounded cross-check is not represented as exhaustive independent review.

Complete small named response vectors and every Boolean coefficient, including zeros and cubic/quartic coordinates, are retained in five_local/named_controls.json:24 vectors and 24576 coefficients, with 2184 nonzero cubic and 2160 nonzero quartic terms. The inherited four-child controls remain pinned. The large root and ambient records are actual full-query values, not a complete named polynomial bank at 29 vertices.

This packet proves no arbitrary-child prime-nine theorem, no arbitrary-ambient scalar-lift theorem, no pairwise injectivity among arbitrary family outputs, no operations on later-created modules, and no unrestricted named-response injectivity. Its source-relative lower bound is finite and conditional on its explicitly reviewed inputs.
