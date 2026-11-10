# Session244 independent adversarial audit

Status: **PASS within stated scopes**. Chan–Pak remains **OPEN**. No all-target existence conclusion follows from these packets. Audited 2026-11-10; producer files were not edited.

**Final audit sealed:** the order-twelve census and independent normalized shallow exclusion are complete. The fixed-order depth-two scalar claim is refuted only at H=97309 in that census; a depth-one order-twenty witness is supplied. The final section supersedes earlier pending-run status notes.

## Scope and independent method

Reviewed `selector/proof.md`, `logic/PROOF.md`, `composition/PROOF.md`, and the outdegree-one structural continuation in `composition/STRUCTURE.md`. The closest inherited mechanisms are the Session199 physical O/Z/J maps and Session243 compact gap compiler; the hostile distinctions are response versus full-count selection and scalar abundance versus physically placed witnesses. The sidecar that must survive is the full rooted profile, not just the difference or endpoint marginal. The audit board was: asymptotic counting; target-bounded pruning; scalar countermodels; physical matrix identities; saturated-root necessity.

`independent_verify.py` imports no producer code. It enumerates directed permutations of every subset of every labelled tournament of orders 2 through 5 (1,098 tournaments), computes total/start/end counts, and independently recomputes the complete seed-bank inverse through 2001. This differs from the producer's subset dynamic programming. All 3/10/42 seed-profile counts, 545/789/985 image counts, expanded-state counts, and holes agree exactly. In particular the order-at-most-five retained-root family omits 105 and 173 at every tail length. It additionally checks the strengthened cover-deficit inequality on every root, exact density enclosures, 10,000 scalar compact witnesses, and 9,000 constant-error anchor witnesses.

The selector, logic, and composition producer verifiers replay byte-for-byte against their saved JSON; optimized Python runs match ordinary runs. The structural verifier also matches its saved output. Reproduction:

    python3 -I -B independent_verify.py > INDEPENDENT.replay.json
    cmp INDEPENDENT.replay.json INDEPENDENT.json

The independent verifier uses explicit exceptions, not removable assertions. Producer replays are stored here for the audit trail. Finite exact checks do not establish the unbounded theorems; those were separately checked by proof.

## 1. Selector audit

**PASS, PROVED:** The valuation estimate counts inputs rather than distinct outputs, so collisions only strengthen the upper bound. For N=2^j m with odd m, g(m)≥(m+1)/2, including m=1. Thus L(N)≤X implies m≤(2X−1)/(2^(j+2)+1). The number of odd m contributes at most X/(2^(j+2)+1)+O(1). Only O(log X) valuations occur since N≤(X−1)/2. Hence the uniform error is O(log X), and passing from the infinite summable series to the upper density is legitimate. The relative omitted fraction among odds is at least 1−2c; this says nothing about the full tournament spectrum.

**PASS, PROVED:** Pair-image redundancy R(N)=L(2N) includes the exceptional N=1 convention. The all-one-apex classification legitimately strips global sources before the first apex-outneighbor and trailing sinks. The residual final Z forces odd U, giving the finite inverse h≤U≤2h−1. The h=13 exclusion is for this whole class, not an order cutoff.

**PASS, PROVED:** The finite-seed inverse does not silently impose a length bound. A final nontrivial Z yielding h has predecessor response U≤h because D≥E. Trailing O moves do not change D. Every earlier response is smaller, so pruning U>h cannot remove a successful branch. U_k≥2^k(U_0−d)+d and U_0≥2d give a genuine finite search for each h. Profile deduplication is sound because subsequent maps depend only on those four coordinates.

For completeness, the two physical inequalities used in the pruning have direct combinatorial proofs. Cut a full Hamilton path just before r: its suffix begins r and its prefix is the other cover path, giving D≤U. Also U_r(T)=d+sum_{a∈N+(r)} U_a(T−r): separate the singleton distinguished path and otherwise remove its first vertex r. Applying U_a(T−r)≥d and using that a root reaching another vertex has an outneighbor gives U≥2d. This proves the inequalities beyond the small census.

## 2. Logic audit

**PASS, PROVED:** The scalar model excludes 7,21 and all primes at least 151. Multiplication cannot create a forbidden composite 21 without a forbidden factor 7; a product of nonunits cannot create an excluded prime. For N≥30 the interval [2N+1,3N] contains exactly N integers and hence a full residue class modulo 30. The chosen endpoints are proper odd multiples of 3 and 5, greater than their divisors and avoiding 7,21. The small cases in the proof, including its different valid witness at N=5, obey b≤5N.

**PASS, PROVED:** Fifteen consecutive odd integers in [m−28,m] contain each odd residue modulo 30 once. For m≥51, the selected multiple of 3 is at least 23 and thus proper; its partner is also a proper multiple of 5. Uniform additive error 28 therefore coexists with infinitely many missing prime targets. No inference about every smaller positive error or physical tournament identities is licensed.

The atom classification and odd-parent conditional forest are valid. Their roles are respectively an inherited arithmetic reduction and an expressly unproved physical closure obligation. Neither asserts the missing physical invariant.

## 3. Composition audit

**PASS, PROVED conditional on inherited transfer laws:** Full four-coordinate multiplication gives the two stated compiler outputs. Final Z closes a root-reaching nontrivial source into a strong graph; the sources used have enough vertices for the nonroot intermediary. The shared E=d and U=D output conditions, gap 2E, self-feeding face maps and two gcd identities are exact. The initially printed inverse equation was wrong; the second audit pass records its correction below. Coprimality requires its stated primitive and 5∤e hypotheses. Trial division certifies the listed three prime pairs, with no infinitude conclusion.

The outdegree-one structural formula h=H(A)H(B)+U_a(A)Q(B), e=H(A)H(B), and its equality boundary h=3e are valid. For |B|≥2, at least one proper-subset term in Q(B) is positive. Small endpoint rigidity e=1 and e=3 follows with the stated inherited count-three classification. A cone-only inverse is not a physical classification.

## 4. Independent pair-cover proof of an inherited necessity bound

The following argument was also checked during the audit. It is valid and supplies an independent proof of the complete nontrivial strong saturated-root class. The classification and a stronger deficit bound are inherited, as documented in the provenance correction below.

**Pair-cover lemma, PROVED.** Given a tournament A and distinct a,b, there is a spanning ordered two-path cover in which one path starts at one mark and the other ends at the other mark. Let C1⇒...⇒Ck be the strong components of A−a, with b∈Ci. Hamilton paths in a strong tournament can start or end at any prescribed vertex (cut a spanning directed cycle; singletons are immediate).

- If i=k, use the singleton (a) and a Hamilton path of A−a ending b.
- If i<k and a beats some c∈C_(i+1), use the path starting a, then a Hamilton path of C_(i+1) starting c, then the remaining later components. Its companion traverses C1 through Ci and ends b.
- Otherwise C_(i+1) beats a. Use a path starting b in Ci followed by C_(i+2) through Ck. Its companion traverses C1 through C_(i−1), then C_(i+1), then a.

All concatenations follow the component order, and empty prefixes/suffixes cause no exception. The use of the standard strong-tournament Hamilton-cycle theorem should be named in the final theorem dependency list.

Cutting a Hamilton path just before r identifies D with exactly those rooted covers whose other path is empty or ends at an inneighbor of r. Thus U−D counts the remaining covers. Write O=N+(r). The singleton-root covers contribute sum_{a∈O} e_a(T−r). Each unordered pair {a,b}⊂O supplies a further cover by the lemma, after prepending r to its marked-start path. Different unordered pairs give different covers, since the first vertex after r and the other path's endpoint recover that pair. Therefore

    U−D ≥ binom(|O|,2) + sum_{a∈O} e_a(T−r).

If |O|=1, every nonsingleton rooted path consumes the sole outneighbor, leaving no outneighbor available as the other path's endpoint. Hence equality simplifies to U−D=e_a(T−r).

For a **nontrivial strong** tournament, outdegree is positive, so U=D forces outdegree exactly one. Then equality holds exactly when a is not in the last strong component of T−r. Strongness forces a into its first component, yielding precisely the A⇒B decomposition with B nonempty and the structural formulas above.

**Required boundary correction:** the singleton tournament is strong under the usual convention and has U=D=E=d=1, outdegree zero. Accordingly the necessity/classification and h≥3e must explicitly exclude it. This was corrected before promotion. No other proof defect was found.

### Parity strengthening, also audited

A subsequent valid strengthening gives the following argument. Let K_ab count the covers starting a and ending b in separate paths. It counts the Hamilton paths using adjacent arc b→a after those paths are joined. K_ab and K_ba are unaffected by flipping arc ab because the marks lie on different paths. Consequently flipping that arc changes the full path count by K_ba−K_ab. Rédei oddness makes this difference even, so K_ab+K_ba is even. The pair-cover lemma makes it positive, hence at least two. Therefore the stronger universal bound is

    U−D ≥ |O|(|O|−1) + sum_{a∈O} e_a(T−r).

The independent permutation verifier checks this stronger bound, not just the weaker binomial bound, on every root of every labelled tournament of orders 2 through 5. The unbounded proof uses the inherited oddness theorem, not extrapolation.

## Consequence and stop boundary

The audited work strengthens the explanation of why even-gap existence is insufficient and gives genuine structural information about a self-feeding physical face. It does not establish a target-solving selector or stable inverse family for every allowed odd integer. Chan–Pak, the required all-target anchor, and any broader physical completeness conclusion remain OPEN.


## Round 2 audit continuation

**Correction to this audit:** The first pass incorrectly endorsed the composition inverse equations while verifying the forward maps. The logic cross-check correctly found a sign reversal. The correct F− inverse is h=(h′−3e′)/2 and e=3h′−10e′. Substitution gives (9,3) from (171,51). The composition proof has now been repaired. Forward maps and all earlier forward numerical results are unaffected; do not rely on the superseded first-pass inverse endorsement.

**Complete saturation STRUCTURE, PASS:** The alternative reachability proof is sound: two outneighbors whose deletions each destroy root reachability lead to disjoint nonempty unreachable sets, each forced to beat the other. The root-reachability/Hamilton-start equivalence correctly uses ordered strong components and Hamilton cycles. The arbitrary-SCC formula U_r(P⇒C⇒B)=H(P)U_r(C)Q(B) is exact: vertices of P must all belong to the complementary path, the C portions supply the rooted cover, and B splits freely between the two paths. The independent permutation verifier now checks this formula for every root of every labelled tournament of orders2–5. Terminal-singleton saturation is explicitly retained. The prime-target gcd restrictions are necessary, not sufficient, as stated.

**Logic ROUND2, PASS:** The difference cone is preserved by every letter; the first J creates positive deletion difference which no later letter can erase. With J excluded, U doubles each step and the full-count coefficient catches up exactly only for all-Z words. The asymmetric universal face-return argument correctly relies on independence of h,e,g and is not extended to target-specific fibres. Same-branch inverse synchronized gaps require divisibility4 for F− and28 for F+. Cross-branch cancellation yields exactly (6(7E−D−U),12E) in the face coordinates, so integral common continuations cannot produce a total gap ±2. These are restrictions on the named alphabet/continuations, not on arbitrary tournaments.

The one-prime countermodel S_p indeed contains the earlier model and all the stated composite witnesses. For an empty finite set F, the cutoff can be written max({151}∪F), avoiding undefined max(F); this is a harmless notation convention rather than a mathematical obstruction.

**Face verifier, PASS:** The C++ subset recurrences for start/end counts, pair-cover coefficients, exact defect decomposition, strengthened bound, and strongness are correctly indexed, including singleton and empty-complement cases. The -DNDEBUG build with C++17, -O2, -Wall, -Wextra, -Werror exactly reproduces `face/CHECK.json`. Its --fail control exits nonzero (134). This independently confirms the normal producer replay is not dependent on assertions. The census is finite-exact through order6; the proof carries the global statement.

**Selector second-round computation, PASS:** `verify_round2.py` replays exactly against ROUND2.json. Its O/Z/J state search legitimately deduplicates full profiles and has no imposed word-depth cutoff. J strictly increases U by D, while the earlier final-Z cutoff applies also to final J because its full count is the same D+U−E. Seed counts and final letter counts are included, so trailing O words lose no count. The new literal eight-vertex H=173 witness is independently reconstructed from the edge, word ZOOZ, reroot2, suffix ZZ, and then checked by all directed subset permutations: profile (173,9,209,9); all eight roots reach the entire graph. This certificate is stronger than the producer's recurrence-only construction and is included in INDEPENDENT.json.

The freely rerooted potential argument is valid with an explicit reaching-root guard and seeds of at least two vertices. Adding a sink preserves exactly the old reaching-root set and doubles every corresponding U, hence doubles μ. After a Z at a reaching root r, D′=D+U_r−E_r≥U_r≥μ(T); every output rooted response is at least D′, so μ cannot decrease even if new reaching roots appear. Every Z increases D by a positive even integer, hence at least2. Before the final Z, μ≤U_r≤h, bounding the number of preceding O moves by log2(h/μ_seed); all intermediate counts are nondecreasing, bounding the number of Z moves by (h−D_seed)/2. Trailing sinks can be removed, leaving a finite graph search. This does **not** prove the same claim for unguarded arbitrary roots: Z at a sink can preserve D, invalidating that length bound. A separate reduction would be required for unguarded roots.


## Provenance correction: proof validity is not novelty

A lineage audit found prior results that the first proof audit did not adequately credit. This does not invalidate the displayed mathematics, but changes its classification from a proposed new structural advance to **INHERITED with independent alternative proofs and exact checks**.

- `05-knowledge/results/session231_overnight_completion_20261109/descent/ROUND2.md`, §2, already proves the complete nonsink zero-defect face, h≥3e, its equality condition, and the e=1 dyadic family. This was verified directly from the tracked HEAD blob because that path is absent in the sparse working tree. Its statement is broader than the strong-only version.
- THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, (A1), explicitly credits THM-4873 and gives the stronger bound σ_r≥E_r+2Σ H(T−{r,a,b}), summed over outneighbor pairs. Since each residual tournament has at least one Hamilton path, the Session244 degree bound is a consequence. Here inherited E_r=H(T−r)−e_r(T) equals the present e_O(T−r), by appending r to a deletion path exactly when its endpoint is an inneighbor. The prior weighted bound must take precedence over the weaker rederivation. The bound was also traced to `two_entry_funnels_20261103/READY.md`, equation(7).
- Targeted searches of that lineage found no explicit e=3 divisibility slice or nonsink saturated 11/13 exclusion. Their safe status is **elementary corollaries of inherited classification**, with no claim of literature priority or exhaustive repository novelty. The prime-gcd restrictions are likewise elementary divisibility corollaries. The closely related primitive-triple argument is already in Session231 `descent/PROOF.md`; those restrictions are not a fundamentally new primitivity mechanism.

The new pair-cover and parity arguments may be retained as independent proofs, not promoted as first proofs of the inherited conclusions. The original audit's mathematical PASS therefore remains, but is not a novelty endorsement.

Selector `ROUND2.md` is now present and was reviewed in full. Its explicit reaching-root restriction, slot-count bijection, k=0 boundary, O/Z/J target cutoff, finite seed-bank scope, literal173 construction, and freely rerooted inverse potential all pass. The complete grammar search is not promoted to universal target coverage.

## Round 3: efficient endpoint source

**PASS, PROVED within the stated construction.** `selector/ROUND3.md` has now been reviewed. The exact profile follows independently by deleting the O-newborn v and applying the J formulas at the final newborn z: d=g+N−ε; endpoint_z=g−ε; hence endpoint_v=N. Since v uniquely beats z, U_v=d+U_z=2g+2N−ε. This proves all coordinates. The stabilizer proof fixes v, then its unique outneighbor z, then z's unique outneighbor old root; the remaining transitive chain is rigid. The N=1 triangle convention is handled separately and agrees with the profile formula.

The cost inequalities, powers-of-two lower equality, all-Z upper equality, order comparison, and synchronized-compiler substitution are correct. In particular the upper synchronized count is 31g+39N−11ε≤70N. This is a cost bound for an endpoint/deletion-synchronized pair, not an improvement on the different compact 5N pair objective.

The fixed-root Z-tail formula H_t=(2^t+1)(g+N)−ε is exact. Its inverse range ceil(A/2)≤N≤A−1 follows from N+1≤g+N≤2N; finite t bounds follow from g+N≥2. The even-N branch is a nontrivial product, and the odd-N branch cannot reach a Mersenne prime because h+1 is a power of two. These are family-specific arithmetic obstructions, not full-spectrum obstructions.

`verify_round3.py` replays byte-for-byte against ROUND3.json. In addition, the audit's permutation verifier now independently reconstructs and checks the entire efficient endpoint profile for N=1,...,16, without importing any producer module. The reported 16 independent source checks are stored in INDEPENDENT.json. The quantitative theorem is supported by its proof, not by this finite range.

Novelty scope remains appropriately modest: endpoint and trivial-stabilizer existence are inherited; the packet states the retained full profile and linear-cost improvement relative to the specifically inspected prior route. This audit does not certify literature priority or exhaustively exclude other previously known efficient constructions.

## Further inheritance correction: compact-image density

Direct inspection of `git show HEAD:session199/hidden-structure/PROOF.md`, §4, confirms that the exact constant c=Σ1/(2^(j+2)+1), the O(log X) error, and the omitted odd fraction 1−2c=0.137667... were already proved there for the same balanced image. Session244's density argument is therefore **INHERITED / independent rederivation and application**, not a new density theorem. Its mathematical validity remains unchanged. Section3 of that source also already proves the chronological caterpillar reversal forest for fixed O/Z/J words. A general forest description for freely changing roots is a different graph-language observation, but must not claim the fixed-word caterpillar result as new. THM-5114 correctly identifies its substantive additional content as the scalar countermodel and logical nonimplication, not the inherited density bound.

## Round 3 logic: six forbidden cores and exact join repair

**PASS, PROVED structural statements plus exhaustive finite core classification.** Inclusion-minimal minimum-outdegree-k tournaments have every vertex beaten by an outdegree-k vertex. Counting arcs from that set gives n≤(k+1)s−s(s−1)/2, whose integer maximum is (k+1)(k+2)/2. For k=2 this bounds all minimal obstructions by6; orders below5 are impossible by the degree sum. The six representative orbits are disjoint and exhaust every candidate through6, so the finite census combined with the size theorem is an all-order forbidden-core characterization.

The unrestricted grammar is precisely outdegree-1-degeneracy by reversing the elimination ordering. The guarded grammar is strictly smaller; its transitive-or-initial-strong-component-plus-sinks invariant follows from each guarded Z on at least two old vertices making the graph strong. The singleton exception is explicitly handled: Z on the singleton produces an edge, not a strong two-vertex tournament. The strong five-vertex obstruction has a unique possible last Z vertex, whose deletion is source⇒C3, so alternate construction histories cannot repair it.

Adding ordinal joins gives exactly the full degenerate graph class. Closure is hereditary and uses a low-outdegree vertex in the last occupied join block. For construction, nonstrong graphs split into proper components; a nontrivial strong degenerate graph has a sole-exit vertex r. A simple path from r to any old vertex proves its unique outneighbor reaches the whole deletion graph, validating the guard. Singletons and edges are covered. No scalar-spectrum strictness follows: all six forbidden-core counts have permitted witnesses, as the packet explicitly states.

`audit/verify_cores.py` independently regenerates the labelled graph classes by **forward insertion at every label position**, rather than the producer's backward cached deletion recursion. It exactly reproduces unrestricted counts1,2,8,64,1000,26720 and guarded counts1,2,8,56,760,18960. It independently checks the six permutation orbits and tests every potentially offending proper induced subset directly. Results are in CORES.json. The producer's optimized round3 replay also matches ROUND3_CHECK.json. This is a distinct computational path, not just a replay endorsement.

## Canon-summary audit

THM-5113, THM-5114 and THM-5115 were read in full against the audited sources. Their mathematical statements and all-order scopes pass. THM-5113 states endpoint realization, not arbitrary full-count realization. THM-5114 explicitly marks the density bound inherited and the countermodels abstract. THM-5115 correctly retains synchronization/continuation restrictions and the corrected inverse sign, while crediting inherited saturation structure. One wording improvement was adopted: say “the gcd of the two branch total counts” rather than “its two output total counts,” since Fminus by itself supplies only one total-count coordinate. No substantive theorem correction is required. These three summaries do not claim Chan–Pak completeness.

## Pure-count compiler and sharp general cores

**Composition ROUND3, PASS.** The seven-letter words have the stated full four-coordinate matrices. Their difference is (2U,0,0,0), so they preserve equality of the entire endpoint/response/deletion sidecar. Both O and Z fix the pure-count direction, proving arbitrary common fixed-port O/Z suffixes preserve the chosen gap exactly. This is stronger sidecar synchronization than the earlier face-direction compiler, but is still not an arbitrary anchored target constructor.

The singleton boundary works: opening O gives enough old vertices before the terminal Z, while J preserves reachability from the marked root. The final marked root has two outneighbors, one inherited from its last J and the O-newborn after it; the final Z points into it. Therefore this channel intentionally lies outside the saturated outdegree-one face. The response selector yields L=17g+40N, R=17g+42N≤59N and the stated order bounds, including the separate eight-vertex N=1 construction. Since common U−d>0, repeated Z grows the anchor unboundedly without changing the gap, but supplies only the stated constrained offset family.

The independent permutation verifier now reconstructs and checks all four literal compiler outputs for N=1,2 (orders8 and9), including full profiles and strongness, without producer imports. Its N=2 outputs certify97/101 with common sidecar(15,121,15). The optimized producer replay matches ROUND3_RESULTS.json. Inherited equal-scalar rebases are explicitly credited; the audit does not claim a global minimality or literature-priority result for this compiler.

**Logic ROUND3 general sharpness addition, PASS.** For k≥2, a strong tournament S of order k+1 exists, and each ownership multiplicity k−d_S(u) is nonnegative. Old vertices have exactly k total outneighbors; a new vertex has k outneighbors in S regardless of its position among the new transitive vertices. On deleting a new vertex its owner becomes peelable; on deleting an old vertex the cascade has already begun. Reverse reachability through strongly connected S then peels every old vertex, since each removed old outneighbor reduces an old vertex's original degree k below k. The remaining new vertices are transitive and peel by sinks. Thus every single deletion is (k−1)-degenerate, which by heredity proves inclusion-minimality for every proper induced subset. The k=1 triangle handles the only unavailable strong-core order2 case. The bound (k+1)(k+2)/2 is genuinely sharp for every k≥1. The updated optimized verifier, including sharp examples through k=8, matches its saved JSON.

The forest lineage statement is now accurate: the fixed-word chronological caterpillar forest is inherited Session199; arbitrary earlier-parent choices in the unrestricted freely rerooted grammar give the broader recursive forest language. No claim that the restricted and unrestricted history classes coincide is made.

## THM-5116 and ordered-forest census implementation

**THM-5116 exact summary: PASS for promotion from RESERVED.** The sharp all-k bound, all six representatives, guarded-plus-join equality, singleton boundary, and separation of graph exclusions from scalar counts agree with the audited full proof. Its inherited fixed-word forest language is correctly credited, and no all-value conclusion is asserted.

**C++ forest census code review: PASS through its explicit limit12.** At birth of vertex v (zero-based) there are v+1 choices: a sink or one old outneighbor. The nested loops enumerate their product n! exactly, including all choices at the final level computed from the exact U recurrence. These are all chronological forest histories, not all tournaments or an isomorphism census. Every outdegree-1-degenerate tournament has at least one such history after relabeling, so the scalar image is complete for that graph class at the specified order.

The incremental start/end DP is sound. A deeper recursion only overwrites masks containing a new higher bit; old-mask values survive unchanged. Restoring out/in arrays restores the parent's graph before each sibling. Root-depth entries for existing vertices are not overwritten by descendants. Values for endpoints not belonging to a mask are never read. The direct witness auditor also safely reuses its DP array because every used entry is recomputed in increasing-mask order for each reconstructed graph.

The uint32 storage is safe for n≤12. Chronological codes plus the sentinel offset are at most12!=479001600. Full path counts are at mostn!. Rooted two-path covers on m vertices are at mostm!: for each distinguished-path length j, ignoring all orientation restrictions leaves (m−1)! possibilities, and there are m lengths. The shortcut's parent has m≤11, so H+U≤2·11! fits comfortably; individual products and partial sums are nonnegative and bounded by the same combinatorial totals. Forest depths are at most11 and fit the byte storage.

Minimum depth is minimized **per scalar value across all enumerated chronological histories of the specified order**. Recounting the min-depth code confirms an attaining witness, while exhaustive enumeration establishes that no smaller depth was missed. Neither alone should be described as a minimum-order theorem for arbitrary tournaments.

A separate audit compiled and ran the census through8 and saved `forest_n8.replay.jsonl`. `verify_small_forests.py` independently enumerates all histories through6 and counts Hamilton paths by full permutations, not subset DP or transfer formulas. It exactly matches both value sets and their minimum depths. Result: SMALL_FORESTS.json. Reproduction: compile the producer census, run it with arguments `8 OUTPUT_DIRECTORY`, then run `python3 -I -B verify_small_forests.py OUTPUT_DIRECTORY`. The ongoing order12 computation is not certified as completed by this code audit; its eventual output and witness checks must still be verified before reporting its results.

## Fresh post-interruption audit, 2026-11-10

A second independent pass reread THM-5113 through THM-5116 and the associated
proofs. **PASS within the stated scopes.** No new false mathematical statement
was found. The all-odd-except-7-and-21 conclusion remains OPEN.

**THM-5115's appended exact no-J criterion passes.** Before the first J, a
pure-count difference q(1,0,0,0), q>0, stays unchanged under O and Z. The first
J sends it to q(1,1,1,1), exactly q times the singleton's physical rooted
state. Every later physical word is the state of an actual nonempty rooted
tournament, whose deletion Hamilton count is positive, including the empty
deletion convention at order one. Therefore no later letter can restore a
zero deletion difference. This proves necessity for arbitrary finite words;
the two displayed O/Z vector identities prove sufficiency. It is not merely
a check that one J fails immediately.

**THM-5116's sharpness and grammar proof pass again.** The quadratic maximum
occurs at both integer values s=k+1 and s=k+2. A strong S of order k+1 exists
for every k>=2, for example a transitive tournament with its long edge
reversed. Each owned vertex has k outgoing arcs into S; its owner's degree
is exactly k. After any deletion, reverse reachability through S causes all
of S to peel, then the transitive remainder peels. Heredity passes this
single-deletion conclusion to every proper subset. The strong-construction
converse uses a *simple* path starting at the sole-exit vertex, so its suffix
really lies in the deletion graph; no path through the removed vertex is
being silently assumed. The singleton and two-vertex boundaries are retained.

`final_independent_verify.py` imports no producer code. It checks all 88,573
O/Z/J words of lengths zero through ten for the no-J criterion; explicitly
constructs sharp cores for k=1,...,20 and peels every single-vertex deletion;
and counts the literal eight-vertex 173 witness by all 8! permutations. Its
outgoing masks are (218,252,57,240,224,193,132,4), and its start counts are
(23,72,5,27,17,14,10,5). Thus the witness does **not** have a unit start. It
cannot be used to refute Session197's rooted unit-start grammar result.
The 173 count was already known in broader settings; the credit here is the
specific freely rerooted construction/certificate, not discovery of 173.

Both ordinary Python and `python3 -I -B -O final_independent_verify.py`
produce the exact `FINAL_INDEPENDENT.json`; the `--fail` control raises under
both modes. Fresh ordinary and optimized replays of the ten selector,
logic, and composition verifiers through Round 3 all match their saved JSON
byte-for-byte. The older independent full-permutation verifier and forward
core generator also reproduce their saved results.

**Census code and integer ranges pass again.** For n<=12, codes plus the
sentinel offset and full path totals are at most 12!=479001600. A rooted
ordered two-path cover on m vertices has at most m! possibilities: for each
rooted-path length, forgetting orientations gives (m-1)! possibilities, and
there are m choices of length. The final-layer parent has m<=11, so even
H+U before subtraction is at most 2*11!=79833600. Every nonnegative DP
product and partial sum is bounded by its combinatorial total. No uint32
wraparound is needed or possible within the stated limit. Forest depth is
at most eleven and fits the byte field. Parent subset tables and depth
entries survive deeper recursion; restoring adjacency arrays suffices.

A fresh C++17 `-O2 -DNDEBUG -Wall -Wextra -Wconversion -Werror -pedantic`
build reproduced the complete order-eight summary exactly. A separately
compiled direct-endpoint-DP auditor recounted both stored codes for all 277
CSV rows, 554 literal witnesses. The independent full-permutation forest
census through order six again matches every spectrum and minimum depth.
An order-thirteen request and an intentionally corrupted CSV path count both
fail nonzero under `-DNDEBUG`. New replay files are
`final_forest8.replay.jsonl` and `final_forest8.witnesses.jsonl`.

At this stage, this fresh audit does **not** certify completion of the ongoing
order-twelve run. A running log or empty/truncated JSON file is not a completed
census. Its final certificates and new Round 4 proofs require separate review.

The fresh checker now also checks 855 literal one-letter graph moves on every
rooted labelled tournament of orders one through four, counting paths by
permutations rather than importing any transfer-law implementation. This
caught a transcription error in the audit's initial J endpoint assignment;
the corrected assignment is D-E. The producer and theorem already used the
correct law. After correction, the full ordinary and optimized replays agree,
and both deliberate failure controls remain effective.

Direct rereading of Session199 `hidden-structure/PROOF.md`, Sections 3 and 4,
confirms the inherited chronological caterpillar theorem and the exact same
density constant and O(log X) bound. The raw `overnight_s3.out` lists
`173: NEVER (forbidden)` inside its order-at-most-seven calculation. That
phrase is a misleading historical bounded-miss label, not an all-order
nonexistence result, and is not used here as one.

The restarted order-twelve census writes to `logic/forest_census/release/`;
the old `N12_DEPTH.*` files are interrupted raw artifacts, not current run
status or completed evidence. Final release certificates remain pending.

## Dyadic root-menu experiment

**PASS, scoped inherited identity plus explicit failure of a proposed
invariant.** `dyadic/PROOF.md` correctly specializes fixed-root Z: when
E=d and U=D+d+c, the new state has D'=2D+c and U'-D'-d=c. Iteration gives
D_t=2^t(D+c)-c. The two explicit children of the order-four H=5 graph do
lose the opposite sign even after checking every output root. The claim is
about these chosen extensions, not about every possible rebase or a global
impossibility of binary induction.

Both Python modes reproduce `dyadic/CHECK.json`, and both `--fail` controls
fail. A fresh C++17 `-O3 -DNDEBUG -Wall -Wextra -Werror -pedantic` enumeration
reproduces `dyadic/PROBE.txt`. The all-labelled scope is orders two through
six. Adding the singleton does not change its scalar sign lists: it has
only the already listed minus sign at H=1. The displayed intersection of
scalar lists is deliberately distinguished from a single graph with both
signs. This packet supports no all-order census claim.

## Order-nine literal depth hostile: independent positive checks

`verify_hostile217.py` independently decodes n=9/code26485 and enumerates all
9! old-vertex permutations. By deleting/reinserting the newborn, it counts
all one-birth child Hamilton paths without using the graft identity or
producer DP. The complete child list, by choices zero through nine, is
(217,2041,2241,653,769,1197,539,437,501,373). The same complete permutation
set checks every chronological labeling of the parent: precisely sixteen
are valid, and minimum forest depth is three. Ordinary and optimized runs
match `HOSTILE217.json`; deliberate failure controls fail in both modes.

These positive checks certify this graph and its lack of a depth-two
labeling. They do not alone rule out a *different* depth-two graph sharing
its scalar count or child-count signature. That stronger finite-fiber claim
awaits independent review of the exhaustive shallow-fiber comparison.

## Order-nine shallow-fiber transition obstruction

**PASS, explicit graph proof and exhaustive finite same-order obstruction.**
The four forced reverse deletions in `release/NORMAL_FORM.md` are correct.
They force depth zero at vertices 2 and 3 under a hypothetical depth-two
presentation. After the only permitted next deletion, both possible remaining
choices contradict that requirement. The direct proof therefore excludes all
chronological labelings, independently of the finite permutation control.

`release/verify_normal_form.cpp` exhausts every shallow history. Its depth
pruning excludes exactly births whose parent already has depth two, not a
potentially successful future reduction. The separate (number of roots,
number of depth-one vertices) recurrence correctly counts 176,850 histories.
Literal endpoint DP then finds 458 H=217 histories, 30 distinct child sets,
and 212 distinct child values in their union. The six displayed deep-child
values are absent from that entire union, not only from one representative.
The direct DP recomputes every relevant state and does not use the census's
last-birth shortcut. All integer counts are well inside uint32 at orders
nine and ten.

A fresh strict `-DNDEBUG` run matches `NORMAL_FORM_CHECK.json` exactly; its
output is saved in `normal_form.replay.json`. A private copy with the first
expected child count changed from 217 to 219 fails nonzero. Code723 was
separately decoded and verified to have depth one. The all-orderings
histogram is {3:8,4:6,5:2}, agreeing with the independent sixteen-ordering
check. The normal-form note correctly distinguishes graph nonflattening,
fixed-order continuation loss, and the still-unrefuted bare scalar image
comparison. Its discovery history is not presented as a global minimal-order
or minimal-scalar theorem.

## Selector Round 4: exact two-index anchors and bounded-bank density

**PASS, with the inherited mechanisms explicitly credited.** The quadratic
binary-tail term satisfies C(A)+(A mod2)=A+1-g(A+1). For an initial zero,
shifting the remaining bits leaves C unchanged; for an initial one, the
cross terms total B+(B mod2). Both parity cases give the claimed induction,
including A=0. Every positive M is realized by its least-significant-first
binary tail, so the formula is an exact scalar image, not a relaxed integer
model. Discarded trailing O letters do change the response, and the text
correctly confines this normal form to the total count.

Substitution gives the stated L(M,N) exactly. Its coefficient bounds imply
L>=40MN+17, including singleton N=1 and empty-tail M=1. Thus the target
inverse covers every order and tail length. The all-order family exclusions
61 and 173 are finite *parameter* inverses backed by this bound, not a
finite graph search extrapolated to arbitrary orders.

The valuation-density proof has the required uniformity: for fixed physical
seed, the additive constant is bounded independently of j, and the number
of relevant valuations is O_seed(log X) because U-d>0. Each valuation
contributes X/[2^(j+1)(U-d)+d]+O_seed(1). Counting inputs overcounts the image.
The coarse beta_N bound is (91/4136)/N; doubling for both endpoints and using
H_32768<11 gives 1001/2068=91/188. The omitted fraction among odd integers is
33/1034=3/94. All three denominators and the factor of two for odd-density
conversion are correct. The bank is fixed before X tends to infinity;
nothing in the proof permits unbounded N or a bank growing with X.

The general tail formula and the stronger parameter-capacity principle are
correctly credited to Sessions 205/206 and THM-4888, while joint residue
access and the four-block language are credited to Session199. This audit
endorses the explicit application, not a new claim of priority for those
inherited mechanisms.

The 27783 fiber enumerates exact four-coordinate states and minimal J costs,
not graph isomorphism classes. Strict increase of the balanced readout under
every letter makes the target cutoff complete. The absence of consecutive U
values refutes the named sufficient invariant, not the four-block coverage
lemma. The consecutive-J matching rebase swaps adjacent newborn pairs only,
fixes every old vertex, and transports the marked root; its restriction to
consecutive J runs is necessary and stated.

Fresh ordinary and optimized `selector/verify_round4.py` outputs match
`ROUND4.json` byte-for-byte. No proof correction was required.

## Higher interfaces: reviewed core statements

The seven-coordinate definitions in `composition/INTERFACES.md` are
consistent with the ordered two-path conventions, including empty paths.
For a root attached to a cut, h+A-T-E0 counts paths by whether they end at
the root or continue from it. The ordinal interface law correctly includes
the exceptional E0*V0 term when the two paths lie entirely in opposite
blocks. The x,y closing gadget therefore has count h+2A-T.

The common deletion-graph description and literal newborn swap are valid.
Added a->r covers decompose uniquely into a rooted source path and the two
ordered source prefixes; removing that rooted block produces the 42 gadget
skeletons. Twelve of those have the first path end at b and none end at a.
Removed r->a covers give six V contributions in total and cancel in the
endpoint difference. Thus q=4(24Q+42W-6V) and Delta G=-12W. Together with
the small-state identities this yields the exact seven-coordinate difference
(0,0,0,-12W,-2U,0,0). The full-path start sums agree by the two displayed
forced gadget orders.

For the unextended module chain, the affine cube coefficients follow because
h,q,A,E0 are common and only G,T vary. A prefix of j modules has h divisible
by 3^j and A divisible by 3^(j-1); the explicit factor six in the second
coefficient term supplies the missing factor three. This proves the stated
2*3^(k-1) divisor, including the first-module empty-prefix case. Strongness
uses the source reaching-root guard, which is stated. These are scoped
composition laws, not an arbitrary-anchor construction.

`verify_interfaces_independent.py` imports no producer code. Full
permutation enumeration of every subset independently recounts four source
cases, including a nonreaching root of a transitive edge, and both deletion
interfaces. It verifies all seven coordinates, the individual T values,
the q/Q/V/W identity, and the literal deletion-graph isomorphism. Ordinary
and optimized outputs match `INTERFACES_INDEPENDENT.json`, and both --fail
controls fail. Producer `verify_interfaces.py` also matches its JSON in both
modes and retains effective failure controls. Further producer extensions
are audited separately when finalized.

**Exploratory search issue found and reported.** The initial
`escape_divisor_probe.py` deduplicated matrices to their first word and
subsequently filtered that representative by whether it ended in Z. This
loses terminal-Z admissibility: already JZ and OJ have the same matrix, but
OJ is retained, dropping the valid terminal-Z matrix. Therefore the original
terminal-Z pair census was not complete as written. This affects the finite
exploratory search, not any displayed interface identity. The producer is
repairing the search to retain terminal-Z outputs before representative
collapse; the original 255 count is not relied upon here pending that replay.

## Final Round 4 additions and corrected search

**Fixed-N modular image: PASS.** Writing g=g(N), the tail image is
A*M+b*g(M)+c, with A=18g+44N, b=3g+6N and c=-4g-10N. Since g is odd,
`gcd(A,2b)=2*gcd(N,3g)`. The offset b+c=-g-4N reduces to -g modulo that
divisor, as does the upper endpoint shift by 2N. The inherited joint-residue
theorem supplies every compatible M,g(M) pair, proving equality with the
stated coset rather than only a divisibility necessity. If gcd(N,g)>1,
both actual endpoints share that proper divisor and exceed it, so they are
composite. The proof does not incorrectly infer compositeness merely from
the extra possible factor three in gcd(N,3g).

**Bounded-response-error compiler: PASS.** Deleting the opening O replaces
the t=>C3 gadget by C3, giving fourteen added skeletons, four ending at b,
and the same six removed V contributions. The seven-coordinate difference
is exactly (0,0,4E,-4W,E-U,0,0). The common h and q formulas, the E=1
response source U=2n-3, the n=2 singleton, strongness, gap2n and upper bound
76n-117 are all consistent. An arbitrary common binary tail has the stated
four-coordinate difference (2n+4a,0,4*2^m,0).

The protected closure's contextual increment follows directly from two
ordinal products. It is correctly distinguished from the original root
extension's gap: in empty contexts, protected closure changes by U+7E,
whereas the unprotected compiler changes by U+3E. The ordered two-flip
coefficient is +16*E_i*W_j times the stated h/q products; no higher choice
terms arise because h and q are constant, A and G are affine, and T contains
only their ordered pair products.

**No exact response cancellation in a fixed chain: PASS.** Appending a
common O/Z tail adds a transitive suffix to the deletion graph, independently
of the selected cut. Consequently h is unchanged and q is multiplied by
2 for each tail letter. With E=1, h is odd, q_base is positive and even,
and the exact A choice weights satisfy
`v2(c_i)-v2(c_(i+1))=t_i+v2(q_(i+1)^base)>=1`. In any signed difference
between two distinct binary choices, the rightmost changed coefficient has
a unique lowest two-adic valuation, so the sum is nonzero. Fixed external
contexts multiply this difference by positive factors. The stated exclusions
of general rerooting, changing the source modules, or scalar collisions are
necessary and retained.

The independent permutation verifier was extended to the six-letter pair;
its updated `INTERFACES_INDEPENDENT.json` confirms the full seven-coordinate
difference and h/q formulas on the same four source/root cases. Updated
selector and interface producer outputs again match their current JSON in
both ordinary and optimized Python modes.

**Terminal-Z search repair: accepted.** The revised code collects Z(a) for
every previous-level matrix before choosing any first-word representative.
Matrix equality preserves every possible next transfer, while terminal
admissibility is now kept separately. This repairs the lost-predicate bug.
An independent fresh replay through length ten matches the corresponding
saved rows. The producer's repaired complete length-fourteen receipt reports
792415 distinct matrices and the same 255 unit pairs, all +2U, with no
deletion-divisor escape. This remains a finite exploratory search, not an
all-word rigidity theorem.

## Shortest channel and final replay runner

**Raw four-letter channel: PASS.** Direct subtraction gives
(E+U,U,E,0), not a constant endpoint discrepancy. With E=1,U=2N-1, the
upper count is 5g(U)+7U-1<=24N-13, the response discrepancy is one, and the
endpoint discrepancy is 2N-1. Singleton input gives 9/11. Strongness follows
from the reaching-root J/Z argument with the explicit singleton boundary.
The singleton cuts make T=0 and force the seven-coordinate difference
(0,0,E,-4W,0,0,-U). Appending O creates a sink and preserves gap2N; a further
Z restores strongness but changes the gap to 2N+2. All these distinctions
are retained in the final text.

The independent permutation interface receipt now covers this raw channel
as well. Final ordinary and optimized producer/interface-audit replays match
their updated receipts, with effective deliberate failure controls.

The release's strengthened CLI parser checks the parsed unsigned-long value
before narrowing; malformed/wrapping inputs are rejected. The exact
running-to-published source difference changes only this validation.
A fresh `verify_controls.py` replay matches all eighteen stored controls,
including reconstruction of the running source from the reverse patch,
verification of its hash, and byte equality of its small-order outputs with
the published wrapper. Saved audit receipt: `release_controls.replay.json`.

The aggregate `verify_all.py` correctly distinguishes explicit development
PARTIAL mode, stored-certificate recounting, and optional full order-twelve
regeneration. Reading its code found no additional verification gap. The
README, SYNTHESIS and expanded THM-5115 preserve the open all-target boundary;
relative density is now stated relative to odd integers. Final completed
census numbers are audited below once available.

## Final completed census, normalization theorem, and sealed verdict

**Final verdict: PASS within all stated scopes. No unresolved mathematical
or implementation defect remains in the reviewed final packet.** Chan–Pak,
the unbounded two-index coverage problem, and the inherited four-block
coverage lemma remain OPEN. Earlier references in this report to a running
order-twelve computation describe intermediate audit stages and are
superseded by this completed section.

The main release enumerated all 12!=479001600 chronological histories and
completed its factorial-count checks. It has 53549 distinct order-twelve
values, maximum115781 and first allowed missing value91081. Minimum-depth
counts are exactly {0:1,1:46092,2:7455,3:1}; the first requiring depth two is
53443. The only depth-three scalar is **97309**. Depth-two scalar equality
holds at each order through eleven and fails at order twelve. These are
class-specific fixed-order statements, not a new global coverage range.

A fresh independently compiled strict `-DNDEBUG` direct endpoint-DP auditor
recounted all **138488** stored literal graphs across the twelve CSV files,
and its row receipts match the release exactly. Results are saved in
`final_all_witnesses.replay.jsonl` and `final_all_witnesses.replay.log`.
The separate Python release auditor and every summary field are consistent
with the final CSVs. All34 files in the frozen release manifest verify;
the manifest SHA256 is
`f8c886844b838c0970388552f14488f44711730e81899466b3b1aa3dd4b72ac9`.
The current **22-control** receipt was freshly replayed and matches exactly,
including normal/optimized controls and original-running-source recovery.

**All-order graph-preserving normalization: PASS.** Swapping consecutive
vertices x,y with y->x changes only their mutual relative order. The
one-parent forest loses precisely that reversal edge; y becomes a root,
x keeps its earlier parent, and every later parent is still earlier. No
vertex depth increases. At most n-1 swaps therefore give a Hamilton-path
birth order with no immediate parent. Each birth with j predecessors then
has j choices, yielding (n-1)! histories. Starting with a minimum-depth
presentation proves that restricting to these histories preserves every
minimum. Root responses are preserved under the vertex identification;
there is no reset of a distinguished root to a new position.

The general bounded-outdegree corollary is also valid: y loses one earlier
outneighbor, x gains none, and other earlier-neighbor sets stay unchanged.
Repeated swaps preserve the bound b and yield the stated product of
binomial sums. The forest depth conclusion is not asserted for the
multiple-parent case. No external priority claim is made.

`verify_normalization.py` independently tests every chronological history
through order seven, recomputing each vertex's depth and the reversal count
after every swap. Ordinary/optimized outputs match `NORMALIZATION.json`;
the intentional failure controls fail. This finite test corroborates the
separately reviewed universal proof.

**Independent shallow exclusion: PASS.** `release/shallow_census.cpp` never
visits a depth-three birth and stores no minimum-depth array. Normalized
parent exclusion is exactly the immediate-predecessor condition. Its
(roots, depth-one vertices, last depth) recurrence was independently
recomputed in Python, giving2490378 order-eleven parents and23670073
order-twelve histories. The completed run produces53548 values and compares
the entire image against the main CSV, excluding exactly97309. The
normalization theorem upgrades this restricted search to an exclusion from
all123489888 unnormalized shallow histories. Fresh normalized and
unnormalized small runs both reproduce all192 order-eight values, and a
private false-depth CSV is rejected. Receipt: `shallow8.replay.jsonl`.

**Literal exceptional graph and variable-order control: PASS.**
`verify_exception97309.py` independently reconstructs code25112126 and
recounts H in both start/end directions. Its forced reverse deletions are
11->1,10->3,9->2,8->6,7->0,6->4,5->0,4->0. Every presentation therefore
contains8->6->4->0, directly forcing depth at least three; the supplied
history attains it. This graph-level proof is distinct from excluding every
other same-count shallow graph. Ordinary and optimized receipts match
`EXCEPTION97309.json`, with effective failure controls.

The independently recounted depth-one factors (n,code,H)=(6,50,31),
(7,296,43),(7,58,73) have product97309. Their twenty-vertex ordinal join
is a depth-one realization, by the exact join concatenation bijection.
Thus the finite refutation is specifically **fixed-order** scalar equality;
97309 is not a variable-order shallow obstruction or a missing tournament
count. The completed NORMAL_FORM, release README, THM-5116, packet README,
SYNTHESIS and S244 preserve that qualification.

The expanded THM-5115 also matches the audited formulas and scopes, including
the shorter response-error-one channel, exact higher interfaces, finite
anchor inverse and fixed-bank density bound. All producer and audit scripts and receipts are frozen for final aggregate verification and publication.
