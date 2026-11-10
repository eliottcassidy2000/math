# Session244 independent adversarial audit

Status: **PASS within stated scopes**. Chan–Pak remains **OPEN**. No all-target existence conclusion follows from these packets. Audited 2026-11-10; producer files were not edited.

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

The parent supplied the following argument during the audit. It is valid and supplies an independent proof of the complete nontrivial strong saturated-root class. The classification and a stronger deficit bound are inherited, as documented in the provenance correction below.

**Pair-cover lemma, PROVED.** Given a tournament A and distinct a,b, there is a spanning ordered two-path cover in which one path starts at one mark and the other ends at the other mark. Let C1⇒...⇒Ck be the strong components of A−a, with b∈Ci. Hamilton paths in a strong tournament can start or end at any prescribed vertex (cut a spanning directed cycle; singletons are immediate).

- If i=k, use the singleton (a) and a Hamilton path of A−a ending b.
- If i<k and a beats some c∈C_(i+1), use the path starting a, then a Hamilton path of C_(i+1) starting c, then the remaining later components. Its companion traverses C1 through Ci and ends b.
- Otherwise C_(i+1) beats a. Use a path starting b in Ci followed by C_(i+2) through Ck. Its companion traverses C1 through C_(i−1), then C_(i+1), then a.

All concatenations follow the component order, and empty prefixes/suffixes cause no exception. The use of the standard strong-tournament Hamilton-cycle theorem should be named in the final theorem dependency list.

Cutting a Hamilton path just before r identifies D with exactly those rooted covers whose other path is empty or ends at an inneighbor of r. Thus U−D counts the remaining covers. Write O=N+(r). The singleton-root covers contribute sum_{a∈O} e_a(T−r). Each unordered pair {a,b}⊂O supplies a further cover by the lemma, after prepending r to its marked-start path. Different unordered pairs give different covers, since the first vertex after r and the other path's endpoint recover that pair. Therefore

    U−D ≥ binom(|O|,2) + sum_{a∈O} e_a(T−r).

If |O|=1, every nonsingleton rooted path consumes the sole outneighbor, leaving no outneighbor available as the other path's endpoint. Hence equality simplifies to U−D=e_a(T−r).

For a **nontrivial strong** tournament, outdegree is positive, so U=D forces outdegree exactly one. Then equality holds exactly when a is not in the last strong component of T−r. Strongness forces a into its first component, yielding precisely the A⇒B decomposition with B nonempty and the structural formulas above.

**Required boundary correction:** the singleton tournament is strong under the usual convention and has U=D=E=d=1, outdegree zero. Accordingly the necessity/classification and h≥3e must explicitly exclude it. This was reported to the parent before promotion. No other proof defect was found.

### Parity strengthening, also audited

The parent subsequently supplied a valid strengthening. Let K_ab count the covers starting a and ending b in separate paths. It counts the Hamilton paths using adjacent arc b→a after those paths are joined. K_ab and K_ba are unaffected by flipping arc ab because the marks lie on different paths. Consequently flipping that arc changes the full path count by K_ba−K_ab. Rédei oddness makes this difference even, so K_ab+K_ba is even. The pair-cover lemma makes it positive, hence at least two. Therefore the stronger universal bound is

    U−D ≥ |O|(|O|−1) + sum_{a∈O} e_a(T−r).

The independent permutation verifier checks this stronger bound, not just the weaker binomial bound, on every root of every labelled tournament of orders 2 through 5. The unbounded proof uses the inherited oddness theorem, not extrapolation.

## Consequence and stop boundary

The audited work strengthens the explanation of why even-gap existence is insufficient and gives genuine structural information about a self-feeding physical face. It does not establish a target-solving selector or stable inverse family for every allowed odd integer. Chan–Pak, the required all-target anchor, and any broader physical completeness conclusion remain OPEN.


## Round 2 audit continuation

**Correction to this audit:** The first pass incorrectly endorsed the composition inverse equations while verifying the forward maps. The logic worker correctly found a sign reversal. The correct F− inverse is h=(h′−3e′)/2 and e=3h′−10e′. Substitution gives (9,3) from (171,51). The composition proof has now been repaired. Forward maps and all earlier forward numerical results are unaffected; do not rely on the superseded first-pass inverse endorsement.

**Complete saturation STRUCTURE, PASS:** The alternative reachability proof is sound: two outneighbors whose deletions each destroy root reachability lead to disjoint nonempty unreachable sets, each forced to beat the other. The root-reachability/Hamilton-start equivalence correctly uses ordered strong components and Hamilton cycles. The arbitrary-SCC formula U_r(P⇒C⇒B)=H(P)U_r(C)Q(B) is exact: vertices of P must all belong to the complementary path, the C portions supply the rooted cover, and B splits freely between the two paths. The independent permutation verifier now checks this formula for every root of every labelled tournament of orders2–5. Terminal-singleton saturation is explicitly retained. The prime-target gcd restrictions are necessary, not sufficient, as stated.

**Logic ROUND2, PASS:** The difference cone is preserved by every letter; the first J creates positive deletion difference which no later letter can erase. With J excluded, U doubles each step and the full-count coefficient catches up exactly only for all-Z words. The asymmetric universal face-return argument correctly relies on independence of h,e,g and is not extended to target-specific fibres. Same-branch inverse synchronized gaps require divisibility4 for F− and28 for F+. Cross-branch cancellation yields exactly (6(7E−D−U),12E) in the face coordinates, so integral common continuations cannot produce a total gap ±2. These are restrictions on the named alphabet/continuations, not on arbitrary tournaments.

The one-prime countermodel S_p indeed contains the earlier model and all the stated composite witnesses. For an empty finite set F, the cutoff can be written max({151}∪F), avoiding undefined max(F); this is a harmless notation convention rather than a mathematical obstruction.

**Face verifier, PASS:** The C++ subset recurrences for start/end counts, pair-cover coefficients, exact defect decomposition, strengthened bound, and strongness are correctly indexed, including singleton and empty-complement cases. The -DNDEBUG build with C++17, -O2, -Wall, -Wextra, -Werror exactly reproduces `face/CHECK.json`. Its --fail control exits nonzero (134). This independently confirms the normal producer replay is not dependent on assertions. The census is finite-exact through order6; the proof carries the global statement.

**Selector second-round computation, PASS:** `verify_round2.py` replays exactly against ROUND2.json. Its O/Z/J state search legitimately deduplicates full profiles and has no imposed word-depth cutoff. J strictly increases U by D, while the earlier final-Z cutoff applies also to final J because its full count is the same D+U−E. Seed counts and final letter counts are included, so trailing O words lose no count. The new literal eight-vertex H=173 witness is independently reconstructed from the edge, word ZOOZ, reroot2, suffix ZZ, and then checked by all directed subset permutations: profile (173,9,209,9); all eight roots reach the entire graph. This certificate is stronger than the producer's recurrence-only construction and is included in INDEPENDENT.json.

The freely rerooted potential argument is valid with an explicit reaching-root guard and seeds of at least two vertices. Adding a sink preserves exactly the old reaching-root set and doubles every corresponding U, hence doubles μ. After a Z at a reaching root r, D′=D+U_r−E_r≥U_r≥μ(T); every output rooted response is at least D′, so μ cannot decrease even if new reaching roots appear. Every Z increases D by a positive even integer, hence at least2. Before the final Z, μ≤U_r≤h, bounding the number of preceding O moves by log2(h/μ_seed); all intermediate counts are nondecreasing, bounding the number of Z moves by (h−D_seed)/2. Trailing sinks can be removed, leaving a finite graph search. This does **not** prove the same claim for unguarded arbitrary roots: Z at a sink can preserve D, invalidating that length bound. A separate reduction would be required for unguarded roots.


## Provenance correction: proof validity is not novelty

A parent-led lineage audit found prior results that the first proof audit did not adequately credit. This does not invalidate the displayed mathematics, but changes its classification from a proposed new structural advance to **INHERITED with independent alternative proofs and exact checks**.

- `05-knowledge/results/session231_overnight_completion_20261109/descent/ROUND2.md`, §2, already proves the complete nonsink zero-defect face, h≥3e, its equality condition, and the e=1 dyadic family. This was verified directly from the tracked HEAD blob because that path is absent in the sparse working tree. Its statement is broader than the strong-only version.
- THM-4875, `01-canon/theorems/THM-4875-old-mark-digit-lower-bound-lorentzian-gadget-adjoint-and-one-apex-joins-realise-every-odd-number-to-2-32-except-7-and-21.md`, (A1), explicitly credits THM-4873 and gives the stronger bound σ_r≥E_r+2Σ H(T−{r,a,b}), summed over outneighbor pairs. Since each residual tournament has at least one Hamilton path, the Session244 degree bound is a consequence. Here inherited E_r=H(T−r)−e_r(T) equals the present e_O(T−r), by appending r to a deletion path exactly when its endpoint is an inneighbor. The prior weighted bound must take precedence over the weaker rederivation. The parent also traced the bound to `two_entry_funnels_20261103/READY.md`, equation(7).
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

THM-5113, THM-5114 and THM-5115 were read in full against the audited sources. Their mathematical statements and all-order scopes pass. THM-5113 states endpoint realization, not arbitrary full-count realization. THM-5114 explicitly marks the density bound inherited and the countermodels abstract. THM-5115 correctly retains synchronization/continuation restrictions and the corrected inverse sign, while crediting inherited saturation structure. One wording improvement was sent to the parent: say “the gcd of the two branch total counts” rather than “its two output total counts,” since Fminus by itself supplies only one total-count coordinate. No substantive theorem correction is required. These three summaries do not claim Chan–Pak completeness.
