# Three block sizes certify a uniform dyadic count law

**PROVED IN NOTE:** M1-M4. **FINITE-EXACT:** the declared search and controls in M5. M6 additionally uses the incoming complete rooted ready fibre at881. Chan-Pak remains OPEN. This packet makes no literature-priority claim.

## Inheritance and scope

Use the inherited ordered nonempty path-cover expansion from [the previous module compiler](../../cycle_opening_six_return_20261112/counts/PROOF.md#3-exact-run-coefficients-and-arbitrary-modules). Fix a tournament context Q with one distinguished vertex v, and replace v by any nonempty tournament M, all comparisons to outside vertices uniform. Let m=|Q|-1 be the number of outside vertices. Let F(M) count either all Hamiltonian paths, or paths with one or both endpoints fixed at specified outside vertices. If both endpoints are specified, they are distinct. The module does not contain a marked endpoint. Write

    F(M)=sum_(j>=1) A_j K_j(M),

where K_j(M) counts ordered j-tuples of nonempty directed paths partitioning M. The nonnegative integer A_j counts skeletons of j module visits and all m outside vertices. This exact expansion is inherited; the no-hidden-higher-layer certificate below is its new use.

Anchor: full-H ready returns. Niche: finite certificates for arbitrary modules. Wildcard: a finite activation threshold which really controls all later sizes. Closest mechanism is the old run compiler; hostile is a two-point fit hiding a third run; corrected near miss is extrapolating a fitted exponential without controlling higher layers; sidecar is the positive run support.

## M1. Positive support descends, and the context imposes a sharp cutoff

**Theorem.** If A_j>0 and j>=3, then A_(j-1)>0. If s is the number of fixed outside endpoints (0,1,2), then

    A_j=0 whenever 2(j-1)+s>m.                    (1)

Also A_j is divisible by (j-1)!. Thus there are no gaps in the positive support above the second coefficient. For unrestricted H, A_1=H(Q)>0; no analogous A_1 assertion is needed for doubly prescribed endpoints.

**Proof.** Put O={outside vertices beaten by v} and I={outside vertices beating v}. A skeleton with j visits has j-1 disjoint nonempty internal outside paths, each starting in O and ending in I. Since O and I are disjoint, each uses at least two vertices. A specified initial endpoint consumes at least one vertex in the prefix before the first module visit; a specified terminal endpoint consumes at least one in the suffix. This proves (1).

Any two disjoint directed paths can be merged into a directed path preserving the order within each. At each step, compare their current heads and take the one which beats the other head. Its next chosen successor is either its own next vertex, along an existing path arc, or the other head, which it beats. Continue until one list is exhausted. The merged path starts at one old start and ends at one old end. In particular, merging two internal O-to-I paths preserves their boundary types.

For j>=3 merge two internal paths, keep all other internal paths, prefix and suffix, and separate the resulting j-2 internal paths by j-1 dummy module visits. The result is a valid skeleton counted by A_(j-1), with the fixed endpoints unchanged. The two selected paths need not have been adjacent: their position is free because every separator is the same uniform module. Finally, freely permuting the j-1 disjoint nonempty internal paths gives a free symmetric-group action on the skeleton set. This proves the factorial divisibility. QED.

The s=0 cutoff is sharp for every j: take Q[M]=C3[M,TT_(j-1),TT_(j-1)] in cyclic block order. The top skeletons have the form M,b,c,M,...,b,c,M. Each of the j-1 vertices in each outside block can be ordered independently, so A_j=((j-1)!)^2>0, while higher coefficients vanish. Prefix a dominating singleton to force a fixed initial endpoint; append a dominated singleton to force a fixed terminal endpoint. These show sharpness also for s=1,2.

## M2. Three sizes are a complete certificate, not a fit

Write F(k)=F(TT_k), k>=1. Since K_j(TT_k)=j!S(k,j),

    F(1)=A_1,
    F(2)=A_1+2A_2,
    F(3)=A_1+6A_2+6A_3.

Therefore

    E=F(3)-3F(2)+2F(1)=6A_3>=0.                 (2)

In fact E is divisible by12, and E>0 implies E>=12, by M1. The following are equivalent:

* E=0;
* A_j=0 for every j>=3;
* for every positive integer k,
  F(k)=A*2^k+B, where A=(F(2)-F(1))/2 and B=2F(1)-F(2);
* for every nonempty tournament M,
  F(M)=A*P_M(2)+B*H(M), where P_M(2) counts ordered two-path covers allowing empty bins.

Indeed M1 makes A_3=0 force every later coefficient to vanish. The identity K_2=P_M(2)-2H(M) gives the arbitrary-module conclusion. Conversely evaluating the dyadic law at1,2,3 gives E=0. A is a nonnegative integer; B can be negative. Positivity is claimed in the K basis, not the changed (H,P2) basis.

For a prescribed target A*2^k+B, checking F at1,2,3 therefore proves an all-module identity whenever the context and endpoints have the declared uniform form. A proposed difference of two arbitrary observations is not itself a positive observation: one must certify each side separately or retain its complete coefficient vector. Selective attachment and moving internal ports are not covered.

More generally, for r>=2, inversion on k=1,...,r+1 recovers A_1,...,A_(r+1). If the last is zero, M1 proves all later ones zero. This is a finite certificate of exponential degree at most r. Nonnegative formal coefficient vectors with a hole do not inherit this implication; the path-merging theorem is essential.

## M3. If the third test fails, the departure has a known scale

Let D(k)=A*2^k+B be the law determined by the first two samples. For every k>=1,

    F(k)-D(k) >= (E/6)*(3^k-3*2^k+3).             (3)

This follows by retaining the nonnegative A_3 K_3 term. At k=1,2 the right side is zero. At k>=3 it is positive exactly when E>0. More precisely, if r is the largest index with A_r>0, then F(k)~A_r r^k. Thus the finite skeleton capacity identifies the leading exponential scale.

A cheap explicit separation bound is available: for k>=5, 3^k-3*2^k+3 >= 3^k/2. Consequently a required absolute accuracy B0>0 fails as soon as k>=5 and (E/12)3^k>B0. This concerns a fixed context and a specified observable; it is not a phase-transition theorem for arbitrary infinite sets.

Two samples cannot detect the change. The actual five-vertex tournament with outgoing rows (10,28,25,16,1), inflated at vertex0, has

    (A_1,A_2,A_3)=(9,10,2),
    F(1),F(2),F(3)=(9,29,81).

Its first-two-sample law is10*2^k-11, which predicts69 at k=3, missing12. Its exact all-k count is2*3^k+4*2^k-5. Four outside vertices are minimal for such an unrestricted hostile by M1. The literal permutation replay at k=3 is independent of the run compiler.

## M4. A ready-return certificate must include its ports

For a strong one-start parent (U,p), append a new root z beaten only by p. Its defect is the old starting-path count s_p(U), hence one. To inflate a vertex v of U, keep v distinct from p and z. Check both the child's full count and the parent's p-start count by M2. If the former has the prescribed dyadic law and the latter equals1 at k=1,2,3, then every transitive inflation is an actual ready source and its full count is correct for every k. Strongness of the fixed quotient persists under nonempty transitive substitution. For a general M the starting count becomes H(M), so readiness one is asserted only for the transitive family, not every module.

This separates three operations: the genuine theorem for all module counts; the special transitive family where H(M)=1; and a physical same-target replacement which still must be exhibited between source and returned graphs. The certificate alone is not such a replacement.

## M5. Exact bounded search, with positive and hostile controls

`module_certificates.py` enumerates all strong unique-start rooted classes through parent order17 using the inherited complete funnel grammar. It appends one indegree-one root, retains targets97,293,881,2645,7937, and tests every vertex except the new root and its sole predecessor as a uniform TT_k module. The target rays have (A,B)=(45,7),(135,23),(405,71),(1215,215),(3645,647).

There are12 candidate rooted sources and100 eligible module choices. Exactly three choices match the requested rays: two for97 and one for293, recovering the inherited positive controls. None matches the next three rays in this bounded universe. This search allows the variable module anywhere in the nested core; it is broader than inserting a run only in the final funnel tail. No relabelled duplicates are added to the unique-start rooted class grammar.

Independent root controls cover every Hamilton-frame presentation through order5, every module vertex, and every external initial endpoint:1,751 observations, each checked at three sizes by direct subset path counting, plus the hostile by literal permutations. The separate [independent audit](../curves/MODULE_AUDIT.md) uses all1,099 labelled contexts with at most four outside vertices and3,293 observers, including fixed pairs of external endpoints; it checks nontransitive modules and the five-core directly. Both engines use explicit exceptions rather than removable assertions and pass normally and under -O.

## M6. Incoming complete881 fibre turns one bounded miss into a global obstruction

The incoming [Session296 ready fibre theorem, Section4](../../session296_multicore_metagraph_20261113/cores/PROOF.md#4-complete-first-target-ready-fibre-not-a-higher-clock-exclusion) proves that the strong rooted defect-one fibre at H=881 consists of exactly two sources, of orders10 and11. Its proof uses the inherited complete ready inverse, strong-order bound and complete unique-start grammar. We inherit that classification, not its discovery.

**Corollary.** There is no fixed tournament context Q with one uniformly substituted TT_k module and a fixed external indegree-one root of defect one for which

    H(Q[TT_k])=405*2^k+71

for all k>=1. The same exclusion holds if the identity and readiness are assumed only for all sufficiently large k in this same direct TT_k parameter.

For the all-k statement, put k=1. The source must be one of the two classified881 graphs. The sole predecessor cannot be the module vertex, since it would become k incoming vertices. Our complete list of the eligible second run coefficients is

* words (11,1,11):881,1529,1038,678,678,657,462,462;
* words (1,empty,00101):770,913,824,617,430,430,501,339,673.

None equals the required405. These are exact skeleton coefficients, with an independently counted total881. A valid dyadic chart must have A_2=405 by M2.

For an eventual identity, finite exponential independence extends it to all positive k: the positive top run coefficient forbids any base above2, and bases1,2 determine the coefficients. The root's indegree-one condition for large k forces the module outside the root and predecessor. Its defect is an external-start count in the root-deleted context. Positivity forces that count to have only K_1 coefficient1 when it equals1 eventually, so readiness also extends to k=1. Strongness of the uniform quotient is independent of positive module size.

Weak ready sources are the classified strong ready source followed by a transitive tail, by the same inherited inverse theorem. A module in the tail changes no full count, and a module in the strong component gives the same listed coefficients. Thus padding does not evade the obstruction.

**Boundary.** This does not exclude two changing modules, selective incidences, a changing quotient/core, a piecewise or reindexed tail parameter with no TT_1 specialization at881, or any non-module construction. It explains why the available finite881 entrances do not simply stretch into the desired next ray. The finite entrances themselves remain valid, including the long ready ladder to71441.

## Reproduction and connection contract

From the repository root run `python -B 05-knowledge/results/mrb_seven_cutoffs_20261113/module/module_certificates.py` and repeat with `-O -B`. The two receipts contain identical mathematical data. The universe, inherited grammar and endpoints are explicit above. The proof of M6 depends additionally on Session296's complete881 classification; the finite search by itself is not an all-order exclusion.

Source -> target: a fixed uniform tournament context -> its positive run coefficients -> three transitive sizes. Preserved predicate: full Hamilton count and specified exterior endpoints. Lost data: the internal module orientation if only sizes are recorded. Sidecar: K_j(M), with M2 restoring the full arbitrary-module formula when E=0. Cheapest hostile: the order-five example in M3. Cheapest positive control: the existing45/135 ready rays. The next genuinely open construction needs a parameterization outside M6, not another one-module inflation of881.


## M7. Two escape menus were tried, with exact stopping reasons

The negative result above suggests breaking uniformity or starting a new block from size zero. `escape_probes.py` tests both on the complete two-source881 fibre.

1. Clone each of the17 eligible vertices once, allowing the old predecessor to see the two clones in any of four directions. This is68 comparisons,36 of which remain ready. None has H1691. The closest ready example has H1681, but no inference from its distance10 is made.
2. Adjoin one completely arbitrary new vertex while preserving the old root and every old arc. The root must beat the newcomer to retain its unique predecessor. The remaining incident arcs give512+1024=1536 possibilities. Exactly38 are ready; none has H1691. Thus a new uniform block beginning empty also cannot make the first881-to1691 step by a single induced vertex extension from either strong source. Weak padded sources, changing old arcs, and adding several vertices are outside this finite menu.
3. Inflate the predecessor in the parent and choose the first or last clone as the new sole predecessor. The first source at k2 has (H,defect)=(2747,2) or(3555,1); the second has(7553,36) or(8217,35). Moving the mark to the last clone preserves readiness in the first example and fails in the second. Hence even this small physical repair needs the parent's endpoint response, not just H881.

The unified receipt records all68 selective outputs and all14 predecessor-clone samples, and exhaustively counts the1536 extension orientations. These are FINITE-EXACT menu exclusions, not all-construction no-go theorems. They justify moving next to old-arc repairs or genuinely changing contexts rather than extending this same one-vertex menu again.
