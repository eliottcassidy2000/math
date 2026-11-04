# Independent audit of the one-child modulo-four anchor

## Verdict

The mathematical argument at immutable commit `97632252` is valid for quotient order **m >= 2** in Theorem B, subject to the cited reviewed original-module reduction inputs. Theorem A already states its required `m >= 3` hypothesis. The all-order contraction, weighted E1, common-end anchor, and fixed-ambient modulo-four lift pass the independent proof audit and controls below.

There is one literal theorem-statement error to correct before importing the result: **Theorem B omits m >= 2**. Its proof invokes n=m+1>=3. For m=1 and unchanged child, the source and target are the identical named tournament and every query contrast is zero. This is a genuine boundary counterexample, not a failure on the intended prime-node domain. The prime-node corollary remains valid after the correction.

The frozen source also overstates the availability of its verification artifacts. At `97632252`, `PREREGISTRATION.md` ends with an empty Outcomes section; `audit/REVIEW.md` returns404 and is absent from the directory listing. The README lists confirmatory output files that are not in that directory. Those source status lines are not treated as completed evidence here. The controls in this audit were independently generated, with no incoming implementation imported.

## Exact source and physical conventions

Reviewed source:

- https://github.com/eliottcassidy2000/math/blob/97632252/05-knowledge/results/reversal_cube_mod4_anchor_20261104/PROOF.md
- Proof blob `51390c8fe6423e7bd5a7326495ae81e01121c3aa`.
- Preregistration blob `859f738c571829202cf6c3155c3f81e01eb76180`.
- Fixed-ambient input: https://github.com/eliottcassidy2000/math/blob/97632252/05-knowledge/results/session188_structural_bridges_20261104/parity_profile_transfer/PROOF.md

The two ports have the fixed physical edge x->y, and a_v=[v->x], c_v=[v->y]. Child sets and vertex names are the same in the compared tournaments. The quotient is changed to its exact whole converse; it is not independently relabeled. Source and target child interiors may otherwise be unrelated, except for the common-end hypothesis on the distinguished child. Equal complete cover profiles are **not** required for the local modulo-four theorem or its fixed-ambient lift.

The core must remain a common module in one fixed exterior for the ambient statement. At every exterior vertex, both ports must copy the original core's incidence. Arbitrary fixed port bits outside the core are not covered by this lifting statement.

## Proof audit

### 1. Odd-cycle contraction

The classical expansion counts nontrivial directed odd cycles, with no loops. Modulo four, products involving two nontrivial cycles vanish. Port-free cycles produce only a query-independent half-contrast constant.

If a port cycle uses a child in r>=2 maximal segments, its distinguished port fixes the traversal origin. Permuting the nonempty disjoint segment contents among the r existing slots preserves all external boundary arcs and the total vertex set. The action is free even when segment lengths coincide, because the physical supports are disjoint. Thus r! kills such cycles modulo two.

For one segment, sum Hamilton paths over nonempty subsets. Each subset contributes an odd number by Redei. A child of size at least two has odd parity of even nonempty subsets and even parity of odd subsets, independent of its internal arcs and independent of its order parity. Therefore its contracted weight is0, not its size modulo2. A singleton contracts with weight1. This justifies arbitrary target interiors and also explains why a query-independent constant must be retained when profiles differ.

### 2. Weighted E1

The coefficient formulas arise from the four possible port placements in a simple cycle. Boolean multilinearity is essential: repeated variable occurrences are reduced idempotently, not set to zero.

The endpoint lemma supplies the required proper-subset cancellation. In particular, for x!=v,

    sum_{S containing v,x; S proper} eps_v(S) = eps_v(V minus x) mod2,

by applying the proper-subset identity to both V and V minus x. In the weighted correction, the term with x=v vanishes by the original identity. The remaining terms are exactly paired with the singleton-complement correction.

For n>=3, the nonempty-subset parity g(v) is w(v), including the case in which all other weights are zero. This is where the order hypothesis matters. An explicit n=2 boundary is the edge0->1 with weights(1,1): both functionals are0, not1. Our finite controls independently confirm this boundary.

### 3. Unchanged-child swap argument

Only cycles containing the reversed u-beta arc can differ between the two local arc orientations. Such a cycle contains both u and beta. Swapping these two names preserves all its other incidences, because they have identical exterior neighbors. It preserves its used vertex set, hence its total weight, even when w(u)!=w(beta). Thus the difference is invariant on coefficient orbit sums under swapping their query variables. The unequal weights do not invalidate the swap argument.

When another contracted child has weight1, its L_z monomial family is swap-invariant and detects both target types. Otherwise all external weights are0. For k>=3, L_u+L_beta plus the missing cubic swap partner equals1+p_j(Q). If that vanishes, p_j=1, and the coarse linear coefficient detects. For k=2, a changed-arc cycle has even old weight and can be odd only with one port, so it changes only pure-a or pure-c coefficients. Every mixed L_u coefficient survives unchanged. Each case is valid, including quotient order m=2.

### 4. Query count

The quoted conservative9(m+1) bound is safe. Put n=m+1. A single L_v uses at most4n actual Boolean queries. In the hard unchanged case, the union supporting L_u, L_beta and the extra cubic has12 internal queries, plus7 for each common outside out-neighbor and3 for each common outside in-neighbor, at most7n-2. The coarse fallback uses two queries already in that union.

The source Battery sentence should mention the singleton branch's L_z as well as L_u,L_beta and the extra cubic. With an outside singleton, L_z alone handles both unchanged and reversed targets within4n queries. This is a documentation omission, not a failed numerical bound.

### 5. Fixed-ambient lift

For every tournament of the same augmented order N>=2, q2=2 mod4, q3=2[N odd] mod4, and qh=0 mod4 for h>=4. Thus the higher cover-coordinate differences of the two augmented cores vanish modulo four without any equal-profile assumption. Exact ambient run substitution leaves only c1 times the Hamilton difference, where c1 is an odd Hamilton count of the contracted exterior tournament. A local residue2 therefore remains exactly residue2 modulo four. No generic integer-gap or modulo-eight/sixteen transfer is inferred.

## Strict scope and effect on Session189

The incoming anchor removes the surviving **unchanged/conversed transitive-child ambiguity** left by Session189's internal shared-name quartics. It works at every child size k>=2 and every quotient order m>=2, independent of quotient primeness, endpoint-parity p=0 or not, other child sizes, other child interiors, or any ordered exterior witness.

This substantially strengthens the earlier Session189 mixed-quartic root-family separation consequence. The older consequence required an odd transitive child and an ordered exterior witness with equal log-concave profiles; the new modulo-four anchor needs none of those extra hypotheses. The older exact mixed identity remains valid and useful: it is an integer local-converse identity with explicit exterior cover skew, not merely a modulo-four existence theorem. The internal quartic distance theorem is also still essential for excluding arbitrary reorderings when k>=5.

For an original-module XOR full-response collision, hereditary equality first keeps a source transitive child transitive. The highest-node reduction gives a prime quotient changed to its converse and permits equal-profile restoration at disjoint changed nodes. Apply the internal-distance theorem for k>=5, then the corrected anchor. For k=2 or3, source and target transitive orders share at least one endpoint; for k=4, avoiding all source endpoints forces the target's endpoint set to be exactly the source's two middle names. Therefore the stated residual frontier follows:

- every highest changed prime quotient has order at least7 and p(Q)=0;
- no singleton children;
- no transitive child of order2,3 or at least5;
- a transitive child of order4 can survive only in a crossed endpoint reordering.

This remains a **necessary condition for a putative original-module full-response collision**. It proves neither unrestricted injectivity nor existence of any collision. It does not settle general nontransitive children or the crossed TT4 targets. The quoted residual C3 exploratory claim was not needed for this audit.

For a source-relative family battery rather than an assumed collision, the scalar query must precede the hereditary-equality reduction: an arbitrary original-module XOR output of a transitive child need not itself be transitive. If the scalar query is equal, the reviewed hereditary theorem supplies the required transitivity and profiles. A TT3 battery should include both source-end anchor choices because the common end can depend on the target.

## Valid stronger common-end lemma

Transitivity is stronger than the actual anchor argument needs. Let m>=2. Suppose the distinguished child has at least two names and one fixed physical vertex u is a source or sink in **both** source and target child tournaments. The remainder beta is then a common module; both interiors of beta may be arbitrary and unrelated. Contraction produces exactly one of the same unchanged/reversed kernels, with w(u)=1 and w(beta)=[|beta|=1]. The entire argument and ambient lift therefore go through.

This extension is a direct reading of the proof, not a claim that arbitrary original-module outputs preserve such a common end. The physical controls include a nontransitive source child consisting of a dominating vertex above C3, compared with TT4 on the same names sharing that dominating vertex.

## Independent controls

`audit.py` independently enumerates all simple port cycles by literal directed-cycle DFS. The two port edges are symbolic Boolean factors; it collects the exact parity polynomial as an affine function of all old vertex weights. It then tests the stated coefficient functionals. It does not import `wmodel.py`, `confirm.py`, or `fcoef.cpp`.

The physical test layer constructs complete named tournaments and computes actual Hamilton counts using a separate subset DP, then compares every Boolean group query against the contracted model including the constant half-contrast. The Q8 controls use a new compiled subset DP; all augmented orders are at most20, so20!<2^64 protects exact integer storage.

Completed normal-mode controls:

- Every labeled weighted tournament of orders3,4,5 and every vertex:168,128 weighted-E1 functional checks.
-80 additional seeded random weighted tournaments of orders6,7.
- Every labeled quotient of orders2..5, every distinguished child, all outside weight vectors, both beta weights, and unchanged/reversed target types:336,288 kernel comparisons. Checks include nonzero contrast, reversed functional, singleton functional, hard-case functional, and swap-orbit identity.
-1,344 complete physical query comparisons across six cases, with arbitrary nonmatching interior profiles, unchanged/reversed TT3, a different TT4 order with a shared endpoint, TT2, singleton-child cases, and the stronger nontransitive common-end case.
- Four literal cyclic fixed-ambient lift controls.
- The m=1 unchanged-child counterexample, with all16 physical port queries zero.
- Every child position of the source prime-eight obstruction, child sizes2..5, both target types:64 kernel comparisons, complete coarse polynomial blindness, and six literal local/ambient witnesses.

The complete physical masks, source/target rows, counts, and kernel records are retained in `controls.normal.json` and `q8_controls.json`. These are corroborating finite-exact controls; the unbounded conclusions rely on the proof above.

### Two concrete controls

1. Q=C3 with children TT3,TT2,TT2; quotient is conversed, distinguished TT3 is unchanged, and one other pair is reversed. On the same physical names, (a,c)=(25,6) gives2715 versus2149, gap566. In a cyclic two-vertex exterior, lifted masks(281,262) give44727 versus37437, gap7290. Both gaps are2 modulo4.

2. Q8=(238,4,240,118,195,146,34,74), all eight children TT2, distinguished child0 unchanged; all other interiors are also unchanged. The complete child-uniform contrast is zero modulo4. The anchored physical masks(a,c)=(13,0) give2,117,386,219 versus1,398,623,265, gap718,762,954. In the cyclic two-vertex exterior, masks(131085,131072) give129,428,799,899 versus88,692,826,785, gap40,735,973,114. Both anchored gaps are2 modulo4.

## Reproduction

Run `python3 audit.py` and `python3 -O audit.py` in this directory. The checks use explicit exceptions and stay active under optimization. Run `python3 q8_checks.py` for the larger physical controls; it uses the installed C++17 compiler and performs no installation. `optimization_match.json` records comparison of all substantive normal/optimized output fields, excluding elapsed time. The directory contains no external write operation.
