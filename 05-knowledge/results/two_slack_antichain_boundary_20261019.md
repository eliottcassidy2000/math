# Two-slack maximum antichains: sparse boundaries and a layered theorem

2026-10-19. PROVED structural reductions and layered-family theorem; root-reviewed and frozen. Structural reductions and a proved layered-family theorem. The unrestricted seven-thin width-six good-pair/very-good-pair question remains OPEN here. No historical novelty claim is made.

## 1. Setup and exact ambient defects

Let A be a maximum antichain of size m in a finite poset P. Suppose every a in A has incomparability degree at most m+1. The other m-1 elements of A leave at most two incomparable neighbors outside A.

Every exterior vertex is comparable with some A vertex, since otherwise A could be enlarged. No exterior vertex can lie above one A vertex and below another, by transitivity. Therefore the exterior partitions into D, vertices below some member of A, and U, vertices above some member of A.

For a in A define its lower and upper defects

    L_a = {d in D : d is incomparable with a},
    R_a = {u in U : u is incomparable with a}.

Its full strict downset is D minus L_a, and its full strict upset is U minus R_a. In particular |L_a|+|R_a| <= 2. No exterior distribution or uniform induced-core measure is being substituted for the ambient poset.

An A vertex is pure lower-double if its defect sizes are (2,0), pure upper-double if (0,2), mixed if (1,1), lower-single if (1,0), upper-single if (0,1), and zero if (0,0). These six types exhaust the possibilities.

The results below use Zaguia's very good pair definition: equal downsets and both directed upset differences chains, or the dual condition. The empty set and a singleton are chains. A very good pair supplies a balanced pair and, by the existing good-pair synthesis, weak Gold. The pair supplied need not belong to A.

## 2. Degree-two boundary theorem

Assume A itself contains no very good pair. Let D* and U* be the exterior vertices with at least one defect to A. Join a in A to each of its defects, retaining the lower/upper color on exterior vertices. Then every vertex of this bipartite graph has degree at most two.

### Proof

The A-side bound is the two-slack hypothesis. For the other side observe:

1. All nonempty singleton sets L_a are distinct. Indeed, if L_a=L_b={d}, then a,b have equal full downsets; both upper defects have size at most one, so both directed upset differences have size at most one. This would be a very good pair. This includes mixed and lower-single vertices together, not just mixed vertices.
2. The lower defect doubletons are pairwise disjoint. Their A vertices all have the same full upset U. If two doubletons intersect, their two directed downset differences have size at most one, giving a very good pair in the dual.
3. A lower-single defect cannot lie in a lower defect doubleton: those A vertices again have equal upsets, and their downset differences have size at most one.

Consequently a fixed d can occur in at most one lower doubleton and at most one lower singleton. Its degree is at most two; if its degree is two, its singleton neighbor must be mixed, by item 3. The upper proof is dual.

Thus the boundary is a union of paths and cycles, with at most 2m edges and at most 2m exterior defect vertices. Every degree-two lower exterior vertex joins a pure lower-double A vertex to a mixed A vertex; degree-two upper vertices join a pure upper-double A vertex to a mixed vertex. A cycle therefore encounters A-types in repetitions of

    pure lower-double, mixed, pure upper-double, mixed.

Each cycle has a multiple of four A vertices. For m=6 the only possible cycle is an eight-edge cycle containing four A vertices; there can be at most one such cycle. These are necessary incidence restrictions, not assertions that every such graph is realizable by an ambient poset without very good pairs.

## 3. The lower and upper boundaries are almost antichains

For d in D put S_d={a in A : d is incomparable with a}; define S_u dually for u in U. Under Section 2 every such set has size at most two, including the empty set for nondefect vertices.

If d<e are in D*, then S_d is a subset of S_e. If u<v are in U*, then S_v is a subset of S_u. These are necessary conditions, not converse implications: actual comparabilities must still be checked.

Suppose two distinct lower defect vertices are comparable. Two distinct two-element supports cannot contain each other. Equal two-element supports would require both corresponding A vertices to have lower degree two, whereas a degree-two exterior vertex must join a pure lower-double vertex and a mixed vertex; the latter has lower degree one. Thus equal two-element supports are impossible. Distinct singleton supports are incomparable by inclusion. The only possibilities are equal singleton supports or a singleton contained in a doubleton. In either case the two exterior vertices are precisely the two lower defects of a single pure lower-double A vertex. Therefore comparabilities within D* can occur only on disjoint pairs indexed by pure lower-double A vertices.

In fact at most ONE of those pairs can be comparable. If the two disjoint lower defect doubletons for two A vertices were both chains, those A vertices would have equal upsets and both directed downset differences chains, giving a very good pair. Thus D* is an antichain with at most one two-element chain inserted. The identical conclusion holds for U* by duality.

If a lower-single or zero A vertex is present, there can be no comparable pair in D*: comparison of that A vertex with a pure lower-double A vertex would have chain differences on both sides. The upper statement is dual. It is not correct to insist unconditionally that every pure doubleton be an antichain: one chain doubleton is allowed by these local tests when no small same-side row is present.

### Every lower vertex precedes every upper vertex when m >= 5

For d in D and u in U, at most four members of A belong to S_d union S_u. If m>=5, choose a outside that union. Then d<a<u. Hence

    D < U.

This applies to all lower and upper exterior vertices, including those comparable with every member of A.

### Honest handling of the remaining exterior

Write D0=D minus D* and U0=U minus U*. A vertex z in D0 is below all A. No d in D* can satisfy d<z, since an A vertex missed by d would then lie above d by transitivity. Thus z<d or z is incomparable with d. Dually, u<z or u is incomparable with z for u in U*, z in U0. The internal orders of D0 and U0 remain arbitrary, and these exterior attachments can change downsets/upsets of boundary vertices. They cannot be discarded or assigned uniform weights.

For a seven-thin poset the exact residual bound on d in D* is

    |I(d) intersect D0| <= 7 - |S_d| - |I(d) intersect D*|,

and d has no incomparable neighbors in U. The upper version is dual. This states the remaining uncertainty explicitly.

In a hypothetical seven-thin poset with no very good pair anywhere, packet141 gives width at most six. Since D* and U* each contain at most one comparable pair, each has at most seven vertices. This last bound is conditional on global absence of a very good pair, not merely its absence inside A.

## 4. A uniform three-layer theorem

**Theorem.** Let A be a maximum antichain of size m>=5 in a finite poset P. Suppose each a in A has incomparability degree at most m+1, and both exterior parts D and U from Section 1 are antichains. Then P has a very good pair. In particular every seven-thin poset with a maximum six-antichain whose lower and upper exterior parts are antichains has a very good pair, balance, and weak Gold.

**Proof.** Assume P has no very good pair. Section 2 bounds all exterior supports by two, and Section 3 gives D<U. Since D is an antichain and no vertex outside D can be below a member of D, every d in D has empty downset and full upset

    U union (A minus S_d).

If some a had two lower defects d,e, the supports S_d,S_e would both contain a and have size at most two. Their directed differences would each have size at most one. Therefore d,e would be a very good pair, contradiction. There are no pure lower-double A vertices. Dually there are no pure upper-double A vertices.

All L_a are now empty or singleton, and their nonempty singletons are distinct by Section 2. Thus every d in D misses at most one A vertex. Any two members of D would have equal empty downsets and singleton-or-empty directed upset differences. Hence |D|<=1. Dually |U|<=1.

Each of these at most two exterior vertices can be a defect of at most one A vertex. At least m-2>=3 members of A therefore have no defects at all. They have identical downsets D and upsets U, and any two are exact twins. This is the final contradiction.

This theorem is NOT a statement about all height-three posets. A height-three poset need not admit a maximum antichain whose lower and upper exterior parts are both antichains. The hypothesis on D and U must be verified, rather than inferred from height alone.

## 5. Descent for the all-mixed pattern

Choose a maximum six-antichain A minimizing the size of its generated downideal. Suppose A contains no very good pair and all six vertices are mixed. Their six lower defects are distinct by Section 2. Each such defect misses only its associated A vertex, so their six supports are distinct singletons. The inclusion rule of Section 3 makes them a six-antichain A' in D.

Every member of A' lies below some member of A, and no member of A lies below any member of A'. Consequently the generated downideal of A' is a proper subset of that of A. This contradicts the extremal choice. More generally the whole lower region D has width at most five: any six-antichain there would have a strictly smaller generated downideal. Hence D* has at most six vertices, because its induced order has at most one comparable pair. If |D*|=6, that exceptional comparable pair must actually occur; then there is no lower-single or zero vertex in A by Section 3. The analogous upper restriction holds for a maximum antichain chosen to minimize its generated upideal, but these two different extremal choices must not be silently combined.

Thus the all-mixed case cannot occur for the minimum-downideal choice of A. This argument works for any m>=3 under the same two-slack hypotheses and absence of a very good pair in A; it does not eliminate the pure-double, singleton, and zero cases.

## 6. Exact obstruction to locating the pair inside an extremal antichain

Use three antichain levels

    D={d1,...,d5}, A={a1,...,a6}, U={u1,...,u5}.

Put every d below every u. Between D and A, omit only d_i<a_i for i=1,...,5. Between A and U, omit a_i<u_i for i=1,...,5, and also omit a6<u1 and a6<u2. All other cross-level relations hold. Take no other relations.

Every A vertex has exactly two exterior incomparable neighbors, so its incomparability degree is seven. The D degrees are five, the u1,u2 degrees are six, and the other U degrees are five. The poset is seven-thin and has height three.

A is its UNIQUE maximum antichain, of size six. Indeed, an antichain meeting D can meet A only in the unique a_i missed by a chosen d_i; if it contains that a_i it cannot contain any other d_j, and D alone has size five. D and U cannot coexist in an antichain. An antichain meeting U can use at most its missed A vertices; each u_i misses at most two. For a6 the only compatible upper vertices are u1,u2, so the largest such mixed antichain has size three. U alone has size five. Thus no antichain other than A has size six.

The six lower defect sets of A are five distinct singletons and the empty set, so its six downsets are all distinct. Its six upper defect sets are five distinct singletons and {u1,u2}, so its six upsets are all distinct. No pair inside A is very good, despite A being the unique maximum antichain and therefore every extremal choice of maximum antichain.

The poset DOES have very good pairs: any two d_i have equal empty downsets and singleton directed upset differences. Thus this example refutes only the proposed claim that an extremal maximum six-antichain must itself contain the pair. It is not a counterexample to the layered theorem, balance, Gold, or the open global width-six assertion.

## 7. Dependencies and scope

The boundary graph, comparability, descent, and layered-family proofs are elementary and do not depend on any morphology scan. The optional global boundary-size bound uses [packet141](saturated_antichain_very_good_pair_20261019.md). The inherited Gold implication uses [packet94](good_pair_gold_synthesis_20261019.md). Zaguia's [author-deposited paper](https://arxiv.org/pdf/1610.00809v3), Definition5 on PDF page3 and Definition1/Theorem2 on page2, supplies the terminology and balanced-pair implication.

The unrestricted claim that every seven-thin width-six poset has a good or very good pair is not proved or refuted here. The unresolved part consists of arbitrary lower-only and upper-only attachments to the sparse boundary. The exact fixture illustrates why moving the desired pair into the chosen maximum antichain is not a valid shortcut. No broad poset census, full frontier exclusion, or historical novelty claim is made.

## 8. Selected exact controls

The self-contained standard-library control `../experiments/verify_two_slack_boundary_20261019.py` checks just three prescribed posets, not a poset census:

- The 16-vertex local obstruction above: width six, exactly one maximum antichain, no very good pair inside it, and 14 very good pairs elsewhere.
- An 18-vertex all-mixed construction with six vertices in each layer and a single matching of omissions on each interface: exactly three maximum antichains, with the lower one having a strictly smaller generated downideal than the middle one.
- An 18-vertex local witness with seven lower defect vertices, of which exactly one pair is comparable, six middle vertices, and five upper defect vertices. The middle six-antichain has no very good pair. This demonstrates why Section 3 allows one chain doubleton and why seven defect vertices are not excluded without an extremal antichain choice. The witness has very good pairs elsewhere and is not a global obstruction.

The receipts give complete predecessor masks, incomparability degrees, maximum-antichain counts and very-good-pair counts. Normal and optimized runs agree on all 30 guards. These controls illustrate and regression-check the proved statements; they do not replace the proofs or claim an exhaustive search for a width-six counterexample.
