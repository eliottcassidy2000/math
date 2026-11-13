# Audit: the complete 881 fibre excludes the next fixed single-module ray

**Status:** PROVED conditional on the explicitly cited inherited complete ready-fibre classification; the implication is audited here. The inherited classification itself is PROVED in its source. Finite profile checks are stated separately. No new Busch theorem or external literature claim is made. This report audits the root's corollary and is not a priority claim.

**Incoming source pin:** [Session296 cores/PROOF, Section4](../../session296_multicore_metagraph_20261113/cores/PROOF.md), as read in `origin/main` at commit `a03df4746338daf016df8819bcce373ec35c95e8`. It classifies every rooted strong indegree-one defect-one source of full H=881 as one of the ready births from the unique-start addresses `(11,1,11)` and `(1,empty,00101)`, of source orders10 and11. It also states that the weak ready sources are precisely those strong sources followed by a transitive tail. The source supplies its ready inverse, complete grammar inversion and inherited strong-order bound; this audit uses that conclusion and does not re-prove the cited outside bound.

## Exact excluded universe

Fix one finite tournament context with a single nonempty uniform module. Replace the module by TT_k. Fix a root z and its proposed sole predecessor q among the **external** vertices, the same physical names for every k. Suppose that for every integer k above some fixed K:

- z has indegree one, with sole predecessor q;
- the defect `H(T_k-z)-s_z(T_k)` is1;
- the full count is `H(T_k)=405*2^k+71`.

This fixed context cannot exist. Strongness of the whole T_k is not required: the inherited weak-source classification accounts for source padding.

## Why a large-k-only identity still reaches k=1

Use the independently audited nonnegative skeleton expansion from [MODULE_AUDIT](MODULE_AUDIT.md):

`H(T_k)=sum_j A_j j!S(k,j)`, `A_j>=0`.

If any A_j with j>=3 were positive, its summand would grow asymptotically like A_j*j^k, contradicting the asserted O(2^k) bound. Therefore all higher coefficients vanish. The remaining formula is

`H(T_k)=A_1+A_2(2^k-2)`.

Its equality to `405*2^k+71` on the infinite tail fixes `A_2=405`, `A_1=881`. Thus it holds for every k>=1, especially `H(T_1)=881`.

The root-star condition also holds at k=1: all external incidences are fixed, and a uniform module either entirely beats z or entirely loses to z. Since its vertices are not q and z has only incoming q at large k, z beats the entire module at every size.

Let U_k=T_k-z. Because z beats every old vertex except q, a Hamilton path starting at z is exactly z prepended to a Hamilton path of U_k whose first vertex is not q. Consequently

`H(U_k)-s_z(T_k)=s_q(U_k)`.

This defect is itself a nonnegative uniform-module path-count observer with a fixed **external** start. If any of its j>=2 coefficients were positive, it would grow without bound. Its eventual constant value1 therefore makes all those coefficients zero and its j=1 coefficient one. Thus the defect at k=1 is also1.

We have obtained an actual ready source T_1 at881, with its original physical root and predecessor. This step does not assume that finite interpolation by itself controls unseen higher coefficients: positivity and the all-tail identity are essential.

## Strong and weak cases

If T_1 is strong, the complete inherited classification puts its rooted source in exactly one of the two listed classes. The single module has become one vertex v distinct from z and q. Replacing TT_1 by TT_k is exactly the inflation of that non-port vertex v in the classified source. Rooted isomorphism causes no loss: enumerate every such vertex, so every possible image of v is included.

If T_1 is weak, the inherited classification writes it as one of those strong ready cores followed by a transitive tail. Root and predecessor lie in the ready core.

- If v lies in the tail, replacing it by TT_k leaves the tail transitive. H stays881, giving A_2=0.
- If v lies in the strong core, the trailing transitive component contributes multiplicative H=1. Its presence therefore leaves the same full-count profile as inflating v in the strong core.

Thus arbitrary finite weak padding creates no additional nonzero A_2 profile. This addresses padding in every fixed context, not merely padding up to a searched size.

For completeness, strongness also behaves correctly under contraction here: the quotient of a strongly connected tournament by a nonempty proper module is strongly connected. Conversely a nontrivial strongly connected quotient with nonempty substituted modules yields a strongly connected tournament. This is not needed to exclude the weak case but prevents an unrecorded strong/weak transition under TT_k collapse.

## Exact vertex profiles required for the contradiction

The root's count engine lists A_2=(F_2-F_1)/2 for every non-port vertex of the two source representatives. Its lists are:

- address `(11,1,11)`: `881,1529,1038,678,678,657,462,462`;
- address `(1,empty,00101)`: `770,913,824,617,430,430,501,339,673`.

None equals405. Some vertices have higher A_j as well; those already violate the growth requirement, so their inclusion in the list is conservative. The list covers all8+9=17 non-port vertices. The independent `verify_881_profiles.py` reconstructs both graphs directly from the funnel words, without importing the root's compiler or counter, then counts their inflated tournaments by a separate subset DP. It verifies the same coefficient multisets (label orders may differ), graph types, sole predecessor, full H881, literal defect1, and parent unique-start count1. Both normal and `-O` runs pass63 explicit checks on all17 vertices. The exact receipt is `881_profiles.json`, including normalized source hashes at the incoming commit. Completeness still comes from the inherited theorem, not the finite recomputation.

The required `A_2=405` is therefore impossible both in the strong case and after arbitrary weak transitive-tail padding. This excludes the entire specified fixed single-module family, including a construction asserted only for all sufficiently large k.

## Boundaries and the next question

The conclusion does not exclude variable contexts, two or more changing modules, nonuniform attachments, root/predecessor inside the module, k-dependent port names, or a construction along a selected finite set of clocks. Nor does the wording automatically cover a differently parameterized shrinking module `TT_(k-c)`: its k=1 extrapolation is not the same881 base problem. A positive fixed shift `TT_(k+c)` can be encoded by c fixed transitive vertices in an enlarged uniform context, but clock conventions must be stated rather than silently exchanged.

The obstruction is to one precise **uniform** return source. The complete fibre is nonempty, and the inherited finite ladder supplies actual ready sources at several later clocks. The classification therefore identifies a missing coordinate in this method: the next ray needs additional changing structure, a different clock, a moving port, or a nonuniform attachment. It does not refute Chan-Pak or the possibility of another repeatable return.
