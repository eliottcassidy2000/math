# Independent adversarial audit 2

**Verdict: PASS, mathematical proof audit.** Both all-order candidates are sound under the stated tournament/observer conventions. The finite replay below is supplementary, not the reason for the verdict. No publication, canon promotion, or source-checkout mutation was performed.

Inputs: `bridge/BRIDGE_SINK_STAR_LEMMA.md`, `triangularity/AVOIDING_WORD_TRIANGULARITY.md`, and the actual grammar/observer definition from pinned commit `d3132eafbfa9c43c6d9f1110566b1e9c3a3286e1`, file `05-knowledge/results/operation_closed_synthesis_20261029/THEORY.md`. Read source AGENTS, active guardrails and research protocol. The shorter path supplied in the assignment omits `05-knowledge/results/`; inspecting the tree resolved that path without changing the pinned revision.

## Exact scope of the combined conclusion

For every finite legal word of the pinned abcde insertion grammar from its specified four-vertex seed, including the empty word, the integer **native ordered-two-path spanning-cover matrix**, with nonempty paths and singleton paths allowed, is singular if and only if the word contains `be` or `ee` as consecutive letters. If it avoids those pairs, it is permutation-equivalent to a unit upper-triangular matrix, hence unimodular, with exactly one perfect matching in its bipartite support. If it contains either pair, a nonzero integer antisymmetric star lies in its right kernel.

This does not assert the same word criterion for arbitrary tournaments, arbitrary insertion conventions, unordered observers, path covers disallowing singletons, or other k-path observers. Turning grammar coverage into a classification of all unique-Hamilton-cycle tournaments uses the separately cited grammar-completeness theorem; the proofs under audit do not prove that classical input anew.

## 1. Bridge lemma: independent sign and exhaustion review

Let `S((u,t),(x,y)) = M_R[(u,t),(x,y)] - M_R[(u,t),(y,x)]`. This convention makes S antisymmetric in both ordered pair arguments, and it is exactly the column action of `e_x wedge e_y` in the native observer (no division by two or lattice identification).

If a's unique in-neighbor in R is s, a two-path cover ending at s and a has the a-path singleton: otherwise its last edge must be s->a, impossible because s is already the other endpoint. For fixed second start t, the positive contribution to `sum_u S((u,t),(s,a))` is `1_(t=a)L`, and the negative contribution is `h(t)`. Thus the proposed identity and its signs are correct. The assumption `|R|>=3` justifies prepending a to a Hamilton path in R-a ending at s: its initial vertex cannot be s. No unproved equality of individual Hamilton-path counts is used; the relevant equality is of their sums.

For a row `(z,t)` with z in P=Q-r and t in R, all of P lies on the first path and its last P-vertex is c. Its final connection to r yields the **negative** term `-F(z)H(t,x)`. If it enters R, the bridge extension is counted by `F(z)sum_u S((u,t),(x,q))`. For q=a, the two values reduce to `-F(z)H(t,a)`; for q=s, they reduce to `-F(z)H(t,s)`. Therefore the sign-sensitive formula really proves equality of the a,r and s,r antisymmetric columns for both bridge choices.

For both starts in P, first fix the two ordered initial P-segments. Their terminal vertices dominate every R vertex. When both continue into R, summing the signed R-cover over all ordered starts (u,v) gives zero by swapping u,v. When only one continues into R, its unrestricted initial vertex produces `sum_u H(u,x)`, equal for a and s. All factors from P and r are independent of x. Empty R segments on one path are therefore included, rather than silently excluded. If a start is r, its path must be singleton r because r is an endpoint; the other path traverses P and then R. Both starts in R are impossible because nonempty P cannot be reached. These exhaust rows, including reversed order. There is no missing start at c or a: each belongs to an already treated set.

This also checks the Q boundary `|Q|=2`: P={c} has a one-vertex Hamilton path, the both-starts-in-P case is simply absent, and the other formulas remain valid.

## 2. Birth and frozen-sink persistence

Useful elementary grammar invariants, which can be stated explicitly to make the candidate self-contained:

* Every maximum-outdegree vertex has indegree one; there are exactly one (A) or two (B).
* The A maximum is the newest vertex.
* In B the maxima a0->a1 are consecutive on the marked Hamilton cycle.

They hold at the seed. On insertion the new vertex has indegree one, the parent retains its indegree, and every other old vertex gains one. The grammar parent choices and cycle insertion immediately preserve the three invariants; in particular only the old indegree-one parent can remain a maximum alongside the new vertex.

At `be`, the A maximum r was created by the preceding a or c, inserted after the predecessor s of an old indegree-one maximum a. The b,e additions b,c leave the triangle r->b->c->r, all b,c->R, and precisely s->r from R to Q. At `ee`, before these letters the adjacent B maxima are a->r. Deleting r leaves a's sole in-neighbor s unchanged, and a->r is the sole R->r edge. These statements use actual parents, not just state names. Minimum legal examples cause no tiny-R exception: `ee` has |R|=3, while `be` requires a prior c (or a later A-producing history).

At a last bad pair choose this partition. Every R vertex receives edges from b,c; r receives edges from q,c. Therefore neither R vertices nor r can be maxima. No R vertex can precede a maximum on the cycle: its only outward edge to Q is q->r and r is not maximum. All future parents consequently remain in Q, inductively preserving the sole R->Q bridge.

The frozen c has outgoing neighbors exactly R union {r}, all of indegree at least two. It cannot precede a maximum. The first suffix letter is c or d, not e; these parents are respectively r and b at the birth state, not c. Thus after the first insertion c has indegree at least two. Before a hypothetical first selection of c as parent, all added vertices dominate c, so its outgoing set stays R union {r}; hence it is neither maximum nor predecessor of one. Contradiction. The empty suffix already satisfies the lemma. Consequently c remains a sink of Q-r and the fixed nonzero star persists. This proves the forbidden direction for arbitrarily long suffixes without claiming monotone nullity or preservation of arbitrary prior kernel vectors.

## 3. DAG and triangularity audit

For a simple loopless digraph with a marked directed Hamilton cycle C, define D by `s(u)->v` for each noncycle edge u->v. A topological weight f for D assigns zero defect to cycle edges and positive defect to every other edge. Summing defects over a spanning two-path cover telescopes to

`f(s(c))+f(s(d))-f(a)-f(b)`.

The cancellation uses every vertex once as a head except starts and once as a tail except ends, and remains valid for singleton paths. Equality means all used edges are cycle edges. Cutting C before starts a,b then gives exactly one cover with ends `(pred(b),pred(a))`, including adjacent-start and singleton cases. Thus reindexing these columns gives the unit diagonal; all other nonzero entries have strictly larger column weight. Tied pair weights cause no problem: equality occurs only at the actual matched column, not just some column with the same endpoint set.

The insertion update is exact. Aside from the new source x and the replaced outgoing row at q=s(p), the old D rows are unchanged. The q-row omits p and q. In the second preservation case this exclusion of p is load-bearing: p remains a source after removal of x, so the only possible predecessor of q cannot lie on a directed cycle. The first case is immediate. For clarity, “old sources other than q” in the update list should mean old **arc-origin vertices** other than q, not only zero-indegree sources; the stronger precise statement is true and used in the proof.

For a,c,d the successor q is an indegree-one maximum, hence a D-source. Before b, the previous a/c inserted the newest p before a then-indegree-one q, so q's unique D-incoming neighbor is p, itself a source. Before an avoiding e, either this is the initial e or the preceding letter is d. Initial e has p=1,q=2 and q's only noncycle incoming tournament edge is 0->2, giving D-edge 1->2; d produces exactly the same local structure as a/c. The seed D consists of 1->2->3, with isolated 0. This closes the induction with no presumption of the conclusion for forbidden words.

## 4. All-k extension: proved by the same potential

The vertex-DAG lemma extends to every ordered k-path observer `M_k`, for `1<=k<=n`, with nonempty paths and singletons allowed. For ordered distinct starts `A=(a_1,...,a_k)`, let `next_A(a_i)` be the next selected start encountered around C (for k=1, it is a_i itself). Define

`E(A)_i = pred(next_A(a_i))`.

Cutting C before the selected starts yields exactly one ordered cover with starts A and ends E(A). This endpoint map is a bijection: given any ordered distinct ends E, cut C after their underlying set and name each resulting component by its prescribed endpoint index; its ordered starts give the unique inverse. This argument is necessary for k>=3; simply rotating the input tuple is wrong unless the tuple was already in cyclic order.

For any ordered k-cover, the summed defect is

`sum_i f(s(e_i)) - sum_i f(a_i)`.

It vanishes precisely for that canonical cycle-only cover. Since the successors of E(A) form the same **set** as A, reindexing by the bijection E and sorting by total start weight gives unit upper triangular M_k. Hence every avoiding grammar tournament has unimodular M_k and unique support perfect matching for every k. This extension is a sufficient DAG consequence only; do not extend the forbidden-word singularity criterion to all k. In particular `M_n` is the identity for every tournament, since all paths are singleton.

## 5. Independent finite checks (corroboration only)

Reproduce: `python3 session103_words_necklaces/audit2/check_independent.py` from the workspace root. The checker imports none of the producer code and enumerates path covers literally by all vertex permutations and cuts.

* Enumerates all tournaments containing the standard marked directed cycle at n=3,4,5; for each with acyclic D, checks every M_k, every canonical matching coefficient, and strict weight increase of every other nonzero entry. DAG tournament counts are 1,4,15; 94 observer matrices checked.
* Enumerates all allowed R on three vertices with prescribed a<-s, both bridge choices, and every Q on 2,3,4 vertices satisfying the rooted-sink hypothesis. Checks the native signed star residual directly for every ordered start pair: 44 cases, including Q singleton-complement boundary.
* Receipt: `receipt.json`, status PASS. Neither these bounds nor the producer's larger census is substituted for the proofs above.

## Portable replay hardening

The four original `assert` checks in audit2/check_independent.py were replaced by explicit raising guards, preserving the same checks under `python -O`. The root `verify_all.py` runs nine selected producers in both normal and optimized modes, each in a fresh explicit-manifest temporary copy, deletes the target receipt before execution, and compares regenerated JSON with its baseline. All 18 runs passed; original input/script/receipt hashes remained unchanged. A separate deliberately corrupted baseline control was correctly rejected after regeneration. See root `verification_receipt.json`.
