# Independent audit of the four-arc spectrum patch

**Verdict: PASS for the stated construction, sharp finite seed claims, and two fixed-parent exclusions.** Global spectrum completeness remains open. This review binds the final producer proof whose SHA-256 is `d1e5debdf607c3eb83f8054c08f94ec1137ef0b4b17aac775973f878429d5479`; `audit_receipt.json` pins every producer release file and the independent checks.

## 1. Literal seed and minimality

The independent counter visits actual directed simple vertex words and records their vertex sets, first vertices, and last vertices. It does not import the producer's code or recurrence. All five physical seven-name tournaments are valid and strong. At the fixed root 2 the triples `(H,s,F)` are exactly

`(29,1,157), (41,2,134), (29,1,89), (25,1,121), (31,1,157)`.

The flips are exactly `{1,2}`, `{0,1}`, `{2,3}`, `{3,5}`. Every step changes one arc.

The original second representative has rows `(38,28,8,1,13,94,31)` and root 3. Relabeling it by `p=(2,0,3,1,5,4,6)`, where new label `i` carries old label `p[i]`, gives the producer's final seed at root 2. Exhausting all 720 compatible identifications gives minimum arc distance four, attained by five identifications. All 120 shortest flip orders were examined, retaining the complete 80 Boolean-cube vertices. Exactly 38 paths remain strong throughout. None retains `s=1` throughout, and the smallest possible peak start count is two. An additional independently checked seed statistic is that the smallest possible peak Hamilton count among these all-strong shortest paths is 41.

The independent small census also reproduces all 32 unit-start `(H,F)` profiles through order six and finds no equal-F pair differing by two. Thus the minimum **maximum component order** is seven, under this exact twin criterion. No statement about arbitrary endpoints, moved roots, converse equivalence, longer paths, or enlarged tournaments follows from this census.

## 2. Exact exterior and nonovershoot proof

Let `Q=Z_k -> T`, where `Z_k` is transitive, and add `w` with sole incoming arc `r->w`. For any rooted tournament `Q`, split Hamilton words at `w`. Words starting at `w` contribute `H(Q)-s_r(Q)`; all others consist of a first nonempty path ending at `r`, then `w`, then a possibly empty second path, contributing `F_Q(r)`. Hence

`H(Q+w) = H(Q) + F_Q(r) - s_r(Q)`.

For `k>=1`, every Q-Hamilton word traverses all of Z before T, so `H(Q)=H(T)` and `s_r(Q)=0`. In each rooted two-bin cover, the transitive Z names may independently choose the first or second bin; their internal orders are forced. Thus `F_Q(r)=2^k F_T(r)` and

`H(W_k(T,r)) = H(T)+2^k F_T(r)`.

The hypothesis `k>=1` is load-bearing: at `k=0` the formula is `H(T)+F_T(r)-s_r(T)`. Literal zero-clock controls give counts `185,173,117,145,187`, but are not used to silently extend the producer's theorem.

For `k>=1`, the exact strongness criterion is that every vertex of T can reach r within T. Necessity follows because r is the only exit from T. For sufficiency, every T vertex reaches r, then w, then all Z vertices; each Z vertex reaches all of T. Strong T therefore suffices at every displayed stage.

Put `x=2^k>=2`. The three internal counts minus the lower endpoint are `12-23x`, `-68x`, and `-4-36x`, all strictly negative. The final endpoint is larger by exactly two. This proves the all-k nonovershooting construction, rather than extrapolating the finite controls. Independent literal counts were computed at k=0,...,4, and all 855 rooted wrapper cases on every labeled tournament of order 1 through 4 with k=0,1,2 agree with the formula and the exact strongness criterion.

## 3. The preserved interface is narrow

The full displayed seven-coordinate endpoint data `(H,s,e,q,F,S2,P2)` are independently reproduced as

- first endpoint: `(29,1,10,15,157,34,580)`;
- last endpoint: `(31,1,9,9,157,31,644)`.

Even the size-refined summands of F differ: including the empty size slot, they are `(0,15,27,27,27,27,24,10)` and `(0,9,25,31,28,30,25,9)`. Equality of F is equality of one scalar two-bin observable, not a full rooted response or attachment equivalence. This exterior reads only H and F. Four literal old-arc flips lift, but their minimality is not claimed after freely relabeling the added vertices.

The lower endpoint has residue 0 or 1 modulo three, so this family cannot itself furnish the ordered ternary sibling pair `(3H+2,3H+4)`. Equal-H unit-start comparisons satisfy the stated even-response difference restrictions. Neither an indefinite successor operation nor a return/reset mechanism has been proved.

## 4. Independent parent-profile exclusion census

The separate C++ audit has no dependency on producer source. It enumerates actual Hamilton words, stopping at 202 only when all that remains relevant is whether H is at most 201. Every candidate marked F profile is then evaluated by enumerating all actual partial directed words, with no cutoff. Strongness uses Floyd-Warshall reachability, rather than the producer's closure loop. Canonicalization minimizes integer upper-edge codes within degree classes, rather than producer row tuples.

The inherited complete extension cover is included locally as `inherited_cover8.json`. The audit rechecks its empty base, every extension choice, every supplied permutation as an actual bijection, and every directed edge under that permutation: 62,423 explicit transports through order eight. Its last bank agrees exactly with the pinned 6,880-representative input. This inductive cover proves completeness without relying merely on the known class count; irredundancy is not needed for exclusion.

The independent census matches all complete target-F sets and multiplicities:

- Every order-eight representative and all eight deletions were checked. Maximum H among the seven-name deletions is 189. Insertion injectivity extends this upper bound to every smaller order.
- All 1,761,280 one-vertex extensions from the order-eight bank were checked. The 5,693 strong low-H extension records reduce to 545 order-nine isomorphism classes. The independently chosen representatives and producer representatives have exactly the same integer canonical codes.
- All 279,040 extensions of those 545 classes were checked at order ten.
- Neither rooted profile `(191,1273)` nor `(201,5303)` occurs, with no unit-start filter imposed on candidate parents.

Completeness of the low-H reduction follows from strong one-vertex deletion and Hamilton-path insertion monotonicity. The latter map is explicit: insert the new vertex before the first old vertex it beats, or after the last old vertex if it beats none. The preceding old vertex, when present, beats the new vertex. Deleting the inserted vertex recovers the original word, proving injectivity.

The final all-order cutoff is directly supported by the publisher's [Busch abstract](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v13i1n3): for strong tournaments, `H>=5^((n-1)/3)`. At `n>=11`, this is strictly greater than 201 since `5^10=9,765,625>201^3=8,120,601`. The full Busch PDF was not accessible, but no inaccessible exact residue-dependent floor is required. The exact deletion implication was independently checked in [Moon's original publisher PDF, Theorem 1, page 297](https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf): its cycle-length theorem supplies a strong induced (n-1)-vertex subtournament for n>=4.

For the two explicit larger twin pairs, the coefficient comparison forces precisely those two parent profiles when one fixed parent and the same varying clock are used. The independent simple-path counter also checks the larger twins themselves. If a nonstrong parent yields a strong wrapper, its root belongs to its final strong component C. Writing the preceding vertices as X gives `F_P=P_X(2)F_C`; swapping the two ordered X-bins is fixed-point-free when X is nonempty, making F even. Both target F values are odd. Thus the exclusion also covers nonstrong parents whose wrappers would be strong. Clock-dependent parents, unequal clocks, and other exterior architectures remain outside this argument.

## 5. Reproduction and bounds

Run `python3 independent_audit.py`, `python3 wrapper_audit.py`, `python3 census_audit.py`, and `python3 bind_audit.py` from this directory, or invoke them by absolute path. The census expects the producer packet as sibling `spectrum/` in publication or `session190-spectrum/` in this workspace. No producer code is imported. The inherited cover needed for completeness is included here.

The Python programs use explicit checks that remain active under `-O`. Local, wrapper, and census normal/optimized outputs agree byte-for-byte; the complete census JSON agrees too. C++ uses unsigned 64-bit counts, far above the possible number of partial directed words at order ten, which is below `e*10!`, and marked covers are bounded by `10!`. The bounded Hamilton counter never infers exact values above 201. The full local receipt retains zero coordinates and nonstrong cube vertices.

The proof and computation support the claims stated above only. No claim of global Hamilton-path spectrum closure, unrestricted response injectivity, or literature-wide novelty is approved.
