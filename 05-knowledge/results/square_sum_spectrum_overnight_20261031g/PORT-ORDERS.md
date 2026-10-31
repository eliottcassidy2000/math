# Rooted order, shared constraints, and the Hamilton polynomial of a closed square-sum carrier

**Status: PROVED for the authenticated cycle families; graph-complete after primitive refinement. Arbitrary square-sum graph orders and arbitrary forest extension remain OPEN.**

Use packet f's fixed-phase edge cube at a legal size N and schedule, with b source bits. A fully refined schedule has primitive leaves in112..2811. It gives a subcubic spanning square-sum graph G with exactly2^b undirected Hamiltonian cycles. The native protected construction retains its original eleven-vertex path, seed and integer labels.

## 1. Normalize the observation before comparing orders

An undirected cycle has neither a start nor a direction. Every cycle here has the invariant anchor edge1--a, where a=3 in the ordinary carrier and a=8 in the protected native carrier. The **rooted order** is its unique cyclic listing beginning1,a. This is the convention throughout; it is not an extra orientation choice counted as another Hamiltonian cycle.

The native template's stored reference happens to run in the opposite direction. The first order consumer exposed that mismatch even though its edge sets were correct. `Orders` rotates the reference to1 and reverses its remaining list when needed, before deriving intervals or returning order claims. No underlying cycle theorem changes.

## 2. Disjoint intervals, with ties explicitly retained

Each primitive reversal or even-block exchange changes a contiguous reference interval, including its two unchanged ordered endpoints. It preserves the interval's vertex set and those endpoints. Different primitive leaf intervals are disjoint: each lies entirely in its own copied leaf, and every recursive map is a fixed label injection and a fixed choice of traversal direction. A split head removes only the invariant1--3 boundary, outside a varying interval. The protected outer cuts likewise avoid all changing intervals.

The same statement holds if several primitive intervals share one coarse source bit. Their motions are then tied; they are not independent regions. Refinement unties them and leaves the union graph unchanged, as packet f proves.

One can recover the primitive intervals directly from the edge atlas. Each connected component of a bit's removed/added edges is a variable4- or6-cycle. Its first and last positions in the normalized reference enclose the changed interval. Apply that component's edge replacement and walk the resulting path between those same ends. This gives its two ordered lists, without enumerating every global cycle. `order_atlas.py` checks this representation against actual full cycle walking, including reversal, exchange, tied sources and the native two-layer carrier.

## 3. Complete order compiler and a sharper port budget

Given a list Q of distinct marked labels in the desired rooted order:

1. Assign each marked label to its primitive interval, or to its own fixed singleton position outside all intervals.
2. Compress consecutive occurrences of the same interval. This region sequence must agree with the reference sequence restricted to the same labels. Otherwise no source choice can satisfy Q: no interval changes its global location.
3. Within each touched interval, the requested label order must match one or both of its two local lists. Matching neither rejects the request. Matching exactly one fixes the controlling source bit.
4. Retain shared-bit constraints across different intervals and across different requests. Two contradictory values reject the conjunction.

These conditions are necessary by interval preservation. They are sufficient because every consistent assignment of the remaining bits is an actual cycle. If r distinct bits are fixed, the extension count is exactly2^(b-r). For a fully refined carrier these are **all** Hamiltonian extensions in G, by packet f's closed-carrier theorem. For a coarse code the answer is family-relative, and may omit cycles of the same union graph.

**Reference-order budget.** Let Q be the reference rooted order of any k distinct marked vertices. Then at least

    2^max(b-floor(k/2),0)

cycles retain that order. Indeed, an interval with zero or one marked label has the same restricted order in both states. Every fixed bit therefore has a witnessing interval with at least two marked vertices. Pick one such interval per fixed bit. These intervals are disjoint, so2r<=k. The reference assignment ensures consistency.

The same bound holds for any collection of reference-order requests involving altogether k distinct labels. A required edge or a complete reference neighbourhood supplies additional single-bit constraints as in packet f. They must be merged, not tested separately. A reference edge request costs at most one additional bit; the combined conservative bound with e such edges is2^max(b-e-floor(k/2),0).

The bound is sharp within these families. Choose an inversion pair in each of r different primitive bits. Their2r labels, ordered as in the reference, fix exactly r bits. A single label per primitive region fixes none, even when many regions share coarse bits.

### Hostiles retain the missing coordinate

At N2812, the port order1438 before1587 is independently feasible, as is requiring edge162--1438. Their conjunction fixes the same bit to opposite values and is impossible.

At N70312, two ordered pairs are

    (6562,33438),    (39088,36537).

They are incompatible in the25-bit coarse A family because the two regions share one bit. They are compatible in the625-bit AA refinement of the identical union graph, leaving623 free bits. The returned actual cycle satisfies both. Thus source ties and graph closure are different predicates even for order observations.

At N270112, all protected prefix labels are constant, while an additional inversion pair fixes one of2401 bits and leaves2400. The native label clock and marked orientation are both checked by the consumer.

## 4. Exact Hamilton polynomial and weighted optimization

Write C for the set of edges common to every cycle. For primitive bit i, let R_i and A_i be its two disjoint variable matchings. All these sets are edge-disjoint, and

    E(C_z)=C union union_i (R_i if z_i=0 else A_i).

For commuting edge indeterminates x_e, the Hamilton-cycle polynomial of a fully refined G is exactly

    Z_G(x) = (product_{e in C}x_e)
             product_i[(product_{e in R_i}x_e)+(product_{e in A_i}x_e)].

This follows by expanding the product and using the bijection between bit strings and Hamiltonian cycles. Coarse tied-bit versions factor in the same way for their chosen family, but do not give the full graph polynomial.

Consistent edge/order queries fix selected factors to one monomial; a contradictory query gives zero. This permits exact weighted counting over nonnegative integer or rational weights. Arbitrary signed additive edge costs are also simple: include the fixed common-edge cost, then choose the cheaper matching independently at each free bit. Zero weights can forbid a matching; both zero means no positive-weight extension.

`joint_ports.py` compares a64-cycle fibre's polynomial specialization with the sum over independently walked actual cycles. It also checks the minimum of signed edge costs against complete enumeration. The optimizer writes each choice into its **numeric bit address**; template dictionary iteration order is not the source address order. The first audit caught that implementation error before release.

Deleting the invariant anchor edge gives the same exact count and polynomial (with its fixed anchor factor removed) for Hamilton paths between the named anchor endpoints. Their rooted path order is obtained by the known global reversal of the cycle listing after its initial vertex. In the protected version deleting1--8 also removes that edge of the native prefix; no claim that the entire original prefix survives that particular deletion is made.

## Limits and next object

These results solve order and weighted extension within a deliberately constructed square-sum subgraph. They do not say every permissible geometric edge or ordering in the full square-sum graph extends. A connector joining independently constructed carriers must preserve its incidences, rooted region sequence and shared bit constraints; mere existence of each local path is insufficient.

The useful next extension is a two-carrier connector with an authenticated endpoint relation. The present compiler supplies the exact local boundary language needed for that test, without replacing it by a graph minor or a count alone.
