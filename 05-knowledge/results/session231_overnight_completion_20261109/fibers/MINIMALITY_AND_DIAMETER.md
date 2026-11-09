# Scoped length minimality and linear gain diameter

Status: PROVED in the stated scope, independently audited in ../fiber_audit/MINIMALITY_PLACEMENT_AUDIT.md. No unrestricted completeness or global minimum claim.

## Unit-loop minimum in the complete core-order≤9 universe

Consider a directed closed carry walk through a PHYSICAL tournament-cut core of order≤9, with singleton step labels and scalar gain±2. Quantization implies h|q. The content filter implies gcd(h,A,G)=1 (h is odd). In this saturated case the recurrence gives

    A_i≡2^i A_0 (modh),  G_i≡2^(−i) G_0 (modh).

Closing the walk, together with gcd(h,A_0,G_0)=1, forces2^length≡1 moduloh. Hence ord_h(2) divides its length. For h=1 the core is transitive; every padded deletion graph is transitive, and its A uniquely determines its entire cut word. No nonzero physical-loop gain is possible there.

The complete saturated banks through core order9 contain302644 distinct seven-coordinate signature rows, of which230568 have content one. Among these, the only banks with ord_h(2)<12 are

    (33,660),(33,1320),(33,2640),(341,5456),(341,10912),

all with order10. The complete physical A/G seed counts are176,704,2112,492,1968, respectively. An exact traversal of ALL ten-edge directed walks from all5452 seeds examines108050 partial walks and finds NO closed walk. Therefore there is no unit-gain moving loop of length<12 through any physical core of order≤9. The explicit eight-vertex SHORT12 witness attains12 in this universe.

This is a finite core-order-bounded minimality theorem for this specific moving-loop mechanism. It does not rule out shorter loops based on larger cores, other construction grammars, or more efficient unrelated interval constructions.

`MINIMALITY_SEEDS.json` is the compact complete filtered seed list with source-order counts and source SHA-256 hashes. `check_shorter_loops.py` checks all walks from that list. With `--regenerate`, it reads the complete `SATURATED3..9.jsonl` outputs of `saturation_bank.cpp`, independently applies content/order filters, rewrites the compact list, and repeats the walk check. The generator and full regeneration route are in README.md. Default and optimized receipts agree.

## Linear diameter for a fixed one-core carrier bank

Fix h,q and any physical core bank with these counts. Let a full-interface fiber be formed by placing one such core among L singleton vertices, allowing all cut labels and core positions, with common total q/h. Any two equal-A/G constructions cancel common prefix/suffix bits and reduce to a moving path of length d≤L. This is the exact decoding argument in PROOF.md.

Every intermediate carry pair obeys−h<a<q and−h<g<q, while q≥2h. The interior cocycle

    w=a y+h p y−x g'−h x x'

therefore satisfies|w|≤2q+h. The first and last physical core T values lie in[0,q−2h]. The endpoint corrections contribute at most4h in total. Consequently the absolute T difference, and hence the star-closure count difference inside a common full interface, is at most

    L(2q+h)+q+2h.                                   (D)

Here is the endpoint calculation, so no formal intermediate is mistaken for a physical core. Put p_0=0 and x_d=g_d mod2. Choose formal intermediate endpoint coordinates S_i=h x_i and E_i=h y_(i−1), with E_0=0; propagate T by the cocycle. Each formal identity C_i⇒s_(y_i)=s_(x_i)⇒C_(i+1) then holds in all seven coordinates. Restore actual initial/final cores C,D and exterior boundary bits u,v. The resulting T difference is

    sum(w)+T_C−T_D+E_C y_0
       −x_(d−1)(S_D−h x_d)
       +u(S_C−h x_0)+v(h y_(d−1)−E_D).

The four correction terms have absolute values at mosth each. Together with|T_C−T_D|≤q−2h this proves(D). For d=0 use the direct core bound instead.

Thus fixed-h,q one-core movement has at most LINEAR scalar diameter in padding length. It cannot by itself create exponentially wide intervals from linearly many singleton positions. In the strong-flank family, scalar size is at least(q/2)2^L, so its attainable interval width is at most O(log target) for fixed h,q. The bound does not obstruct a baseline selector using many different contexts, and says nothing uniform when the carrier itself varies with the target.

## Compatible positions in the short witness

The complete same-(h,q)=(45,1440) bank through order9 was checked along the repeated SHORT12 word. The only physical intermediate positions are the multiples of12. For k≤12, the complete fiber along that path has exactlyk+1 values and maximal runk+1; no additional fractional-period core positions enlarge it. The twelve-step period can also be checked once and then repeated algebraically. This is relative to the complete order≤9 bank, not a claim that no larger core with the same h,q can occur.
