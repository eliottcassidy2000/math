# Exact two-address compensation and the selection bottleneck

**New result.** A fixed asymmetric tournament core has an infinite family of exact ±1 centers even though neither of its continuation sides can be varied alone. The two addresses pay each other's signed drift. This does not give universal absolute placement: on a particularly simple synchronized branch, the distinct center count through X is only O(X^log_9(2)). The general inverse below exposes where this kind of repair spends the variation needed for coverage.

## 1. The continuation state of an arbitrary two-port core

Take any nonempty tournament N of nonports and fix two masks I_x,I_a of vertices beating each port. The variants differ only in x→a (minus) versus a→x (plus). Write D=H(plus)−H(minus), and f=I_x−I_a as a signed vector. Let S,E be the Hamiltonian start/end vectors of N. Let M(u,v) count ordered two-path spanning covers with the first ending at u and the second starting at v. Define

    B_L = <f,S>,                      B_R = <f,E>,
    U_L = −<f,E + row sums of M>,     U_R = −<f,S + column sums of M>.

The bare gap, if needed without computing a whole tournament count, is
D=I_x^T(M−M^T)I_a−<f,S+E+row sums of M>. This counts only paths using the switched arc: contract the adjacent ports, then split into zero or one nonempty path on each side. The cases with an empty path give the separate S and E terms. The remaining case uses M.

Now extend along the nonport order by any shared word at the left and another shared word at the right. At the right the neutral letter is a global sink, and the return letter beats both ports while losing to every nonport. At the left the neutral letter is a global source, and a return loses to both ports while beating every nonport. Let l,r be the binary addresses, with a bit for a return at its chronological birth step, starting at bit zero separately on both sides. Then

    D_final = D + U_L l + U_R r + B_L(l mod 2) + B_R(r mod 2).       (1)

**Proof.** The [one-sided calculation](README.md) works with arbitrary N already assembled: the right-nonempty signed mass is −<f,S+col M>, since it checks the arc leaving the second port and allows the left side to be empty. The right-empty mass is zero before checking the incoming arc, because it merely reads H(N). The incompatible part of that mass is <f,E>. Thus the right correction is exactly U_R r+B_R(r mod 2). Reverse every arc and swap names x,a to keep the two switch branches in their original order; the same argument gives the left correction. Prepending a source relative to N preserves the **right** data: Hamiltonian ends and weighted starts at the old vertices. Indeed a path starting at an old vertex must leave the new source to its complementary path, where it attaches uniquely; for a full path ending at an old vertex prepend the source. Its own contribution is zero since its port contacts are shared. Similarly a right append leaves the left data unchanged. This proves independent composition.

The empty read and the four probes consisting of one or two consecutive returns at either end recover the whole state, since e.g. at right their increments are U_R+B_R and 3U_R+B_R. Hence these five scalars are not just sufficient for every future *gap* read: equality of them is exactly equality of all future gap reads. They need not determine the *center*.

Consequences:

- For right-only continuations, after appending a neutral sink (no count change), exact stability is equivalent to equality of one weighted start sum on the masks, U_R=0. Four separate endpoint matches from the signed product rule are sufficient but not necessary. For two-sided continuations, after neutral padding at both ends stability is equivalent to the two weighted equalities U_L=U_R=0.
- For a fixed signed gap goal d, choose parities p,q and set l=2i+p, r=2j+q. The decoder is the exact lattice line 2U_L i+2U_R j=d−D−(U_L+B_L)p−(U_R+B_R)q with i,j≥0. If the two slopes are opposite and nonzero, any one nonnegative integer solution extends along a positive integer lattice direction to infinitely many. This restores the exact charge without requiring either individual drift to vanish. No assertion about its absolute midpoint follows.

## 2. Why one asymmetric triangle tail cannot be reset

For the original triangle with a transitive exterior and arbitrary port incidences there, the weighted start at a singleton before the triangle is 12 times 2 to the number of later exterior singletons. After the triangle the multiplier is 3. Explanation: a split cover starting at that singleton can include each later singleton or leave it to the complement. If the triangle is later, summing over its two-way assignment gives 2·3+6=12; if it is earlier, it belongs wholly to the complement and contributes 3. At the triangle vertices all weights coincide by rotation, and the two port contacts inside this triangle differ by an equal-sized exchange. Thus, in the contact digit convention of README, if Q_l,Q_r are the signed binary values (difference of the port bits at each position, least significant last) on the two exterior tapes, with |right tape|=b,

    U_R = 12·2^b Q_l + 3 Q_r.

Since |Q_r|≤2^b−1, the two summands cannot cancel unless Q_l=Q_r=0. In one signed binary word with digits −1,0,1 at distinct powers of two, the largest nonzero position outweighs all lower positions together. Consequently U_R=0 iff **every** exterior contact is shared. This alone replaces the small-order “no asymmetric stable tail” evidence with an all-order explanation. One side with a fixed asymmetric seed cannot support infinitely many independently adjustable exact unit suffixes.

For a changed nonport core, those weights need no longer be binary place values. There is a literal order-six core with masks 46,23 and rows (40,1,3,6,15,30). Its ordinary H is 19 and its weighted starts are (21,20,29,35,73,59). Both masks read 143. Their four individual menu sums do not all match: they give (11,17,170,132) and (11,15,160,132) in the previous note's order. Its switch realizes (509,511), and a single neutral sink stabilizes this pair under every subsequent right shared continuation. Thus changing the quotient changes the arithmetic equality that matters, not just the number of tries. Adding k≥1 neutral sinks and then one right return gives counts 2610·2^(k−1)+345 and that value +2. The sink block can occur in at most two runs (leaving it for an old vertex requires the unique return), so its counts are affine in 2^k, determined by the first two literal values.

## 3. One fixed core pays both signs, but still loses absolute coverage

Use the original triangle with prefix 21 and no suffix. Its nonport rows are (30,28,8,16,4), masks I_x=9,I_a=18, and the state (D,U_L,U_R,B_L,B_R) is (−2,3,−12,3,0). The original pair is (69,67). Formula (1) yields exact gap −2 precisely at

    (l,r) = (4t,t) or (4t−1,t),  t≥0 and l≥0.

In particular t=1 yields centers 3634 or 3442, and t=2 yields 13182 or 20346; all four have been independently counted. A direction that cannot be varied freely in isolation can be varied without limit against the other side. This is an integer compensation rule, not a modular approximation.

The first branch l=4t replicates the right address on the left with two low neutral bits. This synchronizes not only the drift but also the forced three-run growth. If t has at least two returns and its last two return positions are i,i+r in chronological order (i,r≥1), the left positions are i+2,i+r+2. Split a Hamiltonian path of the fixed nonport core into three nonempty pieces. The [physical two-sided trace](../signed_ordinal_unit_centers_20261117/TWO_SIDED.md) now gives ≥3^(2i)·2^(2r−2) paths (one independent three-run filling on each side). The number of t for fixed i,r is 2^(i−1). Summing the geometric series over i,r gives at most O(X^log_9(2)) distinct centers at most X on this branch. With only one return, the right tape alone provides ≥2^(k−1) paths for a return at position k, since earlier vertices can choose two runs around the return and port pair; hence only O(log X) such addresses. The sharpened exponent is for l=4t only. The other branch and the whole fixed template still obey the more general O(X^log_3(2) log(2X)) two-sided bound, which does not require the gap to stay two while taking that upper bound.

## 4. What carries over from the archive

| Earlier mechanism | Exact connection | Where the transfer stops |
|---|---|---|
| [Curried births](../curried_births_tooth_diagram_20261117/README.md) | Define a state by equality under all physical reads; the five probes above separate gap states just as the rooted probe words separate continuation states. | A gap state deliberately forgets the common unsigned center. More gap queries cannot recover it. |
| [Centered Christoffel endpoint clocks](../../../01-canon/theorems/THM-778-centered-christoffel-endpoint-skew-product.md) | A reduced rational slope is not the original address: the clocks retain common scale and tie/parity. Equation (1) similarly needs the two parity charges in addition to the lattice slope. | A scheduled pair of directions is an address, not a guarantee of a desired output integer. |
| [Signed sibling ladders](../collatz_complement_ladders_20261007c.md) and [marked child compression](../collatz_child_compression_synthesis_20261007d.md) | Match ordered carries on their actual guarded phase; a longer physical move can narrow, not widen, the set on which the inverse works. The synchronized branch pays a second growth cost for the same free bits. | An identity cancelling drift does not supply a smaller ready source for induction. |
| [The signed Lonely Runner final rung](../../../07-reflections/the-final-rung-is-signed-not-absolute-kps-S124.md) | Taking absolute values would discard cancellation. Here |3l| and |12r| are unbounded on exact unit moves. | This comparison concerns how to retain signs, not a reduction between the problems. |
| [Prime-power lift clocks](../unit_lift_clock_20261015a.md) | States of the same reduced order can have different exact lift delay. Analogously equal current gap does not guarantee equal future behavior, and equal full gap behavior does not fix center. | A larger clock guarantees neither an exact terminal address nor density. |
| [Signed ordinal product](../signed_ordinal_unit_centers_20261117/README.md) | Signed Bézout cancels different tiles; (1) instead cancels two *addresses* on the same tile. Both are linear in the difference coordinate, with different mechanisms for weight selection. | A fixed finite bank of cores and ordinary shared front/back words remains sparse. A universal argument has to change the physical source or admit a richer family of contacts. |

An exact target selector needs a lift of **two** conditions: first choose a physical source with legal mask weights (or a payable two-address line), then prove a bound/recursion for the **absolute center on that line**. The present work solves the first condition explicitly and proves why one profitable-looking way of iterating it still fails the second. A concrete next hypothesis is a reusable neutral rewrite that changes the nonport graph and its center while retaining the two weighted equalities, rather than copying one side's bits to repair the other. Chan–Pak remains open.

## Reproduction

Run two_address.py for 2,160 literal gap comparisons over independent subset DP; the certificate for the sparse branch is the explicit path trace and geometric sum above. Run model.py for the previous finite triangle tape census. All numerical examples in this note can be reconstructed using tri, attach, host, profile and count in two_address.py.
