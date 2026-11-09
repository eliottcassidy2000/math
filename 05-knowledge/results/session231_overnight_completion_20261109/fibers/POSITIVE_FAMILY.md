# A fixed-step unbounded physical fiber from a movable core

Status: independently audited PASS for the all-order construction; exact certificate verified. Chan–Pak completeness and consecutive-odd coverage remain OPEN.

## Theorem

For every k≥0 there is one physical rooted full interface with at least k+1 distinct Hamilton counts forming an arithmetic progression with FIXED difference190. All witnesses are strong tournaments on312k+12 vertices. The deletion Hamilton count stays95; the deletion two-path count is1276·2^(312k+3), coprime to95. Thus spacing does not grow as a power of the number of choices, unlike independent repeated H95 modules.

Let C=[34,20,9,35,13,86,31], cut17. Its seven-coordinate state is

    (h,q,A,G,T,S,E)=(95,1276,270,372,40,20,30).

Let X,Y be the length312 bit words in `LOOP_WITNESS.json`. Both X's last bit and Y's first bit are0. The k+1 deletion objects are

    L_(k,j) = s1 ⇒ s0 ⇒ X^j ⇒ C ⇒ Y^(k−j) ⇒ s0,  0≤j≤k.

Here bit words denote ordinal singleton blocks with those cut bits. Apply star closure L→v→r, with root r beating precisely the cut in L. The initial selected singleton ensures strongness. The total vertex order is312k+12.

Writing b=2^312, R=value(Y), V=value(reverse(Y)), the common full interface is

    h=95,
    q=1276·2^(312k+3),
    A=2822·2^(312k+1)+190R(b^k−1)/(b−1),
    G=1583+5104V(b^k−1)/(b−1),
    S=95, E=0.

The exact scalar is D_(k,j)=D_(k,0)+190j. Baseline D_(k,0) is specified without ambiguity by the displayed physical construction or the exact seven-coordinate recurrence; no closed-form baseline is required for the progression theorem.

## Proof using only one finite identity

Literal permutation enumeration verifies C's state. Repeated use of the inherited ordinal composition law verifies

    s0 ⇒ C ⇒ Y ⇒ s0,
    s0 ⇒ X ⇒ C ⇒ s0

have identical (h,q,A,G,S,E), and replacing the first by the second changes T by−190. This is an exact length312 word identity, not an extrapolation from testing several k.

For the unflanked comparison C⇒Y versus X⇒C, h,q,A,G agree. Under arbitrary nonempty singleton left/right contexts, its T difference is

    ΔT + last_bit(left)·ΔS + first_bit(right)·ΔE.

This follows by applying the composition law twice. For every adjacent pair of positions j,j+1 in L_(k,j), the left context ends0 and the right context starts0, because X ends0, Y starts0, and the external adjacent singletons are0. Therefore each position move has precisely the same T difference−190 as the verified s0-flanked identity. Induction gives T_(k,j)=T_(k,0)−190j. Star closure has count2A+h−T, proving the scalar progression.

The common A,G formulas follow by computing j=0 and using equality after each move. The leading word10 has A2,G1; the trailing wordY^k0 has A2R(b^k−1)/(b−1), G=V(b^k−1)/(b−1). The two full endpoint counts are fixed by the first and last singleton. All divisions are geometric sums of integers.

## How the physical carrier was recovered

The original B95 topology's movable carry graph contains a useful2232-state weighted SCC with nonzero gains, but none of its own cut interfaces lies in that SCC. The complete same-(h,q) query across all classes through order9 found FIVE order7 topologies with h95,q1276. Two other topologies have four cut interfaces inside that SCC. C,cut17 is one of them. Thus the original absence was bank-relative, not a physical impossibility. This recovery of an omitted physical sidecar is essential to the theorem.

`PHYSICAL_BRIDGE.json` lists all matching classes and cuts. `physical_cycle.py` extracts a genuine directed closed walk through C's interface. The theorem's standalone verifier depends only on `LOOP_WITNESS.json`, the composition law, and literal core counting, not on completeness of the census or the discovery SCC implementation.

## Cycle-gain quantization theorem

For any fixed odd h and even q, a moving-core edge with bits x,y sends

    a'=2a+hy−qx,  g'=(g+qy−hx)/2,  x=g mod2.

Augment the state with p=previous y, and let x'=g' mod2. The exact interior scalar cocycle is

    w=a y+h p y−x g'−h x x'.

It is obtained by assigning the formal intermediate core endpoints S_i=h x_i, E_i=h y_(i−1) and comparing C_i⇒s_y with s_x⇒C_(i+1). Those intermediate states are bookkeeping, not assumed physical; boundary corrections restore the real S,E at a physical start/end. Direct multiplication gives

    a'g'−ag ≡ q w (mod h).

Around any closed walk, h divides q times the total gain, so h/gcd(h,q) divides the gain. Also a≡p (mod2) at every cyclic state, because a'≡y. Thus w≡y(a+p)−x(g'+x')≡0 (mod2). Since h is odd, every closed-walk gain is divisible by

    2h/gcd(h,q).

For this unsaturated carrier the minimum possible nonzero magnitude is190; the displayed loop attains it. Hence a unit gap-two drift is impossible in any unsaturated fixed-h carrier with h>1. The theorem does not prohibit other overlap mechanisms, changing h,q, or saturated carriers with h|q.

## Controls, limits, and next question

`verify_family.py` independently enumerates the seven-vertex core by literal permutations, checks every edge of the closed carry word and its modular identities, verifies the finite flanked word identity, and checks the explicit full-interface formulas and all positions through k10. The proof above supplies all k. Normal and optimized outputs agree.

This is an unbounded arithmetic progression, NOT an unbounded consecutive-odd interval. Width grows linearly in vertex order but only logarithmically in scalar magnitude, since q grows exponentially. The constant190 still leaves94 missing odd residue positions between adjacent displayed values. The next positive target is a saturated physical carrier with h|q and nonzero gain2, or a genuine interpolation mechanism filling the190-step gaps while preserving the interface. The complete endpoint-erased census through core order9 still has no run longer than3 (unsaturated maximum2).
