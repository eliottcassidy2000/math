# Where the physical two-tape carrier folds

Use the fixed carrier of [two-tape steering](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md) and [local density](../chanpak_viability_shape_20261115/LOCAL_DENSITY.md):
```
F(x,y)=2855+5760xy+4034x+3914y+45[p(x)+p(y)],
p(2u)=p(u)+u+(u mod 2),   p(2u+1)=p(u),   p(0)=0.
```
The earlier local-density result is about presentations and leaves open exact coverage. A way to see the lost information is to view a row as an integer coordinate chart: even with a physical carrier and fixed x, projecting onto the count coordinate eventually folds. We can locate the sharp lower edge of these folds exactly.

**Theorem.** For each x≥0 there are two addresses y and y+45 with the same count: take y=16384x+11158. Conversely, put A=5760x+3914. For x≥1 and any address collision y<z in row x,
```
87y ≥ 128A−130.
```
This is sharp, with equality for x=2+29t, t≥0, at y=128(177+1920t)+50, z=y+45. Thus the minimum possible onset ratio along unbounded rows is 128/87 relative to A, or 737280/87 relative to x. These collisions concern the map on addresses; no claim about physical tournament nonisomorphism or complement of the attainable spectrum is made.

**Proof.** We first record the exact digit lemma. For q=2^k and 0≤b<q, write r_k(b) for bit reversal using exactly k bits. There is an ε_b∈{0,1}, equal to one iff the most significant bit of b is zero, such that
```
p(qu+b)=p(u)+(q−1−r_k(b))u+p(b)+ε_b(u mod 2).
```
The slope in a fixed low word is the sum of place values reflected from all its zero bits. Only the highest low bit sees the parity of u; the other parity terms depend only on b and are included in p(b). This also follows by iterating the defining recurrence. Induction gives 0≤p(y)≤y: in the even step only an odd parent u needs slack, and then p(u)=p(floor(u/2))≤u−1.

The collision equation is A(z−y)=45[p(y)−p(z)]. Since A is a unit modulo 45, z−y=45j, j≥1. Then y≥jA. If j≥2 this already proves the claimed bound. Consider j=1, so p(y)−p(y+45)=A. Write y=64u+b. First suppose b+45<64 or u is even. Then adding 45 causes no carry past the low seven bits. Write y=128v+b', 0≤b'≤82. The difference is
```
c v+V_e,   c=r_7(b'+45)−r_7(b'),   e=v mod 2,
V_e=p(128e+b')−p(128e+b'+45)−ce.
```
The finite 7-bit calculation, certified by replay_row_folds.py, gives:
- c=90 is impossible: all such V_e are different from 44 modulo 90, whereas A≡44 modulo 90.
- c=87 leaves only (b',e,V_e)=(32,0,8),(34,1,23),(48,0,20),(50,1,35) after the modulo-3 condition. These give 87y−128A=1760,14,1616,−130, respectively.
- All other c are at most 78, with V_e≤92. If 87y<128A−130 then v<A/87 and c v+V_e≤78v+92<A, as A≥3914.

It remains to consider b≥19 and u odd. Adding 45 switches the high 6-bit prefix from u to u+1. For odd u=2w+1,
p(u)−p(u+1)=p(w)−p(w+1)−w−1−((w+1) mod 2)≤−1.
The maximum of r_6(b−19)−r_6(b) over b≥19 is 43. Applying the digit lemma to the two low words gives p(y)−p(y+45)≤43u+63, using p(b)+ε_b≤64, and the negative prefix contribution just established. So from 87y<128A−130 we would deduce u<2A/87 and p(y)−p(y+45)<86A/87+63<A, because x≥1 gives A≥9674. This proves the lower bound. For b'=50,e=1, the parameter v=177+1920t with x=2+29t makes 87v+35=A, giving equality.

For the collision in every row, use the seven-bit suffixes 22 and 67 with odd v. Their surcharge difference is 45v−1, exactly A at v=128x+87. This gives the displayed y. ∎

## Target-first consequences and a construction requirement

Once x and N are given, only one residue y mod 45 need be tested. On that residue, the portion y<(128A−130)/87 is injective: it admits at most one target address. But no selector that globally treats the row as a graph of its count can retain both physical presentations. The fold is not a numerical accident: the leading slope competes with the reversed-bit slope on a common high prefix, and the first available winning low word is dictated by phase. The raw largest difference, c=90, is arithmetically forbidden; the first physical contact is c=87.

A promising positive use of a fold is to retain both its addresses when a further attachment reads an extra interface. If the next reader changes the carriage coefficient from A to A+δ but retains the surcharge coefficient 45 and the same baseline across the pair, their output difference is exactly 45δ. This follows by subtracting the collision equation and is exact at every height. It is a conditional compiler specification, not a tournament construction of that reader. A proof of placement would still need a physically verified paired reader, a controlled baseline, and enough complementary phase classes; merely increasing the local multiplicity or replacing 45δ by its residue is insufficient.

The required change of reader is physical at the difference level. Keep the same 8-vertex core and both guards of the original carrier. The one-apex two-tape transfer has the same Q and h for all painted masks, since these are ordinary unmarked-spine path and two-slot counts. For each mask m, its constant B_m and its y coefficient S_m are the two literal Hamiltonian counts at (x,y)=(0,0),(0,1), taking the difference for S_m (as p(1)=0). Subset-endpoint DP gives:
```
mask     B_m   S_m    S_m−3914   separation on the old pair
 35     2855  3914         0                0
152     3485  3896       -18             -810
 38     3141  3934        20              900
```
This applies to every old pair separated by 45 in row x. The signs and step sizes are literal, with gcd 90. Every mask makes an even coefficient difference, since it is a difference of two odd Hamiltonian counts, so no fixed-spine mask reader of this type gives a finer common gcd. The guarded spine also stays strong after changing the mask.

There is also a physical reader with no mask change: read F(y,x) in place of F(x,y). The surcharge is symmetric, and the coefficient shift is T−S=120, so the old equal pair now separates by 5400. Holding the mask fixed but moving the left address from x to x+d instead gives a separation of 259200 d. These comparisons preserve the phase of the address difference, not the absolute height: in particular they do not supply independent composition of the -810 and 900 steps. The mask differences are two actual calibrations for a future changing-reader return; the outstanding combination and absolute placement are not claimed.

In fact this repaint restriction has a large uniform absolute-height cost. On all 256 masks the same three-count DP gives 2131≤B_m≤4671 and 2880≤T_m≤5670 (T_m reads the left unit tape). Thus, for any repaint that separates the old pair,
```
|F_m(x,y)−F_35(x,y)| ≥ 2y−1636x−1816 ≥ 9884x+6012,
```
using y≥A from the basic collision argument, valid on every row including x=0. The same bound holds at the larger address. The 90 common divisor of the small differences therefore cannot be used as a 90-scale repair of the old absolute target merely by selecting another mask: even the nearest nontrivial repaint displaces the chart by an amount linear in x. Recenter with another changing coordinate before asking for a composition law.
