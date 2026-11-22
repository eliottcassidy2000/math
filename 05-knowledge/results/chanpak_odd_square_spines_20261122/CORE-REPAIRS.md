# The two near witnesses are isolated under all 4+6+6 core replacements

**FINITE-EXACT.** Fix either exterior address pair in [the order-118 bracket](NEAR-WITNESSES.md), with $E=2^{60}(2^{61}-1)$ and $N=E-1$. Replace the three interior blocks by arbitrary tournaments of orders four, six and six, in their given ordinal order. Assign the apex an arbitrary incidence at every interior vertex. Among *all* resulting tournaments, the only ordinary Hamilton path counts within $10^{10}$ of $N$ are the two counts of the original marked arc for that address pair.

| Fixed exterior pair | Complete set of $H-N$ in the band |
|---|---|
| $(735435556549430,1575929833923500)$ | $-844710,-844708$ |
| $(4814468855604222,240731609083994)$ | $184814588,184814590$ |

In particular no single new tournament in this universe hits $E-1$ or $E+1$. The unit-flip search allows any of the fifteen middle arcs to be reversed. It leaves the exterior words themselves and core size changes as open ways to correct the residual.

## Why the product of slot counts can be frozen first

Let $q$ count ways to assign the vertices to two ordered path slots, both possibly empty. The inherited tagged chart is

$$
 H=qxy+tx+cy+h(f(x)+f(y))+r[x\text{ odd}]+e[y\text{ odd}]+k.
$$

The ordinary block count $h$ is at most $q$. The nonnegative endpoint counts $c,t$ are each at most $q$ by restriction of the assignments; the two exterior parity corrections satisfy $|e|,|r|\le q$. The constant $k$ is an ordinary count on the 17-vertex interior with apex, so $k\le17!$. Across an ordinal join the two-slot count multiplies, regardless of the tags. The maximum local $q$ among the [four types at order four](../../../session212/cycle-levels/classes4.rows) is 32, and among the [56 types at order six](../../../session212/cycle-levels/classes6.rows) it is 432. Therefore every chart in this search satisfies

$$
 |H-qxy|\le17!+2(32\cdot432^2)(x+y+1).
$$

Both witness pairs have $xy$ larger than twice the displayed bound by more than $10^{29}+|H_{\rm old}-N|+2$, where $H_{\rm old}$ is the lower reading. Changing the leading integral $q$ thus moves an endpoint more than $10^{29}$ away from $N$. Both endpoints of any competing near unit flip must have the original $q_*=2293760$.

All 12,544 ordered triples of isomorphism types yield only sixteen with $q=q_*$. The possible three factors in their original order are $(32,224,320)$ and $(32,320,224)$. All $2^{16}$ apex tag masks on each of these sixteen were used to compute the table. Because every relabeling transports both the apex tags and any distinguished edge, this also covers every labeling. The single-tournament test includes 1,048,576 surviving tagged presentations. These counts of presentations need not be counts of isomorphism classes after adding the apex.

## Competing arcs and the missing constant correction

For each original word pair, retagging the *original* core in all 65,536 ways while reversing the original middle arc gives 64 tagged unit flips, with 64 distinct midpoint heights. Across replacement types, allow every middle edge, retaining the types with $q=q_*$ in both orientations. This leaves 30 triple/arc choices out of 188,160. Each of their 65,536 tagged presentations was evaluated at each word pair. Of these 1,966,080 presentations per pair, 320 have an actual height difference of two, yielding 160 distinct midpoint heights. These counts are restricted to the surviving leading coefficient. A flip can be tested even when it preserves a difference of two only at the specified pair of words.

| Fixed pair | Nearest midpoint minus $E$ | Next distinct midpoint minus $E$, ordered by absolute distance |
|---|---:|---:|
| First | −844,710 | +10,420,831,809,831,200,142 |
| Second | +184,814,588 | +6,135,804,521,871,654,884 |

To test a pure additive repair, retain a unit flip whose lower orientation preserves the three growing coefficients $(h,c,t)=(2465,1344270,712944)$ as well as $q_*$. Both physical words are even, so parity terms vanish. The only attained additive shift of its lower endpoint, compared to the original lower endpoint, is zero at either address pair. Recognition of the leading product alone did not imply this uniqueness; it follows from the actual tagged middle and endpoint readings.

Run the following replay from this directory. It reads the two atlas files, derives the local slot counts directly from subset paths, and constructs every retained tagged chart from its blocks.

~~~sh
python3 -B core_repairs.py
~~~
