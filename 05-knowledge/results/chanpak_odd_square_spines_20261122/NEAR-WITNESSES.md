# Two order-118 marked arcs bracket the exponent-61 perfect number

**FINITE-EXACT.** Put $E=2^{60}(2^{61}-1)$. There are two canonical spine-four tournaments on 118 vertices whose marked arcs reverse from $E+d-1$ to $E+d+1$, with the following midpoint offsets.

| Left address $x$ | Right address $y$ | Midpoint minus $E$ |
|---:|---:|---:|
| 735435556549430 | 1575929833923500 | −844,710 |
| 4814468855604222 | 240731609083994 | +184,814,588 |

The lower tournament in the first row has exactly **2,658,455,991,569,831,744,654,692,615,952,997,465** Hamilton paths. Reversing its marked arc adds two. Both new offsets improve their respective sides of the earlier order-118 spine-four bracket $(−5,378,156,638,+464,575,890)$ and the order-123 one-receiver bracket $(−2,735,453,618,+380,266,414)$. The rows give two different marked arcs. An arc joining $E-1$ to $E+1$ remains open. [Earlier spine-four pair](../chanpak_decimal_kernel_20261122/FAREY-BOUNDARIES.md) · [Earlier single-receiver pair](../chanpak_e61_fibonacci_target_20261122/README.md#two-physical-edges-bracketing-the-target)

## Reconstruct the pairs from two chosen sums

The positive search used the [weighted sum inverse](SPINE-SEARCH.md) with independent low bits extracted first. The direction and weighted-sum displacement were selected by searching. For these two records both addresses are even, so write $x=2X$, $y=2Y$. With the inherited lower chart parameters, let $Q=4q$, $T=2t+h$, $C=2c+h$ and $K=k$. Then the exact chart is

$$
 K+QXY+TX+CY+h\big(f(X)+f(Y)+[X\text{ odd}]+[Y\text{ odd}]\big).
$$

Thus $S=aX+bY$ fixes the same concave quadratic and its centered lattice coordinate as in the weighted inverse. Let $N=E-1$. Use the following two integers to select a sum near tangency in each direction,

$$
 S=\left\lfloor\sqrt{4ab(N-K)/Q}\right\rfloor+s.
$$

| Weights $(a,b)$ | Displacement $s$ | Tolerance in the lower height |
|---|---:|---:|
| $(15,7)$ | 6,803 | $10^9$ |
| $(1,20)$ | 692,570 | $10^9$ |

For each entire sum slice the replay returns **exactly one** pair of these physical parities within tolerance. It enumerates both signs of the tangent coordinate with integer square roots. The digital terms lie between zero and $h(2+\lfloor S/\min(a,b)\rfloor)$, so enlarging both sides of the quadratic interval by the tolerance keeps all eligible pairs. The two full height computations are performed in separately constructed lower and upper charts with the original marked middle arc in opposite orientations. The order formula is $17+\operatorname{bitlength}(x)+\operatorname{bitlength}(y)$.

An independent recount builds the 118 outgoing adjacency rows of each of the four tournaments. It checks strong connectivity and verifies that the two rows of each differing arc toggle precisely the reciprocal entries. The marked vertex labels, starting at zero, are $(55,58)$ in the first row and $(58,61)$ in the second. The block subset dynamic program splits each Hamilton path at the apex, allocating every block vertex to its side of the apex, and recounts all four endpoints without the digital height formula. [Literal path counter](../chanpak_decimal_kernel_20261122/e15_pair.py)

The displacements are finite search data. The same construction accepts any odd target and candidate direction, while selection of a sum yielding zero displacement in height is still open. Run the standard-library replay from this directory.

~~~sh
python3 -B near_witnesses.py
~~~
