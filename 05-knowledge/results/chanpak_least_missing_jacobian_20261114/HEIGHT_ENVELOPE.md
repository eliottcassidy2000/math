# Exact height bounds for a transitive source with one apex

Normalize the painted transitive spine by deleting neutral leading
zeros and trailing ones. At active length n>=2 the word begins
1 and ends 0. This is exactly the strong part of the chart used
in [modular escape](MODULAR_ESCAPE.md). Its count is
$$
  H=1+\sum_{i<j\,:,w_i=1,,w_j=0}
                \bigl(2^{j-i-1}+[j=i+1]\bigr).
$$
To see it directly, each odd cycle through the apex has fixed
initial and final spine vertices and an even subset in between.
Summing its packing weight two gives the displayed coefficient,
with the adjacent case handled by the empty interior. All cycles
meet the apex, so none can be packed together.

## Maximum and minimum at each active width

A swap 01→10 of an adjacent pair strictly increases H. It adds
two from the selected pair. An old 1 to its left reads a more
distant zero after the swap; a later zero reads a more distant
1. The distance kernel is nondecreasing, including its first
two equal values, so no other contribution falls. Consequently
every maximizer is a sorted word 1^a0^b with a+b=n and both
positive. Summing the kernel, its count is
2+(2^a−1)(2^b−1). Thus
$$
  H_{\max}(n)=2^n-2^{\lfloor n/2\rfloor}-2^{\lceil n/2\rceil}+3.
$$
One needs the two block lengths as equal as possible.

The minimum works the other way. Put r=n−2. First count the
extreme-to-extreme term (2^r for r>=1). For each interior
position i=1,...,r, it contributes a forced edge involving
one extreme. Its cost is a_i if it is zero and a_(r+1−i)
if it is one, where a_i=2^(i−1)+[i=1]. These choices
are independent lower bounds. Since a_i is nondecreasing,
they are all attained by choosing cheap zeroes before
cheap ones, with a monotone tie choice. There are then
no internal inversions. Thus
$$
 H_{\min}(n)
 =1+2^r+\sum_{i=1}^r\min(a_i,a_{r+1-i})\quad(r\ge1),
$$
with n=2 giving three. An interior ordering is minimal iff
its choices meet these local lower bounds and have no
interior inversion. More explicitly, for r≥2 the additive
sum is 2^(t+1) when r=2t, and 3·2^t when r=2t+1.

This shows an exact exchange of geometries: to maximize, put
all the ones before the zeros; to minimize, keep just the forced
outer inversion and sort the interior in the reverse way.
At width 4 the maximum is 11 and at width 5 the minimum is
15. Every larger width has minimum greater than 13 by the
extreme pair alone. Thus 13 is literally between two width
layers. The local compiler sees neither that layer gap nor
its ordering of configurations.

Run `python height_replay.py` to check the two envelope formulas
against every active painted word through width 14.
