# Exact ternary residue counts in the two-tape repair alphabet

Let L be the signed loss of a matched same-core chart, with
width 12k, in the [all-width transport census](../chanpak_viability_shape_20261115/ALL_WIDTH_TRANSPORTS.md).
The half-step L/2 has polynomial
```
D_k(z) = z^[8(k−2)] (1+z^(-16))^(2k−1)
                    (1+z^25)^(2k−1) (1+z^(-34))^(k−1).
```
At a primitive cube root omega, 1+omega=-omega² and
1+omega²=-omega. Hence D_k(omega)=(-1)^(k−1)omega.
By the three-term roots-of-unity filter the exact number in
half-step class r mod 3 is
```
(4*32^(k−1) + (-1)^(k−1)(3 [r=1]−1))/3.
```
It differs from exact uniformity by an O(1) boundary term.

For primitive charts, the first polynomial is D_1 and the period
factor is the full period factor minus D_1, as in the free-monoid
proof. At omega, the two factors are D_1(omega)=omega and
−1−omega=omega² (the all-chart period evaluates to −1).
Thus primitive charts have exactly
```
(4*28^(k−1) + 3 [r=2k−1 mod 3]−1)/3
```
in class r. In each width these three counts differ by at most
one, including at width 12. In particular imposing the
endpoint-congruence restriction on **every** physical repair
would wrongly exclude two entire lanes of this literal alphabet.
This is not an asymptotic equidistribution claim or a claim that
every absolute count is reached.

The distinction is which information is retained. The endpoint
reflection forces a multiple of three half-steps when a rooted
start and the relevant signed data stay fixed and only its
paired end is read. A suffix transport trades two address charts
across the core; it is not that one-coordinate rooted move.
Their exact residue formulas show why a proof needs a menu of
different readouts rather than a universal scalar divisor.
