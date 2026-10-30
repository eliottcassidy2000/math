# The same C5², two different action carriers

**FINITE-EXACT action classification, with exact consequences.** The abstract group alone is not the carrier.

## Recover the actual24-bit representation

In the modulus25 ramified bank, the Sylow5 subgroup of the effective normalizer E_200=C20×C10 is C5². It has commuting representatives:
- multiply by6;
- send δ to6δ.

In coordinates x+yμ these maps are (x,y)↦(6x,6y) and (x,y)↦(x+10y,6y), modulo25. Both have order5. The24 bank bits are antipodal pairs of actual g-orbits.

[action_carriers.js](action_carriers.js) constructs all those orbits and all25 signed coordinate actions. Their joint coordinate orbits have sizes
1,1,1,1,5,5,5,5.
Two five-point orbits have kernel {(0,b)} and two have kernel {(a,−a)} in the displayed generator coordinates. The two kernels are distinct lines.

The action can be made unsigned by choosing a fixed bit assignment as origin. For an affine action of an odd group on a binary vector space, sum all translates of zero: the resulting vector is fixed, since the group order is odd. The checker constructs and verifies this gauge. Hamming weight below means distance from that fixed origin, not weight in an arbitrary original orientation gauge.

Thus the representation consists of four fixed bits, two parallel regular C5 copies for one independent quotient coordinate, and two parallel copies for the other. Its weighted orbit polynomial is exactly
F(t)=(1+t)^4 [((1+t)^10+4(1+t^5)^2)/5]^2.
In particular
F(1)=16·208²=692224.

Equivalently, the fixed-mask histogram is
identity:2^24;
eight elements:2^16;
sixteen elements:2^8.
Burnside gives the same692224 orbits. These are bank-mask orbits under the named subgroup, not graph classes.

## Compare with the old25-point regular action

The very same abstract group C5² acts regularly on the25 gaps of the old50-vertex tournament. Its polynomial is
P(t)=[(1+t)^25+24(1+t^5)^5]/25.

Their cancellation depths differ:
| action | coordinates | rational | characteristic5 | characteristic2 |
|---|---:|---:|---:|---:|
| regular action on old gaps |25|5|21|25|
| actual normalizer action on bank bits |24|8|16|24|

The second row follows either by exact expansion or by applying the general cancellation theorem to each10-coordinate factor, then adding the four fixed coordinates. Same group, different stabilizer kernels, different cancellation law.

This explains why retaining subgroup names and common quotient coordinates is necessary. Counting the abstract symmetry group, or even its order25, cannot recover the action.

## A shared four-phase factor

The clone-doubling recurrence has its exceptional phase when2^j=1 in F5*. The ramified Fibonacci unit reduces to μ=1/2=3 moduloδ. Since3=2^(−1) in F5*, inversion identifies these two four-state phase quotients:
2^j ↔ μ^(−j) moduloδ.

This is an explicit common quotient clock, stronger than merely observing two periods of length4. It does not identify the full clone-dilation dynamics with Fibonacci multiplication: their state spaces, higher congruence information and invariants differ.
