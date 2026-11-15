# Many exact checksums at one prescribed high quotient

The [paired checksum locks](ANCHOR_DECODER.md) change a host and second core simultaneously while fixing an exact high quotient. The following all-order construction gives a growing set of distinct checksums at the **same** quotient, with an inverse that decodes the actual core and host. It still does not give every checksum.

## Equal-Q core dictionary with two independently variable factors in h

Four choices suffice:
```
TT5:           (h,Q)=(1,32)
apex ABB:      (h,Q)=(5,32)
TT8:           (h,Q)=(1,256)
apex ABABB:    (h,Q)=(19,256)
```
The apex words are strong; the profile formulas are from [the weighted area](GAP_HEISENBERG.md), and both are checked by literal path DP below. Ordinal sums multiply h and Q. Fix an integer e≥5. Choose 2e copies from the first pair and e from the second, independently. If i and j are the cyclic-copy counts, the physical middle core W has
```
Q=2^(18e),  t=h(W)=5^i 19^j,  0≤i≤2e, 0≤j≤e.
```
The other core U is a singleton. This dictionary changes the ratio substantially where the fixed-equal-profile norm was too rigid. The coupling resource is *equal Q*, not equal h or equal Q/h.

Set δ=475^e, k=ceil(log2 δ), and take any g≥k+1+18e. Every t divides δ. For each dictionary choice use the binary anchor word 1 0^g v, where v is k bits, its reversed value being r=2^k−δ/t. Each r is a valid k-bit integer. Let b=val(v) and A(v) be the Hamilton count of the apex over v. The exact split identity gives
```
N = 2^g [Q·2^k + δ] + Q(b+1) + t A(v).
```
The leading quotient is identical for all dictionary choices because t(2^k−r)=δ. Their low parts are distinct: A(v)≤b+1≤2^k, so t A(v)≤δ2^k<2δ²<Q. The last inequality holds for e≥5 as 2(225625/262144)^5<1. Thus Q(b+1) occupies a different interval of length Q for each t, using injectivity of width-k bit reversal. Also the low parts are <Q(2^k+1)<2^g. This proves an uninterrupted increase in **the number of distinct available checksums** to (2e+1)(e+1) at one exact high quotient, not an uninterrupted interval of checksum values.

For a prescribed N with this high quotient, let ℓ=N mod 2^g. The only possible host suffix is b=floor(ℓ/Q)−1. Reject if it is not a k-bit value. Reverse its bits to obtain r, and test that t=δ/(2^k−r) is an integer of the specified form. Finally check ℓ mod Q=t A(v) exactly. When it passes, the exponent pair gives an actual W and the suffix gives an actual host. For this particular family there is no search over representations at the chosen e,g. This directly moves beyond a single low checksum per high position, but does not make the remaining low equality automatic.

Replay: `python replay_quotient.py` checks the two nontransitive profiles independently, then reconstructs all checksums for eight successive dictionaries and checks the full seam expression and the inverse. The all-order claim is the symbolic argument above.
