# A variable-size physical radix from the same sampled carry idea

The [base-128 tiles](SPIGOT.md) can stretch independently, not just repeat at one geometric rate. After either digit e=0 or 1, append r≥0 transitive two-vertex pairs internally, each with the first vertex tagged and the second untagged. The resulting tile has h=5, p=128·4^r, z=(p−5)/3 and e_full=0. To see this, the extra pair has h=1,P₂=4, weighted marked start 2, weighted marked end 1, and K=1 (since its joint U is zero). The ordinal transfer gives w'=4w+10, z'=z+p, p'=4p, and K'=K−w. Hence the local quantity w−3(h−K) begins at 10 and increases by 10 per extra pair, for either digit.

For a full word of m stretched digits let s be the sum of the r's. The same initial marked prefix gives J=5 at depth zero and J'=5J. Put V=3W+2J. At each digit:

    V' = 4^r (128V + 72(3e−2)J).

The length-dependent carry is exact: 3K+W=5^m(27+10m+10s). Therefore the same strong two-vertex return reads

    3N = 11W + 5^m(27+10m+10s).

This is a physical identity, not an arithmetic rational-approximation code. The outer root still reaches the tagged leading singleton regardless of all added stretches, and all vertices return via the last block.

## Absolute inverse by valuating the target

Set d=7m+2s, the combined geometric cost. Every count in this layer satisfies

    594·2^d ≤ N < 625·2^d.

Indeed W starts at 162 and multiplies by p at each step before a positive correction. Also w<110·4^r, so at every depth 162≤W/2^d<162+550/123<167, by the geometric ratio 5/128 (stretches can only lower its contributions). The carry divided by 2^d is at most 27: increasing m or s raises its linear factor by 10 but divides by 128/5 or 4; explicitly the factor in m is 5/128 and that in s is 1/4. Thus the displayed layer bound holds and the prescribed N determines d uniquely.

For each m≤d/7 having s=(d−7m)/2 a nonnegative integer, calculate W from the exact read and form V. In the last update, J is odd and 128V_old is divisible by 128. The valuation of V' is 2r+4 for e=0 and 2r+3 for e=1. Thus the parity of the valuation selects the last physical digit, and its value selects the last stretch. Divide by 4^r, subtract its known offset, divide exactly by 128, and repeat, checking total stretch and the initial V=496. This yields an O(log²N) arithmetic inverse and selects the actual sizes from the target, not just the residues. The depth-one stretched word (0,1) gives 306883.

This extension is still provably sparse. At exponent 3/8 the sum of all digit costs to the negative exponent is

    Σ_(r≥0) 2(128·4^r)^(−3/8)
      = 2 / (128^(3/8)(1−4^(−3/8))) < 1.

The strict inequality follows already from 128^(3/8)>6 and 4^(−3/8)<3/5. The geometric sum over all words then implies O(X^(3/8)) possible counts up to X. This is a bound on the number of outputs, not on distance to the nearest output. Variable size makes the inverse substantially richer, but a universal selector must escape this stretched alphabet too.

Replay: `python replay_stretch.py` checks literal local profiles and the valuation inverse.
