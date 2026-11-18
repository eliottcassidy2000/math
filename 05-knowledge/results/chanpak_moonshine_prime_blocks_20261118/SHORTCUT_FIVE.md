# Five signed steps to the two central offset sheets

This extends the [three-step paired window](PAIRED_SHORTCUT_WINDOWS.md) using the same two physically different readout mechanisms (two-run wrapper and a transitive-spine apex).

For every m≥17 and every odd beta from −33 through 31, the count x=2^m+beta and each of its first five strictly shrinking signed shortcuts have literal tournament witnesses, **both for the ordinary values and for their septuples**. Each step chooses the sign that makes 3x±1 divisible by four and removes all its dyadic factors. The replay verifies 148 exact coefficient rows and 296 independent literal graph counts. These are re-encoded witnesses, not a claim of a direct orientation edge for each pair.

The arithmetic part can be read without the graph tables. The offset value sets at depths 0 through 5 have cardinalities 33,17,11,7,4,2 and end at sign(beta). At any later step with sufficient dyadic width w, the count on this sheet has form 3^j·2^w±1 and the next value is 3^(j+1)·2^(w−2)±1. At depth five, the total stripped exponent is at most 14. Thus m≥17 reaches the known positive ordinary two-parameter clone rail with w≥3 for every positive beta. It has further **ordinary** witnessed values for an unbounded number of steps as m grows, until w<3. It does not supply a fixed terminal root or upper septuples indefinitely.

One literal continuation on the minus sheet uses a different mechanism. If beta<0 and w=m−V_5−2 is even and at least 12, the sixth value factors as
```
3^6·2^w−1 = (27·2^(w/2)−1)(27·2^(w/2)+1).
```
The depth-three ±1 rows give both factors and a sevenfold witness of the minus factor at those widths. Ordered strong components multiply H, so joining the appropriate factors realizes the sixth value on both axes. More generally
```
X^(2^a)−1 = (X−1) product_(0≤i<a) (X^(2^i)+1),
X=27·2^u, u≥6.
```
The first factor has the paired minus reader; the other factors have the all-order ordinary positive clone construction (their ternary exponents are at least 3 and their widths at least 6). This gives paired counts along an unbounded family of ternary exponents on the minus sheet, but only for the specified divisible-width subfamilies. No prime claim is made for these products.

This supplements rather than replaces the [physical inverse atlas](INVERSE_METAGRAPH.md): matched endpoint counts, literal arc edits, and genuine merged-metagraph lines are three different claims.

Replay: `python3 paired_collatz33_depth5_replay.py paired_collatz33_depth5.tsv 5`.
