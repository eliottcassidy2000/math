# Disjoint physical repairs: a sparse absolute-count packing

The [determinant family](DETERMINANTS.md) has many distinct determinant values, not necessarily distinct counts. A restricted subfamily, however, gives genuinely disjoint +2 count intervals, not merely a packing in the determinant projection.

Use only internal nibbles ε∈{0,1} and allow every fourth nibble to vary, starting from the high end of the variable nibbles. All intervening nibbles are fixed at zero. The sign of the mutual-arc edit stays +2. We prove that the absolute count in either orientation strictly increases in lexicographic order on the free nibbles, with gaps larger than 2 between successive lower anchors.

## A quantitative endpoint estimate

Let ρ=2/3. Number positions by i=0,…,n−1. Fix both first sites to 3 and both last ones to 0, as in DETERMINANTS.md. For a single early flip 0→3 at index i≥2, any gained path can be compared via its assignment to the three increasing runs and its apex order. To produce guaranteed gained paths, choose apex order A then B, fix site 0 as run 2's start, site i as run 1's start, and the two final sites as runs 0 and 1's ends. Of the earlier sites other than 0, any may go in runs 0 or 2; any free later site may go in any run. These assignments all fail before the flip and succeed after it, giving at least
```
2^(i−1)·3^(n−i−3) = 3^n ρ^i/54
```
gains in either mutual-arc orientation, since the middle run is nonempty. A lost path must have used the flipped site as the end of its run 0 or 1. There are at most 4·3^(n−1)ρ^(n−1−i) such assignments, including both apex orders. If n−1−2i≥13, losses are at most half the guaranteed gains (144·ρ^13<1). Thus the count increases by at least 3^n ρ^i/108. By taking the converse and reversing the spine, the same lower bound holds for a late flip 3→0 at distance i from the right end.

For an arbitrary diagonal flip at end-distance h=min(i,n−1−i), only first-of-run or last-of-run status can change validity. By allowing 2 first runs and 2 last runs per apex order, its absolute effect is at most
```
4·3^(n−1)(ρ^i+ρ^(n−1−i)) ≤ (8/3)·3^n ρ^h.
```
Both estimates hold uniformly in all unfixed letters.

## Separated nibble comparison

In two different parameter words, keep higher agreeing nibbles fixed and compare the highest differing free nibble. A change ε=0→1 adds 0101 in T and 0100 in S without carries. Thus there are two early 0→3 flips at distances i,i+2, and one late 3→0 at right-distance i. All variable positions have n−1−2i≥13 in their respective end coordinates (even the closest nibble to the central patch has i≤m−5 and n=2m+7). The total increase attributable to this nibble is at least
```
(2+ρ²)·3^n ρ^i/108.
```
The same-type positions for subsequent free nibbles are deeper in by multiples of 16, giving a total possible opposing effect at most
```
(2+ρ²)·3^n ρ^i · (8/3) ρ^16/(1−ρ^16).
```
The last factor is less than 1/216. Thus the net increase is at least (2+ρ²)·3^n ρ^i/216, greater than 2 uniformly for sufficiently large j since i≤m−5. The +2 intervals really are disjoint. There are 2^ceil((j−1)/4) free words, distinguishing as many determinants by the previous theorem. Using n=8j+21 gives Ω(X^((log_3 2)/32)) disjoint actual +2 intervals with endpoints ≤X. This is an absolute packing theorem, not an absolute covering theorem.
