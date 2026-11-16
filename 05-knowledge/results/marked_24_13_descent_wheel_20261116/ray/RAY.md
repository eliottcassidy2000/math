# The full 24:13 marked ray is a finite carry language

Classify all spine words whose forward and reversed complement
binary address vectors lie on the common marked ray (13,24).
A letter z=2a+b has debt d=24a−13b. The two ray equations are

    Σ d_i 2^(n−1−i)=0,     Σ (11−d_i) 2^i=0.

Scan left to right with integer states (F,C), initial (0,0), and

    F' = 2F+d,     C' = (C+11−d)/2.

Reject noninteger C' and accept exactly on return to (0,0). The
forward carry has bound −24<F<13 because a remaining suffix
must cancel it and all suffix digits lie in [−13,24].
The second carry is the partial reversed-complement difference
divided by 2^j at width j, so any accepting word has integer
carries. That recurrence stays in [−13,24] by averaging.

## The full grammar

There are 66 reachable states. Trim states not reaching the origin:
58 remain. Compress chains with one outgoing edge. The complete
compressed transition table, replayable with check_ray.py, is

| Start | Tape label | End |
|---|---|---|
| (0,0) | 133220 | (−8,1) |
| (0,0) | 311103311002 | (0,0) |
| (−8,1) | 032220 | (0,0) |
| (−8,1) | 2101012022331 | (−11,11) |
| (−11,11) | 022230022331 | (−11,11) |
| (−11,11) | 200113 | (−3,10) |
| (−3,10) | 301113 | (−11,11) |
| (−3,10) | 1232321311002 | (0,0) |

Write P=133220032220, Q=311103311002, U=1332202101012022331,
V=022230022331, W=200113301113, Z=2001131232321311002.
First returns are precisely P, Q, and U{V,W}*Z. By first returns,
this is a uniquely decipherable code. Every atom begins with
1 or 3 and ends with 0 or 2, so B beats the first spine site and
the last spine site beats B. A has in- and out-neighbors in each
atom. Thus both mutual-arc orientations are strong.

## Gap by a telescoping finite-state readout

For signed gap Δ=H(A→B)−H(B→A), the weighted pair term at (a,b) is
L_b(prefix)(1−a)−L_a(prefix)(1−b). See the
[contraction proof](../../unit_shields_wheel_followup_20261116/apex/SHIELD.md).
Multiplying by 13 and writing X=L_a(prefix), we get
X(11−d)−F(1−a). Put X'=2X+a. Then

    13 · weighted increment = X'C'−XC−aC'−F(1−a).

The first two terms telescope on the whole return. Add the local
adjacency increment 13[(1−a_prev)b−(1−b_prev)a]. The resulting
reduced weighted plus internal-adjacent scores on U,V,W,Z are
−2,−26,+26,+2 in units of 13Δ. The seams in U{V,W}*Z have prior
letter 1 or 3 and next letter 0 or 2, so their adjacency increments
vanish. The two exterior endpoint terms, for 1 and 2, cancel.
Thus Δ=2(#W−#V) on a long first return. On P and Q the gaps are
+2 and −2. Across entire first returns, the determinant is zero on
the common ray, and the endpoint correction is zero for every join
because the four allowed last/first pairs are (0,1),(0,3),(2,1),(2,3).
Therefore all gaps add exactly.

## Width and signed gap generating function

Let x index length and y signed half-gap Δ/2. With
A=(y+y⁻¹)x^12 the first-return series is A+x^38/(1−A). The complete
series, including the empty word, is

    G(x,y) = (1−A)/((1−A)^2−x^38).

Widths of a given half-gap d are exactly 12|d|+24k+38j for k,j≥0:
each 12-letter choice contributes ±1, and each long atom contributes
a 38-letter zero-gap overhead. Conversely, take |d| signed P or Q
atoms, k cancelling PQ pairs, and j copies of UZ. The empty word
is counted only in the formal language, not as a positive ray
carrier. No injectivity on counts or isomorphism classes is asserted.

## Width coverage cannot become an absolute wheel on this ray

Each nonempty ray word of length n has B→first and last→B.
Choose any even-size subset of the n−2 interior spine sites between
those endpoints, then close the resulting directed odd cycle through
B. There are 2^(n−3) such cycles, so the positive odd-cycle packing
formula gives H≥1+2^(n−2) in either orientation. Therefore H≤X
requires n≤log₂ X+2.

Put y=1 in G. At x=23/25, the sum converges because
2x^12+x^19<1 (an integer inequality after multiplying by 25^19).
The number of tapes of length at most n is at most G(x,1)x^(−n).
Taking both apex orientations, the number of absolute counts at most X
is therefore O(X^log₂(25/23)), with log₂(25/23)<1/8
(compare 25^8 with 2·23^8). In particular this fixed marked ray
cannot supply even one cofinite arithmetic lane of positive density.
A width wheel can schedule when a gap fits without providing the
absolute landing; this quantitatively separates the two objects.
