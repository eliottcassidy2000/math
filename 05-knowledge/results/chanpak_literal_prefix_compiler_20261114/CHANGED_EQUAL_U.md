# Four changed physical source pairs with exactly equal retained coordinate

For the rank-13 reset, write \(U(a,b)=47a+56b\).
A pair differing by \((112,-94)\) has exactly equal \(U\),
and an appended length-13 step makes the whole states identical.
For a replacement terminal length 10, the ordinary complete
counts differ by \(-2\); for length 14 they differ by \(+2\).
In both cases the difference is independent of the next positive
width, since
\[
 (56,-47)C_{10}=(-1,0),\qquad (56,-47)C_{14}=(1,0).
\]

## A target-derived changed pair, not a general state scan

Demand that both prefixes end in a genuine zero-length reset
\(C_0\) followed by a positive width. They have states
\((A,2^r A)\) and \((B,2^s B)\), where \(A,B\) are actual
ordinary counts of their predecessor tournaments. Require
\(A=B+112\) and \(2^r A-2^sB=-94\). Since both counts are
positive, \(s>r\). Then the valuation of the difference is \(r\)
(both \(A,B\) odd), giving \(r=1\). It follows that
\[
                         (2^s-2)B=318.
\]
As \(s\ge2\) and \(318=2\cdot3\cdot53\), the only possibilities
are \((A,B;s)=(271,159;2)\) and \((165,53;3)\).
The second solution is immediately physical: the gap words
\((1,0,1,0,2)\) and \((1,0,0,2)\) have ordinary counts 165
and 53 respectively. Appending no zero gaps, then the required
positive terminal widths, gives the two real prefixes
\[
   (1,0,1,0,2,1)\longmapsto(165,330),\qquad
   (1,0,0,2,3)\longmapsto(53,424).
\]
Their common rank-13 coordinate is \(26235\).
The claim only uses the physical second alternative; it does
not need a nonattainment argument for the other candidate.
Both graphs and both seed counts are checked by literal subset
DP. Both repaired rails here remain composite because the
reset source counts 165 and 53 survive as respective factors.

## Three more changed pairs without a reset inside either prefix

Here are explicit physical witnesses, with no equal-state
assumption beyond the displayed transfer calculation. Tuples
are the normalized gap widths, and each is interpreted as a
tournament by reversing the arcs of its consecutive anchors
in the otherwise forward order.

| Common \(U\) | First word | First state | Second word | Second state |
|---:|---|---|---|---|
| 1031563 | (1,0,1,0,4,0,1,0,1) | (5701,13636) | (4,0,0,2,0,1,0,0,1) | (5589,13730) |
| 3338481 | (1,0,3,0,0,1,0,1,0,0,0,1) | (18687,43932) | (3,0,0,1,0,3,0,0,0,0,1) | (18575,44026) |
| 6491169 | (4,0,0,3,0,1,0,0,0,0,1) | (35935,85754) | (1,0,3,0,0,3,0,0,0,0,0,1) | (35823,85848) |

For example, first row with length 10 replacement and terminal
width 1 gives \(26360017\) and \(26360019\), exactly two apart,
with the first prime by trial division. In that row, changing
the terminal width arbitrarily gives the two exact rails
\[
 7791883+9284067\cdot2^r,\qquad
 7791885+9284067\cdot2^r,\quad r\ge1.
\]
The second has a fixed factor 9. For any row and either replacement,
the inverse on its own rail is one division and one power-of-two
test: if the row reads \(A+B\,2^r\), require \((N-A)/B\)
to be a power of two at least 2. No claim is made that the
physical pairs are exhaustive or that these rails cover a general
target. They change the actual source and give more than the
original single retained baseline.

Run replay_changed_equal_u.py after replay_source_switch.py.
The small reset example is checked by literal Hamilton DP;
the longer examples use independent interval packing, and every
word also comes with a concrete adjacency-mask reconstruction.
