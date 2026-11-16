# The first 25 integers and primes, viewed in equal fifths

This builds on the repo's *Ogg's fifteen primes along levels 1|2|6,
the square-sum graph between 12 and 25, and what {7,21} can and
cannot share with the Monster* (monster_levels_square_sum_7_21_20261022.md).
That packet provides the Fi22/B/M prime layers and the square-sum
first-path windows 15,23,25. Here the comparison uses a **marked
prime rank**, not an assertion that p_i preserves square sums.

## Two complete rail matchings occur first at 25

Split [1,5m] into consecutive blocks B1,...,B5 of width m.
A fixed square sum gives a complete matching B1→B3 only when
its reflection constant is 3m+1; similarly B3→B5 requires 7m+1.
Those two integers are squares iff \(a^2=3m+1,b^2=7m+1\),
equivalently \(7a^2-3b^2=4\). The positive solutions are exactly
the forward orbit of \((1,1)\) under
\[
(a,b)\mapsto ((5a+3b)/2,(7a+5b)/2).
\]
Indeed the equation implies the same parity and \(a^2\equiv1\pmod3\).
Its inverse has positive coordinates for \(a>1\), because
\(7/5<b/a<5/3\), and its first coordinate strictly decreases
since \(b>a\). The unique \(a=1\) solution is \(b=1\). The forward
map preserves the equation and is integral, closing the proof.
Thus the positive widths are \(5,120,2760,63365,\ldots\); [1,25]
is the smallest nontrivial equal-fifth rail frame. At each width
the paired reflections are the m disjoint paths
\[
 i\longleftrightarrow 3m+1-i\longleftrightarrow4m+i.
\]
The two reflections compose to the honest translation by 4m.

A square S=5m+r between B1 and B5 contacts rail i and rail m+r−i
when both indices are in [1,m]. It supplies exactly
max(0,m−|r−1|) edge contacts before contracting rails; after
contraction an exchanged pair gives two contacts, and a fixed
rail one loop. In particular m=5, S=25 links rails 1↔4 and
2↔3 twice and leaves rail 5 unpaired. The same square links
B2 and B4 four times (6↔19,7↔18,8↔17,9↔16).
S=5m is not automatically a square at later Pell widths; it
requires the additional simultaneous Pell equations
\(a^2=15t^2+1,b^2=35t^2+1\) with m=5t². A nearby square
always gives a nearly full contact (shift O(sqrt(m))).

## Prime marking of the five rails

The prime quintiles have (Fi22, B-but-not-Fi22, M-but-not-B)
mark sums (5,0,0), (1,3,1), (0,2,1), (0,0,2), (0,0,0).
The third quintile's marked ranks are exactly 11,13,15.
Each doubled reconnection 1↔4 and 2↔3 therefore pairs a
rail with Ogg-marked middle point to a rail without one,
leaving the third marked rail isolated from those reconnections.
The absent even-block square-25 rail skips ranks 10 and 20,
29 and 71, both in M but not B. This is a finite property
of the chosen labels. It does not lift to a square-sum map:
of the ten rail edges using 16 and 36, only (1,15) survives
replacing each rank i by p_i (since 2+47=49).

## The first 25 primes as an exponent sieve

For k≥1, the prime p either never divides \(99\cdot2^k+23\)
or reads a single class k≡a mod d, with d the order of 2
mod p (except p=2, which never divides the odd target).
Compute each class directly inside one order window.
The exact surviving fractions after each block of five are
\[
3/4,\quad 41/72,\quad 41/72,\quad
131893/248820,\quad 5392213/10201620.
\]
They are calculated by CRT inclusion-exclusion, not an estimate,
and describe small-prime survival, not primality. The equality
at blocks two and three has an exact reason: only 37 has a
class in the third block, k≡23 mod36, which was already removed
by 5's k≡3 mod4. Similarly 29's k≡27 mod28 is redundant.
All other displayed active classes are needed. After removing
the two redundancies, the nine active prime tags are
5,13,17,19,53,59,61,67,83; their joint period 61,209,720 is
minimal. For example the reduced survivor denominator is
10,201,620 (ratio six); k=5 distinguishes the half-period
shift and k=21 distinguishes the one-third-period shift.
Thus no proper divisor of the proposed period is a period.
The final quintile only adds 83's class k≡15 mod82; its
first new rejection after the earlier 20 primes is k=917.

Replay with python3 verify.py. The rail theorem is all-order;
the marking and exponent sieve make no claim of any transport
to Hamilton counts or equality at a prescribed integer.
