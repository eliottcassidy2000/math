# The minimum right-readable carry and the limits of a left reader

A background follow-up distinguishes right-only matching from the
all-context switch [whose minimum puncture count is five](MINIMAL_SWITCH.md).
The following examples and classifications were replayed independently.

Write the puncture profile as (h,f,P,Q,S,T), with S,T the nonnegative
endpoint counts (negatives of signed responses), and augmented count K.
For a fixed puncture C followed by any ordinal right puncture D, with
any contact set in D, the apex readout is
\[
 K_{CD}=h_DK_C+h_C(K_D-h_D)+P_CQ_D+S_CT_D .
\]
Thus if two C options match (h,P,S), any right continuation multiplies
their count difference by h_D. Matching f as well keeps the entire
(h,f,P,S) fixed after the continuation, by the ordinal start rule and
the fact that full paths end in the last block. For all transitive
right readers the gain therefore remains literally two.

## Two exact one-sided switches

Rows are out-masks in the puncture, U is the apex out-mask.

| puncture rows | U | (h,f,P,Q,S,T) | K |
|---|---:|---|---:|
| (4,1,2,7) | 7 | (3,24,9,3,3,3) | 3 |
| (10,0,3,6) | 3 | (3,24,9,6,3,2) | 5 |
| (20,1,2,7,14) | 22 | (9,84,41,22,7,4) | 23 |
| (16,1,3,7,14) | 22 | (9,84,41,25,7,5) | 25 |

The first pair is minimum in h: a puncture with h=1 is transitive.
Matching f fixes the width, and the positive response P is its contact
word in binary (leading zeroes retained by that width). Thus it fixes
the whole contact pattern and cannot change K. The first h=3 graph
has an ordinal singleton source, which is neutral for ordinary counts;
minimum h here does not mean every augmented choice is strong.

The second pair changes only the old arc {0,2}. For a transitive
right reader of reverse-complement binary address y, the outputs are
\[
 23+9(g(y)-1)+41y+7(y\bmod2),\quad
 25+9(g(y)-1)+41y+7(y\bmod2),
\]
where g is the usual transitive-spine apex count. At y=925 (contact
word 0100011000) this gives two strong 16-vertex graphs with exact
counts 40,241 prime and 40,243=7·5,749 (prime cofactor). Both were
checked by literal counts. Neither target occurs at **any** order
with a transitive puncture and one apex: normalize the word and
test all binary addresses below 4(N−1), using its mandatory first-
to-last inversion as the order bound. The replay checks this entire
range for both outputs.

## Exactly when an arbitrary left reader leaves the unit gain intact

Prepend any fixed ordinal puncture L, with any contacts, and append
any fixed right puncture D. Empty sides are permitted (h=1, P=S=0).
For either pair above, the exact total high-minus-low difference is
\[
 h_D\bigl(2h_L+3P_L+\epsilon S_L\bigr),
 \qquad \epsilon=-1\text{ for the first pair},\
              +1\text{ for the second}.
\]
This follows by reading the two-block law first with L as prefix and
then with D as suffix. The intermediate (h,f,P,S) still agree.
Since the absolute signed response S_L is at most P_L, this quantity
is at least h_D(2h_L+2P_L). It equals two if and only if both L and D
are transitive (or empty), and L makes no apex-out contacts. The latter
is exactly P_L=0, since every contacted singleton contributes positively
by itself. In particular a genuinely contacted left reader can never
restore the fully matched unit switch, nor make these pairs neutral.
The one-sided h=3 minimum and two-sided h=5 minimum are different
phenomena, not contradictory counts.

Run **python3 replay_right_only.py**. This supplies a useful additional
exact repair chart but not absolute coverage of every odd target.

## Four probes identify the entire right-readable count interface

This distinction is exact, not just a sufficient condition for the two
examples. Read a core with a transitive right word at y=0,1,3,2,
in that order. Since g(0)=g(1)=g(3)=1 and g(2)=3, the four
values are
\[
          K,\qquad K+P+S,\qquad K+3P+S,\qquad K+2h+2P.
\]
They uniquely recover (K,h,P,S). Consequently a K change propagates
as h_D times the same change under **every** right tournament reader
if and only if it does so at these four transitive probes. The count
interface does not need f: the formula for P after appending D uses
f_D, not f_C. If one also asks to preserve the old two-slot count,
f is a separate requirement. This fits the pointed-support lens:
one must keep the finite directed boundary probes that future gluing
actually sees, not every coordinate of every phase observer.

Likewise seven probes recover the full two-sided profile (K,h,f,P,Q,S,T):
use (x,y)=(0,0),(0,1),(0,3),(0,2),(1,0),(3,0),(1,1).
The four right probes give (K,h,P,S), the two left-only probes
give Q,T, and the last probe exposes f. Thus the all-context
fixed-profile switch criterion is both necessary and sufficient
for the change to be independent of arbitrary transitive wings.
