# A tournament countermodel: local factors and logarithmic mean do not give coverage

Consider the transitive spine with one apex. A spine is a binary
word w in which a 1 denotes an apex outneighbor. Restrict to words
starting in 1 and ending in 0, so the resulting tournament is strong.
This simpler model was already known to omit a positive fraction
of odd counts; the point here is to put that omission beside the
same two mathematical inputs one might hope to use for the more
complicated cyclic carrier: all local factors and an unbounded mean.

## Mark an adjacent 10 and count it as two independent addresses

Choose and mark an occurrence of 10 in w. Choose the ordinary
binary address x for the left part, and the reverse unmark address
y for the right part, with 0 for an empty tape. This is an exact
bijection between nonnegative ordered pairs and strong words with
a marked adjacent 10: the word boundaries make each nonempty tape
canonical, and every admissible word has a 10 somewhere.

The [cut transfer](../chanpak_selective_target_20261114/CUT_STATE.md)
with middle word 10 gives
```
H(x,y) = 3+4xy+2x+2y+p(x)+p(y).
```
Indeed the middle has path count 3 with apex, ordinary count 1,
two-slot cover count 4, two unsigned end responses 2, and both
signed readouts zero. The address lemma then gives the displayed
polynomial. Apply the [local and real presentation count](LOCAL_DENSITY.md)
with Q=4,T=S=2,h=1: the same proof only needs the two linear
half-coefficients to be units modulo gcd(h,M), automatic here.
The number of marked presentations through X grows as
(X/4)log X, and each compatible class mod 2M gets exactly 1/M
of that leading term for every fixed M. In particular the
normalized Hardy–Littlewood local factors are all one again.

## Forget the mark and count the actual words

Every strong word begins uniquely with a copies of 1 followed
by 0, a>=1. For this prefix the apex count is 1+2^a and the
unsigned start response is A_a=2(2^a−1) with zero signed bias.
Appending a right address t gives exactly
```
1+2^a + A_a t + p(t).
```
The row estimate from [local density](LOCAL_DENSITY.md), with h=1,
gives a word count X log(1+1/A_a) in this row, at total error
O(sqrt X+log X) summed over a<=log_2 X. In detail, replacing
Y=X−1−2^a by X costs O(1) per row, the row estimates have errors
O(1+sqrt(X)/A_a+log(X+2)/A_a), and the excluded tail contributes
O(1) because the A_a grow geometrically. The logarithmic sum
telescopes:
```
prod_{a=1}^L (1+1/[2(2^a−1)]) = (2^(L+1)−1)/2^L -> 2.
```
So the number of unmarked words with count at most X is
X log 2+O(sqrt X), not of order X log X.

Reversing the word and complementing its bits corresponds to
taking the converse and preserves H. A word of length r has
H>=1+2^(r−2): use the apex, first and last spine vertices, and
an even subset of the interior in the odd-cycle expansion of
[Grinberg and Stanley](https://arxiv.org/abs/2307.05569)
(the r=2 case is immediate). A self-reverse-complement word
has even length r and only r/2 independent bits. Hence there
are O(sqrt X) fixed words through X. Quotient by the converse
involution to obtain an upper bound
```
number of distinct admitted counts through X
    <= (X/2) log 2 + O(sqrt X).
```
At least a relative 1−log 2 of positive odds are therefore
missing in lower density, despite the unbounded first moment
and every compatible local density. Different unmarked words
can still give the same tournament or the same count, so the
display is an upper bound, not an equality.

## Two different generating functions and why they matter

Let P(Z,t)=sum_w Z^H(w) t^a10(w), where a10(w) is the number of
adjacent 10 occurrences. The unrestricted two-address elimination
of the 10 carrier gives (partial P/partial t)(Z,1), while unmarked
word enumeration gives P(Z,1). Both have exactly the same
positive support before or after the converse quotient. But one
has logarithmic inflation of its first moment, entirely due to
the number of marked positions in the same word. A coefficient
cannot be treated as many independent chances to hit a target.

For the larger cyclic carrier this countermodel does not prove
an omission. It gives a test for a proposed major-arc argument:
its witness weight must not grow only from redundant marked
charts. One can use a canonical window rank and count first
witnesses (or bound their dependency), rather than regard each
arc/chart as independent. The new generating-function collisions
in the cyclic carrier make this a live issue, not an analogy
limited to the transitive spine.
