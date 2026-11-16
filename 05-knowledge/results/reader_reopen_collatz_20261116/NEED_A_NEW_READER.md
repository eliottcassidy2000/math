# Changing only the cut is insufficient

At the first next-digit count 7747 and again at r=11, v=2
(next-digit count 16974403), there is no realization by a
tournament with a transitive deletion of one vertex, at any order.
This is a finite exact obstruction to the one-apex family,
not a further global exception.

Such a tournament is a binary spine word w. Initial zeros of w
and final ones are neutral singleton components, and can be
trimmed without changing H. With H>1 the reduced word starts
with one and ends with zero. At reduced length n, the standard
tape update gives 1+2^(n−2) ≤ H ≤ 2^n. The lower bound comes
from the first one closing at the final zero, with every even
subset of intermediate sites. The upper bound follows from
H+p≤A+1 in the forward readout: at a zero H' = H+A+p≤2A+1;
at a one H'=H, p'=1 and A'=2A+1. Thus cutting or changing the
length cannot evade the small set of candidate bit lengths.
The complete reduced-word loops in no_one_apex.cpp find zero
hits at n=13,14 for 7747 and n=25,26 for 16974403.

So a genuinely different reader is necessary at least at those
heights. This does not require adding vertices at the first
height if arbitrary small internal rewiring is allowed. On the
original 15-vertex 7745 reset graph, label tape and sentinel
0,…,13 and the head 14; its marks are 1,3,10,11,13. Reverse the
four arcs on pairs (3,11), (0,14), (1,14), (13,14). Literal subset
DP gives 7747. No minimum fixed-order edit claim is made here.
The three-arc and two-new-vertex conversion in BANK.md is the
one proved uniform in height and preserving every spacer contact.
