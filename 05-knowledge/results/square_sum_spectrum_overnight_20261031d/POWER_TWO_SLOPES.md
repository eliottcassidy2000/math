# Exact power-of-two slopes in the singleton/D grammar

**PROVED scoped; no new forbidden raw path count.** Use actual states(F,q,g,h)
from the initial(1,0,0,1), with singleton generators and the strong four-vertex
D. Then h=5^r. If q=2^a>0, precisely the following states occur:

* a>=5r: g=0 and F=q+h, except r=a=0 gives(F,q,g,h)=(1,1,1,1).
* r>=1 and a>=5r-1: g=0 and F=q+h/5, except r=1,a=4 gives(15,16,2,5).

For every such state delta=F-q+g is therefore **positive**: h or h/5.

## Exhaustion by the last non-padding generator

A singleton1 gives q'=2q+h, which is odd. It can be a power of two only when
q=0,h=1, giving the stated q=1 base. A singleton0 doubles q and preserves
delta, setting g to0; it does not change h.

For a D insertion, q'=32q+h*Q, where Q ranges over the16 actual cut digits

    0,5,6,11,9,14,15,20,7,12,13,18,16,21,22,27.

If q'=2^a with a>=5, reducing modulo32 and using odd h forces Q=0. This is
the unique all-zero cut, and it sends q to32q, h to5h, g to0, delta to5delta.
For a<=4, necessarily q=0 and h*Q is a positive power of two. Since h is a
power of5, this forces h=1 and Q=16. The unique cut is mask12 and the source
is the initial state up to its neutral singleton0 loops. Its output is
(15,16,2,5). These are all possible last generators.

Repeatedly peel singleton0 and D mask0. This proves the displayed exhaustive
classification. Every listed state is attained: start from singleton1 or
D mask12, insert the required singleton0 steps, and then the required D mask0
steps. The necessary slope-exponent thresholds are respectively5r and5r-1.

## What this closes, and what it does not

For any fixed actual seed and m>=1, a pure singleton0 suffix has

    F_m=q*2^m+delta.

It cannot realize infinitely many values2^k-c for any fixed c>0. If it did,
the equality would force delta=-c by divisibility by arbitrarily large powers
of two, and then q would have to be a power of two. The classification makes
delta positive, a contradiction. The q=0 case is constant and also cannot
supply infinitely many distinct values.

This obstruction is **not** a numerical spectrum exclusion. The earlier
Fibonacci fans realize2^k-1 by changing their internal word with k. The digital
reflection identity in the main packet gives exactly those addresses. It
explains why a fixed-source pure-zero continuation cannot replace their
reparsing, while leaving mixed suffixes and changing sources available.

`power_two_slopes.py` independently checks the exact action over the complete
actual-word state set q<=256,h<=25. Monotonicity of q,h makes the filters
hereditary. The universal assertion is the predecessor proof, not extrapolation
from that finite set.
