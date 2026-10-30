# Complete weight frontier after the four-deletion AP2 lemma

**PROVED algebraic classification, conditional on the separately audited structural row.**
No circuit lower-bound record or literature priority is asserted. The old row table
and its unchanged coverage/affine-state assumptions are pinned in the Session 117
source audit. This document classifies weights, not the validity of every source claim.

For W=g s+i I+p q, normalize g=1 and require i,p>=0. With four deletions in AP2,
the nontrivial row constraints become

    i<=2+p,
    t<=2+2i-p,
    t<=4+2p,
    t<=3+i+p.

They represent paired cleanup, one-child GP, strengthened AP2, and AP1 respectively.
All other displayed rows are nonincreasing for nonnegative weights; additional GP
children only improve its charge. Let c=t-i and b=t, and restrict c>=0.

## Exact classification

A pair (c,b) is feasible if and only if

    0<=c<=4,   b>=max(c,2c-2,4c-8).

All corresponding normalized weights are given by i=b-c,t=b and

    max(0,b-c-2,(b-4)/2,c-3) <= p <= b-2c+2.

Proof: substitution in the row constraints gives this interval and b>=c.
Its four lower endpoints are each at most its upper endpoint exactly when
b>=2c-2, c<=4, b>=4c-8, and b>=3c-5. The last is redundant: for c<=3 it is
implied by b>=2c-2, and for c>=3 by b>=4c-8. This proves necessity, sufficiency
and the complete fiber of weights, not just one optimized choice.

The sharp lower envelope has three pieces:

    b_min(c)=c                 for 0<=c<=2,
    b_min(c)=2c-2              for 2<=c<=3,
    b_min(c)=4c-8              for 3<=c<=4.

Boundary witnesses (i,p,t) are respectively

    (0,0,c), (c-2,0,2c-2), (3c-8,2c-6,4c-8).

At c=4,b=8 the weights are uniquely (g,i,p,t)=(1,4,2,8).
Thus the strengthened row reduces the coefficient-four dimension penalty from
nine to eight. It does not increase the leading coefficient beyond four.
The obstruction c<=4 follows from cleanup and GP alone, untouched by AP2.

## Conditional telescope, with endpoint information retained

If an elimination sequence with the stated invariant and coverage exists, starts
at a=n,q=0,I<=n,s=|C|, and ends at a_T<=d+1, each reducing row obeys
Delta W<=8 Delta a for W=s+4I+2q. Hence

    |C| >= 4n-8a_T+s_T+4I_T+2q_T+4(n-I_0)
         >= 4n-8(d+1).

This is a rigorously derived implication from the named sequence assumptions.
It is not an independent certification of the provisional disperser construction,
all-length language assembly, or unrestricted leading-coefficient record. The
new local lemma and exact global bookkeeping can be used without promoting the
entire external candidate to an established theorem.

The sequence accounting companion retains additional surplus and proves that
AP1/AP2 events themselves consume a terminal budget. A real improvement beyond
four must address the GP/cleanup cycle, enforce extra deletion or usage loss, or
use circuit minimality/semantic hardness in a new way. Changing these weights alone
is now completely classified and cannot accomplish that improvement.
