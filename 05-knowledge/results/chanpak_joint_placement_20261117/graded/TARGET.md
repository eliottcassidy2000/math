# An exact joint decoder without prescribing the residual first

The [general fixed-shift inverse](../inverse/PROOF.md) takes both the
absolute count and displacement as inputs. On the original
near-balanced branch, there is also a direct inverse for a **single
prescribed endpoint**: the signed displacement, if any, is recovered
from that target. This does not make every endpoint feasible.

Fix Z>=1, and use 1<=z=x-y<=Z, |Delta|<=G², all complementary
suffixes, and both endpoints of the physical repair (and both
reflected partners). For all sufficiently large dyadic G, for
each fixed choice of source/destination and lower/upper, a prescribed
target has at most one triple (x,y,r) at that G, and one can find it
in O(log N) integer comparisons and O(log N) evaluations of F.
Enumerating G uses O(log N) levels. Finitely many small G are an
explicit bounded enumeration since z and Delta are bounded.

First subtract 2 when the prescribed endpoint is upper. Put
S=u+w, d=u-w=16Gz+2r+1. The source count is

    H = 64 S²+142 S+81 - 64 d² + 5(F(u)+F(w)).

For the destination replace d by d+8G, retaining S. Under the
bounds of the branch, S=Theta_Z(G³), by the displacement identity;
both d and d+8G are O_Z(G). For sufficiently large G, the
difference of either endpoint from T(S)=64S²+142S+81 is at
most 10S, using F(u)+F(w)<=S+2. Both T(S)+10S and T(S)-10S
grow strictly for S>=1, and the intervals for successive integer
S are disjoint since T(S+1)-T(S)>10(2S+1). Binary search
determines the only possible S, without specifying Delta.

The sum now determines K=x+y via S+1=G(16K+10). The repair
equation is

    Delta = (15K - (2A+15)z)/2 - B_r.

Across suffixes B_r is in [4096G²,6200G²+1100G] for G>=4,
and its spread is at most 2048G²+40G (vary r in [1,G)).
Thus the uncertainty in B_r plus the allowed Delta is less
than 4200G², whereas stepping integer z changes the right
hand side by (2A+15)/2>16000G². An exact enclosing interval
using the displayed B bounds leaves at most one integer z,
and hence at most one candidate pair x,y.

It remains to recover r. For fixed x,y the product of the
two source addresses decreases strictly in r, since their sum
is fixed and d>0. Incrementing odd r by 2 decreases its
contribution to H by 64[(d+4)²-d²]>=512d>=8192Gz.
At the destination use d+8G instead, which only increases
this saving. The F-sum can never compensate. To see this
directly, write the two high prefixes as a,b (with a-b=16z+1
at the source and 16z+9 at the destination). Their suffix
decomposition makes its variable part

    RC(r) a + RC(s) b + F(r)+F(s)
       + g(r)[a odd]+g(s)[b odd],

where g is the leading-zero flag. Complementary suffixes
imply RC(r)+RC(s)=G-1. Across a step, the first two terms
therefore change by at most (16z+9)G; F(r)+F(s) varies by
at most G, and the parity term by at most 2. Even times five
this is strictly less than 8192Gz for z>=1. Hence actual endpoint
counts, not just approximations, decrease strictly with odd r.
Binary search the exact F-based count and then check the
address parities and exact Delta. If no r hits, the branch
does not attain that target at this scale.

The sum/imbalance filtration thus decodes this channel without
requesting a separate correction. It can report that the
starting-count fiber is empty, consistent with the recurring
octal holes. The z=0 chamber and additional contact grammars
are outside this near-balanced uniqueness statement.
