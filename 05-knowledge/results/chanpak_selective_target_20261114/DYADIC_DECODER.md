# Run a two-height selective cut backward from the prescribed count

**Exact decoder, not universal coverage.** A common way to invert a
two-height window is to factor a large translated target. That is not
necessary here. The low binary valuation of the target limits the
location of the two cut heights to a constant-size candidate list.

## The literal window

Let B=TT_k → R → TT_m, where R is the cyclic tournament on residues
0,...,4 with i→j when j−i mod 5 is 3 or 4. Attach apex v beating
the whole left slab, none of the right slab, and either U={0,4}
or U={0,3} in R. Take k,m≥1. The selective cut sums internal to
R are 13 and 14, respectively. Ordered two-path cover weight of R
is Q=120, and h=H(R)=15. These finite computations are rechecked
by literal path DP in the replay.

Here is the physical calculation, also valid with any nonempty
vertex-transitive middle tournament R on n vertices. Let r=|U|.
For paths starting in U, weighted by the Hamilton count of the
complement, the total is A=r(Q−h)/n by transitivity. The weighted
odd minus even quantity is a=rh/n: the sum over *all* starting
vertices equals h because the signed nonempty color convolution is
h, by assigning every odd packing block independently to one of
two colors. Similarly, the ending aggregates for the complement of
U are D=(n−r)(Q−h)/n and d=(n−r)h/n. (These are integers.)
Write C for the internal even cut weight and X=2^k,Y=2^m.
On the left full-win spine the nonempty start weights are X−1
in total and 1 in odd-minus-even weight. On the right full-loss
spine the analogous ending weights are Y−1 and 1. Across an
interior R, both color weights total Q and their signed
difference is zero. Summing the contributions of odd cycles
through v, grouped by their initial and final blocks, gives
[
 H=h+2C+(X-1)D+d+A(Y-1)+a+Q(X-1)(Y-1).
]
Consequently
[
 H=QXY+(D-Q)X+(A-Q)Y+3h+2C.
]
This computes an actual cut, not independent freely assigned
coefficients. For the two named C5 masks it gives
[
 H=120,2^{k+m}-57,2^k-78,2^m+{71,73}.
]
The first pair is 281,283. The source is the same in the two
arms at fixed k,m. Both are strong: the left beats the right,
the right beats the apex, and the apex beats the left; the middle
is accessible from the left and leads to the right.

## Choose the actual source and mask in constant many tests

For each offset c=71,73, let M=(N−c)/3. Reject unless M is a positive
integer. Write t=v2(M). The equation is
[
 M=40,2^{k+m}-19,2^k-26,2^m.
]
The product term has strictly higher 2-adic valuation than both
linear terms. If k ≠ m+1, their valuations do not tie and the
smaller must equal t. Thus test:
* k=t≥1 and 2^m=(M+19·2^k)/(40·2^k−26), with m≥k;
* m=t−1≥1 and 2^k=(M+26·2^m)/(40·2^m−19), with k≥m+2.

In each case the quotient must be a positive integral power of two
with the stated guard. If k=m+1, then
(M=2^{m+4}(5·2^m−4)). At m=1 and 2 its valuations are 6 and 10;
at m≥3 it has valuation m+6. Therefore test only m=1 when t=6,
m=2 when t=10, and m=t−6 when t≥9. Check the original equation in
the tie cases. These are necessary **and** sufficient in the
literal family; a successful branch gives the nonnegative orders
and chooses U before building the tournament. Arithmetic uses
big integers, not a search over possible tournaments or divisors.
For example N=13897 is decoded as c=73, (k,m)=(4,3).

## A genuine changing-core triple, not just two masks on one core

The three marked cores in [SWITCH](SWITCH.md) give
$N=432\,2^{k+m}-303\,2^k-174\,2^m+d$ with $d=225,227,229$,
starting at 999, 1001, 1003. Each offset identifies an actual
core and mask. The two unequal-valuation branches are again
$k=t$ or $m=t-1$ for $t=v_2(N-d)$, after which the other height
follows by exact division. In the tie $k=m+1$, one has
$N-d=2^{m+2}(216\,2^m-195)$, so $m=t-2$ and only that pair can
work. The C5 pair and changing-core decoders are both included
in the literal replay; the latter needs at most three height
tests per offset.

## General bounded-branch principle

For fixed nonzero a,b,c, consider
(N-d=a2^{k+m}+b2^k+c2^m). Restrict first to a fixed sufficiently
deep quadrant; the finitely many horizontal/vertical strips
outside it each have a one-height linear inverse, including the
case of a constant function. The depth can be chosen so the
product term dominates the sum of the linear terms in ordinary
magnitude and has strictly higher binary valuation than either
of them. Let u=v2(b), v=v2(c), w=v2(a). For nonzero M=N−d, put
t=v2(M). In the nontie case one of k=t−u or m=t−v is forced,
after which the other power of two is determined by a linear
equation.

In the tie case put s=k+u=m+v and
q=b/2^u+c/2^v. Factoring 2^s, the two internal terms have
valuations z=v2(q) and s+w−u−v. If q=0 they leave
t=2s+w−u−v. Otherwise:
* if the second valuation is lower, t=2s+w−u−v;
* if the first is lower, t=s+z;
* if they are equal, s=z+u+v−w is fixed and the full
  equation is tested.
This gives a constant-size candidate list for a fixed readout,
regardless of the target's size. It is especially helpful for
comparing many physically marked middle modules, where a scalar
modular hit alone is not a certificate.

A fixed finite bank of such windows is still sparse. This decoder
solves exact absolute placement only when the chosen physical
family has a match; it does not build a family for every N.
