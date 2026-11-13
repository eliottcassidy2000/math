# Exact reservoir switches and the synchronization tax

Use the [terminal atlas](../session317_exact_anchor_atlas_20261113/TERMINAL_ATLAS.md).
A terminal word with physical parameters (h,E,s,c) evaluates at prefix B to
F(B)=c+(h+E)B-hK(B)-s(B mod 2). A reservoir change below means an actual
change in h, not merely a different realization of the same scalar profile.

## A literal switch and an interlocking triple at every scale

For every k>=4, the two charts
| terminal | (h,E,s,c) | binary prefix value |
|---|---|---|
| 2 | (1,2,1,3) | 2^k+3 |
| 30 | (3,3,2,5) | 2^(k-1)+1 |

both realize **9*2^(k-2)+9**, with strong literal tournaments. For every
k>=4 the following three charts realize **21*2^k+15, +17, +19** respectively
(the +17,+19 replace the +15 intercept):

| terminal | (h,E,s,c) | binary prefix value |
|---|---|---|
| 320 | (5,16,3,15) | 2^k |
| 30 | (3,3,2,5) | 2^(k+2)+2 |
| 2 | (1,2,1,3) | 2^(k+3)+6 |

In the triple the reservoirs have three different counts 5,3,1. In the
first pair the exact output remains unchanged even though the reservoir
changes. The charts are not formal signed combinations; replay_switch.py
checks their literal parameters, strongness and direct path counts at
small order. All k follow from the theorem below, not from that replay.

Write B=2^m X+y with 0<=y<2^m and m>=1, v=rev_m(y), and eta its highest
bit. The exact inverse identity gives
F(B)=C+AX-hK(X)-h eta (X mod 2),
where A=2^m(h+E)-hv and C=c+(h+E)y-hK(y)-s(y mod 2).
Taking X to be any power of two >=2 kills both remaining correction terms.
For the pair take m=4; the coefficients A are 36,72, both with C=9.
For the triple take m=4; the coefficients are 336,84,42, with C=15,17,19.
The stated changes of exponent match the slopes. Strongness follows
from the atlas reachability criterion once the first prefix bit is 1.

## A changing reservoir cannot synchronize on generic prefixes with fixed seams

Compare any two such charts with a common variable high part X, possibly
different fixed endings, tail lengths and appended low words. Let their
terminal h differ. Fix an integer tolerance L>=0. Among X<2^r there
are O(r^2) values where the two Hamiltonian totals differ by at most L.
The constants may depend on the charts and L. In particular exact equality,
and adjacency of two fixed-width triples, have zero density in the
common prefix. A finite menu of fixed chart pairs does not help this density.

**Proof.** Subtract the two chart formulas. Write the result as
aX-dK(X)+b+e(X mod 2), with d nonzero. Put
V(X)=sum over j>=0 of b_j*2^(-j-1), where b_j are the bits of X from
the right. Unwind K as sum_j b_j(floor(X/2^(j+1))+b_(j+1)).
Consequently |K(X)-X V(X)|<=1+floor(log2 X). On the shell
2^u<=X<2^(u+1), bounded discrepancy forces
|V(X)-a/d| <= [u+1+(L+|b|+|e|)/|d|]/2^u.
But V on that shell is an exact translated grid with spacing 2^(-u).
So there are O(u+1) candidates per shell, O(r^2) in all. If a/d lies
outside [0,1], there are no sufficiently large candidates. The same argument
yields an explicit finite candidate list on each shell: intersect the
rational interval with the dyadic mesh, reverse the digits, and test the
exact count. This is not merely the uniform limiting law; it is an exact
sparse inverse sieve for close or identical outputs. If h is equal, the
K term vanishes and the subtraction is affine on each parity class.

What the literal switch does: it sends two differently cut end reservoirs
to the same absolute count, along an infinite ray. What it does not do:
return a reusable endpoint state, select such a switch at a prescribed
arbitrary odd target, or prove that the many isolated triples overlap.
To escape this obstruction with the atlas, the fixed seams must change
with the target, or the high prefixes must be allowed a more general relation
than sharing the same X.

## Bounded retiming does not remove the tax

The obstruction survives comparing chart 1 at X to chart 2 at X+t for
any fixed integer t: at most O(r^3) near-matches below 2^r, uniformly for
t in any fixed finite set. It is enough to take t>=0 after exchanging
charts and translating X if necessary. Choose fixed p with 2^p>t and
write X=2^p Z+y. For each of the finitely many y, addition either
does not carry to Z or increments Z by one. In the latter case partition
by the length j of the final run of ones of Z. In each branch X and X+t
have identical bits above a fixed low word (of length p in the first
case and p+j+1 in the second), so V(X+t)-V(X) is a branch constant.
The length j on the shell is O(r). Substituting K(n)=nV(n)+O(log n)
in the chart difference and absorbing the bounded t V(X+t) term shows
that V(X) must again lie within O(u/2^u) of a branch-specific rational
on shell u. The constant in that bound is uniform in the branch:
only the rational center changes. There are O(u) branches and O(u)
mesh points per branch, hence O(r^3) in total. Allowing bounded shifts
of the high prefix therefore still cannot yield a dense return rule.
