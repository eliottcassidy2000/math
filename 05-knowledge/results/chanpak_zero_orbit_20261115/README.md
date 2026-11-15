# Genuine collisions and a support-preserving zero-return quotient

This follow-up to the [viability packet](../chanpak_viability_shape_20261115/README.md)
separates three objects: literal graph isomorphism classes, repeated
visits to a count in one transport chain, and unrelated collisions
between different chains. It gives an exact quotient for the second.
It does not assert positivity for all odd targets.

## Addresses are distinct actual tournaments

Use the fixed [cyclic core carrier](../chanpak_absolute_steering_20261115/TWO_TAPE_STEERING.md).
The apex is the unique vertex whose deletion leaves all strongly
connected components of order at most 8. Deleting it leaves the
strong order-8 core in a line of singleton SCCs. If a noncore spine
vertex is deleted, the core together with the apex survives as a
strong order-9 subgraph, because the strong core has contacts of both
directions at the apex. If a core vertex is deleted, the first and last
spine vertices survive as singletons before and after the core. The
apex reaches the former, and the latter reaches the apex; every surviving
vertex reaches and is reached from the apex through these singletons.
This SCC has order at least 10. The identification is intrinsic to
the tournament. Recover its apex, read the ordered SCCs of its
deletion, and remove the two guards directly next to the core. The
marks of the remaining tapes recover the ordered addresses x,y.
Thus equal F values at different ordered addresses are collisions of
**nonisomorphic tournaments**, not just duplicate names for one graph.

This is different from the [marked-word countermodel](../chanpak_viability_shape_20261115/MARKED_WORD_COUNTERMODEL.md),
where a single tournament can generate a large number of aliases.
It does not imply a similar uniqueness claim for other carriers.

## An optimal logarithmic example inside one primitive chain

The first two widths of the
[matched transport spectrum](../chanpak_viability_shape_20261115/ALL_WIDTH_TRANSPORTS.md)
give zero neutral moves at width 12 and ten at width 24. The ten are
all primitive because none of the four width-12 losses
2, −48, 34, −16 sum in pairs to zero. With A=2^24, one of the ten is
U=11418043,V=11767483, giving the identity
`F(x,A z+V)=F(A x+U,z)` for every outer x,z. Its two matched
linear coefficients and its zero constant difference are checked
directly in the replay.

Let R_j=(A^j−1)/(A−1). The k+1 ordered pairs
`(U R_j,V R_(k−j))`, 0<=j<=k, give the same exact
count N_k=F(0,V R_k); their graphs are pairwise nonisomorphic.
For example N_1=46129409497 and N_2=773923063178445579.
Since p(y) lies between 0 and y, N_k has order A^k.
These give honest logarithmic multiplicity in individual
target fibers. A pointwise uniform constant bound on the
number of representations is therefore false.

The example is asymptotically sharp **among equal-count visits
along a single primitive transport chain**. Between any successive
same-count visits, loss is zero, so the total width of its block is
at least 24. With m consecutive such returns, the original right
address is at least 2^(24m)/2, because each block radix is at least
2^24 and its highest consumed block starts with 1. But
y<=(N−2855)/3914. Thus the number of visits is at most
1+log_(2^24)(2(N−2855)/3914) when there is a repeated visit.
There may be additional equal-count visits in other chains.

## First return parses and counts all within-chain collisions

Take the alphabet of **all** primitive matched transports,
with all signed losses, not only unit moves. The primitive V codes
and the primitive U codes are each suffix-free, and each
chart's U or V identifies the other. At any pair there is at most one
forward and one backward successor; both terminate, since the consumed
address strictly shrinks. The component is a path with a signed loss
assigned to every edge. Counts need not be monotone.

A **first zero return** is a nonempty word of primitive transports
whose total loss is zero and whose proper nonempty prefixes do not
have zero loss. It is a physical composite chart, so can be moved
in one step without changing the count. Such words are prefix-free
in chronological order, and their reversals are prefix-free too
(a proper zero suffix of a zero word would leave a proper zero
prefix). Suffix-freeness of the primitive binary codes therefore
makes the composed V codes and composed U codes separately
suffix-free. Along a primitive chain the first zero return joins
one visit to the next visit at the same count. Hence the oriented
zero-return rewrite has exactly one sink for each pair
(primitive chain, attained count). This quotients all identified
within-chain equalities, not unrelated collisions between chains.

Write R(Z)=sum_(x,y>=0) Z^F(x,y). If the zero-return V codes are
indexed by b, of radices A_b, the exact support-preserving series
is
```
R_canon(Z) = R(Z) − sum_b sum_(x,z>=0) Z^F(x,A_b z+V_b).
```
The subtracted cylinder sources are disjoint and each nonsink
is counted once. Every coefficient remains nonnegative, and
positive iff the corresponding raw coefficient is positive.
Unlike replacing a series by its residue image, this subtraction
is exact at the prescribed integer.

## The whole identified quotient removes only a tiny real factor

Put t=1/4096 in the bivariate transport series and let
E=[Z^0]T(t,Z), where T is the series of all nonempty charts,
Z marks displacement, and t marks units of 12 bits. Each word's
weight is the reciprocal of its radix. Parse any zero-total word
uniquely into a sequence of first zero returns. If rho0 is the
sum of their weights, the exact renewal equation is
```
1+E = 1/(1−rho0),   rho0 = E/(1+E).
```
The first eight zero-chart coefficients by width are
(0,10,0,5250,0,3653650,0,2834954850); the tail after these
is below 1.37e−20 by comparison with the full chart counts
4·32^(k−1). This puts rho0 between 0.0000005960647449810
and 0.0000005960647449811. This is the limiting fraction
of raw presentations removed, less than one per million.

For one fixed discarded suffix of radix A, the presentations
through X have leading term X log X/(5760 A), by the
[translated hyperbola calculation](../chanpak_viability_shape_20261115/LOCAL_DENSITY.md).
The same fraction holds in every fixed compatible odd class:
approximate by the positive shifted product count at error O_A(X),
then vary the unrestricted x. The fixed-suffix row coefficient of
x after division by two is a unit modulo gcd(45,M), so the
prefix mixing argument distributes each x interval to its
required class with sublinear error. Summing those errors over
z up to O_A(X) is O_(A,M)(X). The tail across arbitrarily long
discarded charts is uniformly controlled because V/A>=1/2 and
F(x,A z+V)>=c A(x+1)(z+1). Together with the full chart count,
this gives a convergent geometric bound on omitted radix weights.
Therefore the canonical first moment is still asymptotic to
(1−rho0)X log X/5760 overall, and to the same quantity divided
by M for each compatible residue modulo 2M.

This shows why exact quotienting is useful but still not a
coverage proof. The explicit correlation causes unbounded true
multiplicities at some targets, yet has too little mass to explain
away the logarithm in the global first moment. The next question
is multiplicity **between different canonical chains** in
unbalanced ranges.

## Reproduction

`python replay_zero_return.py` also works with `-O`.
The background collision example and the full first-return
quotient are distinguished above.
