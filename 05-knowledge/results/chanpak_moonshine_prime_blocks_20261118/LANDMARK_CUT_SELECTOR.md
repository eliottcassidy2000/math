# Hypercube landmarks and a height-addressed one-contact selector

## The four hypercube questions: current inherited status

The linked [arXiv:2607.28813](https://arxiv.org/abs/2607.28813)
is Allikvere's *Multiset Dimension of Graphs: Extremal Values and
King Grids*; the four questions quoted are from his separate
[edge-multiset hypercube paper, arXiv:2608.09983](https://arxiv.org/abs/2608.09983).
They concern distance histograms from an **edge** to a vertex set,
not ordinary labeled landmark distances or vertex multiset dimension.

1. **Q6 is exactly 15.** [THM-4525](../../../01-canon/theorems/THM-4525-edge-multiset-dimension-of-q6-is-15.md)
   records three independent exhaustive symmetry-reduced searches
   excluding every size through 14. There are exactly 229 orbits
   of resolving 15-sets, all with trivial stabilizer. Failure of
   annealing is not the lower-bound argument. The new replay
   independently recounts the positive masks (including all 229
   given the archived list); it does not redo the exclusion searches.
2. **Growth order is known.** [THM-4534](../../../01-canon/theorems/THM-4534-edge-multiset-dimension-grows-like-exp-theta-cube-root-d.md)
   proves log edim_m(Q_d)=Theta(d^(1/3)), with certified liminf
   lower constant about 0.8146 and limsup upper constant about
   2.0526. The existence and value of a limiting constant are open.
3. **The original table is no longer needed to prove finiteness.**
   The same theorem gives a closed-form probability estimate from
   d=17 onward and explicit resolving sets at the smaller finite
   dimensions. The stricter request for one analytic probability
   estimate from d=11 is still open. The previous Fourier–Hölder
   computation from d=10 is an exact numerical certificate, not
   that analytic estimate.
4. **Half density is far from optimal.** The audited sparse random
   bound achieves expected size exp((2.0526…+o(1))d^(1/3))
   with the usual positive slack in the exponent; sampling at
   density one half instead costs about 2^(d−1). These results
   establish order of growth, not an optimal random threshold.

### A uniform analytic partial bound from d=11

Here is a simple monotonic argument disposing of some types without
the per-d rational forest optimization. It is intentionally scoped
to parallel pairs at Hamming offset h=1,...,5 or h=d−1.
Take independent fair-bit landmark inclusion. By
[Allikvere's forest lemma](https://arxiv.org/html/2608.09983v1),
a level edge of multiplicity 4C gives collision atom at most
sqrt(2/(4*pi*C)) < 1/sqrt(6C). Let n=d−1 and g=n−h.

For h≤5 choose p=floor((h−1)/2) in the opposed coordinates.
As the common coordinate count r ranges from 0 through g,
the level edges {p+r,h−p+r} form a forest: they connect
integers at a fixed positive distance. Their multiplicities
are 4 binom(h,p) binom(g,r). The number of parallel pairs of
this type is d·2^(d−2)binom(n,h). Keep this h fixed while
increasing d by one. Pascal increases every old shell factor,
including the old last factor by g+1; the new endpoint
contributes 6 binom(h,p) in the denominator squared.
Thus the bound for the whole orbit shrinks by at least the
factor
```
[2(d+1)/(d−h)] / sqrt(6 binom(h,p)(d−h)) <= 4/sqrt(36) < 1.
```
For h=n, the level graph is the antipodal matching, with
multiplicities 4 binom(n,p) for p<n/2 (not a connected graph).
Increasing d by two increases every old shell and introduces
one more with binomial coefficient ≥d+1. The orbit bound
therefore shrinks by a factor at most
```
4(d+2)/(d sqrt(6(d+1))) < 1        (d>=11).
```
At the initial d=11 the five fixed-offset terms are bounded
by 0.00001, 0.0008, 0.0013, 0.044, 0.132. The antipodal term
is <0.019 at d=11 and d=12. The replay checks these base
inequalities by squaring positive rational comparisons.
Consequently the sum of the six orbit bounds is **<1/5
uniformly at every d≥11**. In particular these parallel types
do not require separate certificates anywhere in 11–50.
Intermediate parallel offsets and cross directions are not
covered, so this does not close the strict version of question 3.

## How a cube coordinate can act as a tournament landmark

The old Q6 audit also tested natural tournament-structured sets
in the five-vertex tiling cube: unions of isomorphism classes,
Hamilton-count levels, or score classes did not resolve all
edge histograms. Forgetting labels can hide exactly the
distinction needed to invert a move. There is a complementary
way to exploit that observation: **mark the one coordinate
where two physical construction words differ, then invert
only the conditioned fibre**. This is not a reduction from
the hypercube dimension problem; it uses the same distinction
between an observed scalar and the data retained to resolve it.

For two rooted sources at the same (h,e,d) with start counts
t<s, use the [binary transfer](ANS_INVERSE.md). Let W be the
word address of the higher-t source, with bit p zero. Give the
lower-t source the address W+2^p. Their outputs agree precisely
when, writing D=s−t and w_ij=2^(j−i−1)+[j=i+1] for i<j,
```
D W + d sum_(i≠p) w_min(i,p),max(i,p) b_i
    = t·2^p − e[p=0].
```
The shared quadratic terms cancel; only contacts with the
marked bit remain. This is an exact equality at absolute
height, not just a congruence. Split W=x+2^(p+1)y. The left
side becomes
```
D x + d rev_p(x) + d b_(p−1)
  + (D·2^(p+1)+d)y + d(y mod2),
```
omitting the low terms when p=0. Thus the upper half y is
bounded by t/(2D), independent of p; for each possible y
and top low bit the remaining problem is an outside-in carry
on x and its reversal. A binary carry for target T in this
low equation (after subtracting the endpoint debit) moves by
c_(i+1)=(c_i+D b_i+d b_(p−1−i)−T_i)/2, bounded between 0
and D+d when integral. Start at c_0=0 and end at
c_p=floor(T/2^p); peeking from both ends gives a finite exact
search for the marked collision. It does not require a
small affine approximation to the original quadratic kernel.

## A full readout theorem on a literal TT3 pair

Take a single transitive tournament with masks (6,4,0).
Root at vertex 1 for (h,e,t,d)=(1,0,2,1) and at vertex 0
for (1,0,4,1). The y above must be zero. The top low bit
must be one, since otherwise 2x+rev_p(x) is too small.
The equation is exactly
```
2x+rev_p(x)=2^(p+1)−1.
```
Equivalently, the reversed p-bit word is the complement
of 2x in p+1 positions. The remaining p−1 bits of x must
pair with their opposite at reflected positions. A middle
bit cannot be its own complement. Thus p=2k+1 is necessary
and sufficient, with exactly 2^k possible shared addresses.
At p=0 the original equation has no solution.

More importantly these are **2^k different attained heights**,
all with an exact target inverse. Put L=2^k, take y=0,...,L−1
now as the *high-block index* (not the earlier vanished high
tail), and put a=rev_k(y). Then
```
x=L²+L(L−1−y)+a,
H_k(y)=8L²−6L+3+6a+y²−(6L−2)y−2K(y)−(y mod2).
```
For k=0 the single count is 5. To derive the formula for k≥1,
split the bits of x into a low block a, high block given by
the complement of y, and a top one. K within the low block
is K(y). For a complemented block
K(complement y)=L−2−y−a+(y mod2)+topbit(y)+K(y),
obtained by expanding each product of two bits. Cross
contributions between low and high blocks give y(L−1−y);
their adjacency correction vanishes (the boundary bits
are complementary). Cross contributions to the top one
give Ly+L−a−(a mod2). Substitution in 1+4x−K(x)
gives the stated height.

At an even y, increasing it by one increases a by L/2
and increases K(y) by y/2+(y/2 mod2); the exact drop in H
is 3L−2−y+2(y/2 mod2), at least 2L. At an odd y the
increase of the reversed address a is nonpositive, and
the change in K(y) is at least −(L−2). Thus the drop in
H is at least 4L−2y≥2L. These inequalities prove strict
ordering of the whole contact fibre, not merely an
empirical monotone chart.

Finally 2^p+3≤H≤2^(p+2)−3 by H=1+4x−K(x), x of
width p, and 0≤K≤2^p−2. Given an absolute positive N, the
two possible p values from its binary length have only
one odd member. Test that chart by binary search in y
and one exact debit comparison. On a hit, append the
two physical words x+2^p and x to the specified TT3
roots. Without enumerating tournaments, this uses
O((log N)^3) arithmetic operations with a direct
quadratic debit evaluator. The family realizes
Theta(sqrt(X)) distinct counts through X, because its disjoint
charts have 2^k heights of order 4^k. It cannot by itself
close Chan–Pak. The first charts contain 5 and 19,23;
a different literal pair at height 5 gives prime 53
with suffix addresses 3 and 2. The latter is not in
the TT3 chart: changing sources can add exact values.

**Next proof target.** Generalize the marked-cut inverse to
several independently routed contacts or change the two
reader states from the target. The single-cut law is a
selector that reads an absolute height and emits two
witnesses; it is not a surjectivity theorem. On the cube
side, the remaining bottleneck for the d=11 analytic
threshold is the unresolved contribution of other pair
orbits, not the finite existence problem.

Replay: `python3 cube_probes.py` (optionally pass the
[229-mask archive](../procgen_edim_20261001_q6_resolving15_orbits.txt))
and `python3 cut_decode.py`. Tests include direct literal
Hamiltonian counts and decoder comparisons.
