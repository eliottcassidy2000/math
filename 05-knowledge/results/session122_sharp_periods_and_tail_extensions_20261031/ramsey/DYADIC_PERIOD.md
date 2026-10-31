# Exact sharpness lattice and least eventual period

2026-10-31. **PROVED, scoped.** No unrestricted Ramsey bound or novelty claim.

Let D be generated from the singleton by binary join, disjoint union, and
full C5 substitution with all five children nonempty. Put F(A,B)=max |G|
for G in D with omega(G)<=A and alpha(G)<=B, for positive integers A,B.
Let w(0)=0,w(1)=1, w(2h)=5w(h)/2 and
w(2h+1)=3w(h)/2+w(h+1). Put L(A)=2^ceil(log2 A), with L(1)=1.

Dependencies: G2 in `ramsey_neighbourhood_lifts_20261030/GRAPH_CAPACITY.md`
proves w superadditive and increasing, F(A,B)<=w(A)B, and sharpness when
L(A)|B. Session121 `ramsey/EVENTUAL_TRANSLATION.md` proves existence and
effectivity of an eventual translation law with period L(A).

## Theorem 1: exact sharpness lattice, for all budgets

For every positive A,B,

    F(A,B)=w(A) B  if and only if  L(A) divides B.

This is stronger than scalar integrality. For example w(3)=4 is integral,
but equality is possible only when 4 divides B. Also w(6)=10 is integral,
but equality is possible only when 8 divides B.

### Strict-superadditivity lemma

For positive x,y, w(x)+w(y)<w(x+y).

Proof by strong induction on n=x+y, using inherited non-strict
superadditivity (which also gives w(t+1)>=w(t)+1).
If x=2u,y=2v are even, apply the strict induction hypothesis to u,v and
multiply by 5/2. If x=2u+1,y=2v+1, write h=u+v+1. Twice the left side is

    2[w(u+1)+w(v)] + 2[w(u)+w(v+1)] + [w(u)+w(v)].

The first two bracket types are <=w(h), while the last is <=w(h-1)<w(h).
The total is therefore <5w(h)=2w(2h), including u=v=0.
If x=2u,y=2v+1, then u>=1 and h=u+v. The left side is

    (3/2)[w(u)+w(v)] + [w(u)+w(v+1)].

The first bracket is <=w(h); the second is strictly <w(h+1) by induction,
since both arguments are positive and h+1<2h+1=n. This proves the lemma.

### Proof of Theorem 1

Sufficiency is inherited. Necessity is by induction on A, with the graph's
construction tree handled by structural induction for that A. A=1 is trivial.
For A>1 a singleton is strictly below the bound. A join of nonempty children
with actual clique numbers x,y has order at most
[w(x)+w(y)]B<w(x+y)B<=w(A)B, so cannot be sharp. At a disjoint union, the
actual independence numbers B_j add to at most B. Equality requires every
component to be sharp at budget A and its actual B_j; structural induction
then gives L(A)|B_j and hence L(A)|B. (Equality also forces their sum to be B.)

At a C5 root let a_i,b_i be the actual child clique and independence numbers.
All are positive. The a_i satisfy adjacent constraints a_i+a_(i+1)<=A,
and b_i satisfy nonadjacent constraints b_i+b_(i+2)<=B. Set c_i=w(a_i).
The five-cycle edge polytope for b/B has vertices the independent-set
indicators of the constraint cycle and the all-half vector. This is exactly
the elementary polytope used in inherited G2.

Every indicator objective is strictly <w(A): a two-element support is an
adjacent pair of the original C5, so strict superadditivity applies; a
one-element support has a_i<A because the neighboring child is nonempty.
The empty indicator has objective zero. The all-half objective is
(1/2)sum c_i<=w(A) by G2. Thus attaining order w(A)B forces uniquely

    b_i=B/2 for every i;   sum_i w(a_i)=2w(A),
    |G_i|=w(a_i)b_i for every child i.                       (1)

In particular B is even. For A=2 this is the required condition. For A>=3,
write k=ceil(log2 A), so k>=2 and A>2^(k-1). If every a_i<=2^(k-2), then

    sum_i w(a_i) <= 5w(2^(k-2))
                 = 2w(2^(k-1)) < 2w(A),

contradicting (1). Hence some a_i>2^(k-2). Each a_i<A, so the outer
induction hypothesis applies to its sharp child. It gives L(a_i)|B/2,
and L(a_i)>=2^(k-1). Because L(a_i) is a power of two, 2^k divides B.
This completes necessity. Equality also forces actual omega=A and alpha=B,
by strict monotonicity of w and positivity of w(A).

## Theorem 2: the least eventual period is exactly L(A)

For every A the least positive eventual translation period of F(A,.) is

    d_A = 2^ceil(log2 A).

Indeed session121 gives the eventual L(A)-translation. Any eventual
translation with step d has increment w(A)d (by the limiting slope), hence
preserves the deficit w(A)B-F(A,B). Choose arbitrarily large multiples B of
L(A), where the deficit is zero. Then B+d also has zero deficit, so Theorem1
forces L(A)|(B+d), hence L(A)|d. The sufficient period L(A) is therefore
minimal. In particular the scalar denominator can be a proper divisor of
the true period, even denominator one; no denominator-only criterion works.

## Theorem 3: complete recursive equality criterion

Within the specified construction grammar, every sharp graph and only those
sharp graphs are described by these rules (with actual A,B):

* A=1: an edgeless graph on B vertices.
* For A>1, a disjoint union whose nonempty constituents are all sharp at the
  same clique budget A and whose independence numbers sum to B.
* For A>1, a full C5 substitution satisfying (1), the adjacent a-constraints,
  and recursively sharp children. All children have independence exactly B/2.

No nonempty join is sharp for A>1. This is a necessary-and-sufficient
recursive criterion for construction trees, not an isomorphism census or a
claim of unique decomposition. The equality clique-weight patterns remain
finite tests at each fixed A; no classification of their isomorphism orbits
is being inferred.

## Remaining boundary

The full period is now uniform and minimal, and the zero-deficit residue is
completely classified. General nonzero residue penalties and useful numerical
stabilization thresholds for all A are not determined by this theorem.
