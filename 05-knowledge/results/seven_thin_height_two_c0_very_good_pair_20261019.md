# Seven-thin height-two c=0 roots force a very good pair

Date: 2026-10-19. Status: a structural theorem for every finite height-at-most-two poset with a c=0 critical pair and maximum incomparability degree at most7. This is separate from the eleven-thin Fano obstruction. No conclusion for higher-height posets or ordered slave sets is asserted.

## Full height-two statement and same-level roots

Let P have height at most two (no three-element chain), have maximum incomparability degree at most7, and contain an ordered critical pair (x,y): x,y are incomparable, D(x)⊆D(y), U(y)⊆U(x). Assume no third vertex is incomparable with both roots. Then P has a very good pair, and therefore has weak Gold.

If both roots are minimal, there can be no other minimal vertex: it would be incomparable with both. Every other vertex is consequently maximal, and criticality implies it is above x. Its predecessor set is either {x} or {x,y}. The vertices with predecessor set {x} are the exclusive upper neighbors. If there are at least two, they are exact twins. If there is exactly one, x,y have equal downsets and their upsets differ by a singleton, so x,y are very good. If there are none, x,y are exact twins. Arbitrarily many common upper vertices are allowed. The case of two maximal roots is the order dual, allowing common lower vertices.

These cases and bottom/top roots exhaust the possibilities. Indeed, any predecessor of x is also a predecessor of y, which in height two makes both roots maximal. Any successor of y is also a successor of x, which makes both roots minimal. Otherwise x is minimal and y is maximal. For bottom/top roots no common lower or common upper vertex can exist: either would form a three-element chain. The c=0 condition then forces precisely the matrix model below.

## Bottom/top matrix class and conclusion

Let P have vertex set {x,y}∪L∪R, with |L|=l, |R|=u. Its only possible strict relations are

    b<y for every b∈L,
    x<a for every a∈R,
    selected relations b<a with b∈L,a∈R.

In particular x,y are incomparable, L and R are antichains, and no other vertices occur. This is the height-two c=0 critical-root model with x minimal and y maximal. Assume **every vertex**, not merely x,y, has incomparability degree at most7.

Then P has a very good pair in Zaguia's sense: in P or its dual, an incomparable pair has equal strict downsets and both strict upset differences are chains. Consequently P has a balanced pair by Zaguia, and has weak Gold by the attributed good-pair synthesis of packet94. [Definitions1,5 and Theorem2](https://arxiv.org/pdf/1610.00809v3); [local Gold synthesis](good_pair_gold_synthesis_20261019.md).

## Missing-incidence formulation

For b∈L define its zero set Zb={a∈R:b is incomparable with a}. Define Za⊆L dually for a∈R. Root degrees give l,u≤6. The full vertex degrees are exactly

    deg_inc(b)=l+|Zb|,
    deg_inc(a)=u+|Za|.

Thus |Zb|≤7−l and |Za|≤7−u. These include incomparability with x or y in the same level.

If some row has |Zb|≤1, the pair (x,b) has equal empty downsets and upset differences Zb and {y}; both are chains. It is very good. If some column has at most one zero, the dual pair (a,y) is very good.

For b,b'∈L, their equal downsets are empty, and their upset differences have sizes |Zb'\Zb| and |Zb\Zb'|. Therefore two rows whose zero sets differ by at most one element in each direction give a very good pair. The corresponding statement holds for columns by duality.

It remains to rule out the possibility that all row and column zero sets have size at least2 and no such pair exists.

## Small set-system argument

If l=0 or u=0 and the other side is nonempty, a column or row has zero zeros, giving the preceding witness. If l=u=0, x,y are exact twins. Otherwise, assuming all zero sets have size at least2 forces l,u≥2.

By duality take l≥u; then 2≤u≤l≤6.

- If l=6, the degree bound gives |Zb|≤1, contradiction.
- If l=5, every Zb has size exactly2. Any two intersecting zero sets differ by at most one element each way. If there were no such pair, the five sets would be pairwise disjoint, requiring ten distinct elements in R, whereas u≤5. Contradiction.
- If l=4, there are four subsets of a universe of at most four elements, each of size at least2. If there is no near pair, all size2 sets must be pairwise disjoint, so there are at most two. There is at most one size3 set, since two size3 subsets of a four-element universe differ by at most one element each way. There is at most one size4 set, and its presence excludes any size3 set for the same reason. Thus at most three sets are possible, contradiction. Repeated subsets already form a near pair.
- If l≤3, every zero set is a subset of a universe of at most three elements and has size at least2. Any two therefore differ by at most one element each way. Since l≥2, a witness exists.

This exhausts all cases and proves the theorem.

## Scope and interpretation

The matrix argument uses actual seven-thin budgets to force a very good pair; it does not derive thinness from root degrees. Combined with the same-level cases, it handles every height-at-most-two poset with a c=0 critical pair. It handles all incidence matrices without enumerating them. It does not cover the six-point mixed(2,2) core with internal chains b1<b2 and a1<a2, nor its actual exterior realization: those are not this height-two model.

A general c=0 theorem for arbitrary exterior A,B and arbitrary slave orders remains open in this work. No historical priority claim is made for this small combinatorial consequence. The Fano construction explains why dropping the all-vertex seven-thin bound can invalidate the structural conclusion.
