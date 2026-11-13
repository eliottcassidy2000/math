# Neutral arc components, strong components, and all weak 881 paddings

**PROVED:** the component factorization and equivariance implications below.
**FINITE-EXACT:** the two explicit 881 neighbour lists and literal graph
receipts. The all-order ready-881 source classification is inherited from
[FAMILIES](../../multicore_metagraph_20261113/FAMILIES.md), not inferred from
these two examples. Native complement-path edges, temporary births/deletions,
and simultaneous multi-arc replacements are different operation menus.

## 1. Cross-component reversal has an exact monotonicity obstruction

Use the proved odd-cycle collection formula
[THM-002, OCF](../../../../01-canon/theorems/THM-002-ocf.md): H(T) is the sum
of 2^|P| over collections P of vertex-disjoint directed odd cycles of length
at least3. An explicit primary statement is
[Irving–Omar, Corollary20](https://arxiv.org/html/2412.10572v3#S3.SS4),
attributed there to Grinberg–Stanley. The converse/complement convention gives
the same Hamiltonian path count by reversing each path. The formula uses
all odd cycles, not only triangles.

Let u->v be an arc between different strongly connected components (SCCs)
of T, and let T' reverse only that arc. No old directed cycle contains the
arc. Hence every old odd-cycle packing survives unchanged in T'. In the
linear SCC order, one of the following applies:

1. There is an intervening component: take any vertex w in it.
2. The component of u is nontrivial: take an out-neighbour w of u in it.
3. The component of v is nontrivial: take an in-neighbour w of v in it.

In each case u->w->v in T, so T' contains the new directed triangle
v->u->w->v. Its single-cycle packing is new. More strongly, adding this
triangle to any packing of T-{u,v,w} proves

```
H(T')-H(T) >= 2 H(T-{u,v,w}) >= 2.                       (N1)
```

The only remaining case has u and v in consecutive singleton SCCs. They
are homogeneous adjacent vertices; reversing their arc merely swaps their
roles, so H(T')=H(T). Thus:

**A cross-SCC reversal preserves H exactly when its endpoints are consecutive
singleton SCCs. Every other cross-SCC reversal strictly increases H.**

At the polynomial level, the same proof says the odd-cycle packing polynomial
increases coefficientwise, with a strict coefficient, in every nonexceptional
case. The H assertion follows by evaluating at2.

## 2. Exact factorization of the full single-arc neutral component

Let N(T) be the connected component of the labelled tournament T in the graph
whose edges reverse one arc and preserve the full H. Separate the SCC order
into its nontrivial strong blocks C_i and maximal consecutive runs of singleton
components, of lengths a_j. These runs induce transitive tournaments.

An H-neutral flip cannot split a strong component: if it did, reversing it
back would be a cross-SCC reversal in the resulting graph. Section1 would
strictly increase H unless it merely swapped consecutive singleton components;
that exceptional swap cannot merge them into a strong component either.

An internal arc flip in C_i changes the full H by the product of the other
positive component counts times its local change, since all Hamilton paths
concatenate the SCCs in order. Thus it is globally neutral exactly when it
is locally neutral. Inside a singleton run, only consecutive-position flips
are neutral; they are precisely the adjacent transpositions of its total
order. All other cross-block flips are excluded by Section1. Consequently

```
N(T) = Cartesian product_i N(C_i)
       × Cartesian product_j Perm(a_j),                (N2)
```

where Perm(a) is the graph on the a! labelled total orders, with adjacent
transpositions as edges. This is an exact labelled graph isomorphism: an edge
changes one factor, every allowed factor move lifts, and no other neutral
edge can leave the product. In particular the strong-block vertex sets and
their order are invariant throughout a neutral route. This is a classification
of how the strong-component pieces compose; it does not classify N(C) for
every strong C.

## 3. From one exact neighbour list to an entire component

Suppose every H-neutral one-arc neighbour of a tournament A reverses a
homogeneous pair. By [THM-4463](../../../../01-canon/theorems/THM-4463-homogeneous-tournament-contacts-and-pinned-scale-floor.md),
its homogeneous-pair graph is a disjoint union of ordered paths, whose vertex
sets are its maximal transitive modules M_i.

Reversing an adjacent homogeneous pair is exactly conjugation of A by the
transposition of those two names. Both H-neutrality and homogeneity are
equivariant under every relabelling. Thus the assertion that every neutral
edge is homogeneous remains true at every graph reached this way. Induction
on route length proves that no nonhomogeneous neutral edge can first appear
later. Adjacent swaps generate every order on each M_i, and preserve its
outside incidences. Therefore the full component has exactly

```
product_i |M_i|!
```

labelled graphs, obtained by independent internal module reorderings. If the
named root is in a singleton module, every such conjugation fixes it, so its
entire marked state (h,s,q,F) is constant on the full component. If the root
belongs to a larger module, the response law in PROOF.md applies instead;
root invariance is not asserted in that case.

This is an equivariance theorem following an exact local audit, not a claim
that all neutral arcs are homogeneous in arbitrary tournaments. The H5
four-vertex counterexample in PROOF.md shows that stronger statement is false.

## 4. The two 881 components and their arbitrary transitive tails

The independently recounted source rows from the selector receipt have exactly
the following H881-neutral one-arc neighbours:

| Source | Named root | Neutral arcs | Full component size | Constant root state |
|---|---:|---|---:|---|
| A, order10 | 9 | {4,5}, {7,8} | 4 | (881,74,75,882) |
| B, order11 | 10 | {5,6} | 2 | (881,218,219,882) |

All listed pairs are homogeneous, and neither root belongs to a nontrivial
pair component. Section3 proves the all-route conclusion from these finite
lists; a separate exhaustive traversal directly returns the same4 and2 graphs.
No arbitrarily long route consisting only of H-neutral single-arc flips at
these fixed orders can alter their marked response.

Every weak unit-ready881 source is one of these strong rooted sources followed
by a dominated transitive tail TT_a, by the inherited complete source theorem.
Section2 then makes its full component the corresponding4- or2-vertex component
times Perm(a). Its size is4 a! or2 a!, and the named-root state remains constant.
This proves the same obstruction for every finite tail length, not merely a
bounded padded census. It does not prohibit temporary H changes, new helpers,
deletions, or a simultaneous endpoint-neutral replacement.

## 5. Tail padding preserves the full rooted state and commutes with the repair

For any nonempty rooted tournament (T,r), append a transitive tail dominated
by every vertex of T. The full state (h,s,q,F) is unchanged:

- Hamilton paths and deletion paths concatenate their T-part and the unique
  ordered tail; the same is true for paths starting at r.
- In a cover whose first path ends at r, no tail vertex can occur in that
  first path, since it cannot reach T. Hence the entire tail lies in the
  second path, after any remaining T vertices. Deleting/appending it is a
  bijection, including an otherwise empty second path.

The selector's physical role exchange between the root and vertex3 touches
only core arcs: both old roles beat every tail vertex. Every R or N birth at
a core root also beats every tail vertex, so it commutes with appending the
tail. Hence the selector's full preparation states, formulas, gcd conditions,
and balanced rooted outputs hold for **all weak ready881 paddings as well**.
These padded outputs remain weak; no strong-output assertion is transferred.

The separately cited native controller is stated for arbitrary nonempty
rooted sources, and the resulting padded graphs carry the same retained
root-child interface and balanced state. Its scope, rather than a new
strongness assumption, governs renewal after entry. The initial endpoint
replacement still has no trace using only H-neutral single-arc moves;
Sections3–4 prove that such a trace is impossible in this restricted menu.

## 6. Independent selector proof audit and reproduction

The return R^k N R formula and the two explicit exponent sieves in
selector/PROOF.md are correct as stated. For A, t=782*2^k+125 is odd and
avoids the primes5 and7 dividing2D=140; native entry occurs exactly at odd k.
For B, t=912*2^k+51 is odd, is divisible by3, and meets the single bad
condition2^k=7 mod11 precisely at k=7 mod10. The corresponding b values
are2346*2^k+515 and2736*2^k+505. Both role exchanges preserve full H881
at their endpoints. This does not give the missing881->1691 ready return.

Run from this directory:

```
python -B neutral_scc.py
python -O -B neutral_scc.py
```

The script imports no selector implementation. On all1,099 labelled
tournaments through order5 it checks3,212 cross-SCC flips:2,566 strictly
increase H and646 are the adjacent-singleton exceptions. All2,398 directed
neutral adjacencies preserve the proved strong-block/singleton-run skeleton.
It separately enumerates both entire fixed-order881 components, recounts
all21 roots and both prepared endpoints, and recounts all six literal return
graphs for k=1,2,3. Explicit exponent predicates are checked for k=1..200.
The infinite arithmetic and padding claims follow from their proofs, not
those finite controls.
