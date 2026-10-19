# Two disconnected color contacts still admit width-two repair

## Theorem

Let F be a fixed rainbow forest with r+1 components and let r free colors remain. If at most two free-color contact graphs on the contracted components of F are disconnected, every two rainbow spanning-tree completions are joined by at-most-two-color moves. The diameter is at most r+1. If at most one contact graph is disconnected, the earlier bound r applies.

As before, these operations preserve labels, colors and F but can change the abstract tree shape. No one-color GT1 or Graceful Tree Conjecture claim follows.

## 1. Reconfiguring two partial color choices

Temporarily ignore all connected colors, and consider only the two possibly disconnected colors a,b. Their chosen edges, together with F, must be acyclic. On the quotient by F each chosen edge has a nonloop unordered component-pair type. Two such edges are acyclic precisely when their types are distinct: with only two edges the only possible cycle is a parallel pair.

Let A,B be the available type sets for colors a,b. A partial state is an ordered unequal pair(x,y)∈A×B. For a target(u,v):

- If u≠y, first change x to u, then y to v.
- If u=y but v≠x, first change y to v, then x to u.
- Otherwise the target reverses the two types: (u,v)=(y,x). If A contains a third type w, use(x,y)→(w,y)→(w,x)→(y,x). If B has a third type, use the symmetric sequence.
- If neither set has a third type, both sets are exactly{x,y}. Replace both color choices at once. Their uncolored quotient forest is unchanged.

Thus at most three one-coordinate partial changes suffice, except for a simultaneous two-coordinate swap which takes one move. Actual target edge realizations can be used at the corresponding target-type changes. When source and target types already agree but their original endpoints differ inside F components, one-coordinate replacements realize the exact edges. The same two/three-step bounds still apply.

## 2. Lifting a partial change to a full completion

Suppose a single bad-color edge e is replaced by g of that color, and the new two-bad-edge selection remains acyclic over F. Insert g into the current completion, obtaining a fundamental cycle C.

If the old same-color edge e lies on C, replace it directly. Otherwise C must contain an edge f of a connected color: if it did not, C would lie entirely in F, the other bad edge, and g, contradicting acyclicity of the proposed partial selection.

Delete e and f, add g, and reconnect the remaining cut with an edge of f's color, using connectedness of its contact graph. This is one width-two move. The other bad edge and all of F remain fixed. Thus every one-coordinate partial change lifts to a single allowed move.

For the exceptional simultaneous type swap, replace both bad edges directly with their target realizations. The full quotient graph by F has the same uncolored edges before and after, so it remains a tree. This is a valid width-two move without any connected-color repair.

After at most three moves, the two bad edges agree exactly with the target completion. Now install each of the other r−2 target edges by the monotone connected-contact argument: a non-target fundamental-cycle edge cannot have a bad color because both bad edges are already target edges. At most r−2 further moves suffice. Total≤3+(r−2)=r+1.

## 3. Four free colors in every integer-difference upper-tail fiber

For every m≥4, freeze any acyclic selection of colors5,…,m in the difference-colored graph on0,…,m. The completion fiber with free colors1,2,3,4 is connected under width-two moves, with diameter at most5.

For m≥5, Q1 is connected through the adjacent-label path and Q2 is connected through the two parity chains plus the frozen odd color-five edge. Only Q3,Q4 can be disconnected, so the theorem applies with r=4. At m=4 the unique empty upper tail has12 graceful completions and width-two diameter2, already verified in the previous boundary check.

The m9 tail F={(0,5),(0,6),(2,9),(0,8),(0,9)} has both Q3 and Q4 disconnected and is a genuine non-boundary use of the new theorem. Its label complement gives the other residual tail at that size. This is a proof for every m, not an extrapolation from finite tails.

## Scope

This argument relies on the complete elementary classification of an acyclic two-edge partial selection. Three disconnected colors would require analyzing three-edge rainbow partial forests, where cycles of length three also arise. The present theorem makes no assertion for that case. General upper-tail fibers with five or more free colors remain outside the proved all-size conclusion.

## Exact replay and inheritance

At m=4, color4 has the unique possible edge{0,4}, so treating it as free or frozen gives precisely the same12 completions as the [previous boundary check](../graceful_connected_color_contacts_20261019/CONNECTED_COLOR_CONTACTS.md).

Run `python verify_two_bad_contacts.py` and `python -O verify_two_bad_contacts.py`. Their identical receipts contain1498891 checks and155470 individually validated moves. The controls exhaust every ordered pair in the five actual two-bad upper tails through m9:57680 source/target pairs across496 completions. They explicitly check actual edge realizations, the bad-phase bound≤3, protected edges, every intermediate tree, and total bound≤5. The universal theorem follows from the proof, not this finite control set.

The [one-bad-color criterion](../graceful_connected_color_contacts_20261019/CONNECTED_COLOR_CONTACTS.md) supplies the monotone good-color phase. Its source/target cut argument is also implemented independently in this packet's self-contained standard-library checker. No solver, downloaded code, or mutable state is needed.
