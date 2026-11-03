# A named native return at 233

**FINITE-EXACT native certificate; PROVED uniform word identity.** Local
labels R233.1--3 refer to this file. Same-agent analytic review, with
independent literal graph and subset-DP checks. This is a particular return,
not uniform ternary closure.

## R233.1. The actual source and endpoint

Use the inherited two-apex convention `old.source(w,u)`: a transitive base,
apices a,b with a->b, and words recording each apex's outgoing base edges.
The one-hole ready source p=2,q=1 from
[the previous ready grid](../lambda_emergent_20261103/READY_REPAIR.md) has

```
w=11110, u=11010,
rows=(30,28,88,16,96,79,11), (a,b,last)=(5,6,4), H=77.
```

Rows are outgoing bitmasks. Birth x=7 beaten only by a, followed by birth
y=8 beaten only by x, gives

```
rows=(30,28,88,16,96,207,11,351,127), root x=7,
(H,s_x,H(T-x),F_x)=(233,77,77,233).
```

Delete old vertex3, then reverse, in order, arcs
`8->6`, `4->5`, `7->4`, `5->7`. The successive Hamilton-path counts are

`233 -> 113 -> 199 -> 71 -> 209 -> 233`.

Each reversal has the two-native-complement implementation below. All old
graph endpoints are strong. The final root is still original name7; the
other marked vertex is original name6; original names5 and4 also survive.
Its profile is `(233,36,37,234)`, so it is ready. Its next literal two-birth
NR output has profile `(701,233,233,701)`.

The certificate uses eight native complements, four temporary helper
births/deletions, and one old-vertex deletion. Peak order is nine. It does
not keep H constant during the repair, and it does not assert that every
temporary complement graph is strong.

**Native edge-reversal proof.** Let P=(v,w,...) be a Hamilton path beginning
with the arc to reverse. Adjoin h which beats v and loses to all other old
vertices. Both P,h and, after the first native complement, v,h,w,... are
Hamilton paths. A native complement reverses precisely the edges off its
chosen path. Composing these two complements and deleting h reverses only
v--w among old edges. The two path edge sets differ there and at edges
incident to h. [The receipt](return_native.json) supplies every named path;
[the program](return_native.py) checks the path predicate before each move.

## R233.2. The exact minimum and the non-closing next step

The target bank is all all-one-dominant two-apex ready words of orders4--9,
with first other-apex letter1 and last0. Exactly one member has H233:
`w=111110,u=101100`, order8. Hold target root, other apex, and last base at
source names7,6,4; require original apex5 to survive; delete one old vertex.
There are exactly600 candidate maps. Their minimum old-arc discrepancy is
four, attained uniquely by target-to-source map

`(8,5,0,1,2,4,7,6)`.

This is a minimum in that explicitly finite assignment universe, not a
global lower bound on native moves or on other source banks.
The repair is a supplied-source program with original names5,6,4 and child
name7 retained in its instructions. It is not claimed to compile into the
two-live-port positive-birth language of JE1.

There is an all-word reason for the repair. For any normalized source
`w=1^(m-1)0`, other word u with first1,last0, NR creates x,y. Reversing
`y--b,a--last,a--x,x--last`, and deleting old base position j<last, gives
the transitive base order `(y,a,old base without j)` and new apex pair x,b.
Its new other-apex word is exactly

`u' = 10 delete_j(u)`.                                      (R1)

Check the four exceptional incidences: y becomes dominated by b; a becomes
dominated by last and x; last becomes dominated by x. All other incidences
already have the asserted base/apex form. This proves (R1) without a count
formula or an all-size search.

Within the finite one-hole range 1<=p<=9,0<=q<=9, every eligible j was
tested; only p2,q1,j3 preserves the desired new count 3H+2. At the returned
u=101100, the five next word repairs have counts

`571,601,673,673,377`, never701.

Thus the source return succeeds once and the same word-repair family is
not closed. [word_return.py](word_return.py) checks the identity on small
literal graphs and the complete stated finite bank with the exact word DP.

## R233.3. Clocks are transported, not assumed invariant

The full eighteen-polynomial endpoint observer and all five Q4 phase
presentations are in [return_native.json](return_native.json). The inherited
golden coordinates and phases change as follows:

| Stage | Golden (H,x,y) | Golden phase in F7^3 |
|---|---|---|
| NR child233 | (233,987,992) | (2,5,4) |
| Ready233 | (233,344,358) | (5,1,1) |
| Next child701 | (701,3634,3648) | (1,1,0) |

These are sourced charts, not a phase-preserving return. The native cost,
temporary versus permanent vertices, physical names, and all chart phases
are separate coordinates. The H35 ready obstruction remains valid; this
local result does not produce a ready representative for every required H.
