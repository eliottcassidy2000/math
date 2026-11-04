# A common-state all-strong second fork with two permanent births per step

**PROVED actual construction / FINITE-EXACT programs — [independent audit](DIRECT_CHILDREN_AUDIT.md) passed.** This
extends the frozen DIRECT programs and leaves every earlier bundle intact.
It constructs a genuine two-level binary tree on retained physical names,
with a common actual state at each first child. Every native intermediate
and every completed output is strong. Only two of the four final leaves
are already positively prepared; no unbounded closure is asserted.

## 1. The actual common first children

Use the original eleven-name H1429 graph and root8 from DIRECT. Its first
two branches are chosen once and for all:

* DIRECT A gives the thirteen-name `(H,F_8)=(4289,7653)` target, gap2 at12.
* DIRECT B gives the thirteen-name `(H,F_8)=(4291,6403)` target, gap1 at12.

Thus the first+4 branch is the same physical DIRECT B graph for both future
choices. We do not select a different first-child graph after seeing the
second digit. Each first branch spends exactly the two permanent names of
the common zero C2 prefix and retains all eleven old names.

## 2. The two continuations from the actual 4289 child

Apply C2 to its actual auxiliary12. Both births are strong and preserve
root8. They add names13,14 and give the common fifteen-name graph

```
(32510,32388,31872,32758,32390,32662,32694,32000,28183,28804,
 31232,29184,16640,4096,8192),
```

with `(H,F_8)=(12871,22959)`. This full response identifies it: it is not
the earlier weak H12871 child with a different F.

For the next digit+2, reverse in order

```
(2,8),(3,8),(4,8),(8,14),(9,10),(9,11),(11,12),(13,14).
```

The literal completed target is

```
(32510,32388,32128,32502,32646,32662,32694,32000,11787,31876,
 30720,24576,18688,20480,256).
```

This repairs H downward by2, giving12869=`3·4289+2`. Auxiliary14 has
root8 as its sole exit. Its exact response is

```
(H,F,S,L,K,G,sigma,tau,b)
=(12869,16599,12870,4311,12288,16599,1,0,3730).
```

For the next digit+4, instead reverse

```
(1,8),(2,10),(4,5),(7,10),(7,13),(8,11),(9,11),(12,14).
```

The literal target is

```
(32510,32644,30848,32758,32422,32646,32694,22784,26133,30852,
 31364,28928,256,4224,12288).
```

Here H remains12871=`3·4289+4`, while auxiliary12 becomes positively
prepared at root8:

```
(H,F,S,L,K,G,sigma,tau,b)
=(12871,22291,12873,5907,16384,22291,2,0,9420).
```

Both targets are strong. They have independent ready-language descriptions:
C3, entry{0}, tail1001100011 for the first; and C4 `(2,4,9,3)`, entries
{0,1}, tail111001001 for the second. The source reconstructs these actual
matrices with explicit physical permutations and separately verifies the
full named native programs. The joint identity G=F follows from the sole
exit and is independently checked through L and K.

## 3. A same-H rearming from the common actual 4291 child

For its next+2 branch, use DIRECT B's existing gap1 immediately. C2 gives
the strong fifteen-name child `(H,F_8)=(12875,19209)`.

For its next+4 branch, first change the **same DIRECT B graph** to the
literal DIRECT C graph by the following count-neutral old-pair sequence:

```
(2,7),(2,9),(2,12),(4,8),(6,8),(7,9),(8,11),(8,12),(9,11),(11,12).
```

H stays4291 at the endpoints, while the actual marked response changes

```
(F_8,auxiliary,sigma,tau)=(6403,12,1,0) -> (8371,11,2,0).
```

No permanent name is added or deleted in this switch. It is implemented
by twenty guarded Hamilton complements and ten temporary helpers; all
thirty native states and all completed adjacency stages are strong.
Then C2 at11 gives the strong fifteen-name child
`(H,F_8)=(12877,25113)`.

The root response is carried through the intermediate same-H state; it is
not replaced by an H-only equivalence. This is the sidecar required to
compose the selected second branch correctly.

## 4. The full actual four-leaf tree and its boundary

Together with DIRECT, the programs give

```
                       H1429, F1563, order11
                     /                       \
       H4289,F7653,order13             H4291,F6403,order13
             /       \                      /       \
        H12869     H12871                H12875     H12877
        F16599     F22291                F19209     F25113
        order15    order15               order15    order15.
```

All numbers refer to old physical root8. Every old physical name survives,
and each count step uses exactly two permanent births. Temporary helpers
are distinctly named and are the only deleted vertices. The two fifteen-
name repairs use sixteen helper-stage vertices at peak; the same-H
thirteen-name switch peaks at14. Thus the whole tree has peak order16.

The three newly checked repair programs have26 old-pair steps,52 guarded
Hamilton complements,26 temporary helper deletions, and78 native states.
Every one is strong. The source uses the independently audited pair-block
compiler and a separate named-set interpreter; it records every matrix,
guard, helper cut, and `(H,F_8,F_8−H)` along the programs. Direct forward
DP independently agrees with the backward observer implementation at the
shared sources and endpoints. Both modes reproduce the frozen receipt; the independent audit also recounts every native observer.

The12869 and12871 final graphs carry positive joint guards1 and2 as
recorded above. The12875 and12877 graphs are ordinary C2 images of strong
thirteen-name inputs; the inherited every-auxiliary C2 obstruction applies
to them. They are **not** already positively prepared. Therefore the result
is an all-strong four-leaf tree, not a four-leaf tree closed under all future
nonzero steps.

This escapes ORDER_OBSTRUCTION because the operations have changed. The
first nonzero digit is a zero C2 followed by count-changing repair; on the
4289 side the second digit uses C2 followed by repair. These are outside
the narrower H-neutral-preparation-before-C2 grammar. The retained-order
obstruction remains valid and explains why the earlier operation failed.

No further search or third-level tree is part of this bundle. The remaining
unbounded obligation is a realizable returning class with explicit full
port-state transport and a selector valid at every subsequent step.

## 5. Reproduction

```
python3 04-computation/nonzero_preparation_20261103/surgery_direct_children.py --output /tmp/direct-children.json
python3 -O 04-computation/nonzero_preparation_20261103/surgery_direct_children.py --output /tmp/direct-children-O.json
```

`DIRECT_CHILDREN_MANIFEST.json` pins the source, proof and complete receipt.
The source verifies the frozen DIRECT supplier hash and constructs the
actual common first children rather than substituting other graphs with
the same H. The independent analytic/source/replay audit passed.
