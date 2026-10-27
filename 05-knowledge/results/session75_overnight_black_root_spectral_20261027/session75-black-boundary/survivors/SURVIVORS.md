# Two surviving exact-distance statements

**Session75, 2026-10-27. PROOF CANDIDATE for independent audit.**
This is an additional, separate packet; the already frozen nine-file
`/tmp/session75-black-cut-proof` packet has not been edited.

## 1. A unique unbalanced minimum shore gives the lambda equality

**Theorem.** Let `T` be strong, and suppose its minimum outgoing cut has a
unique shore `A`, of cardinality different from half the order. Then its
black path-reversal distance to nonstrong is exactly `lambda(T)`.

**Proof.** Use the inherited all-colour construction, always choosing a
Hamilton path beginning in `A` and ending outside `A`. At step `j`, the cut
of `A` equals `lambda-j`. Every other initial cut was at least `lambda+1`,
and can decrease by at most `j`, so every other current cut is at least
`lambda+1-j`. Thus `A` stays the unique minimum shore. Before step `lambda`,
all cuts remain positive, so a Hamilton cycle supplies the next directed
path with the required endpoints.

At any such intermediate tournament, an anti-automorphism `theta` would
carry the unique minimum shore `A` to the equally minimizing outgoing shore
`V\theta(A)`. Uniqueness would force `A=V\theta(A)`, hence `2|A|=|V|`, contrary
to the hypothesis. Therefore every intermediate tournament is non-self-
converse, and every directed Hamilton word chosen there is black. After
`lambda` moves the cut of `A` is zero. The general cut-current lower bound
proves equality. No predetermined-path blackness is inferred merely from
its endpoint scores. QED.

This repair was suggested by the parent and checked independently here.
It explains both features of the new obstruction: the even family has one
balanced minimum shore; the odd family needs two minimum shores. The theorem
is not an all-order equality when those hypotheses are absent.

## 2. The obstruction families nevertheless attain the core bound

Use all labels and definitions of the frozen packets
`ROOTED_CONVERSE_COUNTEREXAMPLE.md` and `ODD_MIDDLE_EXTENSION.md`.
For `a>=1`, write `m=a+4`, `r=a+3`, with `V` labels `0,...,m-1`, `U` labels
`m,...,2m-1`, the cyclic `B` labels `(a,a+1,a+2)`, and middle vertex `w=2m`
only in the odd case.

**Theorem.** Every even or odd member of those infinite families has black
distance2 to a singleton plus a strong core. Thus all have

`d_black(nonstrong)=d_black(singleton + strong core)=delta=2 > lambda=1`.

This is a statement about these families, not the unresolved all-order
core-distance equality.

Keep the first black path `P` from the frozen packet (insert `w` in its
previous location for odd order). Replace the second path by

`Q_core_even = (r,a+1, 2m-1,2m-2,...,m, 0,1,...,a-1,a,a+2)`.

At odd order insert `w` between `a+1` and `2m-1`:

`Q_core_odd = (r,a+1,w, 2m-1,2m-2,...,m, 0,1,...,a-1,a,a+2)`.

Both paths end at the same named vertex `v=a+2` in `V`, as does the first
path `P`. The initial indegree of `v` is2, so these two steps make `v` a
source. Its minimum-semidegree lower bound is therefore tight once we prove
directedness, blackness, and strongness of the remaining core.

### Directedness and actual blackness

After `P`, the induced blocks are the same `U',V'` described in the frozen
proof. Inside `V'` we have `r -> (a+1)`; the latter vertex beats `r_U` across
the blocks in even order. In odd order it beats `w`, which beats `r_U`.
The complete reversed-alpha word traverses `U'` ending at `m`; its current
cross arc is `m -> 0`. The increasing `A` word in `V'` then continues through
`a -> (a+2)`. These are exactly the listed arcs of `Q_core`.

At the intermediate tournament, `v` has score `n-2`. The new start `r_V` has
score `a+6` in even order and `a+7` in odd order: it was interior in the first
path, so its full score did not change. Their sum is strictly greater than
`n-1`. Thus `Q_core` is black by a valid necessary condition for the actual
anti-automorphism predicate. The checker also tests the complete predicate.
The first path is black and its intermediate strong by the frozen proof.

### The actual terminal core is strong

Both paths reverse the same full alpha path inside `U`, in opposite orders.
Hence the terminal induced `U` is restored to the original strong `U_a`.

The induced tournament on `V\{v}` has these exact arcs:

- every vertex in its nonempty `A` block beats `r`;
- `r -> a -> (a+1) -> r` is a directed triangle;
- both `a` and `a+1` beat every vertex of `A`.

Therefore it is strong: the triangle is strong, and each `A` vertex is reached
from `a` and reaches `r`. The internal order on `A` is irrelevant here.

At even order the only remaining cross arc from `U` to this `V` core is

`r_U -> (a+1)_V`,

while all other cross arcs point `V_core -> U`. Thus the union of the two
strong blocks is strong.

At odd order all direct cross arcs point `V_core -> U`, and the middle vertex
has orientation `V_core -> w -> U` except for

`r_U -> w -> (a+1)_V`.

These two exceptional arcs join the strong blocks and `w` into a directed
cycle of parts. Hence the entire deleted-source core is strong here as well.

The terminal tournament is exactly `v -> strong_core`, not merely nonstrong
or a source plus several SCCs. The endpoint-current lower bound is
`delta=2`, already established for every family member, so its black core
distance is exactly2. QED.

## 3. Portable exact controls

Run from this directory:

```sh
python3 -B check_survivors.py > normal.json
python3 -B -O check_survivors.py > optimized.json
cmp normal.json optimized.json
```

The script checks both parities for `a=1,...,20` (forty cases, orders10–49).
It retains the full named `R` and literal `K` routes, checks directedness and
actual blackness at both ends, verifies the intermediate is strong, verifies
that the named target is a source, and checks strongness of the actual
remaining core. It separately verifies restoration of `U` and strongness of
the induced `V` minus its source. The co-delivered helper files in `inputs/`
are exact copies of the frozen family checkers. All reject conditions are
explicit exceptions, and the normal/optimized receipts are identical.

These finite controls support the all-order structural proof; no finite
census is substituted for its arbitrary-order quantifier.
