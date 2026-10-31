# A concrete meeting point: 60, divisor balance, and icosahedral rotations

**Classical / inherited mathematics, with exact controls; no novelty claim.**
This is a bounded addendum to `README.md`, not another proposed five-way
equivalence. The original lane's files are unchanged. Definitions and the
balanced-integer classification are recorded in the frozen
`dependencies/balance_dynamics/README.md` and its pinned original sources;
the all-time PRIMEGAME conclusion comes from section6a of that same report.
The canonical macro proof is frozen in
`dependencies/fractran/FRACTRAN_PRIME_EXPONENT_REPORT.md`.
These dependency reports are included unchanged; their references to sibling
scripts describe their original audited packets, not required local imports.

## 1. A real shared value, with its different meanings retained

For `N>=2`, use proper nontrivial divisors `D(N)={d:1<d<N,d|N}`. The counts
`F=|D(N)|`, `U=#prime elements of D(N)`, and `S=#squarefree elements of D(N)`
satisfy `F=U+S` exactly for the inherited shapes `p,p^3,p^2qr`, with the primes
distinct in the mixed shape. The prime subset lies inside the squarefree
subset: this is a count identity, **not** a disjoint partition.

The least mixed balanced integer is

`60 = 2^2*3*5 = |A5|`.

Indeed the product of three distinct primes is at least30 and the repeated
prime is at least2. Equality is attained by2,3,5 with2 repeated. The same60
counts the rotations of an icosahedron or dodecahedron, as already sourced in
the main note's [A5 geometry construction](https://math.ucr.edu/home/baez/mathematical/six.bac.html).
This is an exact **meeting of two constructions at one value**. Their defining
operations remain distinct.

| Rotational group | Order | F | U | S | F−U−S |
|---|---:|---:|---:|---:|---:|
| tetrahedral `A4` |12|4|2|3|−1|
| cube/octahedron `S4` |24|6|2|3|+1|
| icosahedron/dodecahedron `A5` |60|10|3|7|0|

Thus only the icosahedral **rotational** order is balanced. The full symmetry
orders are24,48,120, whose defects are respectively1,3,4; none is balanced.
In particular the full icosahedral group is `A5 x C2`, not `S5`, even though
both have order120. Their centres have different sizes. The word “rotational”
cannot be silently dropped from the60 claim.

## 2. Ten proper divisors, seven realized subgroup orders—but not the same seven

At60,

```
D(60) = {2,3,4,5,6,10,12,15,20,30},
squarefree elements of D(60) = {2,3,5,6,10,15,30},
proper nontrivial subgroup orders in A5 = {2,3,4,5,6,10,12}.
```

Both latter sets have seven elements, but they differ: the subgroup-order set
contains4,12 in place of15,30. The missing proper divisor orders for subgroups
are `{15,20,30}`. Consequently the subgroup split `10=7+3` is not the
arithmetic split `F=U+S=3+7`: the prime-divisor set is `{2,3,5}`, not the three
missing subgroup orders.

Existence of the realized orders is concrete. In the natural permutation action
on `{0,1,2,3,4}`, take cyclic groups of orders2,3,5, the Klein four-group fixing
one letter, the point-stabilizer `A4`, and the groups

```
S3 = < (012), (01)(34) >,                         order6,
D10 = < (01234), (14)(23) >,                     order10.
```

The displayed involution conjugates the displayed cycle to its inverse in
each case. Thus the order6 subgroup really is `S3`, not an assumed subgroup
inferred from divisibility. Subgroups of orders15,20,30 would have indices4,3,2.
Since `A5` is simple, each nontrivial coset action would be faithful, embedding
the60-element group into `S4`, `S3`, or `S2`, which is impossible. For a
self-contained simplicity control, its conjugacy classes have sizes
`1,12,12,15,20`; the only unions containing identity whose sizes divide60 have
sizes1 and60.

The stdlib script independently enumerates every subgroup by adjoining
generators until closure. The complete order histogram is

```
order:   1   2   3   4   5   6  10  12  60
number:  1  15  10   5   6  10   6   5   1
```

This counts59 subgroups total, not seven subgroups. Seven counts their distinct
proper nontrivial **order values**.

## 3. Ten axes are ten subgroups, not ten divisor values

The ten unoriented opposite-face axes of the icosahedron correspond to its ten
Sylow3 subgroups. For `P=< (012) >`, the normalizer is precisely the displayed
`S3`, of order6, so its conjugacy orbit has size `60/6=10`.
This is the geometric realization of the main note's equivariant map from the
ten2-subsets to3-cycle subgroups supported on their complementary triples.

By contrast, taking a subgroup's order sends **all ten** Sylow3 subgroups to
the single divisor value3. Conjugation acts trivially on numerical order values,
while it acts transitively and nontrivially on the ten axes. Therefore this
order-forgetting map cannot identify the axes with the ten proper divisors
equivariantly. Having `A5` as a common parent does not repair the destroyed
subgroup identity, nor give a map to the ten tournament nodes.

## 4. The program's last balanced seed2 state is also60

The separately audited balance packet proves that the standard seed2 orbit
has exactly16 balanced visits, the last at **60 on update70** (seed2 is
update0). Its proof combines a complete finite prefix with a token lower bound
excluding the infinite tail; this addendum imports that argument, not a
finite-search extrapolation. Exact local replay confirms

`update69: 8 -> update70: 60 -> update71: 450`.

The first two states are balanced;450 is not. Thus balance is not preserved,
and the program neither computes an icosahedral symmetry group nor explains
its order. The legitimate connection is the common integer60 with all three
roles explicitly typed: program state / divisor-balance input / group order.

## 5. The user's five inputs split in different ways

All five are square-sum path orders, as established in the main note.

| Integer | Factorization | Constructible regular polygon? | Joins standard PRIMEGAME orbit? | F,U,S | Balanced? |
|---:|---|:---:|:---:|---|:---:|
|15|`3*5`|yes|yes|2,2,2|no|
|16|`2^4`|yes|yes|3,1,1|no|
|17|prime|yes|no|0,0,0|yes|
|23|prime|no|yes|0,0,0|yes|
|25|`5^2`|no|no|1,1,1|no|

The productive-seed claims are exact coalescences, not extrapolations from
three observed outputs:15 is seed2's next state;23 reaches2 after11 updates;
16 reaches97540625 after14 updates, the same canonical state reached by seed2
at update145. Hence the seed2 prime-generation theorem applies to their common
tails. Starting16 does not retroactively emit earlier primes; its subsequent
prime-exponent outputs begin5,7,11. In contrast,
`17 -> 1 -> 55 <-> 65` and `25 -> 1375 <-> 1625` enter exact2-cycles.

The three positive sets are respectively `{15,16,17}`, `{15,16,23}`, and
`{17,23}`. These3-versus2 splits are not the same partition and do not define
a canonical matching to solids, Fermat primes, or pure-blue nodes.

## Verification

```
python3 check_a5_balance_meeting.py > a5_balance_normal.json
python3 -O check_a5_balance_meeting.py > a5_balance_optimized.json
cmp a5_balance_normal.json a5_balance_optimized.json
```

The receipts record all group/arithmetic values and exact source hashes.
All checks use explicit exceptions and survive optimized mode. No original
core file, source worktree, or remote was changed by this addendum.
