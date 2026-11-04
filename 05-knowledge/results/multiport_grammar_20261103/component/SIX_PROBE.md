# Complete strong-six atoms: a shifted 44 supplier, while 24 remains blocked

2026-11-03. **PROVED construction and bounds; FINITE-EXACT exhaustive decisions.**
Independent analytic/source/reconstruction audit passed.
The finite atom universe and the two bounded inverse decisions are exhaustive.
The all-word conclusions use the bounded-intercept theorem in
[SHIFT_PROOF.md](SHIFT_PROOF.md). The positive return also has
an explicit old-name-preserving eleven-arc schedule. Compilation into a more
restrictive native operation set and the complete return tree have separate
certificates; global spectrum coverage is not asserted here.

## 1. Precise boundary and the actual positive object

The inherited seven-coordinate signature and ordered-component product are
those of [PROOF.md](PROOF.md), originally the
[S182 induction transport](../../session182_physical_returns_20261104/induction/PROOF.md).
Write a signature as `(h,p,f,g,j,e,s)` and its intercept as `c=h-j`.
A marked suffix attached after the nonempty transitive clock prefix has count
`H_t=c+2^t f`. A clock shift is a fixed integer, legal when the resulting
prefix remains nonempty. The full order-five all-shift theorem leaves both
addresses `24` and `44` unaccepted. Their required `(c,odd(f))` values are
`(-35,63)` and `(-29,63)` respectively.

**New exact conclusion.** In the same scalar replacement architecture:

* `44` has a suffix with SCC orders `(6,1)` and clock shift `-2`.
* `24` has no suffix all of whose SCCs have order at most six, at any word
  length or fixed integer clock shift.
* No single strong order-six atom supplies either target.

The actual strong-six supplier is

```text
rows = (2,4,57,3,11,27), ports = 16 = {4},
sigma = (17,224,28,44,0,2,4).
```

Rows are out-neighbour bitmasks on labels `0,...,5`. Append one marked
singleton after this component, orienting every cross edge toward it. The
resulting literal seven-vertex suffix is

```text
rows = (66,68,121,67,75,91,0), ports = 80 = {4,6},
sigma = (17,448,252,105,46,17,4).
```

For example, its joint coordinate is `0+0+2*1+44*1=46`, and its left
marginal is `28+224=252`. Thus for every `t>=3`,

```text
H_target(t-2) = -29+252*2^(t-2)
              = -29+63*2^t = 3*(-11+21*2^t)+4.
```

This is the `+4` child of the original `4` sibling. Its terminal marked
singleton guarantees strongness of the complete root attachment by the
inherited SCC guard. The source independently counts this entire suffix and
its attachments at three legal clocks. The scalar equality is then an
identity for all legal clocks, not an extrapolation of those controls.

The recorded word is unique only in the source's seven-signature representative
bank and finite prefix quotient. It is not claimed to be the unique labelled
suffix. Unmarked singleton prefixes alone create further representatives with
the same odd marginal and an adjusted clock shift.

## 2. Physical reblocking and eleven strong arc reversals

The original `4` sibling is the root attachment with clock `TT_t`, suffix
`TT5` and suffix ports `{0,2,4}`. For `t>=3`, retain the first `t-2` clock
vertices as the new clock. The last two clock vertices, followed by all five
old suffix vertices, form a physical `TT7` block. No vertex or edge has been
changed by this reinterpretation. Number this block `0,...,6`; its ports are
`{2,4,6}`, and its signature is `(1,128,84,21,12,1,0)`.

Reverse the following ten old internal arcs, written in lexicographic order:

```text
(0,2), (0,3), (0,4), (0,5), (1,3), (1,4), (1,5),
(3,4), (3,5), (4,5).
```

These are exactly the inversions converting the first six transitive vertices
into `(2,4,57,3,11,27)`. Then reverse the old edge from block vertex2 to the
physical root, removing that port. The resulting block is precisely the
supplier in Section1. All old vertex names and the root are retained; the
total tournament order is still `t+6`. This is not a newly born suffix whose
counts merely match the old family.

At every intermediate step, vertex6 is still a sink in the suffix and a port.
It is therefore the final SCC, and the inherited terminal-SCC attachment
criterion proves the complete tournament strong for every `t>=3`. This also
proves strongness for any ordering of the eleven independent reversals: none
changes an arc incident with vertex6 or removes its port. The displayed order
is used for the reproducible path and observer table. The source records all
twelve complete seven-signatures, including entry and exit, and literal
adjacency at original clocks `t=3,4,5`. It verifies that each step changes
exactly its declared physical arc, that the graph remains strong, and that
the named root readout obeys `b_root=F_root-H=j_suffix` at every step.

In particular the actual root defect changes from12 to46. These operations
are also within the named-flip semantics of the separately audited coloured
[transport grammar](../transport/PROOF.md), if the seven block names and root
are retained. No historical anonymous-name resurrection is used. This does
not, by itself, compile the schedule into any more restrictive native guard
language or supply a return from the newly produced state.

## 3. Exhaustive class generation without an external table

Start with the single one-vertex tournament. Given one representative of
every isomorphism class of order `n-1`, append a fresh vertex in all
`2^(n-1)` orientations. For every augmented tournament, inspect all `n!`
permutations and retain the lexicographically smallest row tuple. Deduplicate
these canonical tuples.

This generates every class of order `n`: delete any vertex from an arbitrary
`n`-vertex tournament, identify the remaining tournament with one of the
retained representatives, and use its transported incidence bitmask for the
new vertex. Exhaustive permutation minimization identifies exactly the
isomorphism classes. No strongness or port filter is used before this step.

The resulting exact census is:

| Order | All classes | Strong classes | All labelled mass | Strong labelled mass |
|---:|---:|---:|---:|---:|
|1|1|1|1|1|
|2|1|0|2|0|
|3|2|1|8|2|
|4|4|1|64|24|
|5|12|6|1024|544|
|6|56|35|32768|22320|

For each representative, the source independently counts its automorphisms
by exhaustive permutation stabilizers and adds `n!/|Aut|`. The total equals
`2^(n choose 2)` at every order. Strongness is checked both by Floyd-style
reachability closure and by the inherited per-source reachability routine.

Every one of the 35 strong-six representatives is assigned every one of its
64 port subsets. This gives 2,240 actual port assignments and 1,950 distinct
seven-signatures. Together with all smaller strong components the quotient
has 2,120 representatives; its order-five part is independently equal to the
inherited bank generated from all 544 labelled strong-five tournaments.

This quotient is safe for scalar ordered-word acceptance because every
product and scalar observer uses only the seven coordinates. Every retained
signature still stores actual adjacency rows and ports. It does not identify
physical names for a native transport problem.

Finite strong-six statistics are `min h=15`, `max h=45`, and `min f=15` over
nonempty port sets. The lower bounds are attained together at

```text
rows=(2,4,9,51,7,23), ports=8,
sigma=(15,216,15,45,0,1,3).
```

These are complete order-six statistics, not assertions at other orders.
The direct single-atom tests for `c in {-35,-29}, odd(f)=63` find no match.

## 4. Why the all-word query is finite

The implementation uses the complete prefix normalization and bounded
intercept inverse of [SHIFT_PROOF.md](SHIFT_PROOF.md), with `M=45`.
For a normalized target word, first and last components are marked and

```text
h <= B(c) = max(45^2,abs(c)).
f <= L(c) = max(max_A f_A,
               max_marked_A floor(((f_A+p_A)B(c)-p_A c)/h_A)).
```

The bound on `h` follows because a normalized word with at least three atoms
has `j>=2h`; a word with at most two atoms has `h<=45^2`. The marginal bound
then follows by retaining the first atom's contribution `g_A f_tail` to `j`.
These are actual-atom inequalities, not an assumed maximum word length.

The initial unmarked prefix contributes the multiplicative pair
`(h_prefix,odd(p_prefix))`. For `24`, the only attainable pairs dividing
`(35,63)` are `(1,1)` and `(5,1)`. For `44`, only `(1,1)` occurs. One actual
representative of each pair suffices for existence after allowing all fixed
integer shifts; lost physical information is not recovered by this quotient.

The exact complete queries are:

| Address | Prefix pair | Normalized c | h bound | f bound | Tested f |
|---|---|---:|---:|---:|---|
|24|(1,1)|-35|2025|56799|63 times 2^a, 0<=a<=9|
|24|(5,1)|-7|2025|56395|63 times 2^a, 0<=a<=9|
|44|(1,1)|-29|2025|56712|63 times 2^a, 0<=a<=9|

For each fixed marginal, enumerate attainable products of atom Hamilton
counts `h<=min(f,B)`. Peel the first atom `(H,P,F,...)` by

```text
tail_h=h/H, tail_f=(f-F*tail_h)/P.
```

Only positive integral tails are retained. Since `P>=2`, the marginal strictly
decreases. A marked one-atom word is the base case. Grouping atoms with the
same `(H,P,F)` only avoids repeating this divisibility computation: all actual
representative indices are expanded afterward, with their distinct remaining
coordinates retained. Thus the recurrence is exactly the inherited inverse.
Words whose first atom is unmarked are removed only after the exhaustive
fixed-marginal recursion, as required for normalization.

The complete results are:

| Address | Fixed-marginal calls | Normalized candidate visits | Accepted representatives |
|---|---:|---:|---:|
|24|20|1610|0|
|44|10|805|1|

The sole representative is the object in Section1. The negative `24` result
is therefore an all-word, all-fixed-shift result for this finite SCC bank.
It implies that a `+4` return from the original `2` sibling needs an SCC of
order at least seven, or a departure from this architecture or uniform clock
notion. It does not say that every operation on the actual original source
has been classified. This core receipt tests no successors of `44`; the
separate [complete tree](SIX_TREE.md) supplies those decisions.

## 5. Independent paths, hostile boundaries and reproduction

The seven-signature producer inherited from `component_inverse.py` uses
endpoint subset dynamic programming. The new source also enumerates every
actual directed vertex permutation on every subset of each strong class.
Summing complementary-subset path counts recovers all seven coordinates,
independently, for every raw port assignment through order six. Empty paths
are included only in the two-slot observables, with the inherited convention.

The new grouped left inverse is compared against the inherited independent
right inverse for every query with `f<=1008`. The whole larger inverse has
its own all-length proof above; these comparisons are controls, not the
reason it terminates. The accepted suffix and its three clock attachments
are also counted directly from adjacency. Singleton and all-marked directed
triangle signatures are explicit positive/hostile regression controls.

The class-generation augmentation bound is structural, not a sampled order-six
bank. Conversely, the six-vertex statistics are not extrapolated to order seven.
The failure of every single-six atom does not imply failure of words: the
`(6,1)` supplier is the exact hostile to that inference. Section2 additionally
pays the old-name-preserving physical path, without inferring a native compiler.
In particular, this result neither promises reusable two-branch closure nor
overrides the adaptive two-branch
density warning in
[session154 closure/NOTE.md](../../session154_gadgets_20261103/closure/NOTE.md).

```bash
python3 04-computation/multiport_grammar_20261103/atom_six_probe.py --shift
python3 -O 04-computation/multiport_grammar_20261103/atom_six_probe.py --shift \
  --output /tmp/multiport_atom_six_optimized.json
```

The source's default invocation performs only the complete atom census and
single-atom queries; `--shift` is required for the all-word decisions. The JSON
receipt records the full bank, all class representatives and automorphism
counts, every normalized query and all accepted actual graphs. Both inherited
source modules are hash-pinned before import. Final replay pins follow below.


Both normal and optimized full `--shift` replays pass **2,999 always-active gates**.
Their **590,751-byte** JSON receipts and stdout are byte-identical.

Frozen source: **16,140 bytes**, SHA256
`ed295d2f7ce1a7815a7c377ede9544f3118efe3736458d1fc921d460c21900ba`.
Frozen receipt: [atom_six_probe.json](../../../../04-computation/multiport_grammar_20261103/atom_six_probe.json), SHA256
`1a57397d99deacf460d6e7c01f80da705b0304240d0fec363108cedb5bc4d951`.

The exact displayed physical schedule has these seven-coordinate readouts:

| Step | Reversed arc (block names) | (h,p,f,g,j,e,s) |
|---:|---|---|
|0|entry|`(1, 128, 84, 21, 12, 1, 0)`|
|1|0 to 2|`(3, 192, 123, 63, 36, 3, 1)`|
|2|0 to 3|`(5, 256, 166, 97, 56, 5, 2)`|
|3|0 to 4|`(9, 336, 203, 149, 80, 9, 5)`|
|4|0 to 5|`(17, 448, 286, 211, 120, 17, 10)`|
|5|1 to 3|`(15, 432, 271, 187, 104, 15, 7)`|
|6|1 to 4|`(15, 432, 260, 185, 96, 15, 6)`|
|7|1 to 5|`(17, 448, 269, 167, 84, 17, 5)`|
|8|3 to 4|`(17, 448, 259, 217, 112, 17, 9)`|
|9|3 to 5|`(17, 448, 259, 217, 112, 17, 9)`|
|10|4 to 5|`(17, 448, 269, 167, 84, 17, 5)`|
|11|2 to root|`(17, 448, 252, 105, 46, 17, 4)`|

The full intermediate adjacency matrices and all root readouts at clocks3,4,5
are retained in the receipt. The [independent audit](SIX_AUDIT.md) rebuilds the
complete bank, uses a larger search envelope and independently checks this
physical schedule. Root promoted only the status/routing after that audit;
the producer source and receipt remain frozen.
