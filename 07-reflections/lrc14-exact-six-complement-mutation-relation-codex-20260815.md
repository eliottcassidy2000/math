# The exact-six complement relation: a rigid first-open boundary, not a recursion

## Status

**FINITE-EXACT STOPPING RESULT.**  This is an unnumbered research reflection,
not a proved dependency and not an LRC(14) proof.  Its exact companion freezes
the full pool-14 six-completion relation on the raw THM-3366 `k=2,3` support
rows.  It explains precisely why extending the complement-clock search from
five clocks to six stops at the conjecture itself.
The capped-minimum terminology is repaired under `MISTAKE-395`.

## 1. Inheritance and live concepts

The closest proved mechanism is THM-3366: if at most five integer clocks cover
the unsupported open cells, a hypothetical quotient cover would make at most
twelve global danger clocks, contradicting cited LRC through twelve speeds.
The canonical hostile is the full-period row `D=L`, where the original body
`F` itself is an unavoidable six-clock completion.  The corrected near miss is
to interpret the first completion returned by a branch-and-bound solver as a
canonical mutation; the lawful object is the **full** relation over all
`C(14,6)=3003` completions.  The least-used sidecar is the next chart's aligned
sector and divisor, which the body relation forgets.

The concept board is:

1. THM-3366's pointwise unsupported-cell target;
2. the first open `7+6=13` clock count;
3. the full body-to-body completion relation;
4. the refined post-THM3366 ledgers;
5. the dyadic half/quarter auxiliary charts;
6. the lost next-sector/divisor/tail sidecar.

## 2. Exact universe and relation

For a literal body `F in C({1,...,14},6)`, let

```text
L=14 lcm(F),
D|L,
U_D(F)=the THM-3366 unsupported open-cell target.          (1)
```

For `k=2,3`, retain exactly the support rows passing the inherited THM-2928
cutoff.  For every retained row, the companion computes the minimum size
**through depth six** of a subset of `{1,...,14}` covering `(1)`, with all
strict endpoints and open atoms checked.  Rows of minimum at most five are
already THM-3366 terminals.  The sentinel `None` means no cover by at most six
pool clocks; it does not distinguish minimum at least seven from a target that
the full pool cannot cover.  On a minimum-six row define the full relation

```text
F -> C  iff C in C({1,...,14},6) covers U_D(F)
          for at least one retained divisor D.             (2)
```

This is not a first-witness graph.  Every six-subset `C` is tested.

## 3. Exact boundary atlas

For `k=2`, the minimum-through-six histogram is

```text
1:12662, 2:2764, 3:998, 4:1106, 5:1668,
6:4814, no cover by at most six (>6 or uncovered):3151.     (3)
```

The `4814` exact-six rows have `4918` completion incidences: `4727` rows have
one completion, `70` have two, and `17` have three.  Aggregating by body gives
`4908` directed edges.  After deleting self-loops, the outdegree histogram is

```text
0:1539, 1:1099, 2:299, 3:57, 4:8, 5:1.                    (4)
```

For `k=3`, the corresponding data are

```text
1:12659, 2:2764, 3:976, 4:1052, 5:1602,
6:4778, no cover by at most six (>6 or uncovered):3139;     (5)

4881 completion incidences with row multiplicities
1:4692, 2:69, 3:17; 4872 body edges; nonself outdegrees
0:1559, 1:1088, 2:295, 3:54, 4:6, 5:1.                    (6)
```

At every one of the `3003` full-period rows `D=L`, the unique six-completion
is the original body `F`.  There are six additional lower-divisor self rows in
each sector.  The full relation—not merely the solver's first choices—has only
two nontrivial strongly connected components, both three-cycles:

```text
{(2,3,5,7,9,12), (3,7,8,9,10,12), (4,5,6,7,9,10)},

{(1,3,7,8,9,10), (2,4,5,7,9,12), (4,5,7,8,9,10)}.        (7)
```

The truncated sentinel is active, not cosmetic.  At

```text
F=(1,2,4,6,9,10), D=1260,
```

no subset of at most six pool clocks covers the target, while
`(1,2,3,5,8,9,10)` is an exact seven-clock cover.  Thus the row has pool-14
minimum seven and must not be called pool-14-uncoverable.

Thus the finite relation is strikingly rigid and almost acyclic.  That fact
is an exact structural signal, but it is not by itself an iteration theorem.

## 4. Why this does not close a row

If the inherited seven quotient/aligned clocks covered their carrier and six
clocks covered `(1)`, THM-3366 would construct a global cover by

```text
7 inherited clocks + 6 complement clocks = 13 clocks.     (8)
```

Thirteen is exactly the first open covering count for LRC(14).  No cited
smaller-runner theorem contradicts `(8)`.  A collision between the two banks
would lower the count, but neither the raw support key nor the current refined
key forces one.

One may try to iterate `(2)`, replacing body `F` by `C`.  This is ill-typed.
The map preserves only the existence of a six-clock cover of the present
unsupported cells.  It destroys:

- the next body's aligned-count sector `k`;
- the next divisor and quotient numerator;
- the labels and distinctness of the seven inherited clocks;
- the physical entry time and endpoint owners.

The reciprocal witness

```text
F=(1,2,3,4,5,7), D=1960,
C=(1,5,6,7,9,12)                                         (9)
```

and the cycles `(7)` are decisive hostiles to a scalar descent potential.
Restoring the next-sector/divisor/tail packet is the required sidecar for any
lawful recursion.

## 5. Refined and dyadic stopping boundary

An independent exact composition with the current post-THM3366 refined keys
finds, for `k=2`, `2968` body-only rows (`199075988191` occurrences), `382`
nontrivial-only rows (`133049669` occurrences), and `706` rows with no
six-cover (`860479343` occurrences); for `k=3` it finds `1823`
(`2547058578`), `20` (`313120`), and `54` (`1522136`) respectively.  There is
no body/nontrivial overlap.  These are classifications of an open boundary,
not terminals: the six-completion still makes thirteen clocks.

The same composition tests the midpoint and quarter-pair dyadic auxiliary
charts.  It finds no unconditional or parity-conditional support
contradiction.  The sharp `k=2` hostile has `D=3822`: its two surviving rows
have unsupported sizes `1530` and `1560`.  The full-order scalar
`ceil(3822/7)=546` is **not** the arbitrary transverse capacity: the exact
maximum over every residue modulo `2D` is `1274`, attained by residue `2548`
of quotient order three.  Both targets still exceed that repaired
order-sensitive maximum.  For `k=3`, no post row even enters a base-free
rank-five, rank-six, or odd rank-seven chart.  Aligned parity is absent from
the current row keys, so it cannot be silently inferred.  See MISTAKE-392.

## 6. Connection ledger

| item | exact content |
|---|---|
| source | raw THM-3366 `k=2,3` body/divisor support rows |
| target | full directed relation on the 3003 six-clock bodies |
| map | exact open-atom set cover by every pool-14 six-subset |
| preserved | present body, present divisor witness, and pointwise unsupported-cell coverage |
| destroyed | next sector/divisor, inherited-clock labels, collisions, endpoint owners, and physical time |
| needed sidecar | a chart-to-chart transport of the seven inherited clocks and their distinctness |
| cheapest decisive tests | the unique full-period self completion, the seven-cover depth sentinel, witness `(9)`, and either three-cycle `(7)` |

No tournament is intrinsic: `(2)` is a directed existence relation without a
pairwise orientation gauge, and its SCCs do not encode runner dominance.

## 7. Reproduction

Run

```bash
python3 04-computation/lrc14_exact_six_complement_mutation_graph_probe_20260815.py
python3 -O 04-computation/lrc14_exact_six_complement_mutation_graph_probe_20260815.py
```

The companion is standard-library, `Fraction`/integer/bitset exact, checks all
strict endpoints and open atoms, and enumerates the full relation rather than
sampling or retaining a single solver witness.  Its LF-normalized source and
stored-output hashes are respectively

```text
a799d77af7d930e6a46ab7f22544d78feda327a2a9221025a2c28cf9ed5c85ea
3627ec4fc55ecbca4a571d76dac1e87c9c9f20c4ff6ba9d1346642ccae5117c4
```

and both ordinary and optimized runs reproduce the same stored transcript.
