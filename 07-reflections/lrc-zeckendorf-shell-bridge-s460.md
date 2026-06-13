# LRC, Zeckendorf Supports, and Tournament Fugacity

**Session:** codex-2026-05-31-S460

The useful old dictionary is still the simple one:

```text
Zeckendorf      = independent sets in the Fibonacci path, fugacity x=1
Tournament OCF  = H(T)=I(Omega(T),2), fugacity x=2
LRC endpoint IP = rows covered by speed columns
```

This session adds a missing coordinate to the LRC side: label every speed
column by its Zeckendorf support.  A speed is then a vertex of the Fibonacci
cube, and the LRC endpoint incidence matrix becomes a row-cover problem over
Fibonacci-cube vertices.

## Repo Threads Re-read

I used four old Zeckendorf/tournament strands as the backbone:

- `04-computation/zeckendorf_tournament.py`: Zeckendorf is `I(P_m,1)`,
  while tournament path-conflict graphs give `I(P_m,2)`.
- `07-reflections/zeckendorf-non-consecutivity-pairing.md`: non-consecutive
  Fibonacci terms are path-independent sets; the tournament fugacity shift
  turns Fibonacci/Lucas data into Jacobsthal/Mersenne data.
- `07-reflections/summand-graph-fermat-zeckendorf.md`: the additive summand
  graph has Fibonacci reset points and triangular tournament-staircase
  thresholds.
- `07-reflections/lucas-summand-graph-zeckendorf-geometry.md`: Lucas numbers
  are minimum-gap Zeckendorf pairs `F_{m-1}+F_{m+1}`.

Then I overlaid the recent LRC strands: S386 labelled endpoint cycles, S411
row/column recursion, S420 integer-program rows, S430 incidence core, and S440
the `n=14` gate fan tax.

## Main Computation

The new script is `04-computation/lrc_zeckendorf_bridge_s460.py`.

It verifies the familiar path-fugacity bridge:

```text
I(P_4,1) = 8
I(P_4,2) = 21 = F7
I(C_3,2) = 7 = F4+F2
```

Those are exactly the old forbidden-H anchors: `7` is the odd-cycle branch,
and `21` is the path-conflict branch.

The LRC denominator frontier has a clean Zeckendorf shell:

```text
13 = F6
14 = F6+F1
15 = F6+F2
16 = F6+F3
18 = F6+F4  (Lucas/min-gap pair)
21 = F7     (reset)
24 = F7+F3
```

This is a better local coordinate than raw adjacency of runner counts.  The
`14,15,16` block is the low-payload walk over the same `F6=13` anchor; `18`
is the next Lucas/minimum-gap payload; `21` is the next reset.

## Exact n=14 Fan Factorization

The strongest new fact is the complete family of minimum lower covers for the
`14`-gate endpoint row.  Earlier S440 had one exact cover and the forced
private columns.  The S460 audit enumerates all minimum covers:

```text
{1,3,5,7,9,11,13} union {2}
{1,3,5,7,9,11,13} union {4}
{1,3,5,7,9,11,13} union {6}
{1,3,5,7,9,11,13} union {8}
{1,3,5,7,9,11,13} union {10}
{1,3,5,7,9,11,13} union {12}
```

So the local `14`-gate row factors as:

```text
forced odd fan + one arbitrary even bridge.
```

In Zeckendorf digits, the forced odd fan is:

```text
1  = F1
3  = F3
5  = F4
7  = F4+F2
9  = F5+F1
11 = F5+F3
13 = F6
```

Its digit load is:

```text
F1:2 F2:1 F3:2 F4:2 F5:2 F6:1
```

The private endpoints force this whole fan.  The even bridge is locally free
inside the `14`-gate endpoint row.  That is exactly the kind of local quotient
one wants before proving a global Hall/Farkas dual: collapse the six even
bridges to one fiber, then show the rest of the LRC rows split the fiber with
positive total cost.

## Comparison With n=15 and n=16

The same local computation over lower columns gives:

```text
n=15: exact size 10, 8 minimum covers
n=16: exact size 9, unique minimum cover
```

That makes `n=14` the cleanest Zeckendorf-shell laboratory.  It has one local
free coordinate.  `n=15` has a two-column free part, and `n=16` is rigid,
matching the existing picture: `n=15` is an odd-prime product-building
transition, while `n=16` is pure dyadic endpoint flow.

## Owner Debt

S440's owner-charge idea becomes sharper under Zeckendorf supports.  In the
S380 gate ladder, the top exposed owners are:

```text
154 = F11+F5+F2
168 = F11+F7+F3
182 = F11+F8+F3+F1
```

So the gate-heavy repair does not merely make speeds larger.  It raises the
top Zeckendorf index and exports endpoint debt to higher Fibonacci shells.
This suggests a second potential beyond p-adic depth:

```text
Zeckendorf height = max Fibonacci index in the owner support.
```

A candidate proof could combine p-adic product-depth with Zeckendorf height:
if a gate repair preserves frontier mass by moving in the `2 x 7` product
tree, it still pays by moving owner labels up the Fibonacci shell.

## Connection Back To Natural Operation Graphs

The user's old natural-number graphs fit this lens well.

Addition `X+Y=Z` has a no-carry region in Zeckendorf coordinates: disjoint,
non-adjacent supports simply union.  Multiplication `X*Y=Z` is sparse and
carry-heavy, because products of Fibonacci sums usually leave the Fibonacci
cube coordinate chart.  Product-sum equations such as

```text
X+Y = X*Y
X+Y+Z = X*Y*Z
```

are carry interfaces between those two shadows.

For LRC, the endpoint IP sees exactly this split.  Small-denominator rows are
multiplicative/divisibility invoices.  Gate fan rows are additive-looking
residue fans.  Zeckendorf supports give a way to say when a repair is no-carry
within the additive chart and when it must export debt through multiplicative
carry.

## Proof Routes

1. **Even-bridge quotient.**  Prove the `n=14` gate row as a lemma:
   all minimum lower covers are odd fan plus one even bridge.  Then build a
   reduced IP whose first branch variable is the bridge fiber.

2. **Bridge-breaking rows.**  Add coarse denominator and primitivity rows to
   the local problem and measure which even bridges survive.  If each bridge
   pays a different unavoidable row, the six-way freedom becomes dual weight.

3. **Zeckendorf height divergence.**  Track exposed endpoint owners by max
   Fibonacci index.  Try to prove gate-heavy repairs cannot keep both p-adic
   product-depth and Zeckendorf height bounded.

4. **Tournament analogy.**  Path conflict graphs give Jacobsthal `I(P_m,2)`.
   The LRC fan row may be the row-cover analogue of a path conflict component:
   a locally path-like independent coordinate whose realization constraints
   force extra rows.

## Caution

Zeckendorf support is not a complete invariant.  It does not know modular
residue by itself, and endpoint protection still depends on exact rational
inequalities.  The value here is as a coordinate system for branching,
feature extraction, and possible dual weights, not as a direct standalone
proof.

The immediate next computation should take the six `n=14` bridge fibers and
add row families one at a time: coarse denominators, primitivity, S440 owner
debt, and Bruhat-Tits product-depth leaves.
