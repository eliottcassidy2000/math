# Prime-even half-twist cap-seven finite atlas

**Status: FINITE-EXACT for target-free odd primes `p<=599`; the all-prime
tail is OPEN.  This is an unnumbered research companion, not proved canon and
not an LRC(14) closure.**

Codex, 2026-08-15.

## Result and exact boundary

Put `Q=2p`, where `p` is an odd prime, and define the literal half-twist mask

```text
B_(Q,r)={ell in Z/QZ: ||r(2ell+1)/(2Q)||<1/14}.
```

Use sign representatives `1<=r<Q`, discard the identically empty coefficient
`r=p`, and require distinct owners.  Remove

```text
p in {5,11,13,23,29},
```

because then `Q` contains one of the already-supported lower bases
`10,11,13,23,29`.  In the resulting universe of `103` target-free odd primes
through `599`, an exact cover by at most seven masks exists precisely at

```text
p=7:  Q=14, (1,3,4,5,9,11,13),
p=19: Q=38, (1,9,17,20,21,29,37).
```

This sentence is **only** a bounded finite statement.  In particular, it does
not assert that `7` and `19` are the only prime-even atoms for arbitrary `p`.
The complete negative list, including every new negative prime from `191`
through `599`, is frozen in the output companion.

The `Q=14` witness is a partition.  Its type profile is `A^6 E`, all seven
blocks have size `2`, the six `A` owners have order `14`, and the `E` owner has
order `7`.  The `Q=38` witness is also `A^6 E`; all blocks have size `6`, the
six `A` owners have order `38`, and the `E` owner has order `19`.  Its sheet
multiplicities are `1` on `34` sheets and `2` on `4` sheets.

## Inheritance pass

The closest proved mechanism is
[THM-3423](../01-canon/theorems/THM-3423-odd-interval-ratio-complement-and-dyadic-clique-law.md),
which converts disjoint ordinary odd intervals into a bounded-rational Cayley
graph.  The closest dyadic formulation is
[THM-3435](../01-canon/theorems/THM-3435-dyadic-fibre-grid-decomposition-for-literal-half-twists.md),
now **PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED** for the fibre identity,
endpoint law, Boolean branch boundary, and target-free necessary gates.  Its
bounded even census is separately independently audited FINITE-EXACT through
`Q=362`.  The computation here imports neither companion: it reconstructs
every literal mask from the strict integer inequality, and its type and sheet
identities independently agree with the proved statement.

The relevant audit repair is MISTAKE-393.  The degree-`d` circle map is
bijective on each inverse-image component, not on their union, and the
two-sheet partner law must be stated on canonical sign representatives.  The
coefficient universe used here, `1<=r<Q` with the empty `r=p` removed, is
already in that repaired scope.

The canonical hostile is the `Q=14` partition.  It refutes the tempting
projection of an odd active owner to an ordinary radius-`1/14` block on `p`:
its actual projected support has radius `1/7` and still carries a selected
sheet.  The corrected near miss is to replace a disjointness graph by an
exactly weighted graph.  Covers in the large-prime profiles can spend up to
eight sheets of overlap, so zero/nonzero adjacency alone forgets the decisive
budget.  The least-used necessary sidecar is the `A` orientation bit, or,
equivalently, the explicit partner `r <-> Q-r` when both branches are present.

## Three exact owner types

Write an even coefficient as `r=2u`.  There are three types:

- `A`: `r` odd.  It selects one oriented sheet over a widened radius-`1/7`
  base support.
- `B`: `r=2u` with `u` odd.  It is a two-sheet pullback of an ordinary odd
  radius-`1/14` block.
- `E`: `r=2u` with `u` even.  It is a two-sheet pullback of an ordinary even
  radius-`1/14` block.  Every `E` block contains the reflection-fixed base
  fibre, while `A` and `B` avoid it.

For an odd `A` owner, deck translation sends its mask to the disjoint partner
with coefficient `Q-r`, and the pair fills both sheets over its widened base
support.  Projecting one `A` mask without recording its orientation is not a
literal-cover operation.

Let `p=14k+s`, let `(o,b,a)` be the numbers of `A,B,E` owners in a seven-owner
profile, and put `Omega=sum |B_i|-2p`.  The exact sizes and half-budget are:

| `s` | `|A|` | `|B|` | `|E|` | `Omega/2` |
|---:|---:|---:|---:|---:|
| `1` | `4k` | `4k` | `4k+2` | `a-1` |
| `3` | `4k` | `4k` | `4k+2` | `a-3` |
| `5` | `4k+2` | `4k` | `4k+2` | `2-b` |
| `9` | `4k+2` | `4k+4` | `4k+2` | `b-2` |
| `11` | `4k+4` | `4k+4` | `4k+2` | `3-a` |
| `13` | `4k+4` | `4k+4` | `4k+2` | `1-a` |

The two sheets over the fixed base point already cost `2(a-1)` overlap when
`a>0`.  Consequently the residue class `s=3` is impossible from the budget
alone.  The same comparison sharply restricts `a` or `b` in the other
classes.

There is also a useful proof-facing reduction of the remaining raw budgets.
The target-free necessary gates in THM-3435 require exactly seven owners and
at least two `A` owners, and at least one `E` owner (`v_2(r)>=2`).  The last
condition is visible directly on the central two-sheet fibre, which every
`A` and `B` mask misses.  The odd-owner gate also has a direct explanation: an
all-even row is the pullback of an ordinary cover on `p`, while a row with one
`A` owner would force at most six even pullbacks to cover the opposite sheet of
every base fibre.  Thus all-even rows route to the proved odd-modulus
classification THM-3434, and one-`A` rows route to the rank-six gate.  In class
`s=1`, the only raw mixed profile above overlap `8` is `E^6 A`.
More explicitly, its entire overlap budget is already spent by the six `E`
blocks on the fixed fibre; away from that fibre the even union is
deck-invariant, but one nonempty `A` mask is not.  After these honest exits,
every unresolved mixed profile has

```text
Omega <= 8.
```

This bound motivates the weighted tail graph below.  The finite DFS itself
does not rely on either exit: it searches every mass-feasible profile.

## Why the bounded DFS is exact

Odd units modulo `2Q` permute the odd sheet coordinates `2ell+1` and preserve
the three coefficient types.  They act transitively within each type after
sign normalization.  Therefore any nonempty type in a profile can have one
chosen owner normalized to

```text
A -> 1,   B -> 2,   E -> 4.
```

The search chooses whichever available normalized root has the smallest
low-weight neighborhood, then enumerates combinations of the remaining
owners.  This is an orbit normalization, not a heuristic deletion.

For a partial family with union `U`, add masks one at a time and let

```text
I=sum |B_new intersect U_old|.
```

The elementary union identity gives

```text
I=sum |B_i|-|U|.
```

For a fixed profile, `sum |B_i|=2p+Omega`.  Since intersections can only
increase during the search, a branch with `I>Omega` is impossible.  At a
terminal family,

```text
U=Z/(2p)Z  iff  |U|=2p  iff  I=Omega.
```

Thus overlap-budget pruning is an iff test for literal coverage, including
all higher multiplicities; it is not merely a pairwise necessary condition.
The companion enumerates every profile of one through seven owners with mass
at least `2p`.  Across the `103` primes it processes `2,083` profiles,
`56,008` recursive states, and `53,925` candidate branches.  No node or time
cap is used.

The exact positives reproduce the witnesses above.  The negative rows through
`p=181` agree with THM-3435's separately audited FINITE-EXACT bounded census,
while the range `191<=p<=599` is new in this artifact.  Agreement is a
cross-check, not a proof dependency.

## Finite weighted-graph telemetry

For labelled owners define the honest symmetric weight

```text
w(r,s)=|B_(2p,r) intersect B_(2p,s)|.
```

One large prime in each nonzero class modulo `14` was audited.  By odd-unit
normalization, the zero-weight degrees at roots `A=1`, `B=2`, `E=4` are
respectively

```text
(19,31,23)
```

on all six controls.  The following table records the minimum positive root
weight, the complete root-weight histogram through weight `8`, and clique
data.  `raw/gated` is the maximum seven-owner profile excess before/after the
proved `A>=2`, `E>=1`, and fixed-fibre gates.  An empty gated entry means the
whole residue class is structurally impossible.

| `p mod 14` (control) | min positive `(A,B,E)` | weights `<=8` beyond zero | mixed / A-only clique | raw/gated `Omega` |
|---|---|---|---|---|
| `1` (`547`) | `(10,8,8)` | `B: 8^1; E: 8^1` | `6 / 4` | `12 / 8` |
| `3` (`521`) | `(8,8,8)` | `A: 8^2; B: 8^1; E: 8^1` | `6 / 4` | `8 / impossible` |
| `5` (`593`) | `(12,8,8)` | `B: 8^1; E: 8^1` | `6 / 4` | `4 / 4` |
| `9` (`541`) | `(10,12,10)` | none | `6 / 4` | `10 / 4` |
| `11` (`599`) | `(12,12,12)` | none | `6 / 4` | `6 / 4` |
| `13` (`587`) | `(12,8,8)` | `B: 8^1; E: 8^1` | `6 / 4` | `2 / 0` |

Each displayed mixed clique of size six is exactly three complementary
branch pairs, with type profile `A^4BE`.  For example at `p=547`, `Q=1094`,

```text
(1,273,274,820,821,1093)
 =(1,1093) union (273,821) union (274,820).
```

This is an exhibited equality structure on the six finite controls, not a
forced-pair theorem for all primes or for positive-weight configurations.

Finally, every one of the `19+31+23=73` normalized zero-neighbor incidences
on every control has a raw coefficient ratio represented modulo `p` by

```text
+-a/b,  gcd(a,b)=1,  a+b<=7.
```

The maximum observed height is exactly `7`.  This is the finite shadow of the
circular-gap/Bezout mechanism in THM-3423.  It does not yet say that every
edge of weight at most `8` belongs to a bounded atlas uniformly in `p`.

## OPEN tail connection and loss ledger

The proposed tail object is a typed weighted Cayley graph, not a tournament.

| field | exact content |
|---|---|
| source | labelled literal masks `B_(2p,r)` with types `A,B,E` |
| target | the graph with vertex labels and edge weight `w(r,s)` |
| map | odd-unit normalization followed by the coefficient ratio modulo `p` |
| preserved | exact pair intersections; with all labels retained, the total pair-weight lower bounds the cover excess |
| destroyed by the bare ratio | sheet locations, higher intersection multiplicity, literal union, dyadic endpoint types, and the `A` orientation |
| required sidecars | ordered endpoint types, the `A` sheet/orientation bit, complementary branch-pair labels, and the `E` fixed-fibre incidence |
| desired uniform lemma | every typed ratio with `w<=8` lies in an explicit bounded rational atlas |
| cheapest decisive test after that lemma | classify seven-vertex typed weighted configurations whose total incremental overlap is at most `Omega<=8` |

THM-3423 already supplies the exact zero-weight ratio and clique law for the
ordinary odd `B/B` sector beyond its stated threshold.  A widened-interval
analogue should control the projected `A/A` sector, but it must be proved with
the orientation sidecar still attached.  Mixed `A/B/E` weights then require a
typed near-complement lemma, not an untyped appeal to the ordinary interval
ratio set.

The precise open boundary is therefore:

> Prove a uniform typed height bound for every ratio of pair weight at most
> eight, then show that no seven-owner weighted configuration within the
> profile budget can retain all literal sheets, except the already-known
> `p=7,19` controls.

Until both clauses are proved, the finite census must not be extrapolated.

## Reproduction and hashes

Run with only the Python standard library:

```bash
PYTHONHASHSEED=0 python3 -B 04-computation/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.py
PYTHONHASHSEED=1 python3 -B -O 04-computation/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.py
```

The normal and optimized outputs are byte-identical and match
[the frozen output](../05-knowledge/results/lrc_prime_even_half_twist_cap7_finite_atlas_20260815.out).

```text
script LF SHA-256   ed54b01f9bf155643b6407c6af8ee6f15c7a099f8ac3b60399dd49b496fb1d12
output LF SHA-256   9751da7e0464f16f92b675ca9c95b118fd3fe95fc11cc99c9a7015d8425adb65
semantic SHA-256    19f846db72803c683f608149aac8c7da2015eaca6e1e3524d9954eaeef0fa826
```

No arbitrary common-centre row, physical time, decrement, all-prime
classification, or LRC(14) consequence is claimed.
