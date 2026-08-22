# Missing and double edges are the symmetric channel of a tetrahedral tournament

**Status: VERIFIED-EXACT finite atlas; not independently audited.**  The
companion enumerates all `4^6=4096` labelled loopless digraphs on four
vertices.  It treats tournaments, partial orientations, and semicomplete
digraphs as slices of one pairwise object without forcing ties into arrows.

## One unordered pair is a four-state cross

For `i<j`, retain both arc bits `(A_ij,A_ji)` and set

```text
s_ij=A_ij+A_ji-1,              o_ij=A_ij-A_ji.          (1)
```

The four possibilities are

| state | arcs | `s` | `o` | XOR |
|---|---:|---:|---:|---:|
| missing | `00` | `-1` | `0` | `0` |
| forward | `10` | `0` | `+1` | `1` |
| reverse | `01` | `0` | `-1` | `1` |
| double | `11` | `+1` | `0` | `0` |

Thus

```text
|s_ij|+|o_ij|=1.                                      (2)
```

XOR detects whether exactly one arrow is present, but deliberately merges
missing and double edges.  The centered symmetric coordinate `s` is the
smallest sidecar that separates those two XOR-zero states.  A tournament is
the equator `s=0`; missing/double information is not an orientation tie.

## The six symmetric and six skew coordinates split differently

On `K4`, the symmetric edge space decomposes as

```text
vertex incidence degrees (rank 4)
  + perfect-matching contrasts (rank 2).                    (3)
```

Four centered undirected degrees and any two of the three matching sums form
a nonsingular six-coordinate system of lattice index `8`.  The third
matching sum is the dependent symmetric check.  This is the ordinary
`1+3+2` decomposition of the six unordered edges.

The skew edge space decomposes as

```text
net score/cut space (rank 3)
  + Haar square holonomy H1 (rank 3).                        (4)
```

Three independent net scores and the three primitive complementary-square
cycles from the Haar/XOR sidecar have lattice index `32`.  Combining (3) and
(4) has continuous lattice index `8*32=256`, and reconstructs every one of
the `4096/4096` discrete directed graphs exactly.

The hostile controls locate both losses sharply:

- all-missing and all-double have identical zero skew channel;
- a tournament and its converse have identical zero symmetric channel;
- net scores alone forget cycle holonomy;
- undirected degrees alone forget the matching-contrast plane.

## Two ternary charts and the tournament equator

Forbidding double edges gives `3^6=729` partial orientations.  Forbidding
missing edges gives `3^6=729` semicomplete digraphs.  Their intersection is
the `2^6=64` tournament equator.  Allowing the choice independently on each
edge gives the full `4^6` object, so it is misleading to call the whole space
one ternary tournament.

Exact relabelling under `S4` gives `218` unlabeled loopless digraphs.  This
orbit census is a hostile against accidentally quotienting the symmetric and
skew channels by different vertex gauges.

## Connection contract

| field | content |
|---|---|
| vertices | four labelled atoms |
| pairwise observable | two directed arc bits per unordered pair |
| skew target | orientation/current and `H1` holonomy |
| symmetric target | missing-versus-double occupancy and matching contrasts |
| XOR preserves | single-arrow support |
| XOR destroys | missing versus double |
| tournament restriction | `s=0`, all `o=+/-1` |
| sidecar for generalized tournaments | retain both `s` and `o` |

This is an intrinsic pairwise decomposition, so Tournament Analysis is
legitimate here.  It does not identify any LRC endpoint relation, JC flux,
physical current, or proof obligation with a digraph edge.
