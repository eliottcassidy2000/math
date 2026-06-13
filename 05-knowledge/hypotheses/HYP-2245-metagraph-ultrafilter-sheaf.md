---
id: HYP-2245
status: OPEN descent/sheaf hypothesis with finite quotient-leak evidence
source: codex-2026-06-05-S670
related:
  - HYP-2244
  - HYP-2243
  - HYP-2240
  - HYP-2187
  - HYP-2167
  - HYP-2164
---

# HYP-2245: Metagraph Ultrafilter Means Cube Filter Plus Quotient Address

## Claim

The useful version of "the tournament metagraph is an ultrafilter" is not that
the tournament isomorphism-class metagraph literally satisfies the filter laws.

It is:

```text
Boolean tiling cube + principal side-choice = literal ultrafilter
isomorphism quotient = observer shadow
mixed quotient fibers = missing address/proof-obligation coordinate
```

Equivalently, the metagraph carries an ultrafilter **sheaf**: each local tiling
address has exact upper/lower membership, but the quotient needs descent data
before the choice is proof-usable.

## Divisor-210 Model

The divisor lattice of `210=2*3*5*7`, ordered by divisibility, is the Boolean
lattice `B_4` via prime support.  S670 brute-forces the finite filters and finds
exactly four ultrafilters:

```text
U_p = {d | p divides d},  p in {2,3,5,7}.
```

For each prime `p`, `U_p` has eight elements, its complement lower ideal has
eight elements, and every pair `d <-> 210/d` has exactly one member in `U_p`.

For every principal choice, the Hasse edges split as:

```text
lower -> lower : 12
lower -> upper : 8
upper -> upper : 12
```

This matches the user's upper/lower picture cleanly: black/blue can be modeled
as the two sides of a principal prime decision.

## Tiling-Cube Model

The fixed-base-path tournament tiling cube `Q_m` has
`m=C(n-1,2)` non-base-path tile bits.  Choosing one tile coordinate `k` gives a
principal upper half:

```text
U_k = {tilings b : b_k = 1}.
```

The complement-tiling operation is `b -> b XOR (2^m-1)`, so every raw
blue/black complement pair has exactly one endpoint in `U_k`.

This is the correct MISTAKE-033 rule: blue/black is complement tiling inside
`Q_m`, not tournament complement `T^op`.

## Quotient-Leak Evidence

S670 checks whether a principal tile-coordinate ultrafilter descends to the
tournament isomorphism quotient.

It does not, already at `n=4`:

| n | tilings | tournament iso classes | degree-even graph iso classes | best mixed iso classes | best leaked complement pairs |
|---|---:|---:|---:|---:|---:|
| 4 | 8 | 4 | 3 | 1 | 4/4 |
| 5 | 64 | 12 | 7 | 5 | 30/32 |
| 6 | 1024 | 56 | 16 | 35 | 512/512 |

So the literal ultrafilter lives on the cube, not on the quotient.  The
quotient metagraph can carry it only after retaining a tile/address coordinate
or a coarser invariant that keeps the upper/lower predicate pure.

Complement-tiling lines also cross isomorphism classes heavily:

```text
n=4: same class 1, different class 3
n=5: same class 4, different class 28
n=6: same class 8, different class 504
```

This supports the "observer-coupled" reading: a quotient observer sees the line
but loses the side membership unless the address is named.

## Equinumerosity Caution

S670 also keeps HYP-2244 in view.  Fixed-path tilings and labelled degree-even
graphs both have `2^C(n-1,2)` elements, but degree-even graph isomorphism
classes are not tournament isomorphism classes:

```text
n=4: tournament iso 4, degree-even iso 3
n=5: tournament iso 12, degree-even iso 7
n=6: tournament iso 56, degree-even iso 16
```

The Royle-even theorem is a different, deeper equinumerosity.  It remains a
count shadow until a retained fiber functor preserves the tournament predicate
being transferred.

## Tournament Analysis

Vertices are interpretations of the ultrafilter/equinumerosity analogy:

1. `LRC_owner_filter_profile`
2. `tiling_cube_principal_ultrafilter`
3. `divisor210_principal_ultrafilter`
4. `iso_metagraph_filter_sheaf`
5. `degree_even_cycle_space`
6. `Royle_even_equinumerosity`
7. `raw_cardinality`

Pairwise observable:

```text
(literal filter law,
 complement decision,
 quotient descent,
 retained address profile,
 LRC actionability)
```

Switch: majority.  Tie Hamiltonian path: listed priority order.

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,5:3}`
- `directed_3cycles=1`
- `scc_sizes=[3,1,1,1,1]`
- `hamiltonian_paths=3`

The top SCC contains the three complementary strengths:

```text
LRC owner-filter profile
tiling cube principal ultrafilter
divisor-210 principal ultrafilter
```

## LRC14 Transfer

For LRC14, do not make the tournament vertices automatically be runners or
arcs.  Candidate vertex/section sets include runners, gaps, fixed circle
sections, section boundaries, wall-crossing events, residues, cover arcs,
Fourier modes, matroid circuits, proof obligations, owner bits, and carry
cocycles.

The quotient should preserve the predicate:

```text
floor vs strict loneliness over the Res_27 / carry-owner lift
```

It may destroy:

```text
which owner/carry/tile/residue is responsible for the side choice.
```

Challenged assumption:

```text
the metagraph quotient itself is the filter.
```

Replacement:

```text
the metagraph quotient is the base space of a filter sheaf; proof work is
descent/purity of the upper-lower side choice.
```

## Next Tests

1. Build the same descent audit over actual LRC14 owner/carry fibers from
   HYP-2241: which principal owner choices survive the quotient?
2. Search for non-principal but quotient-pure filters on the tiling cube
   generated by isomorphism-invariant tile families.
3. Test whether Royle-even graph data admits a side-choice functor preserving
   `H`, `beta1`, and disjoint `c3` packets.
