---
id: THM-838
title: Centered-CF transport preserves the n=9 rank-four defect but exposes two dimensions to n=10 raw S2
status: PROVED RECURSIVE COORDINATE COPY + FINITE-EXACT DEFECT/Q/ENDPOINT-COUPLED TRANSPORT
source: codex-2026-07-15-S13/engine
depends_on: [THM-778, THM-812, THM-813, THM-828, THM-832]
related: [THM-825, THM-829, THM-834, THM-846, THM-853, HYP-6880]
verification:
  - 04-computation/continued_fraction_n9_defect_transport_codex_S13.py
  - 05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.out
  - 05-knowledge/results/continued_fraction_n9_defect_transport_codex_S13.json
---

# THM-838 — centered-CF transport exposes half the defect

Continue THM-812/813's centered Christoffel coordinate-copy recursion along
the consecutive pairs

```text
(p,q)=(2,3),(3,4),(4,5),(5,6),(6,7).
```

Starting with phase `s=1`, THM-778's update `s'=(1-s) mod 2` gives phases
`1,0,1,0,1` and increment words

```text
5->6:  (1,2)
6->7:  (1,2,1)
7->8:  (1,1,1,2)
8->9:  (1,1,2,1,1)
9->10: (1,1,1,1,1,2).                                  (1)
```

Let `Phi:X_9->X_10` be the final coordinate copy.  It is injective, commutes
with complement, and intertwines staircase reflection.  On THM-828's
four-dimensional defect space `V`, its linear part has rank four, but

```text
dim(Phi(V) intersect S_10)=2,                            (2)
```

where `S_10` is the target raw-S2 parity-syndrome code.  The surviving
subspace is exactly `span(b0,b1)`.  Thus centered-CF transport preserves the
literal four-dimensional object while making half of it visible to raw
layer parity.

On the 58 literal false-palindrome pairs, all 58 target reflection orbits and
all 58 projected complement-line node-pair cells remain distinct.  Yet only
ten pairs retain raw-S2 equality.  Hence the continued-fraction action is
more faithful on literal `Q` and endpoint-coupled `bar P` carriers than on the
static raw histogram codec.

## The derived coordinate map

Index the 28 source and 36 target staircase coordinates in repository tile
order.  The target-to-source copy word is

```text
rho_9,10 =
(0,1,2,3,4,5,6,6,7,8,9,10,11,12,12,13,14,15,16,17,
 17,18,19,20,15,21,22,23,24,21,25,26,24,27,26,27).      (3)
```

Every source coordinate occurs, so copying by (3) is injective.  The verifier
also reconstructs the preceding maps recursively rather than accepting (3)
as input.  At every size `5->6` through `9->10`, `rho` is surjective and

```text
rho sigma_target = sigma_source rho.                    (4)
```

The `6->7` reconstruction exactly reproduces THM-813's independently stored
map, providing the base regression for the recursion.

## Linear action on the four-dimensional core

For THM-832's ordered basis

```text
b0=0192486, b1=08c2c0c, b2=11b4600, b3=4483414,
```

the target images are

```text
Phi(b0)=000c48906, Phi(b1)=00860980c,
Phi(b2)=110dd0c00, Phi(b3)=48440e814.                    (5)
```

Both lists have binary rank four.  Direct evaluation of the target parity
checks gives

```text
Phi(b0),Phi(b1): no parity defect,
Phi(b2): tau=7 low and high parity defects,
Phi(b3): tau=5 low and high parity defects.              (6)
```

Consequently the only coefficient coordinates `z`, including zero, whose
images have zero target syndrome are

```text
0000,0001,0010,0011,                                    (7)
```

where the binary coordinate convention is THM-832's `(c3 c2 c1 c0)`.
Equations (5)--(7) prove (2).  This is a statement about the linear parity
syndrome.  It does not by itself assert equality of the two particular target
S2 words.

## Literal layer genealogy

Apply `Phi` to both endpoints of every one of THM-828's 58 pairs.  The number
remaining equal through successive target raw-S2 layers
`tau=3,...,9,fixed-10` is

```text
58,58,20,20,10,10,10,10.                               (8)
```

Thus 38 pairs first split at `tau=5`, ten more at `tau=7`, and ten survive
raw S2.  The survivor sectors are exactly

```text
D=0192486:6 pairs, 08c2c0c:2 pairs, 095088a:2 pairs.     (9)
```

The target first positional moments separate all ten: four first split at
`tau=7` and six at `tau=8`.  This is an actual-bank statement; THM-825 remains
the stronger unconditional positional reconstruction theorem.

The source chirality bit is not functorial.  Its sign is preserved on 50
pairs and reversed on eight, precisely the four pairs in each of sectors
`18e4e8a` and `1976a0c`.  Therefore a static one-bit repair at `n=9` cannot be
transported without the layerwise moment word or an explicit action on its
orientation.

## Literal `Q` and projected endpoint-coupled transport

Reflection equivariance gives

```text
Phi(sigma_9 u)=sigma_10 Phi(u).
```

For a line endpoint define `e(u)={u,u xor FULL}` and THM-813's literal carrier
`Q(u)=[e(u)]_sigma`.  All 58 target `Q` cells are distinct.  For the canonical
apex-zero representative, distinguish the following ordered and unordered
coupled node coordinates, retaining multiplicity and the common endpoint
swap:

```text
P^->(u)=(nu(u),nu(u xor FULL)),
bar P(u)={nu(u),nu(u xor FULL)}.                          (10)
```

Exact tournament isomorphism on the 232 source and 232 target endpoint masks
finds

```text
source bar-P cells       58 singleton nonloops
target bar-P cells       58 singleton nonloops
bar-P descent failures    0
reflection failures   0.                                 (11)
```

There is a useful convention audit.  On the 58 chosen source presentations,
the direct-chart node projection has 53 values and the complement-partner
projection 54; the local classifier's bit convention names these in the
opposite order.  Their fibre histograms are

```text
53 projection: {1:48,2:5},
54 projection: {1:51,2:2,3:1}.                          (12)
```

Both `P^->` and `bar P` have 58 singleton values.  This
independently matches THM-834 and proves that endpoint coupling restores the
information lost by either bare node projection on this bank.

THM-813 shows that projected-cell descent fails globally already at `X_6 ->
X_7`.  Equation (11) is therefore a special purity theorem for the 58-orbit
defect bank, not a global functoriality claim at `n=10`.

## Tournament Analysis and boundary

Use seven information carriers as vertices:

```text
source sector, target difference, target syndrome, target S2 orbit,
source bar P, target bar P, target literal Q.
```

The pairwise observable is the number of unordered source witness pairs
separated.  Switching from raw retention to retention per logarithmic cell
cost flips 14 edges.  Ties use the displayed carrier order.  Both gauges are
transitive with score histogram `0,...,6`, no directed triangle, singleton
SCCs, and one Hamiltonian path.

The challenged assumption is that a continued-fraction copy should act on a
bare node or raw histogram.  Here it acts exactly on literal coordinates and
`Q` orbits; the endpoint-coupled `bar P` carrier also happens to descend on this
finite image.  The action preserves no LRC gap, owner, wall, or loneliness
predicate.  It is one centered consecutive CF word, not an arbitrary
`GL_2(Z)` word, and it does not classify the global `n=10` metagraph. ∎
