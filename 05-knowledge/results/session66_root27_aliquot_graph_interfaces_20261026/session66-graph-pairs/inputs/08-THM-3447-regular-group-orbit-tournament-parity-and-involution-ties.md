---
id: THM-3447
title: "Regular group-orbit tournaments exist exactly at odd order"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.
  A regular finite-group orbit admits an invariant tournament exactly when
  the group order is odd.  More generally, every nonidentity involution
  forces one disjoint perfect matching to be missing or bidirected in any
  invariant directed relation; all generalized, tournament, partial,
  semicomplete, symmetric, and maximal-partial counts follow from
  inverse-pair orbits.  For THM-3444's Hensel lattices this separates
  odd-prime half-set gauges from two-adic forced ties.  No canonical
  orientation or LRC current follows.
source: codex2 regular-orbit inverse-pair synthesis, 2026-08-15
audit: independent proof reconstruction; size-four/six invoice; XOR and harmonic-weight boundary audit; normal/optimized/stored exact replay
depends_on:
  - THM-3444-commuting-smooth-hensel-vector-field-lattice-action
related:
  - THM-2504-endpoint-tournament-no-go-and-root-chart-holonomy
  - THM-2595-modular-v4-affine-lift-dichotomy-and-six-vertex-tournament-no-go
  - LEM-004-tournaments-are-odd-functions
script: 04-computation/regular_group_orbit_tournament_parity_thm3447.py
output: 05-knowledge/results/regular_group_orbit_tournament_parity_thm3447.out
script_sha256: 8868aa940cbde11b7b79b5a5c4f805622dca3b1486c4d82e87a0566917330099
output_sha256: d772b7947f59c071d22035ff20c19075a3093fab6b708b8caf5f62339086fae5
semantic_sha256: 0672eab70672e58ff4f44298e7547bb44b9843f0f8f47987374fd03895b142c9
hash_basis: LF-normalized bytes
---

# THM-3447 -- regular group-orbit tournaments exist exactly at odd order

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

The inverse-orbit proof and all finite invoices below have passed an
independent reconstruction and replay audit.

## 1. Exact classification

Let a finite group `G` act regularly on a finite set `Omega`.  A loopless
directed relation may have, on each unordered pair, one of four states:

```text
missing,       forward only,       backward only,       bidirected.
```

Call it invariant when the diagonal `G`-action preserves all arcs.  Put

```text
N=|G|,
iota(G)=#{h in G\{1}:h^2=1},
q(G)=(N-1-iota(G))/2.                                  (1)
```

Then:

1. invariant loopless directed relations are counted by

   ```text
   4^q 2^iota=2^(N-1);                                  (2)
   ```

2. invariant partial orientations (missing or one-way, never bidirected)
   and invariant semicomplete relations (one-way or bidirected, never
   missing) are each counted by

   ```text
   3^q;                                                  (3)
   ```

3. invariant symmetric relations and maximal invariant partial orientations
   are counted respectively by

   ```text
   2^(q+iota),             2^q;                         (4)
   ```

4. invariant tournaments are counted by

   ```text
   2^((N-1)/2),       if N is odd,
   0,                 if N is even.                     (5)
   ```

Thus a regular group orbit admits an invariant tournament **if and only if**
the group has odd order.  This is valid for nonabelian groups as well as
Cayley circulants.

For every involution `h`, the pairs

```text
{x,xh},              x in G,                            (6)
```

form a perfect matching of `Omega`.  In an invariant directed relation every
edge of `(6)` is forced to be either missing or bidirected; it cannot be
singly oriented.  Matchings from distinct involutions are edge-disjoint, so
exactly

```text
N iota(G)/2                                             (7)
```

unordered pairs are forced non-tournament pairs.  If `N>1` is odd, every
invariant tournament contains a directed triangle.

## 2. Proof by inverse-pair orbits

Choose `o in Omega` and identify `g in G` with `g o`.  Diagonal left
translation preserves the difference

```text
d(x,y)=x^(-1)y.                                         (8)
```

Consequently an invariant directed relation is uniquely determined by a
connection set

```text
S subseteq G\{1},       x->y iff d(x,y) in S.           (9)
```

Reversing an ordered pair sends `d` to `d^(-1)`.  Inversion partitions
`G\{1}` into `q` two-element orbits and `iota(G)` one-element orbits.  A
two-element orbit permits all four states, giving the factor `4^q` in `(2)`;
an involution permits only missing or bidirected, giving `2^iota`.  For a
partial orientation the two-element orbit permits missing or either direction,
while an involution is forced missing.  For a semicomplete relation the same
three choices remain and an involution is forced bidirected.  This proves
`(2)--(3)`.

For a symmetric relation, each two-element inverse orbit is either wholly in
`S` or wholly outside it, while every involution is independently in or out;
this gives `2^(q+iota)`.  A partial orientation is maximal among invariant
partial orientations exactly when every two-element inverse orbit contributes
one of its two directions and every involution is absent, giving `2^q`.
This proves `(4)`.

A tournament requires exactly one element from every inverse pair and none
from a singleton inverse orbit.  It therefore exists exactly when
`iota(G)=0`, and then has the count in `(5)`.  By Cauchy's theorem a finite
group has a nonidentity involution exactly when its order is even.  This proves
the parity criterion, including the nonabelian case.

For fixed `h=h^(-1)`, right multiplication by `h` is fixed-point-free and
involutive, proving `(6)`.  An unordered edge determines its inverse orbit, so
different involutions give disjoint matchings and `(7)` follows.  Finally, a
transitive tournament has distinct outdegrees `0,...,N-1` and hence trivial
automorphism group.  A nontrivial regular `G`-action is an automorphism group,
so an invariant tournament with `N>1` is not transitive and therefore contains
a directed triangle.

Changing the chosen basepoint conjugates the difference coordinates and hence
conjugates `S`.  Thus the classification and counts are intrinsic, but a
particular half-set is not basepoint-independent unless the required conjugacy
condition is separately imposed.  Requiring both left- and right-translation
invariance is stronger still: `S` must be a union of conjugacy classes.

## 3. Hensel-lattice specialization

THM-3444 gives, on each of its free orbits at depth gap `n=a-c>0`, a regular
action of

```text
G_(a,c)=(Z/p^n Z)^r,          N=p^(rn).                (10)
```

The present theorem does not identify Hensel points with LRC endpoints; it
only classifies invariant binary relations on one such group orbit.

- If `p` is odd, `G_(a,c)` has odd order.  Each orbit admits exactly

  ```text
  2^((p^(rn)-1)/2)                                     (11)
  ```

  invariant tournaments.  Choosing one requires an inverse half-set—an
  orientation gauge not supplied by the vector fields.
- If `p=2`, the group has exactly `2^r-1` nonidentity involutions.  Every
  invariant generalized relation has

  ```text
  p^(rn)(2^r-1)/2 = 2^(rn-1)(2^r-1)                   (12)
  ```

  forced missing/bidirected pairs, so no invariant tournament exists.
- At the first two-adic lift `n=1`, every nonzero exponent is an involution.
  Formula `(12)` becomes `binom(2^r,2)`: **every pair is forced**.  This is
  precisely the XOR swap `x -> x xor (x xor y)` behind THM-2504's affine-cube
  no-go.

When THM-3444 has `r=d`, one orbit is the whole Hensel fibre.  For `r<d`, the
theorem applies separately to each of the `p^((d-r)n)` orbits; cross-orbit
relations require the orbit-bank labels that THM-3444 deliberately retains.

## 4. Sizes four and six

The regular group structure matters quantitatively, but parity already kills
all tournaments:

| regular group | involutions | forced non-single pairs | free inverse classes | max singly oriented pairs | generalized relations |
|---|---:|---:|---:|---:|---:|
| `C4` | 1 | 2 | 1 | 4 | 8 |
| `V4=F2^2` | 3 | 6 (every pair) | 0 | 0 | 8 |
| `C6` | 1 | 3 | 2 | 12 | 32 |
| `S3` | 3 | 9 | 1 | 6 | 32 |

Ordinary tournaments certainly exist on four and six labeled vertices.  What
cannot exist is one invariant under a regular group action on those vertices.
Allowing missing or both-way edges does not erase the obstruction; it records
it exactly as the involution matchings.  At odd sizes, the obstruction
disappears but a half-set gauge remains: for example `C3` has two invariant
tournaments, the two directed 3-cycles.

The power-set realization makes the even obstruction maximally transparent.
Under symmetric difference, `P([r])` is `F_2^r`; for distinct subsets `A,B`,
translation by `H=A Delta B` swaps them.  Hence every unordered pair is an
involution pair.  The same argument holds on `P(N)`, viewed as a Boolean group.
Representing a subset by its harmonic weighted indicator

```text
(1_A(n)/n)_(n>=1)
```

adds information but not an invariant orientation.  For every finite cutoff
`M`, its scalar partial-sum weight obeys

```text
w_M(A Delta H)=w_M(A)+w_M(H)-2w_M(A intersect H).       (13)
```

The intersection term breaks XOR-translation invariance, and equal weighted
sums can still tie.  Thus harmonic weighting is a useful sidecar, not a
canonical repair of the Boolean tournament obstruction; for infinite subsets
the total scalar sum may converge or diverge, while the weighted-indicator
sequence always remains defined.

## 5. Connection and loss ledger

| field | exact content |
|---|---|
| vertices | one regular finite-group orbit |
| pairwise observable | the group difference `x^(-1)y` |
| orientation gauge | one representative from each inverse pair |
| forced ties | missing/bidirected matchings indexed by involutions |
| preserved | diagonal group action, converse, and the four pair states |
| destroyed by tournament completion | involution data and the distinction between missing and bidirected |
| Hensel sidecar | orbit-bank label when `r<d`; exponent basis/half-set for an actual orientation |
| cheapest odd positive | `C3`, giving the two cyclic orientations |
| cheapest even hostile | `C2`; `V4` makes every pair hostile simultaneously |

This extends LEM-004's cyclic odd-function dictionary and isolates the general
mechanism in THM-2504's XOR no-go.  It does not make arbitrary XOR translations
physical LRC symmetries, turn a Hensel orbit into an owner/current, or supply a
canonical tournament.  No LRC row, `JC(2)`, or D5 bridge follows.

## 6. Exact companion

The deterministic companion validates finite group tables, inverse-orbit and
matching invoices, exhausts every connection set for `C3,C4,V4,C5,S3,C6,C7,
F2^3,C9`, and supplies the nonabelian odd control `Heis(F3)` of order `27`.
It checks the size-four/six rows and two-adic/odd-prime Hensel formulas.  Normal
and optimized runs reproduce the stored transcript byte-for-byte after LF
normalization:

```text
python -B 04-computation/regular_group_orbit_tournament_parity_thm3447.py
python -B -O 04-computation/regular_group_orbit_tournament_parity_thm3447.py
```

The independent audit reconstructed the inverse-orbit proof, checked the
size-four/six and Boolean boundaries, and replayed the exact companion in both
modes against the stored transcript.  This completes the proof.  QED.
