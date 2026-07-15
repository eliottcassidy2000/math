---
id: THM-781
title: Hamiltonian-path inversion of merged tournament-tiling fibres
status: PROVED (general forward quotient, direct inverse, and H/Aut fibre formula; exact bidirectional audit of all 33,866 explorer tilings at n=3..7)
source: codex-2026-07-14-S8
depends_on:
  - HYP-6825
related: [THM-550, THM-646, THM-773, THM-778]
verification: 04-computation/tournament_tiling_metagraph_fibre_inverse_codex_S8.py
  (+ 05-knowledge/results/tournament_tiling_metagraph_fibre_inverse_codex_S8.out and .json)
---

# THM-781 — Hamiltonian-path inversion of merged tournament-tiling fibres

## 1. The two functions

Fix `n`.  Put `m=binom(n-1,2)` and let

```text
Q_n={0,1}^m
```

be the explorer's staircase tilings.  Its declared observer path is

```text
P_*=(0,1,...,n-1).
```

For `t in Q_n`, let `tau_n(t)` be the labelled tournament obtained by fixing
the arcs of `P_*` forward and reading the other `m` arcs from the tile bits.
Let `kappa(T)` be a canonical code for the isomorphism class of `T`, let
`T^op` reverse every arc, and let `A_n` send a converse orbit of canonical
classes to its objective HYP-6825 merged-node address.

The requested forward function is

```text
f_n(t)
 = A_n({kappa(tau_n(t)), kappa(tau_n(t)^op)}).            (1)
```

Equivalently, with the committed total ordering on canonical codes, the inner
merged-orbit key is

```text
min(kappa(tau_n(t)), kappa(tau_n(t)^op)).                 (2)
```

The node rank in (1) is not the raw minimum in (2).  HYP-6825 orders merged
nodes first by local flip depth/path, then by rooted weighted blue/black line
structure, combined refinement, deletion ancestry, and only finally by the
canonical orbit code.  Formula (2) identifies the node; `A_n` gives its
structural address `n-aXXX`.

For the inverse, let a merged node `u` carry the one or two unmerged classes

```text
C(u)={c}                         if c=c^op,
C(u)={c,c^op}                    otherwise.              (3)
```

Choose a representative tournament `T_c` for each `c`.  For every directed
Hamiltonian path

```text
p=(p_0,...,p_(n-1)) in HP(T_c),
```

relabel `p_i` to the fixed position `i` and read every nonconsecutive arc as a
tile bit.  Call the resulting tiling `N_c(p)`.  The requested inverse function
is the set-valued map

```text
g_n(u)=union_(c in C(u)) {N_c(p):p in HP(T_c)}.           (4)
```

It is independent of the chosen representatives and is exactly the fibre:

```text
g_n(u)=f_n^(-1)({u}).                                    (5)
```

Thus tiling-to-node is a function, while node-to-tiling is necessarily a
set-valued inverse.  The merged nodes partition `Q_n` into these fibres.

## 2. Direct inverse and exact fibre cardinality

`Aut(T_c)` acts freely on `HP(T_c)` by applying the automorphism to every path
vertex.  Two paths normalize to the same fixed-path tiling precisely when they
lie in the same automorphism orbit.  Hence

```text
HP(T_c)/Aut(T_c)  <->  {fixed-path tilings in class c},   (6)
```

and every unmerged class contributes exactly

```text
|F_c|=H(T_c)/|Aut(T_c)|.                                 (7)
```

Consequently

```text
|g_n(u)| = sum_(c in C(u)) H(T_c)/|Aut(T_c)|.             (8)
```

Converse classes have the same Hamiltonian-path and automorphism counts.  A
non-self-converse merged node therefore has twice the contribution of either
class, while a self-converse node has only one contribution:

```text
|g_n(u)| = H(T)/|Aut(T)|       for a self-converse class,
|g_n(u)| = 2H(T)/|Aut(T)|      for a converse pair.       (9)
```

This recovers the known odd/even checksum: tournament Hamiltonian-path counts
and tournament automorphism groups have odd order, so every unmerged class
fibre is odd, self-converse merged fibres are odd, and non-self-converse merged
fibres are even.

## 3. Proof

Formula (1) is well defined because canonicalization is constant on labelled
isomorphism classes and converse commutes with relabelling.  Converse is an
involution, so (2) is a canonical key for the unordered orbit (3).

For a representative `T=T_c` and `p in HP(T)`, relabelling the vertices in
path order makes every consecutive arc point forward.  The remaining
`binom(n,2)-(n-1)=binom(n-1,2)` arcs therefore give one and only one explorer
tiling `N_c(p)`.  Its tournament is isomorphic to `T`, so (4) is contained in
the fibre in (5).

Conversely, suppose a fixed-path tiling `t` has class `c`.  Choose an
isomorphism from `T_c` to `tau_n(t)`.  The inverse image of `P_*` is a directed
Hamiltonian path `p` of `T_c`, and normalization along `p` recovers `t`.
Thus the unmerged version of (4) is onto the complete class fibre; taking the
one or two classes in (3) proves (5).

If paths `p,q in HP(T)` normalize to the same tiling, their two path-order
relabelled tournaments are equal.  The map `p_i -> q_i` therefore preserves
every arc and is an automorphism of `T`.  Conversely every automorphism gives
the same normalized tiling.  The action is free because an automorphism fixing
the `n` vertices of one Hamiltonian word is the identity.  Every normalization
fibre consequently has size `|Aut(T)|`, proving (6)--(9). ∎

## 4. Executable API

The companion module exposes the literal maps:

```python
# Run from 04-computation/, or place that directory on PYTHONPATH.
from tournament_tiling_metagraph_fibre_inverse_codex_S8 import MetagraphFibreAtlas

atlas = MetagraphFibreAtlas()

atlas.tiling_to_node(7, 615)
# node_id = "n7-a267", fibre_index = 0

atlas.node_to_tilings(7, "n7-a267")
# the complete ordered list of 25 tiling records
```

`tiling_to_node` computes `tau`, canonicalizes both orientations, merges the
converse orbit, and checks the result against the stored O(1) lookup.
`node_to_tilings` decodes only the node's one or two canonical tournament
representatives and enumerates their Hamiltonian paths.  It does **not** scan
the `2^m` tiling cube.

The inverse ordering remains HYP-6825's declared fibre order

```text
(Hamming weight, unmerged canonical class code, explorer mask).             (10)
```

The browser explorer now exposes the corresponding JavaScript functions
`tilingToMergedNode(mask,n)` and `mergedNodeToTilings(node,n)` for its
`n<=6` in-browser atlas.  The Python API continues exactly through `n=7`.

## 5. Exact audit through n=7

The independent forward formula and direct inverse were checked against every
committed atlas entry:

| n | tilings | classes | merged nodes | class fibre range | node fibre range |
|---:|---:|---:|---:|---:|---:|
| 3 | 2 | 2 | 2 | 1--1 | 1--1 |
| 4 | 8 | 4 | 3 | 1--5 | 1--5 |
| 5 | 64 | 12 | 10 | 1--13 | 1--13 |
| 6 | 1,024 | 56 | 34 | 1--43 | 1--86 |
| 7 | 32,768 | 456 | 272 | 1--159 | 1--306 |

Across all `33,866` tilings, `530` unmerged classes, and `321` merged nodes:

```text
forward-map mismatches                         0
direct class-inverse mismatches                0
Hamiltonian-path/automorphism multiplicity     0
H/Aut class-size formula mismatches            0
merged sum-formula mismatches                  0
node inverse round-trip mismatches             0
```

The prime-seven heptagon node receives a conceptual explanation:

```text
n7-a267: H(T)=175, |Aut(T)|=7, |g_7(n7-a267)|=175/7=25.   (11)
```

These are exactly the 25 masks reached by all `5,040` owner-to-sheet
assignments in THM-773.  They are the automorphism orbits of Hamiltonian
observer cuts on the cyclic heptagon tournament, not an arbitrary collection
created by the explorer.

## 6. Preservation boundary

Equations (1) and (4) solve the exact combinatorial correspondence for the
explorer's binary fixed-path tilings.  They do not make the node a complete LRC
state.

The forward quotient forgets which Hamiltonian-path orbit/tiling was chosen.
The inverse repairs exactly that fixed-cut ambiguity by returning the whole
fibre.  It still does not recover:

- an owner-to-vertex or owner-to-sheet assignment;
- the metric gaps or wall position that realized the tournament;
- inverse winding steps, endpoint ties, and centered-CF phase;
- the chronological next-wall owner and collision-hop target;
- the global sheet carry or core-safe component.

Moreover, a fixed-path tiling is not the same object as a particular
blue/black complement **line orbit**: several simultaneous-isomorphism line
orbits can project to the same endpoint-node incidence.  The correct LRC
carrier therefore remains

```text
merged node base
  + Hamiltonian-path/tiling fibre
  + owner-labelled metric transport stalk.               (12)
```

## Tournament Analysis and assumption challenge

- **Vertices:** directed Hamiltonian paths of a class representative.  Runner,
  gap, wall-event, and proof-obligation vertices were also considered; paths
  are the exact carrier for this inverse problem.
- **Pairwise observable:** whether two paths normalize to the same explorer
  tiling.
- **Switch/gauge:** quotient by the free `Aut(T)` action; equivalently, choose
  one lexicographically least path in each orbit.
- **Tie Hamiltonian path:** lexicographic order on the path words.
- **Fingerprint:** `H(T)` path vertices, uniform equivalence blocks of size
  `|Aut(T)|`, and `H(T)/|Aut(T)|` quotient vertices.
- **Destroyed by the isomorphism node:** the path orbit, fixed observer cut,
  and tiling mask; restored set-valuedly by (4).

The challenged assumption is that the inverse should enumerate all tilings or
that a node should select one preferred tiling.  The intrinsic inverse object
is instead `HP(T)/Aut(T)`, and its union over the converse orbit.
