---
id: THM-3515
title: "U_full all-role endpoint weighted-tree spectral closure in B1"
status: >
  PROVED STRUCTURAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  In the fixed
  THM-3514 U_full owner-boundary atom bank, the eight relation roles occupy
  five distinct refined character classes.  Every class is nonzero on all
  52 active chamber/drift addresses, and the resulting 5x52 response matrix
  has rank five.  For every one of the 72 lawful private-support charts, the
  forced bridge, both K4 weighted spanning-tree factors, and their product
  are nonzero separately at every address, every drift marginal, every
  drift-Fourier mode, and every chamber-resolved Fourier fibre.  The edge
  responses are differences of endpoint potentials and therefore lie in
  graph B1; this is not ancestry, physical current, absolute H1, a grouped
  coefficient, scalar-row exclusion, or LRC(14).
source: codex/ufull-all-role-independent-audit/2026-08-16
audit: >
  Independent five-character contraction, 72-chart generator, modular
  matrix-tree evaluator, four exact bank reconstructions, pointwise and
  modewise rank tests, flat-potential and killed-bridge hostiles, and
  normal/optimized/stored replay agreement.  The audit imports only the
  pinned independent THM-3514 atom-table engine and does not import the
  all-role candidate script.
depends_on:
  - THM-3514-ufull-owner-boundary-k4xf13-endpoint-factorization-and-walsh-spectrum
related:
  - THM-3479-literal-half-twist-relation-current-two-transplant-certificate
  - THM-3482-private-count-gradient-weighted-spectral-closure-without-absolute-h1-flux
script: 04-computation/lrc_ufull_all_role_weighted_tree_spectral_closure_independent_audit_20260816.py
output: 05-knowledge/results/lrc_ufull_all_role_weighted_tree_spectral_closure_independent_audit_20260816.out
script_sha256: 7b2a477845062bf2f40cc6570063edc9b4220ba1945ce71f4f1714a6ea471dba
output_sha256: 3d7eabca4ab7e175eca9118a0e8dd0eabf3d867dcd8e188477b14c5c9934004e
semantic_sha256: 61db05458acb26cd8d6a9a89991d3bbeed91f95f6bdf035720a3a93a5b4f886b
hash_basis: LF-normalized UTF-8 bytes; exact finite-field semantic ledger
---

# THM-3515 -- all lawful endpoint charts close in every tested address and drift mode

**PROVED STRUCTURAL + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

## 1. Fixed universe and five role classes

Use exactly the THM-3514 endpoint atom bank over its certified split field

```text
F_P,  P=572252886246508880869,
zeta=505438565698892403012,       ord(zeta)=13.        (1)
```

The active address carrier is

```text
X={L,R} x {L,R} x F_13,             |X|=52,           (2)
```

where the two bits are the left and right owner-boundary chambers and the
field coordinate is relative guard-sheet drift.  The eight relation roles
have the following refined character classes:

```text
c1          -> (0, 0, 0),
c3          -> (0, 1, 0),
c2,q3,q4,q5 -> (1, 0, 0),
H           -> (1, 0, 1),
q2          -> (1,12, 0).                             (3)
```

These are five distinct classes, with inverse-transform values in `F_P`

```text
(0, 0, 0) -> 405336876493642499425,
(0, 1, 0) -> 518539850465495448196,
(1, 0, 0) -> 503604956476841920373,
(1, 0, 1) -> 320618948602619577408,
(1,12, 0) ->  15703541686881447885.                  (4)
```

Let `B_q(C,D,d)` denote the THM-3514 atomwise pair contraction for class
`q`, chambers `C,D`, and drift `d`.  The independent audit reconstructs this
quantity directly from the pinned owner-boundary atom tables and the two
guard kernels.  For every one of the five classes,

```text
B_q(C,D,d) != 0       for all (C,D,d) in X.           (5)
```

Moreover the five rows of the `5 x 52` matrix `(B_q(x))` have rank five.
Thus no refined role class disappears and no nontrivial linear combination
of all five response rows vanishes on the complete active address carrier.

## 2. The 72 lawful charts and the weighted-tree factorization

The private-support graph has vertices `0,...,7`, hub `4`, leaf `6`, and
edge set

```text
03 04 05 12 14 17 24 27 34 35 45 46 47.             (6)
```

Its two four-vertex blocks are

```text
K_L={0,3,5,4},              K_R={1,2,7,4},            (7)
```

and edge `46` is the forced bridge from the hub role `H` to the leaf role
`q5`.  A lawful chart fixes `H` at vertex `4` and `q5` at vertex `6`, places
the three blocker roles `c1,c2,c3` on one block, places `q2,q3,q4` on the
other block, and permits either block assignment and arbitrary permutations.
Consequently there are exactly

```text
2 * 3! * 3! = 72                                    (8)
```

lawful labelled charts.

At an address or transformed mode, let `p_v` be the role response assigned
to vertex `v`, orient every graph edge in the fixed owner order, and put

```text
w_(uv)=p_u-p_v.                                      (9)
```

For a chart `chi`, let `Tau_L(chi)` and `Tau_R(chi)` be the weighted
spanning-tree polynomials of the two labelled `K4` blocks, equivalently any
reduced `3 x 3` Laplacian determinant.  The reduced determinant of the full
two-block graph factors structurally as

```text
Delta_chi=(H-q5) Tau_L(chi) Tau_R(chi).              (10)
```

This follows either by the bridge form of the matrix-tree theorem or by
splitting every spanning tree uniquely into the forced bridge and one
spanning tree in each block.

## 3. Maximal pointwise and modewise nonvanishing

For each representation below and every lawful chart, the audit evaluates
the ordered quadruple

```text
(H-q5, Tau_L, Tau_R, Delta_chi).                     (11)
```

All four entries are nonzero in every case:

```text
representation                    fibres       charts    zero counts
---------------------------------------------------------------------
individual address (C,D,d)          52            72      (0,0,0,0)
sum over four chambers at fixed d    13            72      (0,0,0,0)
Fourier transform in d               13            72      (0,0,0,0)
fixed (C,D), Fourier transform in d  52            72      (0,0,0,0). (12)
```

At each of the thirteen physical drift slices, the `5 x 4` matrix of role
class against chamber pair has rank four.  The same rank-four statement
holds at every drift-Fourier frequency.  Hence neither the drift domain nor
its character domain contains a rank-deficient chamber lane.

The four exact banks have the independently reproduced digests

```text
pointwise:
a1d5183f37cc39deed976876ed91132f662f91137f41a10a79fd3e974ced2dfc
drift marginal:
6bcdd9fa616ba2dc4883d8334dca7358d06375d7342f43c4aac8dea45deb5027
drift Fourier:
e96a3859e6838fa91d04a55ca55fa4e44cafa45bd86b93c7f6b612b1e65f4dc1
chamber-resolved Fourier:
c9ca8c9b7ed7000b00e39b496523ad75a44379d84ea35e49163d6b64e4196e73. (13)
```

Nonzero reduction in this certified split field proves that the
corresponding characteristic-zero cyclotomic quantities are nonzero.  Rank
five or four after reduction proves that a corresponding characteristic-zero
minor is nonzero.

## 4. Why the closure is in `B^1`, not absolute `H^1`

For every address or Fourier mode, equation `(9)` is the graph coboundary of
the eight vertex potentials:

```text
w=delta p in B^1(G;F_P).                             (14)
```

Therefore its pairing with every graph cycle is zero and its image in
absolute graph `H^1` is the zero class.  The nonzero weighted tree
determinants in `(10)` depend on the actual coboundary representative; they
do not descend to absolute cohomology.  This is the same structural boundary
isolated abstractly by THM-3482, now verified on the complete U_full endpoint
address and frequency banks.

The chamber-pair states form an undirected `K4` carrier.  No intrinsic
pairwise observable or orientation is supplied, so neither that carrier nor
the two private-support `K4` blocks are tournaments.

## 5. Hostiles and exact failure boundary

Two hostile controls identify the source of the signal.

1. If all eight vertex potentials are equal, every edge response, bridge,
   tree factor, and full product is zero.
2. If the `q5` response is replaced by the `H` response, the forced bridge
   and every full product are zero, independently of the two `K4` factors.

Thus the nonvanishing in `(12)` is neither a constant-potential artifact nor
a determinant convention that silently omits the bridge.

The audit also preserves the THM-3514 root-label boundary.  The abstract
guard sheets and collision-root labels are isomorphic regular `F_13` torsors,
but no support relation on a common Boolean ancestry base is constructed.
Cartesian endpoint pairing loses the physical root, common base, owner,
word, source sheet, horizon, chronology, and arrival data.

Accordingly this theorem proves endpoint weighted-tree spectral closure in
`B^1` only.  It proves none of

```text
THM-2471 ancestry or a common-stalk support map;
a physical or lawful Boolean current;
a nonzero absolute H^1 class;
a grouped exact-address coefficient or all-unit projector;
a bispectrum theorem or scalar-row exclusion;
LRC(14).                                             (15)
```

## 6. Independent verification record

The audit imports only the LF-hash-pinned independent THM-3514 atom engine.
It does not import the candidate all-role probe.  It separately implements
the five character contractions, 72-chart enumeration, modular Gaussian
determinants of the two reduced Laplacians, all four bank transforms, and the
rank tests.  It reproduces all four candidate bank digests while producing
an independent semantic ledger.

Run from the repository root:

```text
python -B 04-computation/lrc_ufull_all_role_weighted_tree_spectral_closure_independent_audit_20260816.py
python -B -O 04-computation/lrc_ufull_all_role_weighted_tree_spectral_closure_independent_audit_20260816.py
```

Both executions reproduce the stored output exactly after LF normalization.
The pinned semantic digest is

```text
61db05458acb26cd8d6a9a89991d3bbeed91f95f6bdf035720a3a93a5b4f886b. (16)
```
