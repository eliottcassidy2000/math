---
id: THM-4065
title: "Rule 30 temporal-cylinder transfer, infinite zero columns, and next zero"
status: >
  PROVED universal temporal-cylinder bijection and infinite isolated-zero
  theorem + FINITE-EXACT physical next-zero certificate + VERIFIED-EXACT +
  INDEPENDENTLY HOSTILE-AUDITED. For every cyclic temporal period P, the
  native spatial transfer T_P(A,B)=(B,X) is a bijection from pairs with
  B nonzero to pairs with first word nonzero. Its finite partial graph is a
  disjoint union of interior cycles and boundary-to-boundary paths. Applied
  to the eventually dyadic-periodic Rule-30 left-front columns, every zero
  column starts such a path and therefore has a later zero, at gap at least
  three and at most J(P)+2 when its successor has least period P, where J(P)
  is the exact count of nonzero word-pairs of joint least period P. Hence
  there are infinitely many isolated eventual-zero columns. The THM-4047
  physical boundary state at depth 87866 follows a period-32 path whose next
  zero is exactly depth 1420878968, with gap 1420791102. Moving-center
  behavior, balance, query lower bounds, and all Rule-30 prizes remain OPEN.
source: codex-frontier-synthesis-creative-20260825b / Rule-30 operation lane
audit: >
  PASS. The primary C path solves each cyclic word from its first reset bit;
  the no-import C audit composes affine monodromy from residue zero and then
  materializes predecessor/successor graphs. Both exhaust every state through
  P=8, replay the 12133-edge path from the last THM-4047 zero to the frozen
  depth-100000 state, and independently walk 1420778968 further transfers.
  Seven checkpoints and three full-orbit digests agree. A Python parent audit
  reconstructs all 100001 THM-4047 tails from the separately audited seven-bit
  physical query packet in 30744208 checks. Normal/optimized streams
  byte-match the frozen outputs; the Python audit has no assert nodes or float
  literals. A third hostile referee used direct packed physical evolution and
  a pointer-doubling cyclic solver and recovered the same parent words,
  endpoint, indexing, gap, terminal predecessor, and two successor phases.
depends_on:
  - THM-4047-rule30-left-front-affine-monodromy-clock
related:
  - THM-3458-rule30-right-edge-2-adic-odometer-and-moving-observer-boundary
  - THM-4048-rule30-periodicity-balance-and-model-firewalls
  - THM-4050-rule30-half-arc-marked-cylinder-and-radius-nine-hostile
  - THM-4064-rule30-cyclotomic-kernel-character-and-c60-alias-obstruction
script: 04-computation/rule30_temporal_cylinder_transfer_thm4065.c
output: 05-knowledge/results/rule30_temporal_cylinder_transfer_thm4065.out
independent_audit_script: 04-computation/rule30_temporal_cylinder_transfer_thm4065_independent_audit.c
independent_audit_output: 05-knowledge/results/rule30_temporal_cylinder_transfer_thm4065_independent_audit.out
parent_audit_script: 04-computation/rule30_temporal_cylinder_transfer_thm4065_parent_audit.py
parent_audit_output: 05-knowledge/results/rule30_temporal_cylinder_transfer_thm4065_parent_audit.out
script_sha256: 4f5b43bc39bc40924d7ac3d0fbdade44083a4cbb251568466e99d3cfd3a25137
output_sha256: e4615c0d0b0dc9203c97a58b47a31d5a88d5d98b472d15f4445378b9c5c50964
independent_audit_script_sha256: d4cff42bd4d4a412353ba37c676d301df1e98e4f7fb28d7db2c0536ab0b1a13e
independent_audit_output_sha256: df236987fe89f88df2c0d365fa57368d4ff6e5ff0b63fc293af3328f61b69de4
parent_audit_script_sha256: 387fa9fae0f51c5af8b57e205a1aaa5c5e54dd04fb598726be811a4c87dfc94f
parent_audit_output_sha256: 6ce62b4807c430233db0577b393a8636377f5e1850280471eab3d66018c9517b
orbit_semantic_sha256: >
  checkpoints plus xor 25469bdaa1c33b27, sum 65b743c483016d28,
  and FNV 46485ada15e2f5b9
hash_basis: raw LF bytes
---

# THM-4065 -- the native spatial operation closes the zero clock

**PROVED in the universal scopes below; FINITE-EXACT for the named physical
Rule-30 endpoint.** THM-4047 compiled each fixed left-front column separately.
The missing operation was to move once in source depth while retaining the
entire temporal word. On a cyclic time cylinder that operation is reversible
away from one boundary, and the boundary cannot feed a cycle. This turns every
eventual-zero column into the source of a finite path ending at another zero.

## 1. The temporal-cylinder transfer

Let

```text
W_P=F_2^(C_P)
```

be the `P`-bit words indexed by absolute time residue. For `A,B in W_P`, seek
`X in W_P` satisfying

```text
X_(t+1)=A_t+B_t+(1+B_t)X_t,             t in C_P.    (1)
```

This is exactly THM-4047's fixed-left-front recurrence, with consecutive
spatial columns `(A,B,X)`.

### Theorem 1.1 (existence, uniqueness, and inverse)

If `B!=0`, equation `(1)` has a unique solution `X`. Consequently

```text
T_P:W_P x (W_P minus {0}) -> (W_P minus {0}) x W_P,
T_P(A,B)=(B,X)                                      (2)
```

is a bijection, with inverse

```text
A_t=X_(t+1)+B_t+(1+B_t)X_t.                         (3)
```

Choose `s` with `B_s=1`. Equation `(1)` at `s` forgets the incoming bit and
forces `X_(s+1)=A_s+1`. Propagating through the other `P-1` residues determines
every bit. The final return is compatible because the update at `s` is
independent of `X_s`. This proves existence and uniqueness; `(3)` proves the
bijection.

The transfer commutes with cyclic time shift. It preserves the **joint** least
period of the pair: a common smaller period of the output pulls back through
`(3)`, and one of the input pushes forward by uniqueness. It need not preserve
either individual word period. The minimal hostile is

```text
P=2,       A=B=(1,0),       X=0.                    (4)
```

All later uses therefore concern an ambient common period or the pair's joint
least period, never an unjustified individual-period equality.

## 2. The finite transfer graph

Make `(2)` a partial directed graph on `W_P x W_P`. A vertex has an outgoing
edge exactly when its second coordinate is nonzero, and an incoming edge
exactly when its first coordinate is nonzero. The bijection gives indegree and
outdegree at most one. There are four types:

```text
(A,B), A,B!=0       interior;
(0,B), B!=0         source boundary;
(A,0), A!=0         sink boundary;
(0,0)               isolated.                       (5)
```

Every finite component is therefore either an interior directed cycle or a
source-to-sink path. A source cannot enter a cycle: pulling a repeated state
backward by its unique predecessor would give the source an incoming edge.
The `2^P-1` disjoint source paths end at all `2^P-1` sinks and define a
permutation of the nonzero words.

If a source word has least period `P`, every state on its path has joint least
period `P`. Mobius inversion therefore gives the sharper interior-state budget

```text
J(P)=sum_(d divides P) mu(P/d)(2^d-1)^2.             (6)
```

The path has at most `J(P)+1` transfer edges. In particular this is at most
`(2^P-1)^2+1`. For dyadic `P>1`, the exact count simplifies to

```text
J(P)=(2^P-1)^2-(2^(P/2)-1)^2.                       (7)
```

A path cannot have one edge: `T_P(0,B)=(B,0)` would put `X=0` in `(1)` and
force `B=0`. Thus every boundary path has at least two edges.

The exhaustive controls are:

| `P` | boundary paths | maximum path edges | interior cycles | cycle vertices | maximum cycle |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 2 | 0 | 0 | 0 |
| 2 | 3 | 4 | 1 | 2 | 2 |
| 4 | 15 | 28 | 2 | 14 | 12 |
| 8 | 255 | 667 | 13 | 5510 | 2816 |

These finite graphs audit all vertex types and the explicit inverse; the
general result is the algebraic bijection, not a finite extrapolation.

## 3. Infinitely many isolated eventual-zero columns

Use THM-4047's single-seed left-front columns

```text
ell_r(t)=a_t(-t+r),
ell_r(t+1)=ell_(r-2)(t)+ell_(r-1)(t)
              +(1+ell_(r-1)(t))ell_r(t).             (8)
```

Every fixed column is eventually dyadic-periodic. Call `r` an eventual-zero
depth if its eventual word is zero.

Adjacent zero depths are impossible. If `ell_r` and `ell_(r+1)` were
eventually zero, `(8)` at depth `r+1` would force `ell_(r-1)` to be zero.
Descending would contradict `ell_0(t)=1`.

Let `z` be a zero depth and let `B` be the eventual word of `ell_(z+1)`. It is
nonzero. Take its least period `P`; the pair begins at the source

```text
(ell_z,ell_(z+1))=(0,B).                             (9)
```

As long as the second word is nonzero, THM-4047's reset case and Theorem 1.1
identify the next eventual word with one application of `T_P`. The source
component must reach a sink `(A,0)`, hence a later zero `z'`. If the path has
`m` edges, spatial indexing gives

```text
z'-z=m+1,
3 <= z'-z <= J(P)+2 <= (2^P-1)^2+2.                (10)
```

The new zero again has a nonzero successor, so the construction repeats.
Beginning with THM-4047's zero at depth two proves that there are infinitely
many eventual-zero columns. The lower bound shows every zero has at least two
intervening nonzero columns.

The endpoint also recovers THM-4047's clock. The source pair has joint least
period `P`; hence the sink pair `(A,0)` does too, so `A` has least period `P`.
The successor `Y` after the new zero obeys

```text
Y_(t+1)=A_t+Y_t.                                    (11)
```

If `A` has odd weight, `Y` has exact least period `2P`. If its weight is even,
there are two complementary `P`-periodic phases; either has least period `P`
because its difference word is `A`. The physical choice is one phase bit.

## 4. The exact next zero after depth 87866

THM-4047 certifies

```text
Z intersect [0,100000]
 =(2,7,28,399,53207,58286,87866).                   (12)
```

Its absolute-phase period-32 boundary word is

```text
(ell_87866,ell_87867)=(00000000,cf3030cf).           (13)
```

Exactly `12133` transfers, with no intervening zero, give

```text
(ell_99999,ell_100000)=(be79924b,90f58380).          (14)
```

Starting from `(12)`, both exact C implementations take `1420778968` further
transfers before the second coordinate is first zero. Therefore

```text
z_8=100000+1420778968=1420878968,                    (15)
z_8-z_7=1420878968-87866=1420791102.                 (16)
```

Equivalently, the full boundary path from `(13)` has `1420791101` edges.
Every intermediate second word is checked nonzero, so `(15)` is the **next**
zero after `87866`, not merely another zero.

The terminal sink is

```text
(dae372fa,00000000).                                 (17)
```

The predecessor has weight `20` and least period `32`. Thus this is an
identity, not a doubling, event. The two possible successor words are

```text
93425cac,       6cbda353,                            (18)
```

which are complements and both have least period `32`. The parent bank does
not determine which physical phase occurs at depth `1420878969`.

Seven frozen orbit checkpoints include steps `1`, `10`, `1000`, `10^6`,
`10^8`, `10^9`, and the terminal step. Three whole-orbit checksums are

```text
xor =25469bdaa1c33b27,
sum =65b743c483016d28,
FNV =46485ada15e2f5b9.                               (19)
```

The primary solver walks from the first reset bit of `B`. The independent
solver instead composes the affine monodromy from time residue zero and makes
a second pass from its forced cyclic initial bit. The parent audit separately
reconstructs the THM-4047 bank from its already audited seven physical query
bits and freezes

```text
ell_87867=cf3030cf,   ell_99998=f687520b,
ell_99999=be79924b,   ell_100000=90f58380.            (20)
```

## 5. Scope firewall

This theorem applies the native spatial operation to entire fixed-column
temporal words. It proves that the zero clock never terminates and computes its
next event. It does **not** transfer fixed-depth periodicity uniformly to the
moving center

```text
c_t=ell_t(t).                                        (21)
```

Therefore Rule-30 non-eventual-periodicity, asymptotic balance, and fixed-seed
query-complexity lower bounds all remain **OPEN**. The physical phase in `(18)`
is open as well. The enormous exact gap `(14)` is a hostile to extrapolating a
low-order scalar recurrence from the first seven zero depths.

## 6. Reproduction

```bash
cc -std=c11 -O2 -Wall -Wextra -Wpedantic \
  04-computation/rule30_temporal_cylinder_transfer_thm4065.c \
  -o /tmp/r30_thm4065_primary
/tmp/r30_thm4065_primary

cc -std=c11 -O3 -Wall -Wextra -Wpedantic \
  04-computation/rule30_temporal_cylinder_transfer_thm4065_independent_audit.c \
  -o /tmp/r30_thm4065_independent
/tmp/r30_thm4065_independent

python3 04-computation/rule30_temporal_cylinder_transfer_thm4065_parent_audit.py
python3 -O 04-computation/rule30_temporal_cylinder_transfer_thm4065_parent_audit.py
```

The C runs are exact walks rather than cycle samples. Both compiler
optimization modes reproduce the frozen semantic transcripts.
