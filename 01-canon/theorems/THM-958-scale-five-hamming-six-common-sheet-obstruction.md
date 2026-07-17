---
id: THM-958
title: Scale-five Hamming-six owner-triple obstruction
status: PROVED STRUCTURAL + FINITE-EXACT — the primitive proper AP-centred c=5 Hamming-six common-sheet bank is empty; an exhaustive 14,414,400-context mask scan and an independent capacity/owner-obligation reduction agree, so no metric recursion is needed
source: codex-2026-07-17-S61 scale-five structural audit
depends_on: [THM-765, THM-810, THM-823, THM-859, THM-860]
related: [THM-861, THM-862, THM-957, HYP-6820]
verification:
  - 04-computation/lrc13_scale_five_hamming_six_common_sheet_codex_S61.py
  - 05-knowledge/results/lrc13_scale_five_hamming_six_common_sheet_codex_S61.out
---

# THM-958 — scale five is killed by an owner-triple obstruction

Let `R subset F_13^*` have six elements, put `P=[12] minus R`, and consider
the primitive proper AP-centred packet

```text
A=5P union {w_r:r in R},
w_r=5r (mod 13),                 w_r>0,        w_r!=5r. (1)
```

Assume for contradiction that

```text
M(A)<=1/13.                                             (2)
```

For each replacement put

```text
D_r=5/gcd(5,w_r) in {1,5},
u_r=D_r w_r/5,
e_r=u_r mod D_r in (Z/D_r Z)^*.                        (3)
```

THM-860 gives both the hereditary leave-one-out law

```text
lcm(D_s:s!=r)=5                                      (4)
```

and common-sheet coverage at every replacement owner.  Equation (4) says
that at least two colours have order five.  This theorem proves that no such
order/unit decoration covers all five sheets.

## Theorem

### A. The complete common-sheet bank is empty

The exhaustive scale-five census is

```text
D=5 colours       tested order/unit contexts       survivors
2                         221,760                       0
3                       1,182,720                       0
4                       3,548,160                       0
5                       5,677,056                       0
6                       3,784,704                       0
total                  14,414,400                       0. (5)
```

Thus there are zero common-sheet presentations, zero unit contexts, zero
metric root contexts, and zero first metric edges.  In particular (2) is
impossible and every primitive proper packet in (1) satisfies

```text
M(A)>1/13.                                             (6)
```

### B. Sheet capacity reduces the scan to 64 six-colour supports

At a fixed owner, the literal masks have this unit-independent support:

```text
D=1:  the owner supplies all five sheets; every off-owner supplies none;
D=5:  every nonzero mask is a singleton;
      provider/owner in Z={4,9,12} supplies no sheet;
      every other ratio supplies one sheet.            (7)
```

Moreover, at a `D=5` owner its self mask is independent of its unit.  Every
nonzero off-owner provider maps its four units bijectively to the four other
sheets.

If there are `k<=4` order-five colours, an order-five owner sees at most `k`
singleton masks and cannot cover five sheets.  If `k=5`, every ordered pair
of order-five labels must be nonzero in (7).  Hence every unordered label
ratio must lie in

```text
N={a:a notin Z and a^(-1) notin Z}={2,5,6,7,8,11}.     (8)
```

The element `2` generates `F_13^*`, and the elements of `N` are exactly its
six odd powers.  Therefore the undirected graph defined by (8) is `K_6,6`,
with clique number two.  It has no five-clique, so `k=5` is impossible.

It remains only to consider `k=6`.  Scalar sheet capacity requires every
owner to have at most one incoming zero-provider arc

```text
p ->_0 o       iff       p/o in {4,9,12}.              (9)
```

Exactly 64 six-label sets obey this condition.  Every one has exactly one
zero predecessor and one zero successor at each vertex; ratio `12` never
occurs.  These 64 supports are exactly the 64 signed-doubling supports of
THM-860(D):

```text
p ->_2 o       iff       o/p in {2,-2}.                (10)
```

Their multiplicative orbit sizes are

```text
12,12,12,4,12,12.                                     (11)
```

### C. The faithful unit object is a six-constraint all-different CSP

On each of the 64 supports, relation (9) is two directed triangles, while
relation (10) is one directed Hamiltonian six-cycle:

```text
relation       edges       SCC sizes       Hamiltonian paths
->_0             6           (3,3)                 0
->_2             6             (6)                 6. (12)
```

At an owner `o`, omit `o` itself and its unique zero predecessor.  The other
four providers each choose one of the four nonself sheets through a
provider/owner-dependent permutation of their unit in `F_5^*`.  The owner is
covered exactly when those four sheet symbols are all different.  This gives
six four-ary all-different obligations on the six unit variables.

Write the doubling cycle in order as

```text
C=(v_0,v_1,v_2,v_3,v_4,v_5).                          (13)
```

For every support, the exact unit-word census is

```text
owner obligations satisfied       0       1      2     3     4    5    6
unit words                      2,332   1,248    504     0    12    0    0.
                                                               (14)
```

Thus no unit word covers more than four owners.  Every pair of owner
obligations is compatible, but exactly eight owner triples are incompatible:

```text
{v_i,v_(i+1),v_(i+2)}        for i in Z/6Z,
{v_0,v_2,v_4},               {v_1,v_3,v_5}.           (15)
```

The last two triples are precisely the two zero-provider SCCs in (12).  The
six others are the consecutive length-three windows of the doubling cycle.
Hence (15) is the exact minimal unsatisfiable owner-triple hypergraph.

The twelve maximally compatible unit words in (14) have an equally rigid
description: their failed owners are one of the three opposite pairs

```text
{v_0,v_3}, {v_1,v_4}, {v_2,v_5},                      (16)
```

with four unit words for each pair.

### D. Tournament completion sees the support but loses the obstruction

For Tournament Analysis, take the directed zero-provider arc (9) as the
pairwise observable.  Its nine missing unordered pairs form `K_3,3`.  Choose
the lexicographically first Hamiltonian path of that tie graph and orient all
ties by its total order.  Across all 64 supports, the completed tournaments
have

```text
score multiset             (1,2,2,3,3,4),
directed triangles         6,
SCC sizes                  (6),
sparse-edge flips          {2:16,3:32,4:16},
Hamiltonian-path counts    {29:16,31:32,33:16},
joint fingerprints         3.                           (17)
```

This completion is only an audit.  A bare tournament forgets which arcs are
literal zero-provider arcs and, more decisively, it has no place for the
three-owner hyperedges in (15).  The faithful finite carrier is either

```text
30 owner-sheet obligations + provider/unit incidence hyperedges,          (18a)
```

or, after the capacity reduction,

```text
six unit variables + six all-different owners + forbidden-triple hypergraph. (18b)
```

## Proof

### 1. Literal finite exhaustion

The possible effective-order/unit states are

```text
(D,e)=(1,0),(5,1),(5,2),(5,3),(5,4).                  (19)
```

Among the `5^6` state words, equation (4) removes the all-order-one word and
the `6*4` words with only one order-five colour.  Thus every label set has

```text
5^6-1-6*4=15,600                                      (20)
```

admissible state words.  There are `binom(12,6)=924` label sets, giving the
`14,414,400` rows in (5).

For state `(D,e)`, provider `r`, and owner `o`, let `u` be the CRT solution

```text
u=Dr (mod 13),                       u=e (mod D).
```

Its exact sheet mask is

```text
E_(r,D,e)(o)={ell in Z/5Z:
 -D < [u(o^(-1)+13ell)]_(13D) <= D}.                  (21)
```

The verifier packs the five bits at each of six owners into a 30-bit vector
and tests whether the union of the six providers is `2^30-1`.  It tests every
row in (5) and finds none.  This proves Part A directly.

### 2. Independent capacity and unit obstruction

Direct reduction of (21) gives (7).  The singleton count proves the `k<=4`
case.  Taking both directions of every order-five pair gives (8), and the odd
exponent description proves that graph is `K_6,6`; this proves the `k=5`
case without unit enumeration.

For `k=6`, exhaust the `924` label sets using only the zero arcs (9).  The 64
capacity survivors have the sparse fingerprints (11)--(12) and agree exactly
with the signed-cycle bank (10).  For each survivor, exhaust its

```text
4^6=4,096                                               (22)
```

unit words using the literal masks (21).  The result is (14).  Testing all
`binom(6,2)` owner pairs and `binom(6,3)` owner triples gives (15), while the
twelve four-owner maxima give (16).  This independently proves the empty-bank
verdict and exposes its minimal obstruction.

Finally, complete (9) by the declared `K_3,3` tie path and enumerate all six
vertex paths and triples.  This gives (17) and proves the tournament audit.
No metric height has entered either proof. ∎

## Assumption challenge and scope guardrail

The decisive vertices are not runners:

- runner vertices expose the recurring 64-support signed-cycle scaffold but
  do not encode units or simultaneous sheet coverage;
- individual sheet residues give the exact 30-vertex incidence problem but
  conceal its small obstruction;
- owner obligations give the sharp quotient: the eight three-owner
  hyperedges in (15) already make the system inconsistent;
- a completed pair tournament erases those hyperedges and is therefore too
  lossy even though its fingerprints are rigid.

THM-958 closes only the primitive proper AP-centred common-scale-five
Hamming-six face.  It does not close `c>=6`, the ramified Hamming-five bank,
non-AP-centred/deep-sheet branches, or global `n=12` sporadic emptiness.  Since
the common-sheet bank is empty, a scale-five metric component recursion would
have zero roots and is unnecessary.

## Verification

```bash
python3 04-computation/lrc13_scale_five_hamming_six_common_sheet_codex_S61.py |
  cmp - 05-knowledge/results/lrc13_scale_five_hamming_six_common_sheet_codex_S61.out
python3 -O 04-computation/lrc13_scale_five_hamming_six_common_sheet_codex_S61.py |
  cmp - 05-knowledge/results/lrc13_scale_five_hamming_six_common_sheet_codex_S61.out
```

Frozen integrity data:

```text
source SHA-256
  f01b08e43bf192ba9e5dc844eeda38a867a905e9b5d8b10feaca5821374eeaf7
output SHA-256
  906d35b1d07192b4f039716230192f8559a4677af3552469e57341a2e8afcb37
local-mask payload FNV-64
  faf78482761e2de4
support/obligation payload SHA-256
  bc5aa59e155d645552b8614786c5b9bcbb91e4f83998308855a0963ba1b90538. (23)
```
