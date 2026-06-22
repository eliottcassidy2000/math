---
id: HYP-2883
status: SYNTHESIS / proof-target; supported by HYP-2872, HYP-2881, HYP-2882, HYP-2617, HYP-2632, HYP-2828
source: codex-2026-06-22-S101
tags: [tournaments, even-graphs, lrc14, cycle-space, parity-dual, coimage, character-kernel, relation-depth, tournament-analysis]
related:
  - HYP-2872
  - HYP-2877
  - HYP-2878
  - HYP-2879
  - HYP-2880
  - HYP-2881
  - HYP-2882
  - HYP-2617
  - HYP-2632
  - HYP-2694
  - HYP-2828
  - OPEN-Q-108
---

# HYP-2883: even graphs are parity-dual cycle-space addresses, not obstruction closures

The current even-graph lesson is sharper than the older equinumerosity slogan.
Even graphs do not preserve the tournament Hamiltonian-path obstruction.  They
preserve a parity-dual cycle-space address that can locate packet incidence
before orientation, order, and scalar weights are restored.

That is exactly the kind of address the LRC(14) residual needs, but only if it
is used at the support-six packet layer rather than as a direct `H` analogue.

## The scar calculation

For any graph `G`,

```text
I(G,2) = number of independent sets of G with weight 2.
```

For a clique,

```text
I(K_r,2) = 1 + 2r.
```

The two selection rules on this clique spine are opposite:

```text
tournament H-spectrum:
  K_r conflict cliques are realizable except the forbidden K_3 packet
  -> clique values appear except 7.

even-graph I(G,2)-spectrum:
  K_r is degree-even iff r is odd
  -> odd clique values 7,11,15,... appear
  -> even clique values 5,9,13,... disappear.
```

So the even-graph projection heals the tournament's `7` scar.  It does not
explain the scalar gap by preserving it.  It explains the layer where the scar
was born: cycle-space support before orientation-realizability.

This strengthens HYP-2872.  Degree-2 smoothing and GF(2) contraction are too
coarse for `H`, and now the clique spine shows why: parity membership can be
perfectly correct while scalar obstruction is inverted.

Mac-mini S38 supplies the compatible positive statement.  With a fixed base
Hamiltonian path, the tournament edge space splits as

```text
edge space of K_n = cut half + cycle half.
```

The base path fixes the cut half, and the non-base tiles give the cycle half,
which is exactly the labelled even-graph address.  Therefore `H` can be read
from the labelled cycle half once the cut is fixed.  But the cycle-only
isomorphism quotient `E_n` forgets the cut and scrambles `H`: S38 reports that
at `n=5`, `64` fixed-path tilings collapse to `7` even-graph iso classes, and
`5` of those `7` classes carry more than one `H` value.  This is the exact
discipline HYP-2883 needs:

```text
cycle address useful; cycle-only scalar quotient unsafe.
```

## LRC translation

The LRC residual has the same danger:

```text
raw speeds
  -> exact-period / phi packets
  -> support-six relation lattice
  -> mod-7 coimage
  -> chi_7 plus affine signed kernel
  -> scalar witness margin.
```

Projecting too early can heal the obstruction, just as even graphs heal
tournament `7`.  But the projection is still useful as an address if it keeps
the correct incidence:

```text
tournament:
  primitive odd-cycle packets -> cycle-space support -> conflict graph

LRC(14):
  exact-period packets -> support-six relation support -> signed coimage graph
```

The LRC analogue of a tournament even graph is not a graph of runners.  It is a
cycle-space packet graph of relation supports after exact-period expansion.

## Packet graph target

Define a support-six signed packet graph `P_7(E)` as follows.

Candidate vertices:

```text
v = (A, R, m)
```

where:

- `A` is a projective speed-residue coimage class from HYP-2617;
- `R` is an exact-period Fourier packet with `A dot R = 0 mod 7`;
- `m` is the exact-period or height-wall label retained before scalarizing.

Edges should be tested in several gauges:

```text
compatibility edge:
  two packets can live in the same low-height wall ledger

conflict edge:
  their signed kernels force opposite affine lanes or Legendre signs

energy edge:
  their supports share enough additive relations to prevent L2 decorrelation
```

The switch/gauge is the HYP-2632 signed kernel:

```text
chi_7 phase plus affine zero lane a+b=2 mod 7.
```

The Hamiltonian path tie for Tournament Analysis is not a runner ordering.  It
is the proof-order spine:

```text
support-six signed packet graph
  > relation-depth separator
  > even-graph cycle support
  > E_7 metagraph odd holes
  > strong-atom ear profile
  > raw fixed-path row
  > raw runner set.
```

Fingerprints to compute:

- score histograms under each gauge;
- SCCs of mutually reinforcing signed packets;
- directed 3-cycles as Schur/additive-triple analogues;
- chordless odd holes compared with the E7 `C_5/C_7` layer;
- edge flips under QR/NQR conjugation and under `a+b=2`.

## Candidate theorem shape

After the finite low-height wall ledger is deleted or certified, every
support-six LRC(14) residual packet profile falls into one of three branches.

```text
1. parity-null branch:
   the packet graph is coimage-null or signed-small under the chi_7 + affine
   kernel, so the signed reciprocal tail is below the margin.

2. bounded-relation branch:
   relation-depth <= 2, hence the profile routes to the finite wall/coimage
   atlas from HYP-2617, HYP-2624, and HYP-2828.

3. spread-cycle branch:
   relation-depth >= 3, hence the packet graph has no coherent low-rank
   additive packet.  Freiman/Plunnecke/Ruzsa structure or L2 Parseval then
   gives the decorrelation tail needed by HYP-2694.
```

This is not yet the LRC(14) proof.  It is a sharper closing stack: the
uncontrolled residual is no longer "all wide clusters" but "packet graphs with
relation-depth >= 3 and no parity-null certificate."  That class is precisely
where high additive energy should force a small Freiman model, and low
additive energy should give decorrelation.

## First finite packet graph

Script:

- `04-computation/lrc14_repeated_packet_graph_codex_s101.py`
- output: `05-knowledge/results/lrc14_repeated_packet_graph_codex_s101.out`

The HYP-2632 repeated-residue kernel already gives a concrete signed graph.
Use residue vertices

```text
V = {0,2,3,4,5,6}.
```

Put a negative loop at `a` for the `4+2` packet

```text
(1,1,1,1,a,a),
```

and an edge between `a,b` for the `4+1+1` packet

```text
(1,1,1,1,a,b).
```

Measured in the HYP-2632 unit `U`, the loop weights are

```text
0:-4, 2:-25, 3:-18, 4:-25, 5:-18, 6:-18.
```

The edge weights are exactly `0`, `1`, or `8`, with zero edges precisely on the
affine lane

```text
a+b = 2 mod 7:
  (0,2), (3,6), (4,5).
```

Off that lane, the high/low split is the Legendre selector

```text
Q(a,b)=ab*(1+3(a+b))-1.
```

The new exact identity is local balance:

```text
loop(a) + sum_{b != a} edge(a,b) = 0
```

for every `a in {0,2,3,4,5,6}`.  The script reports no balance failures:

```text
v  loop  incident  balance
0    -4         4        0
2   -25        25        0
3   -18        18        0
4   -25        25        0
5   -18        18        0
6   -18        18        0
```

Equivalently,

```text
sum loops = -108,
sum edges = 54,
sum loops + 2*sum edges = 0.
```

This is stronger than the previous signed ledger.  The scalar sum still has
net `-54 U`, because an edge is counted once.  But the incidence form is
balanced exactly.  The analytic proof target should therefore lift a local
signed-current identity through reciprocal hyperplane sums, rather than
bounding the repeated packets by absolute mass.

## Why this helps the genuine-wide branch

HYP-2694 asks for arbitrary bounded cluster shape compression and a proof that
single coherent blocks are extremal.  HYP-2828 says the hard wide rows seem to
be shallow after two peels.  The even-graph parity-dual lens supplies the
missing discipline:

```text
Do not compare clusters by raw speed shape.
Compare their support-six packet incidence after exact-period expansion.
```

A multi-block cluster can look complicated in raw speeds while becoming a
disconnected or parity-null packet graph.  A single coherent block should be
the extremal case exactly when its packet graph keeps the most low-rank signed
cycles alive.

Thus the next analytic target is a packet compression lemma:

```text
For fixed small-part P and bounded wall ledger W, among all clusters with the
same packet graph cycle-support profile, the maximal p0/decorrelation error is
attained by a single interval-like block.
```

If this lemma is true, then the current finite atlases and spectrum-sum
certificates apply to the representative block rather than every raw cluster.

## Assumption challenge

Candidate vertex sets considered: runners, fixed-path free arcs, even-graph
edges, even-graph iso-classes, E7 metagraph vertices, odd-cycle conflict
packets, strong components, ear insertion profiles, exact-denominator packets,
support-six relation vectors, coimage classes, Fourier residue shells, affine
zero lanes, relation-depth buckets, and proof obligations.

Chosen quotient: support-six signed packet graph.  It preserves the LRC
predicate relevant to the residual: exact-period multiplicity, relation
support, mod-7 coimage, signed kernel, and additive-energy incidence.  It
destroys raw runner identity, raw speed order, and scalar `H` obstruction.

The challenged assumption is that even graphs should mirror the tournament
forbidden values.  They do not.  Their value for LRC is subtler: they warn us
which quotient layers invert scalar scars, and they point to the cycle-space
packet incidence that must be retained before signs and margins are evaluated.

## Immediate work items

1. Build `P_7(E)` for the HYP-2617 named coimage classes and label vertices by
   the HYP-2632 affine/Legendre kernel.
2. Lift the exact local balance of
   `lrc14_repeated_packet_graph_codex_s101.py` from the finite kernel to the
   reciprocal hyperplane sums after low-height wall deletion.
3. Compare odd-hole incidence in `P_7(E)` with the E7 metagraph `C_5/C_7`
   incidence profiles, not with scalar `H` values.
4. Test whether the HYP-2828 depth-3 branch has low additive energy after
   quotienting by exact-period packets.
5. Add a finite packet-compression audit: same packet cycle-support profile,
   raw cluster varied, measure whether the interval-like representative
   maximizes p0 or the decorrelation error.
