# Exact two-port-set coordinates and the finite detour controls

October 20, 2026. **PROVED**: the source parameterization below and the
uniform open end-pair ledger in [probe_new_ledgers.md](probe_new_ledgers.md).
**FINITE-EXACT**: all strong mixed endpoint presentations at orders 6--9
have literal three-black detours. There is no order-ten census here.
The all-order endpoint-three theorem is the separate proof partition in
the main packet; these finite controls are not its source of quantifiers.

## 1. Inheritance, marks, and completeness

The closest mechanism is the [endpoint-two theorem](../endpoint_two_repairs_20261020g/README.md).
Its hostile is a reversed block with no path having the prescribed
endpoints. Its corrected near miss is exporting a singleton reversal
identity to a longer block. We retain both actual reflected port sets,
their overlap, and each intermediate word's actual anti-involution.

The board consists of the port sets; the initial-arc guards; path-edge
XOR; anti-reversal witnesses; and the distinction between failed local
words and failed reachability. The anchor is the next endpoint layer,
the niche is the reusable affine compiler, and the wildcard is extracting
a uniform block word from the smallest residual certificate.

Fix P=(0,...,n-1), with a=0 and z=n-1. Every closed blue source with
endpoint scores (n-4,3) is obtained uniquely as follows:

1. Choose an arbitrary standard-blue middle tournament on positions
   1,...,n-2.
2. Choose a two-element incoming set I contained in {2,...,n-2}.
3. Put J={n-1-u:u in I}. The endpoint incidences are
   N^-(a)={z} union I and N^+(z)={a} union J.

The path condition excludes 1 from I. Blue reversal of P forces the
reflection relation; it imposes no disjointness condition on I,J.
Conversely these incidences and a blue middle give the required closed
blue source, so this is a bijection of marked labelled presentations.
In particular overlaps 0,1,2 are retained. Overlap1 can occur only in
odd order and requires the central vertex in I; overlap2 means I=J.

There are `2^floor((n-3)^2/4)` blue middle presentations: the nonpath
pairs are freely oriented by their reversal orbits. Hence the full
marked universe has size

```
binom(n-3,2) * 2^floor((n-3)^2/4).
```

The exact universes at n=6,7,8,9 therefore have sizes 12,96,640,7680.
This is not an unlabelled tournament count.

## 2. What the inherited bank misses

[probe.py](probe.py) reuses the frozen endpoint-two bit operations and
compiler but retests every guard and every actual anti equation. It
does not inherit the endpoint-two degree-based blackness conclusions.
The old bank includes the cyclic ledgers, all A/B consecutive-block
cuts, the two-block last-port word, the old C word, and their reflected
copies. The compact receipt records illegal joins by step and by type
(endpoint incidence, middle path, or middle chord), and distinguishes
legality from actual blackness.

| n | Full universe | Strong mixed endpoints | Old-bank covered | Final-bank covered | Individual certificates |
|---:|---:|---:|---:|---:|---:|
|6|12|11|5|8|3|
|7|96|96|77|93|3|
|8|640|640|589|638|2|
|9|7,680|7,680|7,357|7,677|3|

The final bank, in [probe_refined.py](probe_refined.py), includes the new
two-outgoing-port C word, movable reflected-port last-c cuts, the
nonmodule TT3-end repair, both end-pair repairs, and all earlier A cuts.
Every candidate is checked on its actual source even when the source
is below the word's all-order proof boundary. A failed template is
never recorded as a bridge.

The old order-nine residuals were mostly disjoint-port closed middles:
304 such sources, versus only three with overlapping ports. The new
block words removed those structural failures, leaving only three
small-order exceptions. This is why retaining the marks was more useful
than another aggregate degree statistic.

## 3. Exact finite completion and the necessary filter

[probe_complete.py](probe_complete.py) finishes the eleven remaining
presentations by labelled three-K search. It enumerates black first and
last Hamiltonian words; their edge XOR with E(P) determines the only
possible middle undirected path, whose orientations are tested. The
search stops at the first certificate and makes no total-count claim.
Each exception was found within 82 first/last candidate pairs.

Every one of the resulting 8,427 retained presentations has its full
selected ledger independently replayed with sets of directed arcs,
rather than the producer's bitmask K implementation. This checks
directedness, an actual failed anti equation at each step, and literal
target equality. All eleven exceptional certificates, the source and
target rows, the marks, and the independent anti witnesses are in
[probe_complete.json](probe_complete.json).

The filter is genuinely necessary at order six. The unique excluded
source is

```
T=(6,4,8,51,39,7),       I={3,4},       J={1,2},
K_P(T)=(58,60,57,16,32,8)=C3 -> C3.
```

Its partner has nine Hamiltonian paths and none is black. Thus no black
detour can end there. This is the established pure-blue boundary, not
a new strong mixed obstruction. It explicitly refutes extending the
present small-order control to *all* sources at n=6 without a filter.
At n=7,8,9 no source is removed by the strong mixed endpoint filter.

## 4. The sharp cyclic-cut hostile is repaired

The root's all-defect cyclic-cut bound is sharp. For endpoint three,
the order-ten source with two cyclic five-vertex parity orbits is

```
(190,44,760,176,994,704,907,770,559,11),    I={6,8}, J={1,3}.
```

Both P and its one-step cyclic cut are blue. This is a hostile to one
blackness argument, not a bridge. A literal three-black repair is

```
(0,1,2,3,4,6,7,8,5,9),
(0,1,2,3,8,5,9,4,6,7),
(0,1,2,3,8,9,4,5,6,7).
```

Its excesses are (0,1,-1); actual anti failure certifies the balanced
first step. The finite checker independently replays this one control.
It does not classify other order-ten presentations.

## 5. Affine certificates and proof extraction

[probe_symbolic.py](probe_symbolic.py) imposes the source's anti equations,
the exact two port sets, all required ledger joins, and a selected
step's proposed anti equations as an affine system over GF(2).
Inconsistency proves that step black within the stated finite guards;
consistency returns a full labelled source and replays it. The receipt
does not turn a bounded list of infeasible systems into an all-order
claim.

The compiler exposed the actual n9 second-step blue source in the
penultimate branch and ruled out several apparent large-order exceptions.
The useful return was the short arc identities extracted from those
contradictions: the root's separate n11/n14 proof and the new end-pair
constant witnesses now prove the relevant all-order statements directly.
The map is therefore from a labelled feasibility system to an explicit
arc contradiction, preserving the exact source guards; forgetting them
would invalidate the conclusion.

Reproduction:

```
python 05-knowledge/results/endpoint_three_repairs_20261020h/probe.py
python 05-knowledge/results/endpoint_three_repairs_20261020h/probe_symbolic.py
python 05-knowledge/results/endpoint_three_repairs_20261020h/probe_refined.py
python 05-knowledge/results/endpoint_three_repairs_20261020h/probe_complete.py
python -O 05-knowledge/results/endpoint_three_repairs_20261020h/probe_complete.py
```

Checks use explicit exceptions rather than assertions. The retained
universe and source hashes are recorded in the compact receipts. The
proof does not claim a Kuratowski/Petersen correspondence, a Pythagorean
invariant, or connectivity of the full black metagraph from these local
detours alone.
