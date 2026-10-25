---
id: THM-4734
title: "Endpoint-three detours at order ten with the order-six pure-blue exception; finite endpoint-four coverage at orders nine and ten"
status: >
  FINITE-EXACT (pre-audit): (A) order_ten_bank.out / .json (86,016 presentations, bank coverage by template, three
  residual rows) and order_ten_residuals.out / .json (the three labelled certificates: words, independent anti-reversal
  failures, endpoint excesses); (B) endpoint_four_bank.out, endpoint_four_n9.json, endpoint_four_n10.json (every
  presentation covered or search-certified with independent replay; pattern table); uniform_ledger_test.py / .out
  (best single template per search-needed pattern covers exactly half or less).
  METHOD: the packet's own engine is reused unchanged -- probe.bank, probe_refined.new_bank (word templates, every guard
  and every anti-reversal re-tested on the actual source), probe_search.find_one (labelled three-K search over black
  Hamiltonian paths of the source and the target, stopping at the first certificate), probe_complete.independent_replay
  (directed-arc-set replay: directedness, actual blackness witness per step, literal labelled target). No mixedness
  filter; the nonstrong-partner filter is recorded (zero cases at n = 9, 10 in both layers).
  CITED: endpoint_three_repairs_20261020h (README theorem for n >= 11, probe.md universe formula C(n-3,2) 2^floor((n-3)^2/4),
  R1 and theory section 2 general-k bounds, movable-cut theorem section 9), endpoint_two_repairs_20261020g (complete
  (n-3,2) theorem), kernels_games_detours_20261020f (K convention), THM-4729 Proposition K (K_P and R_P), THM-4731
  (black graph has exactly 1 + pure-blue components for n <= 8).
  NOT PROVED: any all-order statement for (n-5,4); black-core connectivity; that the three order-ten residual rows
  have no template detour (the search is first-success, not a template census).
source: mac-mini-2026-10-20-debtcones (owner's prompt: keep going on the remaining order-10 case and the (n-5,4) layer; movable-cut repairs toward black-core connectivity)
depends_on:
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
  - 05-knowledge/results/endpoint_two_repairs_20261020g/README.md
related:
  - 01-canon/theorems/THM-4729-closed-and-open-tilings-hamiltonian-cycle-bipartition-of-the-metagraph-lines-and-the-forced-arrangement.md
  - 01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md
  - 05-knowledge/hypotheses/HYP-9348-endpoint-layers-black-core-connectivity-and-the-next-zero-heavy-layers.md
scripts:
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_bank.py
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_residuals.py
  - 04-computation/experiments/endpoint_layers_20261020/endpoint_four_bank.py
  - 04-computation/experiments/endpoint_layers_20261020/uniform_ledger_test.py
outputs:
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_bank.out
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_bank.json
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_residuals.out
  - 04-computation/experiments/endpoint_layers_20261020/order_ten_residuals.json
  - 04-computation/experiments/endpoint_layers_20261020/endpoint_four_bank.out
  - 04-computation/experiments/endpoint_layers_20261020/endpoint_four_n9.json
  - 04-computation/experiments/endpoint_layers_20261020/endpoint_four_n10.json
  - 04-computation/experiments/endpoint_layers_20261020/uniform_ledger_test.out
---

# THM-4734 -- order ten closed, endpoint four opened

**FINITE-EXACT (pre-audit)** as itemized in the status.

**Scope correction, 2026-10-25:** MISTAKE-646 already corrected the
universal endpoint-three threshold to `n>=7`. The incoming
[ports/owners packet](../../05-knowledge/results/ports_owners_f1_20261025/README.md)
repaired the stale title and body. The independent
[complete low-order audit](../../05-knowledge/results/shape_ports_20261025/tournaments.md)
now records all 12 order-six presentations, exactly 11 detours, the two
order-five pure-blue targets, and the odd-order overlap-two witness.
The `C3 -> C3` target has no incident black move. Original finite receipts
and their pre-audit status are unchanged. See also the separately
[audited order-ten packet](../../05-knowledge/results/endpoint_three_order10_20261021/README.md)
and the [guarded endpoint-four continuation](../../05-knowledge/results/ports_owners_f1_20261025/endpoint.md).

## 0. Setting (endpoint_three_repairs_20261020h)

`P = (0, 1, ..., n-1)` is a closed blue Hamiltonian word of a tournament `T` (closed:
`n-1 -> 0`; blue: the positional reversal is an anti-automorphism). A K move along a directed
Hamiltonian word `Q` is `K_Q(T) = T^op xor E(Q)` (complement the non-path pairs); it is black
when `Q` is not blue. The endpoint-`k` layer has `d(0) = n-1-k`, `d(n-1) = k`; there is an
incoming port set `I` of `k-1` middle positions (`1 notin I`) with `N^-(0) = {n-1} u I`,
`N^+(n-1) = {0} u r_P(I)`; the marked universe has `C(n-3, k-1) 2^floor((n-3)^2/4)`
presentations (arbitrary blue middle). A **three-black detour** is a triple of successive
directed black words `Q1, Q2, Q3` with `E(Q1) xor E(Q2) xor E(Q3) = E(P)`, so that the
composite is the literal labelled target `K_P(T)`: the blue line `P` then joins two vertices
of the same black component and is not a bridge.

## 1. Order ten, layer (6,3)

The packet proved the layer for `n >= 11` and certified the **strong-mixed**
`n = 6..9` presentations by its probe; order ten was
"not claimed". Here all `21 x 4096 = 86,016` presentations were run through the inherited bank
(templates: release/cyclic ledgers, closed_no_chord, closed_chord, L2, C, A/B cuts, last_c_x,
C_two_ports, closed_last_c cuts, end-pair ledgers, mirrored TT3). Result:

```text
universe 86,016   bank-covered 86,013   residual 3   nonstrong partner 0
residual rows:  closed, I = {3,6} = J (overlap 2)      -> certificate after 10 pairs
                closed, I = {5,6}, J = {3,4}             -> certificate after 21 pairs
                closed, I = {7,8}, J = {1,2}             -> certificate after 2 pairs
```

The three certificates (order_ten_residuals.json), with independent replay (directedness,
one actual anti-reversal failure per step, literal target):

```text
I={3,6}: (0,1,2,3,4,5,8,6,7,9), (0,1,2,3,5,8,6,7,9,4), (0,1,2,3,5,6,7,8,9,4)   excesses 0,0,0
I={5,6}: (0,1,2,3,4,5,8,6,7,9), (0,1,2,3,9,8,7,4,6,5), (0,1,2,3,9,7,4,6,8,5)   excesses 0,-1,1
I={7,8}: (0,1,2,3,4,6,7,8,5,9), (0,1,2,3,8,5,9,4,6,7), (0,1,2,3,8,9,4,5,6,7)   excesses 0,-1,1
```

(The first residual is the `I = J` overlap-two case. This overlap is possible
in either parity; at order seven, `I=J={2,4}` is an explicit example. All
three residuals have endpoint-excess-zero words, so blackness is witnessed
by actual arcs, never by scores.)
Hence the inherited conclusion is **every closed blue `(n-4,3)` line for
`n>=7`**, and at `n=6` exactly the eleven presentations outside the
pure-blue-target exception. Orders 7--9 use the packet probe, order ten
uses this computation, and orders >=11 use the packet theorem. The complete
unfiltered order-six universe has 12 presentations: source
`(6,4,8,51,39,7)` maps to `(58,60,57,16,32,8)=C3 -> C3`, whose nine Hamilton
paths are all blue, making any black detour impossible. At order five,
both retained endpoint-three presentations likewise have pure-blue targets.
The displayed order-ten cyclic-cut hostile of R1
(blue one-step cut) is among the 86,013 bank-covered presentations.

## 2. Orders nine and ten, layer (5,4) / (6,4)

Same engine, incoming triples `I` from positions `2..n-2` (`C(6,3) = 20` at `n = 9`,
`C(7,3) = 35` at `n = 10`), blue middles `2^9`, `2^12`:

```text
n = 9 : universe 10,240   bank 9,915   search-certified 325    unresolved 0   nonstrong partner 0
n = 10: universe 143,360  bank 141,050 search-certified 2,310  unresolved 0   nonstrong partner 0
```

Every search certificate was replayed independently. The presentations that escape the
inherited bank (counts "searched" per pattern; each pattern has `2^floor((n-3)^2/4)`
presentations: 256 at `n = 9`, 2,048 at `n = 10`):

```text
pattern                       n = 9   n = 10
closed I = {2,3,4}            128     1024
closed I = {2,n-3,n-2}        128     1024
closed I = {n-4,n-3,n-2}       33      128
open   I = {n-4,n-3,n-2}       32      128
closed I = {2,3,5}/{2,3,6}      1        1
closed I = {4,6,7}/{3,4,6}      1        1
closed {2,5,6}, {4,5,6}, {4,7,8}; open {2,3,4}, {3,4,7}     1 each
```

So the inherited (k <= 3) ledgers already cover 96.8% (n = 9) and 98.4% (n = 10) of the
`k = 4` layer, and the misses sit on three consecutive ports at the start or end of the middle,
exactly where the k = 3 proofs used "the other port is not adjacent". For each such pattern the
best single word template among the found certificates covers exactly half of the pattern
(uniform_ledger_test.out: 128/256 and 1024/2048 for closed `{2,3,4}`; 64/256 and 512/2048 for
closed `{2,n-3,n-2}`; 64/256 and 256/2048 for the end-triples). The structural `(n-5,4)`
theorem therefore needs at least two complementary ledgers per consecutive-port pattern, split
by one middle arc, in the style of the packet's `epsilon` branches. The half-covering template
for closed `{2,3,4}` at `n = 10` is

```text
Q1 = (0,1,2,3,4,5,6,8,9,7),  Q2 = (1,2,3,4,5,6,8,0,9,7),  Q3 = (1,2,3,4,5,6,7,8,0,9),
```

a movable cut placed at the last port with the source relocated: the analogue of the
packet's section-9 words with `L = (1..6)`, `R = (8)`, `c = 7`.

## 3. What is already proved for k = 4 and what is not

From the packet, verbatim for general `k`: the one-step cyclic cut of a closed blue
endpoint-`k` word is black for `n > 4k-2` (sharp), and the common final word `(a,z,b,W,c)` is
black at `B = K_P(T)` for `n > 4k-2`; the movable-cut theorem needs only `z -> b`, `a -> y`
and one further middle out-neighbour of `z`, "with any number of additional exceptional
ports". For `k = 4` these start at `n = 15`. The orders `11..14` and the exhaustive case
partition (which rows of the `k = 3` table survive with a third port, and the new
consecutive-port rows) are open; the data above say which rows need new words. In the
trijection coordinates (graphs on the ordered set `[n-1]`, THM-4733 Proposition A) a K move
along `Q` is a change of coordinate system followed by complementation, so a detour is a
closed walk of three coordinate changes whose complement bits cancel to the base path; the
movable cut is the coordinate change that moves one cut vertex to the port, and it is the
natural candidate for a uniform `k`-independent repair of the closed non-adjacent-port rows.

## 4. Files

`04-computation/experiments/endpoint_layers_20261020/`: `order_ten_bank.py` (`.out`,
`.json`), `order_ten_residuals.py` (`.out`, `.json`), `endpoint_four_bank.py` (`.out`,
`endpoint_four_n9.json`, `endpoint_four_n10.json`), `uniform_ledger_test.py` (`.out`).
All import the packet's modules by path and change nothing in them.
