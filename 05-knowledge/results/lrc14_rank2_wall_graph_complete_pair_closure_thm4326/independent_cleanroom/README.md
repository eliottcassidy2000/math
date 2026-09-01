# Independent rank-two wall audit (2026-09-01)

Status: **FINITE-EXACT**.  Scope: the rank-at-most-two functional wall
certificate on the stated finite pair ledgers.  This packet makes **no physical
entry claim and no LRC(14) claim** by itself.

## Result

The clean-room event sweep gives:

| Pair universe | Rows | Coarse-positive | Coarse-nonpositive | Exact-positive among coarse failures |
|---|---:|---:|---:|---:|
| endpoint-588 typed universe | 22,647 | 22,540 | 107 | 107 / 107 |
| THM-4231 finite remainder | 181,194 | 181,087 | 107 | 107 / 107 |

The hostile set is identical in the two universes.  The broad screen's data-row
FNV-1a is `d7f639e550399b4c`; the 107-row exact ledger's data-row FNV-1a is
`f042132698af02b7`.  GCC `-O2` and `-O3 -DNDEBUG` runs are byte-identical.

The sole broad grid outside signed 64-bit range is asserted exactly:

`(q,r)=(713,719)`, `D=9,351,275,651,380,222,560`.

All grid positions, widths, graph weights, objectives, and ticks in the broad
implementation use signed 128-bit integers.

The smallest exact ticks among the 107 coarse failures are
`245,428,469,244` at `(50,70)`, body `031c7400`.  The full coarse ledger is
needed before calling this global: one coarse-positive row has a smaller raw
lower bound, and three have a smaller grid-normalized lower bound.  Independent
exact replay gives:

| Pair | Grid | Exact minimum ticks | Body |
|---|---:|---:|---|
| `(50,212)` | 4,833,907,245,367,200 | 4,207,251,055,549,752 | `053c6400` |
| `(50,274)` | 12,495,194,200,288,800 | 11,006,069,470,557,714 | `11187401` |
| `(100,110)` | 91,205,797,082,400 | 63,178,284,254,904 | `04f06408` |

Consequently, within this 181,194-row rank-two functional universe, `(50,70)`
is both the true raw-tick minimizer and the true normalized minimizer:

`245428469244 / 91205797082400 = 973922497 / 361927766200`.

This global statement is a minimum of the **rank-two functional surplus**, not
a physical-entry minimum.

## Exact mechanism

Fix a pair `(q,r)` and exact wall grid

`D = lcm(D0,14q,14r)`, where `D0=18,241,159,416,480`.

On every pair-safe open wall cell `C`, let `F(C)` be the subset of the fixed
30-speed pool that fails there.  Retain cells with `|F(C)|<=2`, and aggregate
their widths as a weighted graph:

- `w0`: rank-zero mass;
- `wi`: rank-one mass at vertex `i`;
- `wij`: rank-two mass on edge `{i,j}`;
- `W = w0 + sum wi + sum wij`;
- `ai = wi + sum_j wij`, the weighted degree including the singleton weight.

For a rank-nine body `B`, the retained survivor mass is exactly

`L2(B) = w0 + sum_(i notin B) wi + sum_({i,j} disjoint B) wij`

and therefore

`L2(B) = W - sum_(i in B) ai + sum_({i,j} subset B) wij`.

Because all edge weights are nonnegative,

`L2(B) >= W - sum(the nine largest ai)`.

The coarse certificate is `63*lower_bound - 4D > 0`.

For each of the 107 nonpositive coarse rows, the C++ optimizer exactly maximizes

`sum_(i in B) ai - sum_({i,j} subset B) wij`, `|B|=9`.

At a partial set `S`, an unprocessed vertex `v` has current marginal
`ai - sum_(u in S) wuv`.  For any completion `T`, its added reward is the sum
of these current marginals minus the nonnegative edge mass internal to `T`.
Thus the sum of the largest required number of current marginals is an
admissible upper bound.  Branches strictly below the incumbent are pruned;
equality is retained to preserve the least-mask tie rule.  This proves exact
optimality without enumerating all `C(30,9)` bodies.  A second direct evaluator
recomputes the winning mass from `w0,wi,wij`.

The independent Python replay then reconstructs a fresh sorted wall set,
classifies open cells by exact midpoint residues, and directly reproduces all
107 hostile winners plus the three coarse-positive challengers (110 distinct
pairs total).  It does not use the C++ event states or degree identity.

## Clean-room separation and controls

The graph builder in `src/rank2_event_sweep_wide_cleanroom.cpp` preaggregates
7,132 pool event coordinates once on `D0`, rescales them exactly, and merges
them with monotone `q` and `r` event streams.  This differs from the primary
sorted-wall/midpoint implementation.  The frozen fieldwise cross-check agrees
with that separately implemented primary screen on all 181,194 rows and with
its separate exact B&B on all 107 hostile minima/bodies.

Input hashes:

- `inputs/universe22647.csv`: `14f9be0d9472bc573e582ec6f4cb92c7def6f583f6afaf0b747f2a9713330317`
- `inputs/thm4231_remainder181194.csv`: `9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1`

Important output hashes are recorded in `results/optimization_parity.txt`; all
packet files are covered by `SHA256SUMS.txt`.

## Reproduction

From this packet directory on Windows PowerShell:

```powershell
powershell -ExecutionPolicy Bypass -File .\reproduce.ps1
python .\verify_packet.py
python -O .\verify_packet.py
```

The reproduction script compiles fresh `-O2` and `-O3` binaries in a temporary
directory, regenerates both universes, checks optimization parity against the
frozen ledgers, reruns the selected exact challengers, the raw-cell replay, the
global filters, and finally the hardened verifier.  The temporary build is
removed afterward.

## Packet map

- `results/thm4231_coarse_screen_O3.csv`: retained full 181,194-row broad ledger.
- `results/thm4231_hostile_exact_O{2,3}.csv`: exact 107-row ledgers.
- `results/coarse_screen_O3.csv`: retained 22,647-row narrow ledger.
- `results/hostile107_exact_O{2,3}.csv`: narrow exact ledgers.
- `results/selected_ratio_exact_O{2,3}.out`: exact raw/normalized challengers.
- `results/rawcell_winner_replay.csv`: independent 110-pair cell replay.
- `results/global_minimum_filter.out`: raw and normalized filter certificate.
- `results/primary_crosscheck.out`: frozen clean-room/primary fieldwise comparison.
- `verify_packet.py`: manifest, algebra, filters, parity, overflow, and raw-cell checks.
