---
id: THM-3033
title: "Projected k3 z246 to z244 descent and z243 high-floor addendum"
status: >
  PROVED-CANDIDATE + VERIFIED-EXACT; AWAITING FINAL INDEPENDENT PROMOTION.
  In the lossless projected k=3 atlas, all 197 body rows at z1=246,245,244
  are empty, as are the three first-at-or-above-high-floor rows at z1=243.
  The projected cap therefore descends to z1<=243 and the necessary-row
  ledger becomes 375,051.  The other 151 z1=243 wall rows remain open.  This
  is not LRC(14).
source: cap246-recover-2026-08-01
audit: >
  A read-level hostile audit checked inheritance, quantifiers, the two
  high-floor branches, translated-band/max-gap direction, row universe, cap,
  and ledger.  Independent reduced replays reproduced all z1=245,244 and the
  three z1=243 addendum rows, including every rational Farkas check and the
  six-case z1=245 terminal.  Two full development executions independently
  reproduced the frozen 29,236 Farkas checks, 1,758 terminal cases, transcript,
  and semantic digest.  Post-freeze ordinary, optimized, and stored
  transcripts are byte-identical and the documentation gate passes.  Final
  independent promotion review remains required by the status above.
depends_on:
  - THM-2981-projected-k3-z270-to-z247-cardinality-torsion-descent
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
related:
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
  - MISTAKE-334
script: 04-computation/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.py
output: 05-knowledge/results/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.out
script_sha256: 128bea25d4720a15dabe824f53e5702f37fadf9d6657c8b1a883f76bbc163c48
output_sha256: d012f84d8e380019b4da8efbe9afd6695079c585f47e4aa2c474828ca62014e4
semantic_sha256: 3f9945d333b13adef8aa8e7b162960bfde97f1e5a77945edcc735013d8ce4231
hash_basis: LF-normalized bytes
---

# THM-3033 -- projected k3 z246 to z244 descent and z243 high-floor addendum

**PROVED-CANDIDATE + VERIFIED-EXACT; AWAITING FINAL INDEPENDENT PROMOTION.**

## Statement

In the lossless projected `k=3` scalar atlas inherited from THM-2981, every
body row at `z_1=246`, `245`, or `244` is empty.  The three `z_1=243` rows
whose first label is already at or above the integer high floor are also
empty.  Consequently the candidate conclusion is

```text
projected k=3 cap: z_1<=243,                              (1)
```

while exactly `151` first-below-floor rows at `z_1=243` remain open.

The closed bank contains `197+3=200` body rows.  Relative to THM-2981 its
candidate necessary-row ledger is therefore

```text
375,251-197-3=375,051.                                   (2)
```

This is a projected necessary-sector result.  It is not LRC(14), closes no
remaining `z_1=243` wall row, and has no consequence for `k<=1` or the final
rung.

## 1. Exact inherited universe

Let `L` be a body period and put

```text
H=floor(13L/132)+1.                                      (3)
```

The hash-pinned THM-2981 engine derives the next three occupied levels
directly from the THM-2941 atlas.  Their body counts and exact attained-state
partition are

```text
z_1   rows   states   crude   common status   residual
246    194   49,578   18,946       29,116        1,516
245      1      119        0          113            6
244      2        7        0            7            0
--------------------------------------------------------
total  197   49,704   18,946       29,236        1,522.  (4)
```

Every one of the `29,236` common-status exclusions is rebuilt from its exact
marginals, allowed capacities, and residue-load histogram.  The companion
then checks the returned Farkas identity and strict direction over the exact
rationals.  A solver-selected dual basis is verification data, not semantic
data.

The integer high-floor split is

```text
z_1=246: 156 below H, 38 at or above H;
z_1=245:   1 below H,  0 at or above H;
z_1=244:   0 below H,  2 at or above H.                  (5)
```

The forty first-at-or-above-floor rows have exact profile

```text
1,490 states = 782 crude + 708 common status + 0 residual. (6)
```

Thus they close before any terminal translated-band argument.  This is the
order branch: strict label order makes every subsequent label high.  It is
logically distinct from the wall branch below.

## 2. The first-below-floor terminal

For a residual first-below-floor row, THM-2941's projected-safe-arc wall
forces at least one subsequent label to be at least `H`.  The exact
duplicate-permitting upper bound for packets with two or more such labels
misses the scalar floor on every residual body.  Its smallest exact gap is

```text
86147/532003290>0.                                       (7)
```

Hence every genuine residual packet has exactly one high label.  Literal
finite enumeration of both signs below `H`, together with the height-free
unit-ray maximum above `H`, yields

```text
27 residual bodies;
1,408 zero-high scalar hostile passes;
1,758 one-high cases.                                    (8)
```

The zero-high rows in `(8)` are stress controls: the inherited wall, not a
new scalar inequality, excludes them.

Fix one of the one-high cases.  Let `d` be the exact denominator of its high
label, and let `S` be the distinct residues modulo `d` of body cells that are
safe on their whole closed cell for the body, the fixed first label, and the
two low labels.  Put

```text
kappa=ceil(d/7).                                         (9)
```

Writing the high label as

```text
z=(L/d)u+hL,             gcd(u,d)=1,                    (10)
```

THM-2984 makes the unit phase independent of `h`.  At each local cell
coordinate the strict high-danger set is an arbitrary translated open band
of length `d/7`, so it contains at most `kappa` residue classes.  In every
one of the `1,758` cases the exact companion proves

```text
|S|>kappa.                                               (11)
```

Thus the translated band cannot contain all of `uS`, for any primitive unit
`u` and any local translate.  The complete-cell projection corollary of
THM-2984 supplies a high-safe selected cell at every local coordinate, so the
projected section is the full circle.  This contradicts an aligned
completion and closes every terminal case.

The actual closure mechanism in this census is exactly the cardinality gate
`(11)`.  No terminal needs the stronger shape-sensitive fallback.  For audit
and future lower frontiers, the companion nevertheless checks the exact
THM-2984 maximum-gap criterion: if `S` does not pass `(11)`, then for a fixed
unit `u` it lies in a cyclic block of `kappa` consecutive residues if and
only if the cyclic gap word of `uS` has maximum at least

```text
d-kappa+1.                                               (12)
```

All subsets and all units for `d=1,...,12` are compared with literal affine
block containment (`8,190` exact checks), including empty, singleton, and
wraparound boundaries.  This verifies the fallback direction but does not
misstate it as a mechanism used in `(11)`.

## 3. The `z_1=243` high-floor addendum

The atlas has `154` rows at the next occupied height.  Exactly three already
lie at or above their integer high floor.  Reusing only the exact inherited
state/status screen gives

```text
E=(1,2,3,8,10,12):  1=1 crude+0 status+0 residual;
E=(1,3,4,8,10,12): 13=13 crude+0 status+0 residual;
E=(2,4,6,8,12,14):  1=0 crude+1 status+0 residual.       (13)
```

The last status exclusion is independently checked over the exact rationals.
The other `151` rows have first label below `H` and are deliberately left
open.  Equations `(4)`--`(13)` therefore prove precisely the candidate cap
and ledger in `(1)`--`(2)`, with no silent claim at the remaining boundary.

## 4. Exact verification and promotion boundary

The compositor pins the LF-normalized THM-2981 source, output, and semantic
digest, rejects bare carriage returns, reconstructs its row universe from the
atlas, and freezes every row, status instance, terminal certificate, boundary
control, cap, ledger, and semantic digest.  The final semantic digest is

```text
3f9945d333b13adef8aa8e7b162960bfde97f1e5a77945edcc735013d8ce4231.
```

The canonical commands are

```text
python 04-computation/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.py
python -O 04-computation/lrc14_j7_k3_z246_to_z244_translated_band_descent_thm2941.py
```

The frozen ordinary and optimized transcripts are byte-identical to the stored
output.  The theorem remains a proof candidate until a final independent
review promotes its status.  No cap, ledger, or LRC consequence may be cited
from this file before that promotion.
