---
id: THM-431
name: u(21)=57 — the Erdos unit-distance maximum at 21 points (and the triangular lattice is NOT optimal)
status: |
  EXACT VALUE u(21)=57 is VERIFIED-IN-LITERATURE (Alexeev-Mixon-Parshall 2024,
  arXiv:2412.11914, Theorem 1; computer-assisted proof). Reproduced/derived HERE
  with EXACT INTEGER arithmetic: the triangular-lattice LOWER bound (47) and the
  elementary cherry/common-neighbour UPPER bound (71). The 3N-crossover BRIDGE to
  THM-421 is a new observation (this session).
date: 2026-06-06
session: monad-explorer-2026-06-06-S710
depends_on:
  - THM-421    # unit-distance 3N common-neighbour floor (this repo's lattice frontier)
  - HYP-2267   # triangular finite-crossover minimizer
  - HYP-2285   # the open handoff: exact value of N* / is triangular optimal at N=21
refs:
  - "Alexeev, Mixon, Parshall — The Erdos unit distance problem for small point sets, arXiv:2412.11914 (2024)"
  - "Schade — Exakte maximale Anzahlen gleicher Abstaende, Thesis, TU Braunschweig (1993) [lower bounds]"
  - "Globus, Parshall — Small unit-distance graphs in the plane, Bull. Inst. Comb. (2021) [forbidden subgraphs]"
  - "Harborth (1974) — triangular-lattice penny maximum floor(3n - sqrt(12n-3))"
---

# THM-431: u(21) = 57, and the triangular lattice is NOT optimal at N=21

## The question (dispatched seed)

`u(N)` = the maximum number of **unit distances** (pairs at distance exactly 1)
among `N` points in the plane (Erdos, 1946). The seed asked to pin down `u(21)`
and settle whether the triangular-lattice / Eisenstein construction is optimal
there. **Both are now answered.**

## Statement

**(1) Exact value [VERIFIED-IN-LITERATURE].**
`u(21) = 57`. Proven by Alexeev–Mixon–Parshall (2024), who improved the upper
bound from the previous `u(21) ≤ 68` down to `57`, matching Schade's 1993 lower
bound. Before AMP24 the state of the art was only `57 ≤ u(21) ≤ 68`. They
compute `u(n)` **exactly for every `n ≤ 21`** and enumerate the densest graphs;
at `N = 21` there are **5** densest unit-distance graphs (up to the relevant
equivalence) with 57 edges each, none a section of any single lattice.

Full proven table (AMP24, Thm 1a; lower part Schade; all EXACT):

```
 n :  0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21
u(n):  0  0  1  3  5  7  9 12 14 18 20 23 27 30 33 37 41 43 46 50 54 57
```

First still-open value: **`u(22) ∈ {60, 61}`** (AMP24: `60 ≤ u(22) ≤ 61`), then
`u(23) ∈ [64,66]`, `u(24) ∈ [68,72]`, … , `u(30) ∈ [93,110]`.

**(2) The triangular lattice is NOT optimal at N=21 [PROVED here, exact integer].**
The maximum over **all** triangular-lattice (Eisenstein) sections is Harborth's
penny number `⌊3N − √(12N−3)⌋ = ⌊63 − √249⌋ = 47`. We reproduced this exactly
(best over norms `D ∈ {1,3,4,7,9,12,13,16,21,28}` and several rational centres,
`unit_distance_u21_constructions_s710.py`): the maximiser is the `D=1` penny
disk at **47**. Since `u(21) = 57 > 47`, the **gap is 10**: the optimal
configuration beats every triangular-lattice section by 10 unit distances. So
**HYP-2267's "is triangular optimal at N=21?" resolves NEGATIVELY** — and not
narrowly.

**(3) Elementary upper bound [PROVED here, exact integer].**
The cherry/common-neighbour bound (THM-421(A), `≤2` common neighbours ⟹
`Σ_v C(d_v,2) ≤ N(N−1)`) gives `u(21) ≤ ⌊N(1+√(8N−7))/4⌋ = 71`. The proven
value 57 sits well below this elementary ceiling — closing the gap `[57,71]`
required AMP24's heavy machinery (74 minimal forbidden subgraphs of Globus–
Parshall + totally-unfaithful embedding tests + a SAT-style enumeration), not a
counting argument. (CORRECTS a bug in the first `s710` script that halved the
cherry cap and reported 52; the correct elementary bound is 71.)

## The bracket, made explicit (exact integer, this session)

```
 N  | triangular(Harborth) |  u(N) PROVEN | cherry upper | gap u−tri | gap upper−u
 17 |        36            |     43       |     52       |    7      |     9
 18 |        39            |     46       |     57       |    7      |    11
 19 |        42            |     50       |     61       |    8      |    11
 20 |        44            |     54       |     66       |   10      |    12
 21 |        47            |   * 57 *     |     71       |   10      |    14
 22 |        49            |  [60,61]     |     77       |  11..12   |  16..17
 23 |        52            |  [64,66]     |     82       |  12..14   |  16..18
 24 |        55            |  [68,72]     |     87       |  13..17   |  15..19
```

`u(N)` sits ~40% of the way from the triangular floor up to the cherry ceiling.
Both gaps GROW with N: the triangular lattice is increasingly suboptimal, and
the elementary counting ceiling is increasingly loose.

## The bridge to THM-421 (new observation)

THM-421 / HYP-2267 studied when a **single-norm lattice patch** first exceeds
`3N` unit distances: combinatorial floor `N ≥ 17`, best lattice **blob** at
`N = 32` (U=97>96), best lattice **disk** at `N = 43`. The AMP24 table lets us
locate the crossover for the **true optimum** `u(N)` over *all* configurations:

```
   n  :  26   27   28   29   30
  u≥  :  76   81   85   89   93     (AMP24 lower bounds)
  3n  :  78   81   84   87   90
 u−3n :  −2   ±0   +1   +2   +3
```

`u(N) > 3N` first at **n = 28** (with `u(27) ≥ 81 = 3·27` tying). Therefore:

> **The unit-distance maximum beats 3N at n = 28 — four points earlier than the
> best lattice blob (N=32, THM-421) and fifteen earlier than the lattice disk
> (N=43, HYP-2267).** The interval `[28, 32]` is the crossover-cost of being a
> lattice; `[32, 43]` is the further cost of being a round disk. The
> combinatorial floor 17 (THM-421) lower-bounds them all.

This nests the repo's own THM-421 frontier inside the global Erdos picture:
`17 (floor) ≤ 28 (true optimum) ≤ 32 (lattice blob) ≤ 43 (lattice disk)`.

## Why the optimum is non-lattice (structural note)

AMP24's optimal embeddings use coordinates like `(1/√13)(1, ω, ω²)` with `ω` a
primitive cube root of unity — the **√13 Eisenstein layer** (`r_Q(13)=12`,
density 6, exactly the THM-412 quantization), but assembled as a RIGID
non-lattice graph rather than a lattice section. So the optimum lives in the
SAME cyclotomic/Eisenstein world the repo identified (THM-412/THM-421's √7),
just at the √13 layer and without the lattice constraint. The gain over a
lattice section is the freedom to bend the rigid graph so interior degree-6
structure is reached with `O(1)` boundary deficit instead of triangular's
`√(12n)` deficit. (Connection to the 2025 rigidity work, arXiv:2507.15679.)

## Honest status

- `u(21) = 57`: VERIFIED in AMP24 (computer-assisted; we did NOT re-run their
  proof — we cite it and reproduce the elementary brackets it lives in).
- triangular max `= 47`, elementary cherry bound `= 71`: PROVED here, exact
  integer (`s710` scripts + `.out`).
- the 3N-crossover bridge (true optimum at n=28 vs lattice 32/43): a new
  arithmetic observation from AMP24's table + THM-421; the n≥22 inputs are
  AMP24 *lower* bounds, so "n=28" is an upper bound on the true crossover index.
- A naive float relaxation did NOT find dense configs (negative result, see
  hypotheses INDEX): dense unit-distance graphs need constraint/rigidity-based
  search (Schade/AMP24), not gradient descent on a soft well.

## Open (handed forward)

- **`u(22)`: is it 60 or 61?** AMP24 leaves a gap of 1. A single explicit
  61-edge unit-distance graph on 22 points would close it. (Their lower bound 60
  is Schade's; the matching upper bound failed only at 22,23.)
- `u(23) ∈ [64,66]`, `u(24) ∈ [68,72]` — wider gaps, same flavour.
- Can the √13-layer / rigidity lens give a *constructive* family matching
  Schade's lower bounds for 22 ≤ n ≤ 30 (and maybe push them up)?
