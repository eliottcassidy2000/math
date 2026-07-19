---
id: THM-1257
title: The D-GRADED GATE TOWER is a PRIMORIAL CASCADE — F_D(N) = {1..N}\{N−1} ∪ {D(N−1)} attains the slack-1 rung D/((N+1)D−1) iff N ≡ 1 (mod L_D) and N ≢ 1 (mod 2D−1), where L_D = 2·3·5···(largest prime < 2D−1) and 2D−1 must be PRIME; verified exactly for D=3 (N ≤ 100, extending HYP-4516's ≤ 37), D=4 (members exactly {31, 61, 91} in [8,100] — including the out-of-sample confirmations 4/247 at N=61 and the new 4/367 at N=91), D=5 (binder 9 = 3² composite: NEVER attains; F_5(211) degrades to the floor 1/212), and D=6 (opens at N=211: 6/1271, in a window of width 10⁻⁵). Each level's rung-killer (N ≡ 1 mod 2D−1) is the next level's enabler: L_{D'} = L_D·(2D−1). The N=61 single-far stratum is COMPLETELY classified (unique member 4/247).
status: >
  VERIFIED-EXACT (the gate law): exact M with witnesses for all 4 × 93 canonical
  families (D = 3..6, N = 8..100) matching the fitted gate with ZERO exceptions, plus
  the N=211 probes; member witnesses independently re-verified by residue checks.
  MECHANISM-DERIVED (branch level, this file §3): the b=2 parity/slot-absorption kill,
  the mod-3 double-duty kill (b=3 via 3|4N−1 and b=6 via 3|2N+1), the mod-5 see-saw
  kill (b=5 via 5|4N+1), and the rung-alive condition 7∤4N+3 — jointly equivalent to
  the fitted gate at D=4; not yet a theorem-grade proof of gate completeness (the
  exhibited-multiplier branch analyses need "no other competitor" sealed per N, which
  the exact sweep supplies empirically on [8,100]).
  PROVED — the N=61 single-far classification (absorption + sweep-first per-instance
  certificates): unique member (i=60, x=240), M = 4/247.
source: death-star-2026-07-19-S59c (HYP-7900; owner: extract the full D=4 gate and test N=61)
depends_on:
  - THM-1256  # the N=31 discovery + the see-saw (this file is its promised extraction)
  - THM-1255  # absorption machinery (the N=61 classification leg)
  - THM-1002  # exactness of the evaluator
related:
  - HYP-4516  # the mod-30 gate = the D=3 level, herewith verified to N=100
  - HYP-7890, HYP-7900
  - THM-1240/HYP-7840  # the N=13 wall: 4/55 is the D=4 rung one gate down (closed at 13 since 13 ≢ 1 mod 5)
scripts:
  - 04-computation/lrc_D_graded_gate_table_deathstar_S59c.py -> 05-knowledge/results/lrc_D_graded_gate_table_deathstar_S59c.out
  - 04-computation/lrc_N61_singlefar_classification_deathstar_S59c.py -> 05-knowledge/results/lrc_N61_singlefar_classification_deathstar_S59c.out
  - (probe) 05-knowledge/results/lrc_N211_D56_probe_deathstar_S59c.out
---

# THM-1257 — the gate tower

## 1. The law

For `D ≥ 3` with `p := 2D−1` PRIME, let `L_D := 2·∏{odd primes < p}` (so `L_3 = 6`,
`L_4 = 30`, `L_6 = 210`, `L_7 = 2310`). Then the canonical-D single-far family
`F_D(N) = {1..N}∖{N−1} ∪ {D(N−1)}` attains the slack-1 rung — and thereby populates
the first gap `W_N` —

```text
M(F_D(N)) = D/((N+1)D − 1)  ∈  W_N    ⟺    N ≡ 1 (mod L_D)  AND  N ≢ 1 (mod p).
```

For `2D−1` composite (first case `D = 5`, binder `9 = 3²`) the family NEVER attains
its rung. Verified exactly, all `(D, N) ∈ {3,4,5,6} × [8,100]` plus `N = 211`:

```text
D=3 (p=5):  members in [8,100] = {13,19,25,37,43,49,55,67,73,79,85,97}
            = the HYP-4516 gate exactly (previously verified ≤ 37; now ≤ 100).
D=4 (p=7):  members in [8,100] = {31, 61, 91}  = {N≡1 mod 30} ∖ {N≡1 mod 7}:
            4/127 (N=31, THM-1256), 4/247 (N=61, predicted then found, t=70/247),
            4/367 (N=91, NEW, t=53/367) — all with active pair (7, 4(N−1)).
D=5 (9=3²): NO members anywhere tested; F_5(211) DEGRADES TO THE FLOOR (M = 1/212,
            a second tight family at N=211 — the F(31)-degrade phenomenon recurs).
D=6 (p=11): no members in [8,100]; OPENS at N=211 = 1 mod 210, 211 ≡ 2 mod 11:
            M(F_6(211)) = 6/1271 exactly (t = 115/1271, pair (11, 1260)) — an
            in-window hit inside a window of width 1/(212·423) ≈ 1.1·10⁻⁵.
```

The N=61 and N=211 rows were PREDICTIONS of the fitted gate before computation
(HYP-7900 stub); both confirmed. Every member witness re-verified independently
(residue check at the stated `t`, active pair as stated).

## 2. The cascade identity

The gate's two clauses interlock across levels: the level-D rung dies iff
`p | (N+1)D − 1 ⟺ N ≡ 1 (mod p)`, and that SAME congruence is a factor of the
next level's enabler:

```text
L_{next(D)} = L_D · (2D−1).
```

So along `N ≡ 1 (mod L_D)`: every gate below D is closed (their rung-killers are
implied), and the level-D gate is open unless `N ≡ 1 (mod p)`, in which case the
cascade hands off to the next prime binder. The populated first gaps ride the
primorial progressions `N ≡ 1 (mod 2·3·5···p_k)` — the LRC first-gap spectrum is
governed by primorial arithmetic. Consequences: (a) attained ORDER `k = D−1` is
UNBOUNDED along the tower (`k = 5` already at N=211) — no cross-N uniform order
bound exists, sharpening why O-korder must be N-specific; (b) the next falsifiable
rung is `D=7` (binder 13) at `N = 2311`: predicted `M(F_7(2311)) = 7/16183`
(a 2311-speed exact computation — needs a pruned evaluator; filed as a lead).

## 3. The branch mechanism at D=4 (derivation, sealed by the sweep)

Competitor branches live at `Q_b = 4(N−1)+b`; clearance 4 at `b ≤ 5` (or 6, edge)
makes the family loose. Exhibited-multiplier analysis:

- **b=2** (`Q = 2(2N−1)`, `a ≡ ±2 mod 2N−1`): the residue ladder puts odd `u` at
  distance `|2u−2N+1|`; the distance-1 slot lands on `u = N` when N is ODD (kill)
  and on the DELETED `u = N−1` when N is even (alive — the same vacated-slot
  mechanism that powers the rung). So b=2 dies ⟺ **N odd**. (This explains
  N=16 ≡ 1 mod 15 failing: even.)
- **b=3** (`Q = 4N−1`): `3a ≡ ±4` unsolvable ⟺ `3 | 4N−1` ⟺ **N ≡ 1 (mod 3)**.
- **b=4** (`Q = 4N`, `a ≡ ±1 mod N`): the `u ≡ 3 (mod 4)` ladder's near-zero slot
  is deleted only when `N ≡ 0 (mod 4)`; members need N odd anyway.
- **b=5** (`Q = 4N+1`): `5a ≡ ±4` unsolvable ⟺ `5 | 4N+1` ⟺ **N ≡ 1 (mod 5)** —
  the THM-1256 see-saw (same congruence closes D=3's rung).
- **b=6** (`Q = 2(2N+1)`, value exactly the window top): `6a ≡ ±4` unsolvable ⟺
  `3 | 2N+1` ⟺ **N ≡ 1 (mod 3)** — the mod-3 condition does DOUBLE DUTY.
- **b=7 = rung**: alive ⟺ `7 ∤ 4N+3` ⟺ **N ≢ 1 (mod 7)**.

Intersection: `N ≡ 1 (mod 30)`, `N ≢ 1 (mod 7)` — the fitted gate. The sweep
seals "no other competitor" on [8,100] empirically; upgrading §3 to a full proof
(each branch analysis + an all-other-moduli bound) is the natural Lean target and
would make every gate member a THEOREM-grade gap witness at its N.

## 4. N=61 fully classified; scope ledger (typed)

- PROVED: the N=61 single-far stratum `{1..61}∖{i} ∪ {x}` (all i, all x) has the
  UNIQUE member `(i=60, x=240)`, `M = 4/247` — absorption + sweep-first
  certificates (safe-interval interiors certify `M(base) > 2/123` directly; no
  open-LRC citation; all 61 bases certified).
- VERIFIED (bounded): N=61 non-single-far census budgets (bordered/two-defect/
  targeted/repair) — see the out-file; and the [8,100] gate table is exact but
  says nothing about NON-canonical species at those N.
- OPEN: gate-completeness proof (§3), the D=7/N=2311 rung, non-canonical species
  along the primorial progressions, and everything at N=12/13 — the tower does
  NOT touch the (C)/(1/14,3/41) walls except as contrast: at N=13, D=4's gate is
  closed (`13 ≢ 1 mod 5` keeps the b=5 competitor alive at Q′=53), which is WHY
  4/55 needs a non-canonical realizer or none.
