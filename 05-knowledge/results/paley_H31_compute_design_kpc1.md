# Design doc: computing H(T_31) exactly — the run that falsifies or confirms HYP-2371

*kind-pasteur-2026-06-10-S1, Thread E. Companion artifacts:
`04-computation/paley_R31_prediction_kpc1.py` (+ `.out`) — the falsifiable prediction;
`04-computation/paley_rotdp_smallp_verify_kpc1.py` (+ `.out`) — the harness logic
validated exactly at p = 11, 19.*

## 0. What the run decides

HYP-2371 predicts (from the proven cluster expansion HYP-2307/THM-438 — functional
form `R(p) = e(1 − C/p − D/p² − …)` proven-justified, coefficients pinned on the 5
exact values, order-laddered holdout at p=23):

```
R(31) = H(T_31)·2^30/31!  ∈  [2.58949, 2.60249]      (central 2.59599 ± 0.00650)
H(T_31) ∈ [19830629617139608462365775, 19930130881568868002912737]   (~1.988e25)
H(T_31) ≡ 465 (mod 930)        (i.e. H = 465 × odd; |Aff(QR_31)| = 31·30/2 = 465)
```

One exact computation of H(T_31) lands inside or outside. Acceptance tests for any
claimed value (run them ALL before believing the number):

1. `H` is odd (Rédei).
2. `465 | H` (freeness of Aff(QR) on Hamiltonian paths, THM-F3 / Thread B).
3. `H/465` is odd (orbit parity, HYP-1713 PROVED; re-verified on all 5 known values).
4. `R(31) = H·2^30/31!` vs the interval above — the falsifiable claim.
5. Same binary must reproduce `H(T_23) = 15760206976379349` first (seconds).

## 1. Why naive Held–Karp is out

Start-fixed Held–Karp stores `dp[S][v]`, `S ⊆ V∖{0}`, i.e. `2^30 × 31 ≈ 3.3e10`
states × 16 B (uint128) ≈ **532 GB**. Dead on a desktop.

## 2. Design A — Z_31-rotation-reduced layered Held–Karp (the brief's design)

### 2.1 States up to rotation

Drop the fixed start; let `dp(S, v)` = number of directed Hamiltonian paths of the
induced subtournament on `S` ending at `v ∈ S` (start anywhere in `S`). The Paley arc
relation `v→w ⟺ χ(w−v)=+1` is invariant under `x ↦ x+r`, so the rotation group `Z_31`
acts on states by `(S,v) ↦ (S+r, v+r)` and `dp` is constant on orbits.

**Freeness (proved, elementary):** `S+r = S` with `r ≠ 0` forces `S ∈ {∅, Z_31}`
(31 prime ⟹ ⟨r⟩ = Z_31); and for `S = Z_31`, `v+r ≠ v` for `r ≠ 0`. So the action is
FREE on every state with `S ≠ ∅`: all orbits have size exactly 31, and one canonical
representative per orbit suffices. Verified mechanically at p = 11, 19
(`paley_rotdp_smallp_verify_kpc1.out`: exact H values reproduced).

- Canonical form: `min over r∈{0..30} of key(rot_r(S), v+r)`, `key = (mask << 5) | v`
  (36-bit key; min over 31 rotations; with the doubled-mask trick `(S | S<<31) >> r`
  this is ~31 shift+compare ≈ 10–30 ns in C).
- Total canonical states: `Σ_k C(31,k)·k/31 = 31·2^30/31 = 2^30 ≈ 1.074e9`.
- Initialization: single orbit `({0}, 0)` with `dp = 1`.
- Answer: `H(T_31) = 31 × dp(Z_31, 0)` (one full-set orbit).

### 2.2 Popcount-layer streaming and the memory arithmetic

Transitions only go from layer `k = |S|` to layer `k+1` (append `w ∉ S` with arc
`v→w`), so only two layers are ever live. Canonical states in layer `k` number
**exactly ≤ C(30, k−1)** (equality up to unreachable states, which never materialize;
at p=19 the deficit was 10 states out of 48620 at peak):

| layer k | states = C(30,k−1) | counts (uint128) |
|---|---|---|
| 14 | 119,759,850 | 1.78 GiB |
| 15 | 145,422,675 | 2.17 GiB |
| **16 (peak)** | **155,117,520** | **2.31 GiB** |
| 17 | 145,422,675 | 2.17 GiB |
| 18 | 119,759,850 | 1.78 GiB |

Peak-layer arithmetic: `C(31,16)·16/31 = 300540195·16/31 = 155117520` states ×16 B =
2,481,880,320 B = **2.31 GiB ≈ the "~2.3 GB per layer" claim**. Two live layers of
counts ≈ 4.5 GiB, plus 5-byte keys (sorted canonical keys per layer: 0.72 + 0.68 GiB)
≈ **5.9 GiB** — fits an 8 GB desktop only with nothing else resident; comfortable on
16 GB. Alternative: run twice mod two 63-bit primes with uint64 counts (halves count
memory to ≤ 2.95 GiB total; CRT at the end — `H < 2.1e25 < m₁·m₂`).

### 2.3 Canonical-form indexing

Per layer, keep a sorted array of 40-bit keys + parallel count array. Successor
generation emits `(canonical key, count)` records; aggregate by radix-bucketing on the
top 12 key bits into RAM-sized buckets, sort/hash-reduce each bucket, then merge into
the next layer's sorted arrays. No pointer hash table needed (cache-friendly,
disk-spillable: peak transition stream from layer 16 ≈ 155M states × ~7.5 successors
× 21 B ≈ 24 GB — stream through NVMe if RAM-tight).

### 2.4 Work and wall time

Total transitions: `Σ_k C(31,k)·k·(31−k)/(2·31) = 31·30·2^29/(2·31) = 8,053,063,680`
(~8.1e9; avg out-degree `(31−k)/2`, exactly verified by the formula in the harness
script). Per transition: successor mask + canonicalization + emit ≈ 50–150 ns
optimized C/C++ → **7–35 min core-time**, plus bucket sort/merge I/O (~10–60 min
depending on RAM vs NVMe). Order of magnitude: **~1 hour single node (8 threads),
range 0.5–3 h**. NOT provably under 30 min on this desktop, and Python is ~100×
off — hence prediction-first, compute later (per the session brief).

### 2.5 Where overflow forces 128-bit

Average dp value at layer k ≈ (k−1)!/2^(k−1) (coin-flip heuristic; Paley runs ~e×
higher near the end):

| k | (k−1)!/2^(k−1) |
|---|---|
| 24 | 3.7e16 |
| 26 | 4.6e17 |
| 27 | 6.0e18 |
| 28 | 8.1e19 > 2^64−1 = 1.84e19 |

Counts blow past uint64 at layer ≈ 27–28 on average, earlier for the largest states
(distribution skew) — switch to uint128 at **layer 24** to be safe, or simplest: use
uint128 throughout (the 2.31 GiB figure already assumes 16 B/state), or CRT with two
uint64 moduli end-to-end. Use checked adds if running uint64 in early layers.

## 3. Design B (recommended cross-check) — Karp inclusion–exclusion, O(KB) memory

`H = Σ_{S⊆V} (−1)^{31−|S|} · W_30(S)`, where `W_30(S)` = #length-30 walks inside the
induced subtournament on `S` (computed by 30 mat-vec iterations, `O(30·|S|²)` int
mults). Rotation classes: orbits of subsets number `(2^31 − 2)/31 + 2 ≈ 6.93e7`;
orbit sizes are 31 except `∅` and `Z_31` (size 1). Work ≈ `6.93e7 × 30·E[|S|²] ≈
6.93e7 × 7.4e3 ≈ 5.2e11` mults per modulus; ×2 CRT moduli ≈ 1e12 — **~0.5–2 h on 8
cores**, embarrassingly parallel, ~zero memory, trivially checkpointable (stride the
class index). Walk counts overflow even uint128 (≤ 31^30), so work mod two 63-bit
primes from the start; CRT the final alternating sum (`0 < H < 2.1e25 < m₁m₂`).
Independent of Design A in both algorithm and failure modes — running both is the
gold standard for a number this size.

## 4. Validation chain for the compute node

1. `paley_rotdp_smallp_verify_kpc1.py` logic at p = 11, 19 — DONE this session (exact).
2. Port to C/C++; reproduce p = 19 (instant) and **p = 23 = 15760206976379349** (~secs).
3. Design B at p = 19, 23 with CRT — must agree with Design A.
4. p = 31 full run, both designs; apply acceptance tests of §0.
5. Report `R(31)` to 12+ digits and the verdict on HYP-2371's interval
   `[2.58949, 2.60249]`; update the hypothesis INDEX (PREDICTION → CONFIRMED/REFUTED),
   OPEN-Q-013 line, and the A038375-adjacent data (new record term).

## 5. p = 43, 47 outlook (honest)

- Design A peak layer at p=43: `C(42,21) ≈ 5.4e11` states × 16 B ≈ **8.6 TB** — dead.
- Design B at p=43: `2^43/43 ≈ 2.0e11` classes × `42·E[|S|²] ≈ 1.9e4` ≈ 3.9e15 mults
  per modulus — **months of core-time**; feasible only as a distributed campaign.
- The asymptotic predictor (same machinery as HYP-2371) gives `R(43) ≈ 2.628 ± 0.01`,
  `R(47) ≈ 2.635 ± 0.01` — but p=31 is the practical falsification frontier; p≥43
  needs genuinely new mathematics (or the sub-leading-coefficient route: prove C
  exactly from the cluster expansion and test it on the p=31 value instead).
