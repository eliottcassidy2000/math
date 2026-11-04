# The sole-exit tower is tribonacci, and the three root-response numerations of the odd spectrum

klein-2026-11-04 (worktree `codex/session-spectrum-20261104`), continuing [`two_branch_rearming_sigma_sweep_20261104`](../two_branch_rearming_sigma_sweep_20261104/README.md) (THM-4862, HYP-9454) and [`spectrum_analogies_nearjoin_20261102`](../spectrum_analogies_nearjoin_20261102/) (HYP-9453), inside the nonzero-preparation lineage ([`nonzero_preparation_20261103`](../nonzero_preparation_20261103/README.md), [`multiport_defect_grammar_20261103/construction/RESULT.md`](../multiport_defect_grammar_20261103/construction/RESULT.md)). **Status: PROVED (the tower law, elementary, [PROOF.md](PROOF.md); exhaustive and random verification with byte-identical receipts), FINITE-EXACT (the root-response censuses over all rooted tournaments of orders 4–9), and a reformulation of the two sweep conjectures as one statement about the realisable triples `(h, e, B)`. Nothing about Chan–Pak completeness is closed.** Canon entry: THM-4865.

## 1. The theorem

For a tournament `T` and a vertex `t` let `h = h(T)` (Hamiltonian paths), `e = e_t(T)` (those ending at `t`) and `B = B_t(T) = Σ_{A ⊆ V∖t} h(T[A]) s_t(T[V∖A])` (the *root response*: pairs of a Hamiltonian path on a subset and one on the complement starting at `t`). The *sole-exit tower* over `(T, t)` is `T_k = T + x_1 + … + x_k`, where `x_1 → t`, `x_j → x_{j−1}`, and every other vertex beats `x_j` (each new vertex has out-degree one, pointing at the previous one). Then

`h(T_1) = h − e + B`,  `h(T_2) = h + 2B`,  `h(T_k) = h(T_{k−1}) + h(T_{k−2}) + h(T_{k−3})`  (`k ≥ 3`),

so `h(T_k) = h·A_k − e·C_k + B·D_k` with `A = A000213 (1,1,1,3,5,9,17,…)`, `C = A001590 (0,1,0,1,2,3,6,…)`, `D_k = C_{k+2}`. The tribonacci recurrence holds for **every** base `(T, t)`; the base enters only through the three numbers `(h, e, B)`. Proof: [PROOF.md](PROOF.md) (three lemmas: the one-exit law THM-4862 (a), `e_x(T+x) = h − e_t`, `B_x(T+x) = B_t + h`, then two lines of algebra). Verification: [`scripts/tower_check.py`](scripts/tower_check.py) — every rooted tournament of order ≤ 4 to height 5 and of order 5 to height 3, 40 seeded random bases of orders 6–7, the closed form at every height, and the tower over one vertex reproducing `A000213` ([receipt](receipts/tower_check.normal.out)).

## 2. What it unifies (all inherited statements are cited, none re-proved as new)

- **The lineage's C2 operator is the tower of height two.** The births `x → a` (all else `→ x`) and `y → x` (all else `→ y`) are `T_2` over `(T, a)`, so the child count is `h + 2B_a(T) = 3h + 2σ` with `σ = B_a(T) − h(T)`, which is the lineage's `σ = S − h` in law (10) of `construction/RESULT.md`, valid for **any** auxiliary vertex `a`. The sole-exit condition of the "positive joint guard" (JOINT_CLASSIFICATION.md) is needed for the second functional `F` (`τ = 0`), not for `H`. The sweep packet (THM-4862 (d)) used sole-exit vertices only; §3 drops that restriction.
- **The insertion-response calculus THM-1900** (opus-2026-07-20-S439; corrected 2026-10-20) treats the general insertion of a vertex beating exactly a set `P` and characterises the `H`-neutral patterns (down-sets of the condensation) by counting insertion positions on old paths plus the orphan paths of THM-012b. The one-exit law is the exact `H`-response of the pattern `P = {t}`: the orphan-plus-insertion count is `B_t(T)`, the lost paths are the `e_t(T)` old paths ending at `t`. Nothing in THM-1900 is re-proved here; the tower law is the iteration of this one pattern.
- **The digit identity THM-4862 (b)** (`σ = e_r(T − a)` at a sole exit `a → r`) is one line from the lemmas.
- **The tribonacci partner THM-4731 (6)** (the full tiling tournament `T*_n`, `H = A000213(n)`, discovered and proved by kind-pasteur-2026-03-05-S11; MISTAKE-645 records a later re-derivation) is the tower over a single vertex, `(h, e, B) = (1, 1, 1)`; its auxiliary recursion `e(n) = h(n−1) − e(n−1)` is the second tower relation.
- **The near-join of HYP-9453** (`T → M_k → z`, `z → t` reversed; `H = h + 2^k B(T, t)`, proved there) reads the same `B`; at `k = 1` it has the value of the height-two tower. Its census bound (all odd ≤ 22,227 from bases of order ≤ 8, except 1, 7, 11, 13, 21, 47, 67) is reproduced exactly by the present `B` ([receipt](receipts/coverage_8.normal.out)), an independent check of both computations.

## 3. The three root-response numerations (FINITE-EXACT)

[`scripts/allpairs.c`](scripts/allpairs.c) computes `(h, e_t, B_t, strong)` for every tournament `T` of order `n` (nauty `gentourng`) and every vertex `t`; [`data/triples_n.txt`](data/) holds the distinct rows with multiplicities for `n = 4..8` (orders 4–9 together: 1,782,468 rooted tournaments, 1,728,779 distinct rows; the order-9 file, 29 MB, is regenerated in ten seconds by [`scripts/dedupe_triples.sh`](scripts/dedupe_triples.sh) and is not committed). [`scripts/coverage.py`](scripts/coverage.py) forms three value sets from the rows and reports the largest `X` such that every odd number in `[69, X]` is attained:

| family (value of the base `(h, e, B)`) | bases of order ≤ 8 | bases of order ≤ 9 | bases of order ≤ 10 (streamed) |
|---|---|---|---|
| two-step `h + 2B` (C2 child at any vertex; tower height 2; near-join `k = 1`) | `[69, 3141]` | `[69, 15449]` | `[69, 86229]` |
| the same, strong bases only | `[69, 3141]` | `[69, 15449]` | — |
| near-join `h + 2^k B`, `1 ≤ k ≤ 8` | `[69, 22227]` | `[69, 382687]` | `[69, 6171561]` |
| tower `h_k = h A_k − e C_k + B D_k`, heights `3 ≤ k ≤ 10` alone | `[69, 77]` | `[69, 77]` | — |
| tower, heights `1 ≤ k ≤ 10` (`k = 1, 2` are the one-step and two-step values) | `[69, 20351]` | `[69, 432485]` | `[69, 7854129]` |
| union of the three | `[69, 36165]` | `[69, 702189]` | `[69, 8079323]` |

The order-10 column streams the 99,113,028 rooted bases of orders 4–10 (1.7 GB of rows, not stored) through [`scripts/coverage_stream.c`](scripts/coverage_stream.c) (bitsets to `2³¹`); the record is [`receipts/coverage_10_stream.out`](receipts/coverage_10_stream.out), a single C run reproducible by `gentourng -q 10 | ./allpairs 10 | ./coverage_stream` (about forty minutes).

Below 69 the two-step family misses exactly `1, 11, 13, 47, 67` (besides the never-attained 7 and 21, THM-1370), the same exceptional set as the near-join; strong bases alone miss `1, 3, 5, 9, 11, 13, 25, 33, 47, 63, 67` there; the tower values of height ≥ 3 alone are thin at the bottom (gaps at 79, 83, 103, …) and matter only through the union, which they extend from 382,687 to 702,189 at order ≤ 9. Receipts: [`coverage_8`](receipts/coverage_8.normal.out) (from the committed rows) and [`coverage_9`](receipts/coverage_9.normal.out) (from the regenerated order-9 rows), each byte-identical under `python3 -O`.

Compared with the sole-exit sweep of THM-4862 (d) — every odd in `[307, 5153]` from sole-exit pairs of order ≤ 9, `[307, 26551]` from order ≤ 10 — using every vertex as the root lowers the threshold to 69 and raises the ceiling threefold at the same order. The digit sets are correspondingly larger: for `h ∈ [300, 1000)` the all-vertex digit set `{B_t(T) − h : h(T) = h}` has about 857 members on average at orders ≤ 9 (sole exits at orders ≤ 10: about 86), with maximal digit about `10h` (sole exits: `0.29h`).

## 4. The reformulated conjecture

Let `𝒯 ⊆ ℕ³` be the set of triples `(h(T), e_t(T), B_t(T))` over all rooted tournaments. Then the odd spectrum contains `{h + 2B, h + 2^k B, h A_k − e C_k + B D_k : (h, e, B) ∈ 𝒯, k}`. The ready spectra of the sweep packet are the `e`-slices of the one-step values: `R_σ = {h − σ + B : (h, e, B) ∈ 𝒯, e = σ}` (the digit of the new vertex is `e_t(T)`, Lemma 2), and the digit-1/digit-2 floors `9, 17, 33, 61, 115, 219` and `9, 23, 51, 93, 177, 341` at orders 5–10 are recovered from the strong rows of the census as the minima of these slices. HYP-9453 (1) and HYP-9454 (1) are the statements that already the first family, with bases of order `O(log v)`, contains every odd `v ≥ 69`; the census above is the evidence through 702,189. The proof problem is now a statement about `𝒯` alone: **which pairs `(h, B)` are realisable?** Known points: transitive bases give `(1, 0, 2^b)` and `(1, 1, 1)`; the tower itself maps `(h, e, B)` to `(h_k, h_{k−1} − e_{k−1}, B + h_0 + … + h_{k−1})`, so `𝒯` is closed under an explicit affine map; the near-join and the tower give two different weightings (`2^k` and tribonacci) of the same `B`. What is missing is a family in which `B` moves by controlled small steps at fixed `h` (HYP-9453 (2) in the old language); the one-exit law moves `h` by `B − e`, which is at least `h − e ≥ 0`; `h`-preserving one-vertex additions are sinks (`B ↦ 2B`) and sources (`B` unchanged) only. A negative control ([`scripts/triple_bank.py`](scripts/triple_bank.py), [receipt](receipts/triple_bank.normal.out)): the affine bank generated from the single vertex `(1, 1, 1)` by the sink map `(h, e, B) ↦ (h, 0, 2B)`, the sole-exit map `(h, e, B) ↦ (h − e + B, h − e, B + h)` (every element is an explicit tournament) and the `C₃` prefix and suffix joins `(3h, 3e, 3B)`, `(3h, 0, 12B)` reaches only 22,403 of the 100,000 odd numbers below 2·10⁵ (first misses 7, 11, 13, 19, 21, 23, 35, …): the numerations need genuinely general bases, and no small closed family of constructions can replace the census.

## 5. Two integer sequences not in OEIS (2026-11-04)

- `5, 13, 31, 77, 193, 489, 1255, 3253, 8493, 22293` (orders 4–13): the Hamiltonian-path count of the digit-1 ready tournament whose entry attachment is the pure backward cycle (`F = ∅` in the forward-chord picture of [`data/chords_all.txt`](data/chords_all.txt): every digit-1 class of orders 8–13 written as the set of forward chords of the unique Hamiltonian cycle of `Φ(T)`, one line per class, columns `order H chords`); it equals the tower over a point of height `n − 2` with one more sole-exit vertex pointing at the base vertex ([`scripts/fempty_check.py`](scripts/fempty_check.py), [receipt](receipts/fempty_check.normal.out)).
- `5, 9, 17, 33, 61, 115, 219, 423, 823, 1615, 3183` (orders 4–14): the digit-1 ready floors (minimal `h(T)` of a strong tournament with a sole-exit vertex of digit 1), from the sweep packet's censuses; the deficits from `2^{n−2} + 1` are `0, 0, 0, 0, 4, 14, 38, 90, 202, 434, 914`.

## 6. Reproduction

```
cd scripts
cc -O2 -o allpairs allpairs.c
for n in 4 5 6 7 8 9; do gentourng -q $n | ./allpairs $n | sort | uniq -c | awk '{print $2,$3,$4,$5,$1}' > ../data/triples_$n.txt; done
python3 tower_check.py; python3 fempty_check.py; python3 coverage.py 8; python3 coverage.py 9
cd .. && shasum -a 256 -c SHA256SUMS
```
