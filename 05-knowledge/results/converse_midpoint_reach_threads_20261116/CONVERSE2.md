# The open half of the converse: what is proved, what cannot be structural, how far the census reaches

Session klein-2026-11-16-conv2 (2026-11-16). Owner: *"the open half of the converse (a near-join rewriting of
h(A)h(D) + B·P_D(2)) can be proven and extended."* Context: THM-5264 (this machine, earlier today) proved the out-degree-1
law and the dichotomy; the concurrent wave 56 of `chanpak_waves_20261113` (mac-mini) states the same tail formula and
reads the open half as *closure of near-join counts under the midpoint* `(h(A ⇒ D) + NJ_1)/2`. This note settles what can
be settled: the extended converse is a theorem, the rewriting cannot be structural (two values are proved to have no
same-core rewriting at any order), and the census says the open half is a density statement. Canon: THM-5265.

Notation. `NJ = {h(T) + 2^k B_t(T) : k ≥ 1}` (near-join values, HYP-9453). A **dead-end append** is `h(T) + B_t(T)` with
`e_t(T) = 0`; `NJ_ext = NJ ∪ {dead-end appends}`. `P_D(2) = Σ_W h(D[W]) h(D[W^c]) = 2 Q_D` (the terms for `W`, `W^c` pair
up). A **midpoint value** is `h(A)h(D) + B_{t_0}(A) P_D(2)` with `D` strong (THM-5264 (D)).

## 1. The extended converse (PROVED)

**Theorem 1.** (a) `NJ_ext` is exactly the set of counts of tournaments with a non-sink digit-0 vertex, and equals
`{counts} · {counts of strong tournaments with a digit-0 vertex}`. (b) `NJ ⊆ NJ_ext`, and `NJ_ext ∖ NJ` consists of
midpoint values only. (c) The odd numbers missing from `NJ_ext` are `1, 7, 21` and the *sterile-only* counts (every strong
realiser sterile); below 125 these are `11, 13, 47, 67` (THM-5264 (F)). (d) HYP-9453 (every odd number except the seven
is in `NJ`) is equivalent to: every odd `v ≥ 23`, `v ≠ 47, 67`, is the count of a strong tournament with a sole exit at a
dead-end port **and** every midpoint value is in `NJ`.

*Proof.* (a) THM-5264 (A): a non-sink digit-0 vertex is a sole exit at a dead-end port, so its count is a dead-end append;
conversely a dead-end append is such a count. Its strong component containing `z` is a strong digit-0 tournament and the
other components multiply (ordinal sums), giving the product form; closure `{counts}·NJ ⊆ NJ` is THM-5264 (E) and the same
argument gives `{counts}·NJ_ext ⊆ NJ_ext`. (b) THM-5264 (I3)/(D): a dead-end append whose base has a sink is a near-join
with `k ≥ 1`; otherwise it is a midpoint value. (c) and (d) restate (a), (b) and THM-5263 (B). ∎

So the "open half" is precisely `{midpoint values} ⊆ NJ`, i.e. `NJ_ext = NJ`. Through 398,129 it holds (§3).

## 2. No structural rewriting (PROVED for 15 and 37; FINITE-EXACT for the natural families)

**Theorem 2 (same-core obstruction).** Let `v = h(D) + P_D(2)` be the one-vertex-core midpoint value (`A = {t_0}`,
`S = t_0 → D → z → t_0`). For `D = C_3` (`v = 15`) and `D` the strong tournament of order 4 (`v = 37`) there is **no**
tournament `D″` of any order and no `k ≥ 1` with `v = h(D″) + 2^k P_{D″}(2)`: the midpoint cannot be written as a near-join
over the same core `{t_0}` with any right source.

*Proof.* `P_{D″}(2) = h(D″) + Σ_x B_x(D″) ≥ (n+1) h(D″) ≥ 2h(D″)` for `n = |D″| ≥ 1` (each `B_x ≥ h`), so
`v ≥ h(D″) + 4h(D″)` and `h(D″) ≤ v/5`. A tournament's count is the product of its strong components' counts, and strong
tournaments of order ≥ 4 have count ≥ 5 (Busch's floor), ≥ 9 at order ≥ 5. Hence `h = 1` means `D″` transitive,
`P = 2^n`; `h = 3` means one component `C_3`, `P = 12·2^j`; `h = 5` means one component the strong order-4 tournament,
`P = 32·2^j`; `h = 7` is impossible (7 is a hole). For `v = 15`: `h ∈ {1, 3}`, and `14 = 2^{k+n}` or `12 = 12·2^{k+j}`
force `k = 0`. For `v = 37`: `h ∈ {1, 3, 5}`, and `36 = 2^{k+n}`, `34 = 12·2^{k+j}`, `32 = 32·2^{k+j}` force `k = 0`. ∎

The rewritings that do exist are arithmetic: `15 = 3 + 4·3` (base `C_3`, `k = 2`, i.e. `3·5`) and `37 = 11 + 2·13` (base
the order-5 tournament of count 11 — the strong order-4 tournament with a sole exit appended at its port of end-count 1 —
with a root of response 13). Neither base is built from `{t_0}` and `D`.

**Identity (PROVED).** `P_D(2) = h(D) + Σ_{x ∈ D} B_x(D)` (pairs with `P_2` empty, plus pairs by the start of `P_2`), and
for `T_X = D + w` with `N⁺(w) = X ⊆ V(D)`: `B_w(T_X) = h(D) + Σ_{x ∈ X} B_x(D)`. So a rewriting of the one-vertex-core value
with base `T_X`, root `w` and `k = 1` is the single Diophantine condition
`h(T_X) = Σ_{x ∉ X} B_x(D) − Σ_{x ∈ X} B_x(D)` — satisfiable only by coincidence.

**The natural families, censused** ([receipts/hyp_w.out](receipts/hyp_w.out), [receipts/inside.out](receipts/inside.out),
[receipts/halved.out](receipts/halved.out)):

| rewriting family | success |
|---|---|
| base `D + w`, any out-neighbourhood of `w`, any root, `k ≥ 1` (strong `D` of orders 3–6) | 36 of 43; fails for seven order-6 sources (values 231, 319, 349, 349, 359, 393, 477) |
| the same with `w` a sole exit | 13 of 43 |
| base an induced subtournament of `S` (the 62 strong-source pairs of orders ≤ 8) | 38 of 62; 24 have none (37, 93 and 22 order-8 cases) |
| a one-vertex deletion `S − u` as base | 32 of 62 |
| a halved partner `(h(T), B_t(T)) = (h(D), P_D(2)/2)` of order ≤ 8 (strong `D` of orders 1–7) | 53 of 397 (by order of `D`: 1/1, 1/1, 0/1, 5/6, 20/35, 26/353); `C_3` has `C_3 → M_1` (the `k = 2` near-join in disguise), the strong order-4 source has none, the `(15, 216)` source has one of order 6 and also `strong4 → C_3` with the root of response 9 (`9·12 = 108`, THM-5260 (B)) |

Every known rewriting of a midpoint value uses a base unrelated to `(A, t_0, D)`. There is no uniform construction.

## 3. The reach of bases ≤ 9 (FINITE-EXACT; [receipts/nj9.out](receipts/nj9.out), `nj9.c`, 5 s)

All 198,948 tournament classes of orders 1–9, every root: **`NJ` realises every odd number ≤ 398,129 except
`1, 7, 11, 13, 21, 47, 67`**; the first miss is 398,131 (HYP-9453 had 22,227 from bases ≤ 8; wave 55's seeded grammar
43,343 with bases ≤ 12). Odd misses up to `2^23`: 1,378,932; by dyadic window `[2^j, 2^{j+1})`: none for `j ≤ 17`, then
8/131072, 1111/262144, 34040/524288, 269818/1048576, 1073955/2097152 for `j = 18, …, 22`. The dead-end appends from bases
≤ 9 (`NJ_ext`) have exactly the same misses: every `k = 0` value from a base of order ≤ 9 lies below the reach. Q-1
confirmed.

## 4. Density, not structure (NUMERICAL; [receipts/da_test.out](receipts/da_test.out))

Midpoint values from cores of order ≤ 7 (2,378 distinct rooted pairs) and strong right sources of order ≤ 8: 944,539
distinct values below `2^23`, of which 752,551 lie above the reach and 31,788 of those are missed by the bases-≤ 9 family.
Window by window, the miss rate of midpoint values against the gap density of all odd numbers:

| window | midpoint values | missed | rate | odd gap density | observed/expected |
|---|---|---|---|---|---|
| `[2^18, 2^19)` | 129,428 | 8 | 0.0001 | 0.0001 | 1.01 |
| `[2^19, 2^20)` | 255,900 | 1,069 | 0.0042 | 0.0042 | 0.99 |
| `[2^20, 2^21)` | 372,623 | 21,072 | 0.0566 | 0.0649 | 0.87 |
| `[2^21, 2^22)` | 61,771 | 9,640 | 0.1561 | 0.2573 | 0.61 |

In the first two windows above the reach, midpoint values are missed exactly at the rate of random odd numbers (Q-2
confirmed): nothing in the bases-≤ 9 family favours them. The deficit at `j ≥ 20` is the arithmetic of large products
(many midpoint values there are multiples of small counts, which the closure of THM-5264 (E) covers), not a rewriting.

## 5. Verdict on "can be proven and extended"

- **Extended — yes, and proved**: Theorem 1 (`NJ_ext`, its exceptions, the product form, the equivalence with HYP-9453);
  THM-5264's closure and sink trick; the identities of §2.
- **Proven — not available as a rewriting**: Theorem 2 shows that for 15 and 37 no same-core rewriting exists at any order;
  the four natural one-step families fail on 7 of 43, 30 of 43, 24 of 62 and 30 of 62 cases; and above the reach the
  midpoint values behave like random odd numbers. The open half is therefore a density statement, equivalent to
  `NJ_ext = NJ`, implied by HYP-9453's completeness and not provable by a construction on `(A, t_0, D)`. Chan–Pak OPEN.
