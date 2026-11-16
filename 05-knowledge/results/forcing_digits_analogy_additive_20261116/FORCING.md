# Forcing the midpoint and the +2 child from an arbitrary prime: what can and cannot be forced

Session klein-2026-11-16-force. Scripts `scripts/digits.c`, `scripts/strongB.c`, `scripts/padding_check.py`,
`scripts/module_blind.py`; receipts `digits_4..9.out`, `strongB_3..9.out`, `padding_check.out`, `module_blind.out`.

## 1. The midpoint is forced; no history-blind gadget forces the +2 (THM-5260 (A))

Let `Q` be a tournament with a distinguished vertex `τ` and `Q[T]` the tournament obtained by substituting an
arbitrary `T` for `τ` (every other vertex of `Q` beats all of `T` or loses to all of `T` as it beat or lost to `τ`).
This is the general "history-blind" extension: the gadget cannot see inside `T`.

**Theorem A.** (a) If `τ` lies on no directed cycle of `Q`, then `Q = Q_1 → τ → Q_2` is an ordinal sum and
`H(Q[T]) = H(Q_1)·H(T)·H(Q_2)` for every `T`. (b) If `τ` lies on a directed cycle of `Q`, then `H(Q[T])` is not a
function of `H(T)`: FINITE-EXACT for all 110 pointed gadgets of order ≤ 4 (`module_blind.out`: "depends only on `H`"
⟺ "`τ` on no cycle" with no exception), and the reason is structural — every directed cycle of `Q[T]` meets `T` in one
contiguous run (a cycle cannot enter `T` from a vertex that beats all of `T` and later leave to a vertex that loses to
all of `T` twice over), so the odd-cycle formula (THM-002) writes `H(Q[T])` as a positive combination of the run
sums `Σ_P H(T − P)` over directed paths `P` of `T` of prescribed parities, which equal-count tournaments do not share
(the 3-cycle and the 3-cycle-with-sink both have `H = 3` and odd-run sums 6 and 12).

So the midpoint `3p` is forced from any realiser of `p` by the ordinal sum `C_3 → T` (and, equally, by the height-two
tower at any digit-0 vertex), while **no uniform gadget forces `3p + 2`**: the `+2` must read the realiser. The
owner's "rather than finding them by search" therefore has to mean a *selection rule* on the realiser, and the rest of
this note finds the exact one for the height-two tower.

## 2. Digits, charges, and what padding does to them (THM-5260 (B), PROVED; `padding_check.out`)

For a vertex `t` of `T` let `h = H(T)`, `e_t` = Hamiltonian paths ending at `t`, `B_t = Σ_{A ⊆ V∖t} h(T[A]) s_t(T[V∖A])`
the root response (THM-4865), `σ_t = B_t − h` the digit (the height-two tower at `t` has count `3h + 2σ_t`), and
`a_t = B_t − e_t − h + 1` the charge (the rooted append at `t` has count `2h − 1 + a_t`; wave 45's charge
`2E_t + 1 − h` agrees because `B_t − e_t = 2E_t` for every rooted tournament of order ≤ 9 — `digits_*.out`).

**Ordinal-sum laws.** For `U = A → B`: a vertex `t ∈ B` has `B_t(U) = h(A) B_t(B)`, `e_t(U) = h(A) e_t(B)`, hence
`a_t(U) = h(A)(a_t(B) − 1) + 1`; a vertex `t ∈ A` has `B_t(U) = B_t(A) P_B(2)`, `e_t(U) = 0`, hence
`a_t(U) = B_t(A) P_B(2) − h(A)h(B) + 1`, where `P_B(2) = Σ_S h(B[S]) h(B[V∖S])` (286 random checks, no failure).
*Proof.* A path starting at `t ∈ B` cannot enter `A`, so the complement of `A` in the sum must be all of `A`; a path
starting at `t ∈ A` runs through its `A`-part and then any Hamiltonian path of its `B`-part. ∎

**Corollary (padding).** A transitive prefix leaves `(h, e_t, B_t)` unchanged; a transitive suffix of length `k`
gives `B_t ↦ 2^k B_t`, `e_t ↦ 0`, so `σ_t ↦ 2^k B_t − h` and `a_t ↦ 2^k B_t − h + 1`. Consequently, for a count `h`
whose realisers are the paddings of strong realisers `S` (every prime, and every `h < 125` with one nontrivial strong
component — the strong floor `f(10) = 125`):

- **digit 1 (the tower gives `3h + 2`) is available at some order iff some strong realiser has a vertex with
  `B_t ∈ {h + 1, (h + 1)/2, (h + 1)/4, …}`;**
- digit 2 is available iff some strong realiser has `B_t = h + 2` (odd numbers are not halved); **digit 0 is always
  available at some order** — the sink `s` of `S → s` has `B_s = h` (MISTAKE-776; the column "digit 0" of §3 refers to the
  vertices of strong realisers only);
- a rooted doubling child `2h − 1` (charge 0) exists at some order iff some strong realiser has `a_t = 0` or
  `B_t = (h − 1)/2^k` for some `k ≥ 1`.

## 3. The exact table for `h < 125` (FINITE-EXACT, strong realisers of orders ≤ 9)

The table below (receipt `forcing_table.out`; strong realisers of orders 3–9 cover every count below the strong floor `f(10) = 125`) gives, for each count `h`, whether a strong realiser has a vertex with `B_t = h` (digit 0 among strong vertices; every count has a digit-0 vertex once a sink is padded, MISTAKE-776), with `2^k B_t = h + 1` for some `k ≥ 0` (digit 1 at some order, i.e. `3h + 2` tower-forced), with `B_t = h + 2` (digit 2), and whether a rooted doubling child `2h − 1` exists at some order (`a_t = 0` or `2^k B_t = h − 1`), with the first `B`-values.

| `h` | digit 0 | digit 1 (any order) | digit 2 | rooted doubling | smallest `B_t` of strong realisers |
|---|---|---|---|---|---|
| 3 | 1 | 0 | 0 | 1 | [3] |
| 5 | 1 | 1 | 1 | 1 | [5, 6, 7, 9] |
| 9 | 1 | 1 | 1 | 1 | [9, 10, 11, 15, 16, 21, 23, 25] |
| 11 | 0 | 1 | 1 | 0 | [12, 13, 17, 20, 27] |
| 13 | 0 | 1 | 1 | 1 | [14, 15, 18, 26] |
| 15 | 1 | 0 | 0 | 1 | [15, 18, 19, 21, 23, 29, 40, 45, 57, 59, 63] |
| 17 | 1 | 1 | 0 | 1 | [17, 18, 20, 21, 28, 31, 33, 44, 55, 61, 69] |
| 19 | 1 | 1 | 1 | 0 | [19, 20, 21, 23, 27, 29, 35, 43, 52, 59, 73] |
| 23 | 1 | 1 | 1 | 1 | [23, 24, 25, 27, 29, 30, 31, 33, 36, 40, 48, 50] ... |
| 25 | 1 | 0 | 1 | 1 | [25, 27, 28, 29, 30, 34, 39, 42, 49, 50, 52, 60] ... |
| 27 | 1 | 0 | 0 | 1 | [27, 30, 33, 34, 35, 45, 47, 48, 53, 68, 69, 75] ... |
| 29 | 1 | 1 | 1 | 0 | [29, 30, 31, 33, 34, 35, 39, 41, 43, 44, 45, 49] ... |
| 31 | 1 | 1 | 1 | 0 | [31, 32, 33, 34, 37, 39, 41, 43, 45, 51, 54, 60] ... |
| 33 | 1 | 1 | 0 | 1 | [33, 34, 36, 37, 38, 39, 40, 41, 42, 45, 46, 48] ... |
| 35 | 1 | 0 | 0 | 0 | [35, 38, 39, 42, 43, 55, 62, 81, 91, 113, 136, 173] ... |
| 37 | 1 | 1 | 1 | 1 | [37, 38, 39, 43, 44, 45, 46, 47, 49, 50, 51, 53] ... |
| 39 | 1 | 1 | 1 | 0 | [39, 40, 41, 47, 59, 60, 83, 89, 104, 107, 156, 159] ... |
| 41 | 1 | 0 | 0 | 1 | [41, 44, 45, 46, 47, 48, 49, 50, 51, 53, 54, 57] ... |
| 43 | 1 | 0 | 0 | 0 | [43, 48, 49, 53, 54, 55, 56, 57, 64, 66, 68, 72] ... |
| 45 | 1 | 1 | 0 | 1 | [45, 46, 49, 50, 53, 54, 55, 57, 58, 59, 63, 65] ... |
| 47 | 0 | 1 | 1 | 0 | [48, 49, 61, 65, 72, 80, 108, 120, 161, 176, 180, 207] ... |
| 49 | 1 | 1 | 1 | 1 | [49, 50, 51, 52, 53, 54, 56, 57, 58, 60, 63, 66] ... |
| 51 | 1 | 1 | 1 | 1 | [51, 52, 53, 54, 57, 59, 60, 61, 62, 63, 67, 70] ... |
| 53 | 1 | 0 | 1 | 0 | [53, 55, 56, 58, 59, 61, 62, 63, 65, 67, 69, 71] ... |
| 55 | 1 | 0 | 1 | 0 | [55, 57, 59, 60, 61, 63, 75, 85, 87, 91, 93, 94] ... |
| 57 | 1 | 0 | 0 | 1 | [57, 60, 63, 67, 68, 69, 71, 79, 81, 83, 85, 87] ... |
| 59 | 1 | 1 | 1 | 1 | [59, 60, 61, 62, 63, 64, 65, 67, 69, 70, 71, 73] ... |
| 61 | 1 | 1 | 0 | 0 | [61, 62, 66, 67, 69, 71, 75, 79, 83, 89, 90, 91] ... |
| 63 | 1 | 1 | 0 | 0 | [63, 64, 74, 77, 89, 103, 116, 121, 179, 220, 331, 344] ... |
| 65 | 1 | 1 | 1 | 1 | [65, 66, 67, 68, 69, 70, 73, 74, 75, 76, 77, 78] ... |
| 67 | 0 | 0 | 0 | 0 | [71, 72, 73, 74, 75, 77, 92, 100, 104, 107, 110, 112] ... |
| 69 | 1 | 1 | 1 | 1 | [69, 70, 71, 72, 73, 74, 75, 76, 77, 78, 79, 81] ... |
| 71 | 1 | 1 | 1 | 0 | [71, 72, 73, 74, 75, 77, 78, 79, 81, 84, 85, 86] ... |
| 73 | 1 | 1 | 1 | 0 | [73, 74, 75, 76, 77, 78, 80, 81, 82, 83, 85, 86] ... |
| 75 | 1 | 0 | 1 | 1 | [75, 77, 78, 79, 81, 83, 84, 85, 87, 90, 91, 93] ... |
| 77 | 1 | 1 | 1 | 0 | [77, 78, 79, 80, 81, 82, 83, 85, 86, 87, 91, 94] ... |
| 79 | 1 | 1 | 1 | 1 | [79, 80, 81, 82, 83, 84, 85, 87, 88, 89, 90, 91] ... |
| 81 | 1 | 1 | 1 | 1 | [81, 82, 83, 85, 86, 87, 89, 90, 91, 92, 93, 96] ... |
| 83 | 1 | 0 | 0 | 1 | [83, 86, 89, 90, 91, 92, 93, 95, 99, 101, 102, 103] ... |
| 85 | 1 | 1 | 0 | 0 | [85, 86, 90, 91, 94, 98, 100, 102, 103, 105, 109, 111] ... |
| 87 | 1 | 1 | 0 | 0 | [87, 88, 91, 92, 94, 97, 101, 102, 103, 105, 107, 109] ... |
| 89 | 1 | 0 | 1 | 1 | [89, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101] ... |
| 91 | 1 | 1 | 0 | 1 | [91, 92, 94, 95, 96, 97, 98, 99, 100, 101, 102, 104] ... |
| 93 | 1 | 0 | 1 | 1 | [93, 95, 97, 98, 101, 102, 103, 105, 107, 111, 112, 114] ... |
| 95 | 1 | 1 | 1 | 1 | [95, 96, 97, 99, 100, 101, 103, 104, 105, 108, 109, 112] ... |
| 97 | 1 | 1 | 1 | 1 | [97, 98, 99, 100, 101, 102, 105, 106, 108, 110, 111, 112] .. |
| 99 | 1 | 0 | 1 | 1 | [99, 101, 102, 105, 111, 113, 114, 117, 118, 119, 121, 123]  |
| 101 | 1 | 0 | 0 | 0 | [101, 104, 106, 107, 110, 111, 112, 113, 114, 116, 117, 118] |
| 103 | 1 | 0 | 1 | 1 | [103, 105, 107, 108, 111, 112, 113, 114, 118, 120, 121, 122] |
| 105 | 1 | 0 | 0 | 1 | [105, 108, 109, 111, 114, 117, 119, 120, 121, 123, 125, 126] |
| 107 | 1 | 1 | 1 | 1 | [107, 108, 109, 110, 112, 113, 115, 116, 117, 118, 119, 121] |
| 109 | 1 | 1 | 1 | 1 | [109, 110, 111, 112, 113, 114, 115, 116, 117, 118, 120, 121] |
| 111 | 1 | 0 | 0 | 1 | [111, 114, 116, 117, 118, 119, 120, 121, 122, 123, 125, 127] |
| 113 | 1 | 0 | 1 | 0 | [113, 115, 116, 119, 123, 124, 126, 127, 129, 132, 134, 135] |
| 115 | 1 | 1 | 1 | 1 | [115, 116, 117, 118, 120, 121, 122, 123, 124, 125, 126, 128] |
| 117 | 1 | 1 | 0 | 1 | [117, 118, 120, 122, 123, 127, 128, 129, 130, 131, 132, 133] |
| 119 | 1 | 1 | 1 | 1 | [119, 120, 121, 122, 123, 124, 125, 127, 128, 129, 130, 131] |
| 121 | 1 | 0 | 0 | 1 | [121, 125, 126, 128, 130, 132, 133, 135, 138, 139, 143, 144] |
| 123 | 1 | 1 | 1 | 1 | [123, 124, 125, 126, 127, 128, 130, 131, 132, 133, 135, 137] |

- **Primes `p < 125` with no tower-forced `3p + 2` at any order: 3, 41, 43, 53, 67, 83, 89, 101, 103, 113** (ten of the
  twenty-eight; pre-registered P-F refuted on the count, confirmed on 41); the other eighteen primes have a strong realiser
  vertex with `2^k B_t = p + 1`, so `3p + 2` is forced by the tower over the padded realiser. Padding rescued no prime:
  the ten are the primes without a direct digit-1 vertex at order ≤ 9 as well.
- **Counts below 125 none of whose strong realisers has a digit-0 vertex: 11, 13, 47, 67** — exactly the near-join
  exceptions of HYP-9453 other than the holes and 1. (Corrected per MISTAKE-776: a padded sink is always a digit-0 vertex, so
  the statement is about strong realisers; the coincidence is explained in klein-2026-11-16-nearjoin: the near-join's new
  vertex `z` is a digit-0 vertex, so these four primes are never near-join values.)
- Counts without digit 2: 3, 15, 17, 27, 33, 35, 41, 43, 45, 57, 61, 63, 67, 83, 85, 87, 91, 101, 105, 111, 117, 121.
- **Counts with no rooted doubling child at any order: 11, 19, 29, 31, 35, 39, 43, 47, 53, 55, 61, 63, 67, 71, 73, 77,
  85, 87, 101, 113** (twenty; P-G refuted on the count). It contains THM-5272's containment failures 11, 29, 35 and the
  sterile forks 11, 29, 47, 101 but not 65, 83, 119; the rooted append is far more restrictive than containment.


## 4. Reading

Forcing `3p + 2` from an *arbitrary* realiser is impossible (Theorem A); forcing it from a *chosen* realiser is the
tower at a vertex with `2^k B_t = p + 1`, and §3 says exactly which primes below 125 have one. The `+2` difference at
a prime is thus a property of the prime's strong realisers' root responses — the same object (`B`) that carries the
near-join rows `h + 2^k B` (THM-4865, HYP-9453): the `+2` child and the doubling child are the `B`-direction and the
`k`-direction of one array, read at `B = h + 1` and at `2^k B = h ± 1`.
