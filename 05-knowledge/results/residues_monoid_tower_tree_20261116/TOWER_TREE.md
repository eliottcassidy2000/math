# The mod-3 tower, the period-12 classes, the accelerated chains, and the pruned Collatz tree with edge lengths 1, 2, 3, 4

Session klein-2026-11-16-monoid. Script `scripts/collatz_tower_tree.py`; receipts `collatz_tower_tree.out`,
`tree_mod27.out`. The owner's description is taken as a sequence of claims and each is made exact.

## 1. The tower (PROVED, elementary; receipt)

Rows: `row(3j) = 2j` (a short row holding one multiple of 3) and `row(3j + 1) = row(3j + 2) = 2j + 1` (a long row
with `3j + 1` on the left and `3j + 2` on the right). Period 12 splits each block of twelve integers into three classes
of four: the multiples of 3 (`0, 3, 6, 9`: two odd, two even), the even non-multiples (`2, 4, 8, 10`: two multiples
of 4, two doubles of odds) and the twin-prime candidates `6k ± 1` (`1, 5, 7, 11`). For the accelerated map
`T(n) = n/2` (`n` even), `(3n + 1)/2` (`n` odd):

- **Halving arrows between short rows.** The multiples of 3 halve to multiples of 3: the chain `3 ← 6 ← 12 ← 24 ← 48`
  sits at rows `2, 4, 8, 16, 32`, moving up 2, 4, 8, 16 rows at a time; the odd multiples of 3 (rows 2, 6, 10, 14, …,
  every 4 rows) each start such a chain.
- **The shoelace.** For `m` not divisible by 3 the halving arrow `2m → m` spans `row(2m) − row(m) = 2(⌊(2m − 1)/3⌋ −
  ⌊(m − 1)/3⌋)` rows: `0, 2, 2, 4, 4, 6, 6, 8, 8, …` for `m = 1, 2, 4, 5, 7, 8, 10, 11, 13, …` (the owner's 3, 5, 5, 9,
  9 count the rows inclusively), growing by exactly 2 when `m` grows by 3, so consecutive arrows overlap like laces:
  `8 → 4` over rows 3–5, `10 → 5` over 3–7, `14 → 7` over 5–9, `22 → 11` over 7–15.
- **Three families of odd arrows.** By `n (mod 6)`: left (`n ≡ 1`), short (`n ≡ 3`), right (`n ≡ 5`). The first
  members are `1 → 2` (row 1 to row 1), `3 → 5` (row 2 to row 3), `5 → 8` (row 3 to row 5); within each family the
  next member (`n + 6`) starts 4 rows higher and is 2 rows longer (`T(n + 6) = T(n) + 9`, three rows more at the head,
  two at the foot). All checked to `10^5`; pre-registered P-8 **confirmed**.

## 2. The pruned accelerated tree (THM-5259, PROVED; receipt)

In the inverse tree of `T` every `n` has the child `2n`, and `n ≡ 2 (mod 3)` also the odd child `(2n − 1)/3`. Along
a ray `n, 2n, 4n, …` the residues alternate `2, 1, 2, 1` (for `3 ∤ n`), so forks (nodes `≡ 2 (mod 3)`) occur every
second node; multiples of 3 are sterile (their rays never fork). **Tree lengths keep the two unary states and split
the fork state in three**: for a fork `n`, the odd child `k = (2n − 1)/3` is `≡ 0, 1, 2 (mod 3)` according to
`n ≡ 5, 2, 8 (mod 9)` — sterile, unary-then-fork, fork.

**Theorem (THM-5259).** Prune the sterile rays and contract the unary nodes. Then every remaining node has exactly
two children (the ray continuation and the branch), every edge length lies in `{1, 2, 3, 4}`, ray edges have length
2 or 4, branch edges 1, 2, 3 or 4, and the lengths are determined by `n (mod 27)`:

| fork `n (mod 27)` | ray | branch |
|---|---|---|
| 2, 20 | 2 | 2 |
| 11 | 2 | 4 |
| 8 | 4 | 3 |
| 17, 26 | 4 | 1 |
| 5, 14, 23 | (sterile branch: the fork is contracted) | |

*Proof.* The ray from `n` reaches the fork `4n` after 2 edges; `4n ≡ 5 (mod 9)` (a sterile fork, contracted) iff
`n ≡ 8 (mod 9)`, and then `16n ≡ 2 (mod 9)` is fertile, so a ray edge has length 2 or 4 and **no two consecutive forks
on a ray are sterile** (`n ≡ 5 ⟹ 4n ≡ 2 (mod 9)`). The branch reaches `k` (a fork, length 1) when `k ≡ 2 (mod 3)`,
i.e. `n ≡ 8 (mod 9)`, or `2k` (length 2) when `k ≡ 1 (mod 3)`, i.e. `n ≡ 2 (mod 9)`; if that first fork is sterile
(`k ≡ 5` or `2k ≡ 5 (mod 9)`, conditions modulo 27) the next fork is 2 edges further, giving 3 or 4, and never more by
the same no-two-sterile fact. With both the ray and the branch fertile after contraction every node is binary. ∎

**Reconstruction and symmetry.** The child labels are explicit functions of the parent label and the edge type and
length (`4n`, `16n`, `k`, `2k`, `4k`, `8k`), so the lengths together with the left/right structure recover every
label and every pruned ray (receipt: all `2^19 − 1` labels to depth 18 distinct); without the lengths the tree is
the complete binary tree — "that rigid tree becomes symmetric". Pre-registered P-7 **confirmed** (depth 18: 524,287
nodes, length counts ray 2: 145,346, ray 4: 116,797, branch 1: 77,955, 2: 96,974, 3: 38,842, 4: 48,372).

## 3. The two pictures are one

The rays are the rows `L_c` of COLLATZ_ROWS.md (the doubling chains); a fork every second node is the statement that
the Syracuse map lands on a row at height `v_2(3c + 1)`; the contracted lengths `{1, 2, 3, 4}` are the residues
modulo 27 of the row index. In the tower the forks of a row are its right-hand entries (`2 mod 3`), the sterile rows
are the short rows, and the three odd-arrow families are the three residues of `n` modulo 6.
