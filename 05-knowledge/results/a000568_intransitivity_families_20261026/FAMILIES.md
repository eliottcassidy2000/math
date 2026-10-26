# A000568 enumerated from its families

klein, 2026-10-26. Every table below is exact. Sources: `scripts/tourn_enum.c` on nauty's
`gentourng` output (receipts `receipts/census-n3-n9.txt`, `receipts/census-n10.txt`,
`receipts/analyze-normal.txt`), the Pólya evaluations `receipts/species-normal.txt` and
`receipts/species10-normal.txt`. "Family" means a partition of the isomorphism classes by an
exact structural coordinate; each row sums to A000568(n).

## A. Root type: linear versus prime (the strong-component family)

`A000568(n) = L(n) + S(n)`, `L(n)` = classes with `>= 2` strong components (linear root),
`S(n)` = strong classes (prime root, or the vertex). `P(n)` = prime tournaments (A280772),
`SC(n)` = self-converse classes (A002785), merged metagraph `|G_n/Z_2| = (A000568 + SC)/2`,
`c(n)` = intransitivity structures (INTRANSITIVITY.md §4).

| n | A000568 | L(n) | S(n) | P(n) | SC(n) | merged | c(n) |
|---|---|---|---|---|---|---|---|
| 3 | 2 | 1 | 1 | 1 | 2 | 2 | 2 |
| 4 | 4 | 3 | 1 | 0 | 2 | 3 | 3 |
| 5 | 12 | 6 | 6 | 3 | 8 | 10 | 9 |
| 6 | 56 | 21 | 35 | 15 | 12 | 34 | 32 |
| 7 | 456 | 103 | 353 | 197 | 88 | 272 | 247 |
| 8 | 6880 | 872 | 6008 | 4008 | 176 | 3528 | 3328 |
| 9 | 191536 | 13403 | 178133 | 136655 | 2752 | 97144 | 93612 |
| 10 | 9733056 | 377107 | 9355949 | 7906161 | 8784 | 4870920 | 4774417 |

`L(n)` by number of strong components `k` (`= [x^n] S(x)^k`, ordered compositions; e.g.
`[x^10] S^2 = 2(178133 + 353 + 35) + 36 = 357078`):

| n | k=2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|
| 8 | 719 | 111 | 30 | 5 | 6 | – | 1 | | |
| 9 | 12098 | 1099 | 152 | 40 | 6 | 7 | – | 1 | |
| 10 | 357078 | 18273 | 1494 | 195 | 51 | 7 | 8 | – | 1 |

## B. The species equation made effective (prime cycle-index census)

Session56's identity (`tournament-hierarchy/REPORT.md`) is `T = SEQ(S)`, `S = X + P ∘ T₊`,
i.e. `A(x) = 1/(1 - S(x))` and `S(x) = x + Σ_{Q prime} Z_{Aut Q}(A(x) - 1, A(x²) - 1, ...)`.
It was left there as an identity ("plain prime counts alone do not provide the necessary
cycle-index data"). The data are now supplied: every prime quotient of every class with
`n <= 9` was canonicalised, its automorphism group enumerated from nauty generators, and the
groups fall into **23 cycle-index types** (`data/prime_cycle_types_n9.tsv`; columns: order,
self-converse flag, |Aut|, cycle types of Aut with multiplicities, cycle types of the
anti-automorphisms, number of prime classes). With the 7,906,161 primes of order 10 (each
contributes `1` at `n = 10`), the two equations reproduce `A000568(n)` and `S(n)` **exactly for
n = 1..10** (`receipts/species-normal.txt`, `receipts/species10-normal.txt`). Primes by order and
symmetry (`|Aut|`: count), self-converse in brackets:

| n | primes | |Aut| = 1 | 3 | 5 | 7 | 9 | 21 |
|---|---|---|---|---|---|---|---|
| 3 | 1 [1] | | 1 [1] | | | | |
| 5 | 3 [3] | 2 [2] | | 1 [1] | | | |
| 6 | 15 [3] | 14 [2] | 1 [1] | | | | |
| 7 | 197 [55] | 192 [50] | 3 [3] | | 1 [1] | | 1 [1] |
| 8 | 4008 [72] | 3996 [72] | 12 [0] | | | | |
| 9 | 136655 [2235] | 136553 [2205] | 100 [28] | | | 2 [2] | |
| 10 | 7906161 [5525] | (cycle types not collected at n = 10) | | | | | |

A280772 ("number of simple tournaments") lists `1, 1, 1, 0, 3, 15, 197, 4008` (n = 1..8, Dammak
2016, Sagemath); this census extends it with `a(9) = 136655`, `a(10) = 7906161`. A003505 carries
the same name with `197 -> 203`, `4008 -> 3785` (Müller–Nešetřil–Pelant 1975, Table 1), and
N. J. A. Sloane asked in 2017 which is right: the independent nauty census here confirms A280772
at `n = 7, 8` (and THM-1960's `1, 1, 1, 0, 3, 15` for `n <= 6`), so A003505 is erroneous at
`n = 7` and `n = 8`.

## C. Tree shapes (number of linear nodes, number of prime nodes)

Classes by `(#linear nodes, #prime nodes)` of the modular tree; `(1,0)` is the transitive
tournament, `(0,1)` the primes.

| n | (0,1) | (0,2) | (0,3) | (0,4) | (1,0) | (1,1) | (1,2) | (1,3) | (1,4) | (2,1) | (2,2) | (2,3) | (3,1) | (3,2) | (3,3) | (4,1) | (4,2) | (5,1) | (5,2) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 5 | 3 | 1 | | | 1 | 4 | | | | 3 | | | | | | | | | |
| 6 | 15 | | | | 1 | 22 | 8 | | | 7 | | | 3 | | | | | | |
| 7 | 197 | 14 | 2 | | 1 | 142 | 12 | | | 56 | 22 | | 8 | | | 2 | | | |
| 8 | 4008 | 101 | | | 1 | 1910 | 151 | 23 | | 504 | 58 | | 84 | 35 | | 5 | | | |
| 9 | 136655 | 1585 | 53 | 5 | 1 | 42158 | 1263 | 36 | | 7817 | 651 | 106 | 972 | 129 | | 66 | 39 | | |
| 10 | 7906161 | 36447 | 432 | | 1 | 1549115 | 19792 | 744 | 69 | 192623 | 6432 | 302 | 17785 | 1590 | 278 | 1056 | 172 | 25 | 32 |

## D. Information families: the fibre spectrum of the intransitivity observable

Classes by `N(T) = 2^{p} ∏ k_N!` (INTRANSITIVITY.md §1), i.e. by how many labelled relations
share the class's intransitivity structure.

- `n = 7`: 2:197, 4:130, 8:46, 12:20, 16:22, 24:12, 32:2, 48:8, 72:4, 96:8, 240:6, 5040:1.
- `n = 8`: 2:4008, 4:1850, 8:521, 12:131, 16:89, 24:99, 32:35, 48:81, 96:22, 144:9, 192:9, 288:9, 480:9, 1440:7, 40320:1.
- `n = 9`: 2:136655, 4:41617, 8:7886, 12:1946, 16:1263, 24:1034, 32:161, 48:402, 64:39, 72:55, 96:239, 144:37, 192:74, 240:26, 288:5, 384:6, 432:4, 480:23, 576:20, 960:10, 1152:5, 1440:10, 2880:10, 10080:8, 362880:1.
- `n = 10`: 2:7906161, 4:1539180, 8:194686, 12:44040, 16:20483, 24:16217, 32:2134, 48:5844, 64:303, 72:473, 96:1710, 128:32, 144:242, 192:468, 240:161, 288:279, 384:158, 480:121, 576:92, 864:5, 960:92, 1152:11, 1440:29, 1728:11, 1920:7, 2304:11, 2880:52, 5760:22, 8640:11, 20160:11, 80640:9, 3628800:1.

The values `N` are the orders of the groups `∏ S_{k_N} × (Z/2)^p` that can occur; `12 = 3!·2`
is the first that is neither a power of 2 nor a factorial.

## E. Levels (`c3`, the atlas's axis), with information and primeness per level

Level widths (classes per `c3`), extending the atlas (`n <= 7`) and opus-S316 (`n = 8`):

- `n = 8`: 1, 6, 11, 14, 29, 44, 84, 108, 182, 220, 406, 396, 614, 628, 892, 740, 929, 622, 607, 262, 85 (`c3 = 0..20`).
- `n = 9`: 1, 7, 16, 23, 46, 83, 149, 234, 405, 569, 1042, 1307, 2200, 2993, 4339, 5380, 8169, 8679, 12526, 13801, 16868, 17119, 21413, 17482, 19089, 15181, 11988, 6645, 3310, 457, 15 (`c3 = 0..30`).
- `n = 10`: 1, 8, 22, 36, 71, 140, 262, 440, 793, 1240, 2264, 3220, 5672, 8080, 13416, 18232, 30218, 38020, 61172, 76420, 116350, 137784, 206344, 230868, 333066, 370112, 490744, 520012, 690172, 665864, 837084, 785312, 874245, 757500, 794222, 582472, 507631, 306180, 190310, 63724, 13333 (`c3 = 0..40`; peak 874245 at `c3 = 32`).

Mean `log₂ N(T)` per level and fraction of prime classes per level (`n = 9`): the mean falls
monotonically `18.5, 13.3, 9.9, 8.2, 7.4, 5.9, 5.3, 4.4, 3.9, 3.3, 3.0, 2.6, 2.4, 2.1, 2.0, 1.8,
1.7, 1.5, 1.5, 1.35, 1.31, 1.26, 1.26, 1.21, 1.23, 1.19, 1.19, 1.18, 1.19, 1.14, 1.20`; the prime
fraction is `0` for `c3 <= 6` and then `.02, .04, .08, .09, .18, .21, .30, .33, .44, .48, .59,
.63, .72, .75, .79, .79, .82, .81, .84, .84, .85, .84, .88, .93`. Primes need `c3 >= n - 2`
(Moon: a strong tournament has at least `n - 2` cyclic triples); the bound is attained by primes
for `5 <= n <= 10`, by `1, 1, 2, 3, 5, 8` prime classes, and by `2, 4, 8, 16, 32, 64` strong
classes in all. **The chain family.** Take the transitive tournament on `1 < ... < n` and reverse
the arcs of a chain `1 = a_0 < a_1 < ... < a_r = n`, `S = {a_1, ..., a_{r-1}} ⊆ {3, ..., n-2}`.
Each such tournament is strong with exactly `n - 2` cyclic triples (`Σ_i (a_{i+1} - a_i - 1)`
triples `{a_i, x, a_{i+1}}` plus the `r - 1` consecutive chain triples; chain elements `2` or
`n - 1` only relabel). FINITE-EXACT (`n = 5..9`, by canonical forms; `n = 10` by the counts):
the `2^{n-4}` chain tournaments are pairwise non-isomorphic and are **all** the strong classes
with `n - 2` cyclic triples, and a chain tournament is prime iff every gap `a_{i+1} - a_i` is at
most 2 (a gap `>= 3` makes the intermediate vertices a module; `3 ∉ S` or `n - 2 ∉ S` likewise),
so the minimum-intransitivity primes are counted by compositions of `n - 5` into parts 1 and 2:
`F_{n-4}` (Fibonacci) — HYP-9413.

## F. Slater layers (minimum feedback arc set), extending `defect_boundary_20261020`

- `n = 8`: 1, 21, 151, 584, 1429, 2141, 1797, 690, 66 (`fas = 0..8`).
- `n = 9`: 1, 28, 295, 1765, 7147, 20305, 40167, 53436, 44432, 19986, 3748, 218, 8 (`fas = 0..12`).

Checks of that packet's closed forms: `a_{n,<=1} = 1 + C(n-1,2)` gives 22, 29;
`a_{n,<=2} = (n⁴ - 6n³ - n² + 62n - 72)/8` gives 173 = 1+21+151 and 324 = 1+28+295. Maximum
Slater index 8 and 12 (A003141: 1, 1, 3, 4, 6, 8, 12 for n = 3..9 is consistent with
`diameter-is-feedback-arc-set.md`'s reading).

## G. Symmetry strata

`|Aut|` histograms: `n = 10`: 1: 9677605, 3: 54370, 5: 274, 7: 10, 9: 738, 15: 32, 21: 12,
25: 1, 27: 10, 63: 2, 81: 2 (new; `n <= 9` agree with the atlas and the stored TSVs). The
order-25 class is the stack of two copies of the circulant `Z_5{1,2}`; 63 = 21·3 and 81 = 3⁴
are stacks with the Paley-7 and `C_3[C_3]` factors (THM-2450's product law).

## H. What this does and does not do

It enumerates A000568 through `n = 10` from (A) the strong-component composition, (B) the prime
cycle-index census, (C) tree shapes, (D) intransitivity fibres, (E) levels, (F) Slater layers,
(G) symmetry, each exactly. It does not give a closed form for any family beyond the species
equations, and it does not reduce the Burnside count (THM-514/516, the `04-computation` engines
to `n = 200`), which remains the only way to reach large `n`. The prime cycle-index census at
`n = 10` (needed for `n = 11` predictions) was not collected; `run.sh` documents how.
