# The aliquot table in the first handfuls, the untouchable primes {2, 5}, the forbidden path counts {7, 21}, the powers-and-doublings chain structure, and abc on the row of 1

**Status: PROVED (elementary; THM-4805 Theorems A-D, proofs below) / FINITE-EXACT (printed outputs) / CITED (S57 untouchables, THM-1370, THM-4790, THM-126, S17/S18, the guarded prime-extension formula, THM-3833) / ANALOGY (the crevasse dictionary between untouchables and forbidden path counts) / NUMEROLOGY (the coincidences 2 + 5 = 7, 4 x 5 + 1 = 21, 63 = 7 x 9 = 21 x 3, recorded so they are not regenerated) / HYP-9415 (pre-registered tests, outcomes in section 8). Session opus-2026-10-26-S60, worktree `/tmp/math-wt-s60`.**

Owner's directive (verbatim): "merge in past work regarding the 1-2-4 decomposition of 7 and the idea that {2,5} are the only untouchable numbers that are primes, and how they relate with {7,21} being the tournament forbidden hamiltonian path values. consider prior work on the idea of the natural numbers being arranged with their powers extending below them in one chain and then continual doubling in the other chain, and notice which natural numbers or their chain values have an already existing occurence in the structure up to that natural number. think of the abc conjecture and compare this with the function that associates 0 with 1, 1 with the primes, 2 with nothing, 4 with 3, 9 with 3, nothing with 5, 25 and 6 with 6, 8 with 7, 49 and 10 with 8, 15 with 9, 14 with 10 and pay particular attention to microcosm/macorocosm patterns within the first few handfuls. note in each, 6 is a first instance, in the former, it is the first instnce whose square is assocaited with a later doubling, 9's. and in the latter it gets both 25 and 6, getting itself making it extra rare."

The function is the aliquot sum `s(n) = sigma(n) - n`, listed as *source -> value* (the owner's list switches orientation after "1 with the primes"; "9 with 3" is read as `s(9) = 4`, since `s(9) = 1 + 3`). Scripts: [`aliquot_chains_abc_20261026.py`](../../04-computation/experiments/aliquot_chains_abc_20261026.py) (+ [`.out`](../../04-computation/experiments/aliquot_chains_abc_20261026.out), 26 s) and the pre-registered [`aliquot_sources_preregistered_20261026.cpp`](../../04-computation/experiments/aliquot_sources_preregistered_20261026.cpp) (+ [`.out`](../../04-computation/experiments/aliquot_sources_preregistered_20261026.out)).

## 0. Inheritance (read before drafting; nothing below re-proves these)

- **Untouchables (S57, [`s57_untouchables_root27_20261026`](s57_untouchables_root27_20261026/README.md)):** exact census to `10^10` (A057978: 2, 5, 89, ..., 1,681,871,718 untouchables below `10^k`), the Pollack-Pomerance constant `0.17865689`, the parity lemma (sources of an even value `k` are even `m <= 2k - 2` or odd squares), section 3: an odd `k >= 9` is `s(pq)` whenever `k - 1` is a sum of two distinct primes, so a second odd untouchable is a distinct-Goldbach exception beyond `4 x 10^18`, the only odd `k` whose `k - 1` has no distinct-prime partition are 3, 5, 7, and below `10^8` exactly five odd numbers have a single source (`3 = s(4)`, `7 = s(8)`, `9 = s(15)`, `11 = s(21)`, `39 = s(217)`); section 4: the `F_1` dictionary (`s^-1(0) = {1}` and `s^-1(1) = primes` agree with `tau - 1`; "2 with nothing" is the first place where prime weights matter; `s(p^a) = |P^(a-1)(F_p)|`, so `s(4) = 3`, `s(9) = 4`, `s(25) = 6`, `s(49) = 8` are projective lines and `s(8) = 7` is the Fano plane); section 5: the fibres 7 and 21 of the graph-pair chart (THM-4790).
- **Classical, through S57:** Erdos 1973 (positive density of untouchables), Pollack-Pomerance 2016, Chen-Zhao and Luca-Pomerance bounds; 5 is the only odd untouchable under the distinct-prime Goldbach conjecture (so `{2, 5}` are the only untouchable primes under it: 2 unconditionally, 5 unconditionally, no prime `p >= 7` since `7 = s(8)` and `p - 1` is a distinct-prime sum for `p >= 11`).
- **Aliquot formulas:** the guarded prime-extension formula `s(Cp) = (sigma(C) - C) p + sigma(C)` for a prime `p` not dividing `C` ([`aliquot_guarded_cycles_20261014e`](aliquot_guarded_cycles_20261014e.md)), the dyadic formula `s(2^a m) = (2^(a+1) - 1) sigma(m) - 2^a m` ([Lehmer five note](collatz_aliquot_lehmer_five_20260927.md)), and the perfect-core bijection `s(6p) = 6(p + 2)` for primes `p >= 5`, which makes twin primes actual aliquot edges `6p -> 6(p+2)` ([perfect-core twin edges](perfect_core_twin_edges_20261019/PERFECT_CORE_TWIN_EDGES.md)). The twin collision of Theorem B is a different twin-prime statement (two sources of one value, not an edge), and is not in those notes.
- **{7, 21}:** [THM-1370](../../01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md): no tournament has 7 or 21 Hamiltonian paths (multiplicativity over strong components, the insertion lemma, Moon's subtournament theorem, the Busch floors `f(n) = 3, 5, 9, 15, 25, 45, 75`), with the 2016 MathOverflow discussion as prior art and completeness of the spectrum OPEN; [THM-4790](../../01-canon/theorems/THM-4790-presentation-count-fibres-7-and-21-are-absent-for-every-tournament-on-at-most-11-vertices.md) (presentation counts `tc = H/|Aut|` never 7 or 21 for `n <= 11`); THM-4741 and HYP-9352 (2).
- **{1, 2, 4} and F_21:** [THM-126](../../01-canon/theorems/THM-126-paley-circulant-uniquely-maximizes-z7.md) (`P_7` with `S = {1, 2, 4}` is the unique `H`-maximizer at `n = 7`, `H = 189`); [S18 section 3](mersenne_switch_parity_f21_compression_20261006.md) (`x -> 2x + 1` and `x -> 4x + 1` generate `F_21 = Aut(P_7)` only at the Mersenne prime 7, the equality case being Catalan's `2^3 + 1 = 3^2`); [S17](ramanujan_heegner7_trivial_cycle_20261006.md) (Ramanujan's split of the primes into `{1, 2, 4}` and `{3, 5, 6}` mod 7; the Heegner point of discriminant -7).
- **abc:** [THM-3833](../../01-canon/theorems/THM-3833-abc-conditional-cube-radical-growth-and-hyperbolic-power-finiteness.md) (ABC-conditional radical floors; ABC itself OPEN). The squares/doubles transport ([S9 foundry](collatz_sqdbl_20260925_squares_doubles_foundry.md): "multiplication relates to the squares the way addition relates to the doubles") is the owner's earlier squares-and-doubles seed; it transports Collatz, and does not contain the two-chain arrangement of section 3.
- **Prior statement of the chain arrangement:** a sweep of `05-knowledge`, `07-reflections`, `01-canon`, `inbox` and `.scratch` for the two-chain arrangement with first-occurrence bookkeeping is reported in section 9 (inheritance check).

## 1. The aliquot table in the first handfuls (THM-4805, Theorem A)

Write `t(n) = tau(n) - 1` for the number of proper divisors. Exact table of sources `s^-1(v)` (a composite `n` has `s(n) >= 1 + sqrt(n)`, since `n/p` is a proper divisor for the least prime factor `p`; so every source of `v` is below `v^2`, and the table is exact):

| `v` | sources | shape of each source |
|---|---|---|
| 0 | 1 | - |
| 1 | the primes | `p` |
| 2 | none | untouchable |
| 3 | 4 | `2^2` (column of 2) |
| 4 | 9 | `3^2` (column of 3) |
| 5 | none | untouchable |
| 6 | 6, 25 | `2 x 3` (row of 3), `5^2` |
| 7 | 8 | `2^3` (row of 1 = column of 2) |
| 8 | 10, 49 | `2 x 5` (row of 5), `7^2` |
| 9 | 15 | `3 x 5` (off-grid: the first composite new number) |
| 10 | 14 | `2 x 7` (row of 7) |
| 11 | 21 | `3 x 7` (off-grid) |
| 12 | 121 | `11^2` |
| 13 | 27, 35 | `3^3` (column of 3), `5 x 7` |
| 14 | 22, 169 | `2 x 11`, `13^2` |
| 15 | 16, 33 | `2^4` (row of 1), `3 x 11` |
| 16 | 12, 26 | `4 x 3` (row of 3), `2 x 13` |
| 17 | 39, 55 | `3 x 13`, `5 x 11` |
| 18 | 289 | `17^2` |
| 19 | 65, 77 | `5 x 13`, `7 x 11` |
| 20 | 34, 361 | `2 x 17`, `19^2` |
| 21 | 18, 51, 91 | `2 x 3^2` (row of 9), `3 x 17`, `7 x 13` |
| 22 | 20, 38 | `4 x 5` (row of 5), `2 x 19` |
| 23 | 57, 85 | `3 x 19`, `5 x 17` |
| 24 | 529 | `23^2` |
| 25 | 95, 119, 143 | `5 x 19`, `7 x 17`, `11 x 13` |
| 28 | 28 | itself (perfect) |

Rows `v = 26..60` are in the output (section 1). Untouchable values below 1000, exact (every `n < 10^6` scanned): 89 of them, `2, 5, 52, 88, 96, 120, 124, 146, 162, 188, 206, 210, ...`; the only primes among them are 2 and 5, and 5 is the only odd one (S57's census extends this to `10^10`).

**Lemma 1 (shapes by number of proper divisors).** `t = 2`: `n = p^2`, `s = p + 1`. `t = 3`: `n = p^3` (`s = 1 + p + p^2`) or `n = pq` with `p < q` (`s = 1 + p + q`). `t = 4`: `n = p^4`, `s >= 15`. `t = 5`: `n = p^5` (`s >= 31`) or `n = p^2 q` (`s = 1 + p + p^2 + q + pq >= 16`). `t >= 6`: `s >= 22` (six distinct proper divisors sum to at least `1 + ... + 6 = 21`, and equality would force `60 | n`, which has eleven proper divisors). The minimal source with `t = 4` is 16 (`s = 15`) and with `t = 5` it is 12 (`s = 16`).

**Theorem A.** (i) `s^-1(2) = s^-1(5) = {}` (classical). (ii) `s^-1(7) = {8}`: the value 7 has the single source `8 = 2^3`, so `7 = 1 + 2 + 4` is forced. (iii) `s^-1(21) = {18, 51, 91}`. (iv) Every value `v <= 14` has only sources of the three shapes `p^2`, `pq`, `p^3`; the first value with a source of four proper divisors is `15 = s(16)` and the first with five is `16 = s(12)`. (v) The Mersenne values `2^k - 1 = s(2^k)` have the single source `2^k` exactly for `k = 2, 3` (`3 <- 4`, `7 <- 8`); `15 <- {16, 33}`, `31 <- {32, 125, 161, 209, 221}`, `63 <- {64, 177, 817}`, `127` has eleven sources.

*Proof.* By Lemma 1, a value `v <= 14` can only come from `p^2` (`v = p + 1`), `pq` (`v - 1 = p + q`, distinct primes) or `p^3` (`v = 1 + p + p^2`), which gives (iv). (i) `v = 2`: every composite has `s >= 3`. `v = 5`: `p + 1 = 5` needs `p = 4`; `p + q = 4` needs `2 + 2`; `1 + p + p^2 = 5` has no solution. (ii) `v = 7`: `p = 6` fails; `p + q = 6` needs `3 + 3`; `1 + p + p^2 = 7` gives `p = 2`, `n = 8`. (iii) `v = 21`: `p = 20` fails; `p + q = 20` gives `{3, 17}` and `{7, 13}`, i.e. 51 and 91; `p^2 + p = 20` fails; `t = 4`: `s(16) = 15`, `s(81) = 40`; `t = 5`: `p^5` gives at least 31, and `s(p^2 q) = (1 + p + p^2) + q(1 + p) = 21` has the solution `p = 3`, `q = 2` (`n = 18`) and no other (`p = 2` gives `3q = 14`; `p >= 5` exceeds 21); `t >= 6` gives at least 22. (v) `s(2^k) = 2^k - 1` is the dyadic formula; the other sources are read off the table (exact, `n < v^2`), and for `k >= 4` the distinct-Goldbach partitions of `2^k - 2` (`14 = 3 + 11`, `30 = 7 + 23`, `62 = 3 + 59`) supply a second source; the five single-source odd values below `10^8` (S57) contain no further Mersenne number. QED

**The first handful (values 0-5) uses only columns:** `1`, the primes, nothing, `2^2`, `3^2`, nothing. The second handful (6-10) brings the first row source, `6 = 2 x 3` on the row of 3, the one row-column coincidence `8` (the row of 1 is the column of 2), and the first off-grid source `15 = 3 x 5`, which is also the first composite *new* number of section 3. Values 11-14 are still in the `t <= 3` regime; the regime changes at 15/16 with the row of 1 (`16`) and the row of 3 (`12`) delivering the first four- and five-divisor sources.

## 2. The twin collision, and why 6 is first in the aliquot table (Theorem B)

Column and row images (CITED formulas): `s(p^2) = p + 1`, `s(2q) = q + 3`, `s(2^a q) = (2^a - 1) q + 2^(a+1) - 1` (the guarded formula with `C = 2^a`), `s(p^k) = (p^k - 1)/(p - 1)`.

**Theorem B.** (i) If `p - 2` and `p` are primes (`p >= 5`), then `s(p^2) = s(2(p - 2)) = p + 1`: the square column of `p` and the doubling row of `p - 2` collide. The collision values are `6, 8, 14, 20, 32, 44, 62, 74, 104, 110, ...`. (ii) For a prime `p`, the sources of `p + 1` of the row form `2^a q` (`q` an odd prime, `a >= 1`) are exactly the twin source `2(p - 2)`: `(2^a - 1) q = p - 2(2^a - 1)` forces `2^a - 1 | p`, so `a = 1`, `q = p - 2`. (iii) The collision value is itself a source only at `p = 5`: `s(2q) = 2q` iff `q = 3`. So `6 = s(6) = s(25)` is the unique twin collision at a fixed point of `s`, which is the owner's "getting itself". (iv) 6 is the only even perfect number with a prime-square source: an even perfect `N = 2^(p-1)(2^p - 1)` with `p` odd is `1 mod 3`, so `N - 1` is a multiple of 3 larger than 3 and not prime. (v) FINITE-EXACT: `s^-1(6) = {6, 25}`, `s^-1(28) = {28}` (28 is its own only source), `s^-1(496) = {496, 652}` with `652 = 4 x 163`, `s^-1(8128) = {8128, 10682, 10828}` with `10828 = 4 x 2707` and `10682 = 2 x 7^2 x 109`; the row source `4q` of an even perfect `N` requires `q = (N - 7)/3` prime, which holds for `N = 28` (where `4 x 7` is `N` itself), 496 and 8128 and fails for the next four perfect numbers. Of the 205 twin collision values `p + 1 <= 10^4`, 115 have exactly the two twin sources (`74 <- {70, 142, 5329}` is the first with a third); in `(10^3, 10^4]` the source counts run from 2 to 6 with mean 2.68.

*Proof of (i)-(iv).* (i) and (ii) are the displayed formulas; in (ii), `2^a - 1 = p` would give `q = -p`. (iii) `q + 3 = 2q`. (iv) `2^(p-1) = 1 mod 3` for `p - 1` even and `2^p - 1 = 1 mod 3` for `p` odd. (v) is the parity lemma plus the printed scans (section 9 of the output). QED

So in the aliquot table 6 is first three times over: the first value with a row source, the first value with two sources, and the one collision that is a fixed point. The pattern it starts is the twin-prime pattern, which repeats at every scale (`8 <- {10, 49}`, `14 <- {22, 169}`, `20 <- {34, 361}`, ...): the microcosm is the pair `(2 x 3, 5^2)`, the macrocosm the pairs `(2(p - 2), p^2)`. The perfect numbers repeat the "getting itself" at a larger scale with fewer companions: 28 gets only itself.

## 3. The chain structure: powers below, doublings across (Theorem C)

For each natural number `n` write its **column** `P(n) = {n, n^2, n^3, ...}` and its **row** `D(n) = {n, 2n, 4n, ...}`. Build the structure by adding `1, 2, 3, ...` in order, each with its two rays. A number is **old** when it already lies on a ray of a smaller number, **new** otherwise; a ray value of `n` "has an already existing occurrence" when it lies on a ray of some `m < n`, and is a **later** occurrence of `m > n` when it lies on a ray of `m` that is not yet built.

**Theorem C.** Write `n = 2^a r` with `r` odd.
1. `n` is old iff `n` is even or a perfect power. The rays of a new number meet no ray of a smaller number: the structure is a disjoint union of fresh ray pairs, one per odd non-perfect-power, together with the old numbers.
2. `n^2` lies on the row `D(r^2)` at position `2a`. It already occurs (on a ray of a smaller number) iff `r^2 < n` or `n` is a perfect power; it is a later doubling of `r^2 > n` iff `a >= 1` and `r > 2^a`. The first such `n` is 6 (`36 = 4 x 9`, and 9 is added after 6); the instances up to 40 are `6, 10, 14, 18, 20, 22, 26, 28, 30, 34, 36, 38`, and the count of instances up to `X` is `X/2 - O(sqrt X)`.
3. The same holds for every power: `n^k = 2^(ak) r^k` lies on the row of `r^k`, a later number iff `r^(k-1) > 2^a`; for `k >= 2` the first instance is again 6 (`216 = 8 x 27`). The whole column of an even number `2^a r` is the column of its odd part shifted along the rows: it is the diagonal of the grid spanned by the row of 1 and the column of `r`.
4. The row of an odd `n` meets an earlier column iff `n` is a perfect power: `9` meets `36 = 6^2`; `25` meets `100 = 10^2` and `400 = 20^2`; `27` meets `216 = 6^3`, `1728 = 12^3`, `13824 = 24^3`.
5. The row of 1 is the only row that is also a column (`D(1) = P(2)`), because 2 is the only prime on a row. Its proper initial segment `1, 2, 4` is the divisor set of 8 and the 1-2-4 decomposition of `7 = s(8)` (Theorem A (ii)).

*Proof.* (1) Even `n` lies on `D(r)` with `r < n` (or on `D(1)` when `r = 1`); a perfect power `b^k` lies on `P(b)`. Conversely let `n` be odd and not a perfect power, and `m < n`. `D(n)` meets `D(m)` only if the odd parts agree, forcing `n = m`. `D(n)` meets `P(m)` only if `2^j n = m^k`; `n` odd gives `n = (odd part of m)^k`, so `k = 1` and `m = 2^j n >= n`. `P(n)` meets `D(m)` only if `n^k = 2^i m` with `i = 0`, `m = n^k >= n`. `P(n)` meets `P(m)` only if `n^k = m^l`, so `n` and `m` are powers of a common base; `n` is not a perfect power, so `n` is the base and `m = n^l >= n`. (2) `n^2 = 2^(2a) r^2`. A ray `D(m)` containing `n^2` has `m = 2^i r^2` for some `i <= 2a`, so `m < n` is possible iff `r^2 < n`; a column `P(m)` with `m < n` contains `n^2` iff `n` is a perfect power. `r^2 > n = 2^a r` iff `r > 2^a`; for `a = 0` the "doubling" is trivial, so `a >= 1`. The even `n <= X` failing `r > 2^a` number `sum_(a>=1) 2^(a-1) = O(sqrt X)` over `2^a <= sqrt X`. (3) Replace `r^2` by `r^k`. (4) `2^j n = m^k` with `n` odd forces `n` to be a `k`-th power, and conversely `n = b^k` gives `2^(jk) n = (2^j b)^k` with `2^j b < n` for small `j`. (5) `D(1) = P(2)`; a row `D(m)` with `m > 1` contains the prime 2 only for `m = 1`. QED

**Microcosm and macrocosm.** The first handful shows each coincidence once: 4 is the first number on two rays (`2 x 2 = 2^2`, the row of 1 meeting the column of 2 at their second entries), 6 is the first even number off the row of 1, 8 the first cube, 9 the first odd old number (the first odd square, whose row then meets the column of 6 at 36), 16 the first number on three rays. The macrocosm is Theorem C (3): every even number's column is a shifted copy of its odd part's column, so the first instance at 6 repeats for every `2^a r` with `r > 2^a`, and the coincidence `36 = 6^2 = 4 x 9` is the template for `n^2 = 2^(2a) r^2`.

## 4. The aliquot function on the chain grid

The grid of the chain structure is `{2^j p^k}` (rows of prime columns). The aliquot function is explicit there: columns go to repunits, `s(p^k) = 1 + p + ... + p^(k-1)` (the hyperplane at infinity, S57 section 4), rows go to affine functions of the prime, `s(2^a p) = (2^a - 1) p + 2^(a+1) - 1`, and the row of 1 goes to the Mersenne numbers `2^k - 1`. Off the grid the first shape is `pq -> p + q + 1`, the Goldbach shape of S57 section 3. The owner's list reads, entry by entry: `3 <- 4` and `4 <- 9` are the second column entries of 2 and 3; `6 <- 25`, `8 <- 49`, `12 <- 121` continue the square column; `6 <- 6`, `8 <- 10`, `10 <- 14`, `14 <- 22` are the first row entries of 3, 5, 7, 11 (`q -> q + 3`); `7 <- 8` is the row of 1; `9 <- 15` and `11 <- 21` are off-grid. The two sequences `p + 1` (squares) and `q + 3` (doublings) coincide exactly at the twin primes (Theorem B), which is why the owner's list has double entries at 6 and 8 and the table has them at 14, 20, 32, 44, ....

## 5. {7, 21}: the crevasse dictionary (ANALOGY) and the coincidences (NUMEROLOGY)

**Exact on both sides.** THM-1370: `H(T) = prod H(S_i)` over strong components, strong `H`-values are at least the Busch floor `f(m) = 3, 5, 9, 15, 25, 45, 75` (`m = 3..9`), and `7 in (f(4), f(5)) = (5, 9)`, `21 in (f(6), f(7)) = (15, 25)` with no rescue by factorisation (`21 = 3 x 7` recurses into the first gap). Theorem A: sources are classified by shape, the three small shapes have the value sequences `p + 1 = 3, 4, 6, 8, 12, ...`, `p + q + 1 = 6, 8, 9, 10, ...` and `1 + p + p^2 = 7, 13, 31, ...`, and `2` lies below all of them while `5` lies in the crevasse `(4, 6)` of the squares and below the first semiprime value 6, with no rescue from any larger shape (Lemma 1). S57 section 3 phrased the same thing as "failure would need both shapes to fail at once". The dictionary *shape of source <-> strong component, floor `s(n) > sqrt n` <-> Busch floor, Goldbach partitions of `v - 1` <-> factorisations of `H`* is an ANALOGY: `s` is not multiplicative and `H` has no Goldbach shape; what the two proofs share is "a small value falls between consecutive terms of every finite family of value sequences and cannot be factored into attained values".

**CITED, exact.** `{1, 2, 4} = QR_7 = <2>` in `(Z/7)^x` is the multiplier group of `F_21 = Aut(P_7)`; `P_7` has 7 vertices, 21 arcs, `|Aut| = 21` and `H = 189 = 9 x 21 = 27 x 7` (THM-126; recomputed here, with `H(P_3) = 3`, `H(P_11) = 95095`); `x -> 2x + 1` and `x -> 4x + 1` generate `F_21` only at 7 (S18 section 3), where the equality `2^3 + 1 = 3^2` is Catalan's. The decomposition `7 = 1 + 2 + 4` is at once the proper divisors of 8 (the unique aliquot source of 7, Theorem A), the row of 1 below 8 (the only row that is a column, Theorem C (5)), and the subgroup `<2>` whose orbits `{0}, {1, 2, 4}, {3, 6, 5}` give `7 = 1 + 3 + 3`.

**NUMEROLOGY (no mechanism; recorded so it is not regenerated).** `7 = 2 + 5` and `21 = 4 x 5 + 1 = 2^4 + 5`; `7 = s(2^3)` and `21 = s(2 x 3^2)` have one and three sources; `63 = 2^6 - 1 = 7 x 9 = 21 x 3` is the first non-squarefree Mersenne number, the first abc hit on the row of 1 (Theorem D), and the first attained `H`-value of the form `7 x 3^k` (THM-1370's corollary, `H = 63` at `n = 8`). `H(P_7) = 189 = 3 x 63`. None of these is a relation between the untouchables and the forbidden path counts; the only relation established is the analogy above.

## 6. abc on the chain grid (Theorem D)

**Theorem D.** (i) `(1, s(n), n)` is an additive triple of pairwise coprime numbers iff `sigma(n) = 2n - 1` (`n` almost perfect; the powers of 2 are the only known ones). (ii) For `n = 2^k` the triple `(1, 2^k - 1, 2^k)` is an abc hit (`rad(abc) < c`) iff `2^k - 1` is not squarefree, iff `k` is a multiple of `ord_(p^2)(2)` for some odd prime `p` (`6, 20, 21, 110, 156, 253, ...` for `p = 3, 5, 7, 11, 13, 23`). The hits with `k <= 64` are `k = 6, 12, 18, 20, 21, 24, 30, 36, 40, 42, 48, 54, 60, 63` (FINITE-EXACT, matching the union of the multiples of 6, 20, 21 below 64; the first is `1 + 63 = 64`). (iii) The first abc hit of all, `1 + 8 = 9` (Catalan), is the row of 1 meeting the column of 3; the eight hits with `c <= 130` are `1 + 8 = 9`, `5 + 27 = 32`, `1 + 48 = 49`, `1 + 63 = 64`, `1 + 80 = 81`, `32 + 49 = 81`, `4 + 121 = 125`, `3 + 125 = 128`, every one of which has an odd perfect power among its terms. Among the 121 hits with `c <= 10^4`: 94 contain an odd perfect power, 48 have a power of 2 as the even term, 30 have `a = 1`; exactly one term of every hit is even (parity), so "every hit has an old term" is empty and the content is in the odd perfect powers.

*Proof.* (i) `1 + s(n) = n` iff `sigma(n) = 2n - 1`, and `gcd(n, n - 1) = 1`. (ii) `rad(1 x (2^k - 1) x 2^k) = 2 rad(2^k - 1)`; if `2^k - 1` is squarefree the radical is `2(2^k - 1) >= 2^k`; if not, `rad(2^k - 1) <= (2^k - 1)/3 < 2^(k-1)`. `p^2 | 2^k - 1` iff `ord_(p^2)(2) | k`. (iii) is the printed scan. QED

The abc conjecture (THM-3833 for its conditional consequences) says the hits are rare and the quality `log c / log rad(abc)` stays bounded; on the chain grid the hits of the row of 1 are governed by the squareful parts of Mersenne numbers, so the classical open question whether every `2^p - 1` (`p` prime) is squarefree is the statement that no prime `k` is a Mersenne abc hit (HYP-9415 (5)). The aliquot function and abc meet only on the row of 1: elsewhere `(s(n), n, sigma(n))` is an additive triple but not of the form `a + b = c` with a unit.

## 7. Typing summary

| claim | type |
|---|---|
| `s^-1(2) = s^-1(5) = {}`, `s^-1(7) = {8}`, `s^-1(21) = {18, 51, 91}`, the `t <= 3` regime below 15 | PROVED (Theorem A; (i) classical) |
| single-source Mersenne values exactly `3, 7` | PROVED for `k <= 26` via S57's `10^8` scan; beyond, conditional on a distinct-Goldbach partition of `2^k - 2` |
| twin collision, row-source uniqueness, fixed point only at `p = 5`, 6 the only even perfect number with a prime-square source | PROVED (Theorem B) |
| `s^-1(28) = {28}`, sources of 496 and 8128, `(N - 7)/3` primality, twin-value source counts to `10^4` | FINITE-EXACT |
| old = even or perfect power; disjoint fresh rays; the square-doubling law with first instance 6; odd rows meet earlier columns iff perfect powers; the row of 1 is the only row-column | PROVED (Theorem C) |
| almost-perfect triples; Mersenne abc hits = non-squarefree Mersenne numbers; the list to `k = 64`; the census to `10^4` | PROVED (Theorem D) + FINITE-EXACT |
| `{2, 5}` the only untouchable primes | 2, 5 PROVED; "only" conditional on distinct-Goldbach (S57), FINITE-EXACT to `10^10` through S57 |
| untouchables vs forbidden path counts | ANALOGY (crevasse dictionary) |
| `7 = 2 + 5`, `21 = 4 x 5 + 1`, `63 = 7 x 9 = 21 x 3` | NUMEROLOGY |

## 8. Pre-registered tests (HYP-9415), written before the C++ output was read

1. Twin collision values `p + 1` in `(10^4, 10^5]`: the share with exactly the two twin sources lies in `[0.40, 0.60]` and the mean number of sources in `[2.5, 3.2]`.
2. The fifth perfect number `33550336` is not its own sole source (28 is the only self-only perfect number): at least one further source exists.
3. abc hits with `c <= 10^5`: the share containing an odd perfect power lies in `[0.65, 0.85]`, the share whose even term is a power of two in `[0.30, 0.50]`, the share with `a = 1` in `[0.15, 0.30]`.

Outcomes: see section 8a below (filled after the run).

## 9. Inheritance check and audit

Filled after the repository sweep and the independent audit (section 10).

## 10. Reproduction

```
python3 04-computation/experiments/aliquot_chains_abc_20261026.py > 04-computation/experiments/aliquot_chains_abc_20261026.out   # numpy, sympy; 26 s
g++ -O2 -std=c++17 -o aliquot_sources_pre 04-computation/experiments/aliquot_sources_preregistered_20261026.cpp && ./aliquot_sources_pre   # about a minute, 600 MB
```
