# Small graphs that encode arithmetic: why the square-sum problem starts at 15 with ends 8 and 9 (a zigzag law), Friedman's fences as plane graphs, and the Collatz-alphabet sum graph

**Status block.**

- **PROVED (hand proofs below; every finite ingredient re-checked by the runner).**
  - *Square sums.* Degree formula (Lemma 1). The complete list of vertices of degree at most 1 of the square-sum graph `Q_n` for every `n` (Theorem 1); in particular `min deg Q_n >= 2` iff `n >= 31`. No Hamiltonian path for `2 <= n <= 14`; exactly one at `n = 15`, `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`, whose ends `8` and `9` are forced (Theorem 2). This chain is the *3-square zigzag*: its sums are `(9,16,25,16)^3,9,16`, and its ends are the two defects of the middle square `16` (`8 = 16/2` has no `16`-partner, `9` has no `9`-partner) (Proposition 2A).
  - *General target sets.* Leaf Lemma and its corollary (Theorem 3): whenever at most two vertices have degree at most 1, each of them is a target, a half-target, or the successor of a target. Zigzag Lemma (Theorem 4): if `2k+1, 4k, 6k+1` are targets, an explicit chain of `{1..4k-1}` runs from `2k` to `2k+1`. **Zigzag Threshold Theorem** (Theorem 5): for the shifted squares `S_c = {j^2 + c}` with `c = k(4-k)` (`k >= 2`), the arrangement problem first becomes possible at `n = 4k-1`, with exactly one chain, from `2k` to `2k+1`. The squares are the member `k = 4`.
  - *Pell and Catalan.* `t^2 - 2^m = 1` only for `(3,3)`, so `(8,9)` is the only pair that is both a Pell pair `(2s^2, t^2)` and a Catalan/Gersonides pair `(2^3, 3^2)` (T1.E).
  - *Fences.* A plane-graph identity: `#fields = n + c - V0` (Proposition F1). An isoperimetric upper bound `A(n) <= ((sqrt(1 + 4n/sqrt(pi)) - 1)/2)^2 < n/sqrt(pi)` (Theorem F2). Superadditivity, so `lambda = lim A(n)/n` exists, with `1/2 <= lambda <= 1/sqrt(pi) = 0.5642` (Proposition F3). Complete rectangular grids of `n` unit fences exist iff `2n+1` is composite, which is Sundaram's sieve (Proposition F4). Polyomino thresholds: `A(2k + ceil(2 sqrt k)) >= k`, via `ceil(sqrt k) + ceil(k/ceil(sqrt k)) = ceil(2 sqrt k)` (Proposition F5).
  - *Small graphs encoding arithmetic (T3).* Target sets with at most one element in every `(M, 2M)` give forests (Theorem L). Sums equal to powers of `p` preserve `v_p`; for `p = 2` the components are exactly the valuation classes (Theorem V). Two consecutive targets `m, m+1` give a chain of `{1..m}` (Lemma Z2). This is how the Catalan pair `(8,9)` produces the first perfect-power chains (`n = 7, 8`), and how the Gersonides pairs produce chains of `C_n` at `n = 1, 2, 3, 7, 8`.
  - *The Collatz-alphabet sum graph `C_n`* (`x ~ y` iff `x+y` is a power of 2 or of 3), a new instance. The largest power of two `<= n` always has degree at most 1, so `C_n` never has a Hamiltonian cycle (Theorem C1). **Window Theorem** (Theorem C2): the `n` for which `C_n` has no isolated vertex and at most two vertices of degree at most 1 are exactly `[2,3]` and one window `W_a` for each `a >= 2`. With `B = 3^(a-1)`:
    - `W_a = [max(3^a - 2^p, 2^(p+1) - 3^(a-1)), 3^a - 1]` if `(B, 3B)` contains two powers of two `2^p < 2^(p+1)`;
    - `W_a = [2*3^(a-1), 3^a]` if it contains one.

    So the windows follow the Beatty word of `log_2 3`, and their left ends are the gaps `|2^q - 3^b|` (or `3^a - 3^(a-1)`).
  - *The parity graph.* For constant signs (Collatz, `3n-1`, `5n+-1`, ...), every node has in-degree and out-degree 2, so degree carries no threshold (T3.A).
- **FINITE-EXACT** (exhaustive; every non-local "no" is decided by an exact solver and cross-checked by a second method, section 4):
  - Square sums: paths exactly for `n in {1,15,16,17,23} u [25,300]` and cycles exactly for `n in [32,300]`. Path counts for `n <= 36` agree with A071983 and cycle counts for `n <= 40` with A071984. The forced endpoints are `{8,9}`, `{8,16}`, `{16,17}` at `n = 15, 16, 17` and `{18}` at `n = 23, 25..30`; there are none for `31 <= n <= 300`.
  - The Zigzag Threshold Theorem, exhaustively for `k <= 40`. Its converse over `-300 <= c <= 100`: the zigzag signature at the first possible `n` occurs **iff** `c = k(4-k)`.
  - Variants (the first possible `n`, its chains and ends, cycles): triangular, pentagonal, perfect powers, `s^2 +- 1`, powers of 2, powers of 3, powers of 2 or 3, and cubes. For cubes, paths below 400 occur exactly at `305, 333, 385..399` (first `305`, reproducing A304120), and cycles up to 600 occur exactly at `473..600`.
  - `C_n`: the Hamiltonian-path set for `n <= 2200`, and the Window Theorem re-checked to `3^13`.
  - The hitting-time table. The requested hybrids.
- **CITED:** OEIS A071983, A071984, A090461 (Gerbicz's theorem that chains exist for all `n >= 25` and cycles for all `n >= 32`), A304120 (cubes 305, Resta; fourth powers 6479, Fredette), A116980, A399783/A399784; Mihailescu (Catalan); THM-4484 (Gersonides); Harary–Harborth (minimal polyomino perimeter; UNCITED-RECOLLECTION, used only for context); E. Friedman's *Maximum Fence Area* page (best known values, fetched 2026-09-26); Bollobás / Komlós–Szemerédi hitting times (analogy only, UNCITED-RECOLLECTION of the exact statement); repository results THM-4474, THM-4481, THM-4484, THM-4486 and the Kuratowski, brackets and square-sum lanes, as referenced inline.
- **EMPIRICAL:** the records fit `A(n) ~ n/2 - 0.641 sqrt(n) + 0.093` (`12 <= n <= 50`). A numerical optimizer reproduces the `n = 6` record `1.475853` in its pentagon-plus-chord topology, and finds at best `2.033 < 2.10306` for `n = 8` in a hexagon-plus-two-chords family. **No fence record is improved.**
- **ANALOGY:**
  - Grids versus irregular records, compared with free versus sporadic cycles (THM-4484).
  - The `C_n` windows versus Collatz cycle denominators `2^p - 3^a`: the same numbers appear, with no map between the two problems.
  - Square-type thresholds versus random-graph hitting times.
- **NUMEROLOGY:** reading the square-sum ends `9 - 8 = 1` as `3^2 - 2^3` (Catalan, the Gersonides identity behind the free cycle `-5,-7,-10`); "area 5 first needs 15 fences". (That `15` is also the first `n` where `Q_n` has a cycle is not numerology but the same edge count as the chain, see the remark in 1.3.)
- **OPEN:**
  - Whether `lambda = 1/2` for fences.
  - Which admissible `n` in a window `W_a` are actually Hamiltonian for `C_n`. They are only partly so, for example `W_5 = [162,243]` has paths exactly on `179-194, 224-243`. Also, whether `C_{3^a - 1}` or `C_{3^a}` always has one.
  - Whether the simultaneous Pell system `2s^2+1 = t^2`, `6s^2+1 = u^2` has only `s = 2`. It is checked here for the first 400 Pell solutions.
  - The converse of the Zigzag Threshold Theorem beyond the computed range of `c`.
- **REFUTED (own intermediate guesses, recorded in section 4):** that every sign-strategy parity graph is 2-regular; that the Catalan-type endpoint pair is a general law for "first possible" arrangements.

Collatz is not addressed and remains OPEN. No novelty is claimed for any individual elementary fact.

Session `collatz-procgen-20260922`, lane `smallgraph` (2026-09-26). Scripts: `04-computation/experiments/procgen_smallgraph_20260926_{lib,t1,t2,t3,fence6,run}.py` and `procgen_smallgraph_20260926_ham.c`. Output: [procgen_smallgraph_20260926.out](procgen_smallgraph_20260926.out) (one runner; every claim there is a `check(...)` that raises on failure). Reproduction: section 5.

## 0. The owner's three prompts, answered in one paragraph each

**Square sums at 15, with 8 at one end and 9 at the other.** The graph `Q_n` (vertices `1..n`, `x ~ y` iff `x + y` is a square) has a Hamiltonian path first at `n = 15`, because every smaller `n >= 2` has an isolated vertex or three vertices of degree 1. The `n = 15` chain is unique. Its ends are the only two vertices of degree 1: `8`, whose partner via `16` would be itself, and `9`, whose only partner is `7`. The chain is a zigzag through the three squares `9, 16, 25`. The pair `(8, 9)` is the `k = 4` case of a law: in the family `S_c = {j^2 + c}` with `c = k(4-k)`, the first chain always appears at `n = 4k - 1`, always unique, always from `2k` to `2k + 1`. It has half-target/target ends at distance 1 (proved). In the computed range, only this family shows that signature. The Pell form `3^2 - 2*2^2 = 1` is exactly "`c = 0` at `k = 4`". The Catalan form `3^2 - 2^3 = 1`, the identity behind the free Collatz cycle `-5, -7, -10`, is a coincidence of the small number `8 = 2*2^2 = 2^3` (NUMEROLOGY). In problems where 8 and 9 are both *targets* (perfect-power sums, sums in `{2^p} u {3^a}`), the same pair produces the first chains through a different, proved mechanism: consecutive targets zigzag.

**Fences.** A fence configuration is a plane graph whose bounded faces are the fields. Euler's formula becomes the exact count `#fields = n + c - V0`. Per-field isoperimetry plus the isoperimetric inequality for the union bounds the area by `U(n) < n/sqrt(pi)`; the records reach 54% to 88% of `U(n)`. The arithmetic content is a small graph: a complete `i x j` grid needs `n = i + j + 2ij` fences, so grids exist exactly when `2n + 1` is composite (Sundaram's sieve). The "trivial" records `n = 4, 7, 12` are the grids with `2n+1 = 9, 15, 25`. The least `n` whose record reaches area `k` is the polyomino count `2k + ceil(2 sqrt k)` for every `k <= 20` except `13` and `17`, where the records save one fence. Like the square-sum threshold, these thresholds are archimedean (square-root, isoperimetric). The Collatz provability thresholds are not: they are Diophantine, cycle-driven.

**Small graphs encoding arithmetic.** The repository's instances sort into four threshold mechanisms (section 3.1). The new instance pushed here is the *Collatz-alphabet sum graph* `C_n`: arrange `1..n` so that adjacent sums are powers of 2 or of 3. It never closes into a cycle (proved), and its admissibility windows are cut out exactly by the powers of 3 and by the gaps `|2^p - 3^a|` (Window Theorem, proved). Hamiltonian paths occur, for `n <= 2200`, exactly on `1-3, 5-8, 18-26, 49-63, 65-66, 68-80, 179-194, 224-243, 473-575, 665-728, 1319-1418, 1536-1616, 1620-1663, 1703-1713, 1792-1802, 1920-2114, 2120-2186`.

## 1. T1: the square-sum problem at 15

### 1.1 Setting and the degree formula

For a set `S` of positive integers (the *targets*) let `G_S(n)` have vertex set `[n] = {1..n}` and an edge `{x,y}` iff `x != y` and `x + y in S`; `Q_n = G_S(n)` for `S` the squares. A *half-target* is an `x` with `2x in S`.

**Lemma 1 (PROVED).** For `1 <= x <= n`, `deg_n(x) = |S cap [x+1, x+n]| - [2x in S]`. For squares, `deg_n(x) = floor(sqrt(x+n)) - floor(sqrt x) - [x = 2s^2]`.

*Proof.* `y -> x + y` maps the neighbours of `x` bijectively onto `S cap [x+1, x+n]` minus `{2x}`, and `2x` lies in `[x+1, x+n]` because `1 <= x <= n`. ∎ (T1.A1 re-checks it for all `x <= n <= 300`; T1.C1 on 13 target families.)

Two elementary consequences are used throughout.

- **Leaves are ends.** In a Hamiltonian path of a graph on `n >= 2` vertices a vertex of degree 1 is an end, and a vertex of degree 0 cannot occur. So three vertices of degree `<= 1` rule out a Hamiltonian path.
- **Minimum degree for cycles.** A Hamiltonian cycle needs minimum degree 2.

### 1.2 Theorem 1: every vertex of degree at most 1 in every `Q_n`

**Theorem 1 (PROVED).** A vertex `x` of `Q_n` (`x <= n`) has degree `<= 1` iff `x <= n <= N(x)`, where

| x | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 16 | 17 | 18 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| N(x) | 7 | 13 | 5 | 11 | 10 | 9 | 8 | 16 | 15 | 14 | 13 | 12 | 19 | 18 | 30 |

and no other `x` is ever a leaf. Degree 0 occurs only for `x = 1` (`n <= 2`), `x = 2` (`n <= 6`), `x = 4` (`n = 4`). In particular `Q_n` has minimum degree `>= 2` iff `n >= 31`.

*Proof.* Let `sigma_1(x) < sigma_2(x)` be the two smallest squares above `x` other than `2x`. By Lemma 1, `deg_n(x) <= 1` iff `sigma_2(x) > x + n`, i.e. `n <= sigma_2(x) - x - 1`. Similarly, degree 0 iff `n <= sigma_1(x) - x - 1`.

- (i) The number of squares in `(x, 2x]` is `floor(sqrt(2x)) - floor(sqrt x) > (sqrt2 - 1) sqrt x - 1`. Since `68^2 = 4624 > 4608 = 2*48^2`, we get `(sqrt2-1)^2 * 24 > 4`, so for `x >= 24` there are at least 2 squares in `(x, 2x]`. Since `150^2 = 22500 > 22472 = 2*106^2`, there are at least 3 for `x >= 53`.
- (ii) Every `n` with `x <= n` has `x + n >= 2x`. So if `x >= 24` is not a half-target, `sigma_2(x) <= 2x <= x + n` and `x` is never a leaf. If `x >= 53` is a half-target, two of the at least three squares in `(x, 2x]` differ from `2x`. Same conclusion.
- (iii) The half-targets in `[24, 52]` are `32` and `50`: `(32, 64]` contains `36, 49` and `(50, 100]` contains `64, 81`, so neither is ever a leaf.
- (iv) For `x <= 23` the table is `sigma_2(x) - x - 1`, read off directly. For example `sigma_2(8) = 25` because `16 = 2*8` is excluded, so `N(8) = 16`; `sigma_2(18) = 49` because `36` is excluded, so `N(18) = 30`; for `19 <= x <= 23`, `sigma_2(x) - x - 1 < x`. ∎

(T1.A2-A5 re-check the inequalities exactly, scan `x <= 60`, and compare predicted and actual leaf sets for all `n <= 300`.)

### 1.3 Theorem 2: why 15, and why the ends are 8 and 9

**Theorem 2 (PROVED).**
- (a) `Q_n` has no Hamiltonian path for `2 <= n <= 14`.
- (b) `Q_15` has exactly one, `8,1,15,10,6,3,13,12,4,5,11,14,2,7,9`, up to reversal, and every Hamiltonian path of `Q_15` ends at 8 and 9.

*Proof.*

(a) For `2 <= n <= 6`, vertex `2` is isolated (its only square partner below `9` is `4 - 2 = 2` itself). For `7 <= n <= 14`, Theorem 1 gives at least three leaves:
- `n = 7`: `{1,2,4,5,6,7}`
- `n = 8`: `{2,4,5,6,7,8}`
- `n = 9`: `{2,4,5,6,8,9}`
- `n = 10`: `{2,4,5,8,9,10}`
- `n = 11`: `{2,4,8,9,10,11}`
- `n = 12`: `{2,8,9,10,11,12}`
- `n = 13`: `{2,8,9,10,11}`
- `n = 14`: `{8,9,10}`

(b) By Lemma 1 at `n = 15`:
- `8` and `9` have degree 1;
- `1` (partners `3, 8, 15`) and `3` (partners `1, 6, 13`) have degree 3;
- the other 11 vertices have degree 2.

So `Q_15` has `(2 + 6 + 22)/2 = 15` edges. A Hamiltonian path has 14 edges, ends `8` and `9`, and path-degree 2 at the other 13 vertices. It therefore omits exactly one edge, and that edge must lower both `deg 1 = 3` and `deg 3 = 3`, so it is `{1,3}`. The graph `Q_15 - {1,3}` has degree sequence `(1,1,2^13)`, so it is a path from `8` to `9` together with disjoint cycles. The displayed sequence is a path of `Q_15 - {1,3}` through all 15 vertices, so there are no cycles and the path is unique. ∎ (T1.B1-B5; the exhaustive count gives 1.)

**Proposition 2A (anatomy; PROVED).** The only squares in `[3, 29]` are `4, 9, 16, 25`, and the edges of `Q_15` are:
- `4`: `{1,3}`;
- `9`: `{1,8},{2,7},{3,6},{4,5}`, all inside `L = {1..8}`;
- `25`: `{10,15},{11,14},{12,13}`, all inside `H = {10..15}`;
- `16`: `{1,15},{2,14},{3,13},{4,12},{5,11},{6,10},{7,9}`, each joining `{1..7}` to `{9..15}`.

Along the chain, which avoids `{1,3}`, the following hold.
- Consecutive sums differ.
- A `9`-edge and a `25`-edge never share a vertex.
- Two `9`-edges cannot be joined by a `16`-edge: two vertices of `L` summing to 16 would both be `8`.
- Two `25`-edges cannot be joined by a `16`-edge either.

So the sum word is forced to be `9,16,25,16,9,16,25,16,...`, the **3-square zigzag** (T1.B6-B7). The ends are the two defects of the middle square `16`:
- `8` is the vertex of `L` whose 16-partner is itself (`2*8 = 4^2`) and whose 25-partner `17` is out of range;
- `9` is the vertex whose 9-partner is `0` and whose 25-partner `16` is out of range.

In the zigzag the odd positions run `8,15,6,13,4,11,2,9` (steps `+7, -9`) and the even positions run `1,10,3,12,5,14,7` (steps `+9, -7`). The drift `+-2` per period is the second difference of the squares.

*The window `{15, 16, 17}` (FINITE-EXACT, exhaustive).* `Q_16` and `Q_17` also have exactly one chain each. They are the zigzag with `16` appended at the `9` end (`9 + 16 = 25`), and then `17` prepended at the `8` end (`17 + 8 = 25`), with forced ends `{8,16}` and `{16,17}`. At `n = 18` the leaves are `{16, 17, 18}` and the window closes (T1.B9-B10).

*Remark.* `Q_n` is a forest for `n <= 14`, and `Q_15` is unicyclic: its unique cycle `(1,3,6,10,15)` has sums `4,9,16,25,16` and is closed by the unused edge `{1,3}` (T3.B2). So "first cycle" and "first chain" coincide at 15. By the edge count above this is the same fact, not a second mechanism.

### 1.4 Theorem 3: the Leaf Lemma for arbitrary target sets

**Theorem 3 (PROVED).** Let `s_1 < s_2 < ...` be the targets, and let `x in [n]` be a non-half-target with `deg_n(x) <= 1`, where `s_i <= x < s_{i+1}` (for `x < s_1` read `s_i` as `1`). Then every `y in [s_i, x]` has `deg_n(y) <= 1`.

*Proof.* `deg_n(x) <= 1` and `2x notin S` give `|S cap [x+1, x+n]| <= 1`. The elements of `S` above `x` are `s_{i+1} < s_{i+2} < ...`, so `s_{i+2} > x + n`. For `s_i <= y <= x`, the targets in `[y+1, y+n]` are among `s_{i+1}, s_{i+2}, ...`, and `s_{i+2} > x + n >= y + n`. So there is at most one, and `deg_n(y) <= 1 - [2y in S] <= 1`. ∎

**Corollary.** If `G_S(n)` has at most two vertices of degree `<= 1`, each of them is a target, a half-target, a successor `s_i + 1` of a target `s_i` that is itself of degree `<= 1`, or one of `1, 2`. Consequently, at any `n` where the leaves force the ends of all Hamiltonian paths, the forced ends are of these types. Checks: T1.C1-C2 on 13 families, `n <= 160`, 2067 graphs.

For squares:
- `n = 15`: ends `8` (half-target) and `9` (target);
- `n = 16`: `8` (half-target) and `16` (target);
- `n = 17`: `16` (target) and `17` (successor of the target leaf `16`).

The *types* of the forced ends are therefore a law. Their *difference* is not: `9 - 8 = 1`, `16 - 8 = 8`, `17 - 16 = 1`.

### 1.5 Theorems 4 and 5: the zigzag law

**Theorem 4 (Zigzag Lemma, PROVED).** Let `k >= 1` and `{2k+1, 4k, 6k+1} subset S`. Define `Z_k` of length `4k - 1` by

```text
z_{4j+1} = 2k - 2j,   z_{4j+2} = 2j + 1,   z_{4j+3} = 4k - 1 - 2j,   z_{4j+4} = 2k + 2 + 2j.
```

Then `Z_k` is a Hamiltonian path of `G_S(4k-1)` from `2k` to `2k+1`, with sums `(2k+1, 4k, 6k+1, 4k)^(k-1), 2k+1, 4k`.

*Proof.* For `0 <= j <= k-1` the four families are:
- the evens `2..2k`;
- the odds `1..2k-1`;
- the odds `2k+1..4k-1`;
- for `j <= k-2`, the evens `2k+2..4k-2`.

That is `4k - 1` distinct values filling `[1, 4k-1]`. The sums are
- `z_{4j+1} + z_{4j+2} = 2k+1`,
- `z_{4j+2} + z_{4j+3} = 4k`,
- `z_{4j+3} + z_{4j+4} = 6k+1`,
- `z_{4j+4} + z_{4j+5} = 4k`.

The first term is `2k` and the last is `z_{4k-1} = 2k+1`. ∎ (T1.D1, `k <= 400`; `Z_4` is the square chain, T1.D1b.)

**Theorem 5 (Zigzag Threshold Theorem, PROVED).** Let `k >= 2`, `c = k(4-k)` and `S_c = {j^2 + c : j >= 1}`. Then `G_{S_c}(n)` has no Hamiltonian path for `2 <= n <= 4k-2`, and exactly one for `n = 4k-1`, namely `Z_k`, from `2k` to `2k+1`.

*Proof.* Write `j = k - 2 + i`; then `j^2 + c = 2ki + (i-2)^2` (with `i >= 3-k`).

**The relevant targets.** For `n <= 4k-1` only targets `<= 8k-3` matter. Among `i >= 1` these are `2k+1, 4k, 6k+1` (`i = 1, 2, 3`), since `i >= 4` gives `>= 8k+4`. For `i = 0` the target is `4`, present iff `k >= 3` (for `k = 2`, `i = 0` means `j = 0`). For `i <= -1` the value is `9 - 2k, 16 - 4k, ...`, which is below 3 when `j >= 1` (for `k <= 3`, `i = -1` means `j = 0`). Also `3 notin S_c`. (T1.D2.)

Let `k >= 3`.

- **`2 <= n <= 2k-2`.** Vertex `2` is isolated: its partners would be `2` (itself), `2k-1`, `4k-2` and `6k-1`, all out of range.
- **`n in {2k-1, 2k}`.** Only the targets `4` and `2k+1` are `<= 2n-1`. The three vertices `2, 4, 5` (all `<= 2k-1`) have at most the partner `2k+1-x` each.
- **`n = 2k+1`.** The three vertices `2, 4, 2k` each have one partner: `2k-1`, `2k-3` and `1` respectively.
- **`2k+2 <= n <= 4k-2`.** The three vertices `2k, 2k+1, 2k+2` each have one partner: `1`, `2k-1` and `2k-2` respectively. Their other candidate partners are `2k` itself, `4k+1`, `4k` and `4k-1`, all out of range.
- **`n = 4k-1`.** The degrees are as follows:
  - `1` has neighbours `{3, 2k, 4k-1}` and `3` has neighbours `{1, 2k-2, 4k-3}`;
  - every other `x <= 2k-1` has neighbours `{2k+1-x, 4k-x}`;
  - `2k` has neighbours `{1}` and `2k+1` has neighbours `{2k-1}`;
  - every `x >= 2k+2` has neighbours `{4k-x, 6k+1-x}`.

  So there are `4k - 1 = n` edges. As in Theorem 2, the ends `2k, 2k+1` are forced and the one omitted edge is `{1,3}`. `G - {1,3}` has degree sequence `(1,1,2,...,2)`, and `Z_k` (in which `1` sits between `2k` and `4k-1`, and `3` between `2k-2` and `4k-3`) spans it, so the path is unique.

For `k = 2` (`c = 4`, targets `5, 8, 13`) the cases `n <= 6` are disconnected or have an isolated vertex. `G(7)` has exactly the six edges of `Z_2 = (4,1,7,6,2,3,5)`. ∎ (T1.D3: exhaustive for `k <= 40`, e.g. `k=9, c=-45`: `n = 35`, ends `(18,19)`; T1.D4: every ingredient of this proof for `k <= 119`.)

**Converse (FINITE-EXACT, T1.D5).** Consider every `c` with `-300 <= c <= 100` and the first `n >= 3` with a chain in `G_{S_c}(n)`. That chain is unique, with ends a half-target `u` and the target `u + 1`, for exactly the 18 values `c = k(4-k)`, `k = 2..19`:

```text
-285, -252, -221, -192, -165, -140, -117, -96, -77, -60, -45, -32, -21, -12, -5, 0, 3, 4
```

For the other 383 values of `c` it is not. Uniqueness is decided exactly: `#paths >= 2` iff `G - e` has a Hamiltonian path for some edge `e` of the first path. This was validated against exhaustive counts on 213 instances (T0.2). Examples of non-members:
- `c = 1` (`s^2 + 1`): first at `16`, unique, ends `(5,13)`;
- `c = 2`: first at `15`, ends `(9,11)`, half-target and target but at distance 2;
- `c = -1` (`s^2 - 1`): first non-trivial at `20`, 2 chains.

### 1.6 The Pell/Catalan question, typed

The owner notes that `8 = 2^3 = 2*2^2` and `9 = 3^2`, the Catalan/Gersonides pair behind the free cycle `-5,-7,-10` (THM-4484) and the first Pell pair `3^2 - 2*2^2 = 1`. The findings split into three parts.

- **Law (PROVED).** Two facts are forced by general statements.
  - The forced ends at a threshold are a target and a half-target or successor (Theorem 3).
  - When the first chain is a zigzag, its ends are consecutive (Theorems 4, 5).
  For squares, the zigzag needs the three squares `2k+1, 4k, 6k+1`. With `4k = (2s)^2`, `h = 2k = 2s^2`, this is the simultaneous Pell system `t^2 - 2s^2 = 1`, `u^2 - 6s^2 = 1`. Among the first 400 Pell solutions only `s = 2` works (`9, 16, 25`) (T1.E1; the general case is OPEN here). So `3^2 - 2*2^2 = 1` *is* the zigzag condition, and it is the `c = 0` member of the family `c = k(4-k)`: `k = 4`.
- **Not a law.** The same numbers do not recur at other thresholds of `Q_n` (`(8,16)` at 16, `(16,17)` at 17, `18` at 23..30). Among the variants, only the squares and the powers of 2 or 3 show the signature (T1.G*); the latter at `n = 3`, through the Gersonides pair `(3,4)` (`Z_1`).
- **NUMEROLOGY.** Reading `9 - 8` as `3^2 - 2^3` imports nothing into the square-sum problem. The leaf property of `8` uses `2*8 = 4^2`, never `8 = 2^3`, and the coincidence `2*2^2 = 2^3` is the elementary fact that `t^2 - 2^m = 1` only for `(3,3)` (T1.E2-E4). The Catalan pair *is* a mechanism in two neighbouring problems where 8 and 9 are both targets. There, consecutive targets `m, m+1` give the zigzag chain `m, 1, m-1, 2, ...` (Lemma Z2, section 3.4). So the first perfect-power chains (`n = 7, 8`; A399783) and the chains of `C_n` at `n = 7, 8` are Catalan zigzags. The map from these to Collatz stops at the shared identity: the Collatz cycle equation needs `2^p` against `3^a`, and the square-sum threshold needs `2s^2` against `t^2`. They meet only at `(8,9)`, by the elementary theorem above.

*Remark on the Pell wall (typed).* The owner's Kuratowski source has a "Pell wall" `1 + sqrt2`. The square-sum degree constant is `sqrt2 - 1 = 1/(1 + sqrt2)`: `deg_n(n) ~ (sqrt2 - 1) sqrt n`. Both come from a factor 2 under a square root: the interval `(x, 2x]` here, the isosceles limit of the Berggren B-ray there. Beyond that shared `sqrt2`, the connection is NUMEROLOGY: no map relates the two objects.

### 1.7 Square-sum thresholds and forced endpoints up to 300 (FINITE-EXACT, T1.F)

| quantity | result | certificate |
|---|---|---|
| Hamiltonian path | `n in {1,15,16,17,23} u [25,300]` | `2..14`, `18`: local (Theorem 2, three leaves); `19..22`, `24`: exact search, three methods agree (forcing solver, CP-SAT, exhaustive DFS) |
| Hamiltonian cycle | `n in [32,300]` | `n <= 30`: a leaf (Theorem 1); `n = 31`: exact search, three methods agree |
| path counts `15..36` | `=` A071983 | exhaustive DFS (C) |
| cycle counts `32..40` | `=` A071984 | exhaustive DFS (C) |
| forced ends | `15:{8,9}`, `16:{8,16}`, `17:{16,17}`, `23, 25..30:{18}`; none for `31..36` | exhaustive endpoint histograms |
| no forced end | every `37 <= n <= 300` | two verified chains with disjoint endpoint pairs per `n` |

The number of vertices that can be an end is `2, 2, 2, 4, 11, 10, 17, 20, 16, 16, 26, 32, 33, 34, 35, 36` for `n = 15, 16, 17, 23, 25..36`. For `32 <= n <= 36` every vertex is an end of some chain; beyond 36 only the absence of forced ends was checked. The all-`n` statements (chains for all `n >= 25`, cycles for all `n >= 32`) are Gerbicz's theorem: CITED via A090461 and the audited [square-sum lane](collatz_mod6_20260922_w6_square_sum_hamiltonicity.md).

### 1.8 Variants (FINITE-EXACT, T1.G)

"First" means the first `n >= 3`. The zigzag signature means one chain whose ends are `u, u+1` with `2u in S` and `u + 1 in S`.

| targets | chains for `n <=` range | cycles | first chain: `n`, #chains, ends, leaves | signature |
|---|---|---|---|---|
| squares | `1, 15-17, 23, 25-60` | `32-60` | 15, 1, `(8,9)`, `{8,9}` | yes (`Z_4`) |
| triangular | `1-2, 9-60` | `12-13, 15-60` | 9, 1, `(3,5)`, `{3,5}` | no (distance 2) |
| pentagonal | `1, 29-34, 36-40, 45-80` | `33, 57-80` | 29, 8, 8 pairs, none | no |
| perfect powers | `1, 7-8, 15-60` | `17-60` | 7, 2, `(1,4),(4,6)`, `{4}` | no (Catalan acts through sums) |
| `s^2 + 1` | `1, 16, 20-80` | `27-28, 32-80` | 16, 1, `(5,13)`, `{5,13}` | no |
| `s^2 - 1` | `1-2, 20-27, 31-80` | `23, 39-80` | 20, 2, `(10,12),(12,19)`, `{12}` | no |
| powers of 2 | `1` | none | never (Theorem V) | - |
| powers of 3 | `1-2` | none | never for `n >= 3` (Theorem V) | - |
| powers of 2 or 3 | `1-3, 5-8, 18-26, 49-63, 65-66, 68-80` | none | 3, 1, `(2,3)`, `{2,3}` | yes (`Z_1` on `3,4`) |
| cubes | `1, 305, 333, 385-399` (`n < 400`) | `473-600` (`n <= 600`) | 305, 3, `(95,256),(186,256),(222,256)`, `{256}` | no |

For cubes, every `n < 295` has an isolated vertex or three leaves. The values `295-304, 306-332, 334-343` are excluded by exact search (forcing solver and CP-SAT agree), and `344-384` by the leaf count. The first chain at 305 reproduces A304120 (CITED: G. Resta via Rivera's Puzzle 311). Its forced end `256 = 512/2` is the half-target of `8^3`, exactly as Theorem 3 predicts. The cube cycle threshold `473` equals the minimum-degree-2 threshold.

## 2. T2: Friedman's fences read against the repository

**Setting.** A configuration is `n` closed segments of length 1, the fences, subject to three rules:
- two fences meet in at most one point, which is an end of at least one of them (no crossings, no overlaps);
- every end lies on another fence;
- every field (bounded face of the complement) has area `<= 1`.

`A(n)` is the supremum of the total field area. The table on Friedman's page (fetched 2026-09-26 with a generic user agent; page sha256 `78504fff3a667821649fdc7cfbb6020d8a098b65c4f2ba46b19aecbcd87afa27`) lists **best known** values for `3 <= n <= 50`. These are lower bounds, not proven optima. Most were found by several contributors in August–September 2026. Our reading of the rules is confirmed by the numerical reproduction of the `n = 6` record (section 2.6).

### 2.1 Euler's formula as an exact field count

**Proposition F1 (PROVED).** Call a junction a *T-point* if some fence passes through it in its interior, and a *V0-point* otherwise (only ends meet there). With `c` components, `T` T-points and `V0` V0-points:
- the number of fence pieces is `n + T`;
- `#fields = n + c - V0`;
- every component that encloses a field has at least 3 V0-points, so `#fields <= n - 2` and `A(n) <= n - 2`.

*Proof.* Every fence end lies on another fence, so the vertices of the plane graph are exactly the junctions: `V = T + V0` (at most one fence passes through a point, since fences do not cross). Each fence is cut into `1 + (number of T-points on it)` pieces, so `E = n + T`. Euler's formula `V - E + F = 1 + c` counts the unbounded face once, and every bounded face is a field, so `#fields = F - 1 = n + c - V0`. An extreme point of the convex hull of a component is a fence end. The fence it lies on must end there too, since otherwise the point is interior to a segment and not extreme. So it is a V0-point, and a component that encloses a field has at least 3 extreme points. ∎

Checks: T2.D1 on 18 exact rational configurations (grids, polyominoes, the `n = 5` record, and a square with two chords); T2.D3 rejects a crossing control.

This is the fence form of the Kuratowski-type counting of the Kuratowski lane: `E <= g(V-2)/(g-2)` becomes an identity, because every vertex is a junction and every bounded face is a field. Fields are created by T-junctions: `#fields = n + c - V0` is largest when few junctions are V0-points.

### 2.2 An honest upper bound

**Theorem F2 (PROVED).** `A(n) <= U(n) := ((sqrt(1 + 4n/sqrt(pi)) - 1)/2)^2`, so `A(n) <= n/sqrt(pi) - sqrt(A(n))` and `U(n)/n -> 1/sqrt(pi) = 0.5642`.

*Proof.*
1. **Per field.** Let the fields have areas `a_i <= 1` and perimeters `p_i`. The isoperimetric inequality gives `4 pi a_i <= p_i^2`, and since `a_i <= 1`, `a_i <= sqrt(a_i) <= p_i/(2 sqrt pi)`.
2. **Total perimeter.** A non-junction point of a fence has two sides, each in one face. So `sum p_i <= integral over the fences of m`, where `m in {0,1,2}` counts the sides lying in fields.
3. **The outer boundary.** A boundary point of the union `U` of the closed fields has one side in the unbounded face (the only face that is not a field), so there `m <= 1`. Hence `sum p_i <= 2n - H^1(dU) <= 2n - 2 sqrt(pi A)`, by the isoperimetric inequality for `U`, which has area `A`.
4. **Conclusion.** Combining, `A <= (2n - 2 sqrt(pi A))/(2 sqrt pi)`. Solve for `sqrt A`. ∎

T2.A1: every record satisfies the bound.

| n | record `A(n)` | `A/n` | `U(n)` | record/`U` | best grid | polyomino bound |
|---|---|---|---|---|---|---|
| 3 | 0.43301 | 0.144 | 0.799 | 0.542 | - | 0 |
| 4 | 1 | 0.250 | 1.173 | 0.852 | 1 | 1 |
| 5 | 1 | 0.200 | 1.569 | 0.638 | - | 1 |
| 6 | 1.47585 | 0.246 | 1.979 | 0.746 | - | 1 |
| 7 | 2 | 0.286 | 2.400 | 0.833 | 2 | 2 |
| 8 | 2.10306 | 0.263 | 2.831 | 0.743 | - | 2 |
| 9 | 2.63630 | 0.293 | 3.270 | 0.806 | - | 2 |
| 10 | 3.04687 | 0.305 | 3.715 | 0.820 | 3 | 3 |
| 12 | 4 | 0.333 | 4.621 | 0.866 | 4 | 4 |
| 13 | 4.16199 | 0.320 | 5.080 | 0.819 | 4 | 4 |
| 15 | 5.06345 | 0.338 | 6.011 | 0.842 | - | 5 |
| 17 | 6.01086 | 0.354 | 6.954 | 0.864 | 6 | 6 |
| 24 | 9.02394 | 0.376 | 10.327 | 0.874 | 9 | 9 |
| 33 | 13.01887 | 0.395 | 14.774 | 0.881 | - | 12 |
| 42 | 17.02275 | 0.405 | 19.303 | 0.882 | 16 | 16 |
| 50 | 20.55829 | 0.411 | 23.375 | 0.880 | - | 20 |

`A(3) = sqrt(3)/4` is optimal (PROVED). By Proposition F1 there is at most one field. Its boundary is a closed polygon whose sides lie on three straight unit fences, so it is a triangle with sides `<= 1`. For a longest side `c <= 1` the area is at most `(sqrt3/4) c^2`.

### 2.3 The density constant

**Proposition F3 (PROVED).** `A(m + k) >= A(m) + A(k)`, by disjoint union. By Fekete's lemma `lambda := lim A(n)/n = sup A(n)/n` exists, and `1/2 <= lambda <= 1/sqrt(pi)`. The lower bound comes from `m x m` grids, `A(2m(m+1)) >= m^2`; the upper bound from F2. ∎

The records are superadditive (T2.E; a violation would have been beaten by a disjoint union). They fit `A(n) ~ n/2 - 0.641 sqrt(n) + 0.093` (EMPIRICAL, `12 <= n <= 50`, residual `<= 0.13`). The polyomino has coefficient `1/sqrt2 = 0.707`, so the records save about `0.07 sqrt n` over squares. **OPEN:** is `lambda = 1/2`? Beating squares asymptotically would need cells with half-perimeter `< 2` per unit area (hexagon-like) built from straight unit fences meeting in T-junctions. A T-junction gives a corner to only two of its three faces. Hales's honeycomb theorem (UNCITED-RECOLLECTION) would cap `lambda` at `12^(-1/4) = 0.537` if it transfers to these configurations; this was not checked.

### 2.4 Grids are Sundaram's sieve

**Proposition F4 (PROVED).** An `i x j` grid of unit squares uses `i(j+1) + j(i+1) = 2ij + i + j` fences, so `2n + 1 = (2i+1)(2j+1)`. A complete rectangular grid with exactly `n` fences exists iff `2n + 1` is an odd composite. This is Sundaram's sieve `n = i + j + 2ij`. ∎ (T2.B1.)

The records equal the best grid exactly at `n = 4, 7, 12` (`2n + 1 = 9, 15, 25`), which are the page's "trivial" values besides `3` and `5` (T2.B2). The first record found by search, `n = 6`, has `2n+1 = 13` prime. Elsewhere the records beat the grid, by margins from `+0.011` (`n = 17`) to `+3.821` (`n = 46`).

### 2.5 Integer thresholds are polyomino thresholds

**Proposition F5 (PROVED).** Let `w = ceil(sqrt k)`. The quasi-square `k`-omino (`w` columns, full rows, one partial row) has perimeter `2(w + ceil(k/w)) = 2 ceil(2 sqrt k)`. It therefore uses `(4k + P)/2 = 2k + ceil(2 sqrt k)` unit fences, is a legal configuration, and has `k` fields of area 1. Hence `A(2k + ceil(2 sqrt k)) >= k`.

*Proof of the identity.* `(w-1)^2 < k <= w^2`.
- If `k <= w(w-1)`: then `ceil(k/w) = w-1`, and `4(w-1)^2 < 4k <= 4w^2 - 4w < (2w-1)^2`, so `ceil(2 sqrt k) = 2w-1`.
- If `k > w(w-1)`: then `ceil(k/w) = w`, and `4k >= 4w^2 - 4w + 4 > (2w-1)^2`, so `ceil(2 sqrt k) = 2w`. ∎

(T2.C1: `k <= 10^5`; T2.C2: the fence rules, re-checked exactly. Minimality of this perimeter is Harary–Harborth (UNCITED-RECOLLECTION) and is not needed here.)

The least `n` whose record reaches area `k` equals `2k + ceil(2 sqrt k)` for every `k <= 20` except `k = 13` (record 33, polyomino 34) and `k = 17` (42 against 43) (T2.C3). EMPIRICAL remark: these two `k` have the largest rounding loss `ceil(2 sqrt k) - 2 sqrt k` among `k <= 20` (0.789 and 0.754). The next two, `k = 7` (0.708) and `k = 10` (0.675), are not undercut. So the fence problem's "first possible" thresholds are square-root (isoperimetric) thresholds, softened by irregular cells where the polyomino wastes most.

### 2.6 Optional: numerical optimizer (EMPIRICAL; no improvement)

A multi-start SLSQP search over "outer `k`-gon with unit sides plus unit chords ending on sides or on earlier chords" has these ingredients:
- fields computed as half-edge faces;
- every field constrained to area `<= 1`;
- a crossing test (T2.H).

**Results.**
- For `n = 6`, the optimum over all side pairs of the pentagon-plus-chord family is `1.4758530`, with fields `1.000000 + 0.475853`. This matches the record `1.47585+` (Lagache, 2026) to the printed digits, with the chord between non-adjacent sides. So our reading of the rules reproduces a record found by others.
- For `n = 8`, 40 random combinatorial types of "hexagon plus two chords" (6 starts each) reach at best `2.03296 < 2.10306`.

The pictures of the larger records show topologies outside these families. For example, in the `n = 8` record the upper side of the left field continues past a reflex corner of the outer boundary into the interior, and ends on the vertical fence. **No record is claimed or improved.** An improvement would have required an exact re-verification: coordinates, unit lengths, ends on fences, no crossings, field areas. None was triggered.

### 2.7 Typed relations to the repository

| fence item | repository item | type |
|---|---|---|
| fields = bounded faces; `#fields = n + c - V0` | Kuratowski/Euler counting (Kuratowski lane) | REAL, elementary |
| complete grids exist iff `2n+1` composite | a small graph encoding a factorization | REAL (Sundaram) |
| grids optimal only at `4, 7, 12`, irregular beyond | free cycles forced by `\|2^p - 3^a\| = 1`, sporadic beyond (THM-4484) | ANALOGY: both are identity-forced structures at small size, but grid optimality is an extremal statement (best known), not a classification |
| thresholds `2k + ceil(2 sqrt k)` (square-root) | square-sum thresholds (degree `~ (sqrt2-1) sqrt n`) | ANALOGY with a shared archimedean mechanism (square-root density) |
| fence thresholds | provability first at level 7 for `5n+1`, never for `q >= 23` (THM-4481, THM-4486) | not analogous: those thresholds are set by one extreme cycle (Diophantine), not by density |
| area 5 first needs 15 fences, and 15 is the square-sum threshold | - | NUMEROLOGY |
| the word "fence" in the LRC14 corridor-fence hypotheses (HYP-3431 and relatives) | lonely-runner certificates | no relation (name only) |

## 3. T3: small graphs that encode arithmetic

### 3.1 The principle

**Definition.** A *small arithmetic graph* is a family `Gamma_N` of finite graphs on integers or residues, with edges given by an arithmetic relation. A global property `Pi` (Hamiltonian, acyclic, connected, planar, "every cycle contracting") encodes an arithmetic statement at scale `N`. Its *threshold* is the first `N` (or the last `N` without `Pi`) after which `Pi` holds for good.

**Four mechanisms occur in the repository:**

1. **Archimedean / degree (local obstruction).** `Pi` fails while a degree obstruction persists. The thresholds move smoothly with the density of the relation and end for good.
   - Square sums: leaves until 18 and 30; chains from 25, cycles from 32.
   - Cubes: 305 and 473.
   - Triangular numbers: 9, then cycles from 15.
   - Fences: polyomino thresholds.
2. **Diophantine / cycle.** `Pi` is decided by one extreme cycle whose existence is an exponential Diophantine accident.
   - Provability of `qn+-1` sign strategies (THM-4474, THM-4486): `5n+1` first at level 7, where `rho* = 3/7 < log_5 2 = 0.4307`; `7n+-1` not up to level 30; `q >= 23` never (THM-4481).
   - Free and sporadic cycles (THM-4484).

   Degree carries no information here. For constant signs every node of the parity graph has in-degree and out-degree 2 (T3.A1, PROVED: halving is a bijection from the evens and `s -> (qs+-1)/2` from the odds, onto `Z/2^(k-1)`). For general signs the in-degrees are `1, 2, 3`, with as many 3s as 1s (T3.A2). This is PROVED by counting. Each class mod `2^(k-1)` has one even preimage and `0, 1` or `2` odd ones, and the odd preimages total `2^(k-1)`, so classes with two odd preimages (THM-4481's merge pairs) and classes with none are equinumerous. The repository already shows that such stationary or degree information cannot settle `q <= 7` (THM-4486 (F)).
3. **Local (congruence or valuation).** `Pi` fails forever because the graph splits along a class.
   - Odd-square sums split mod 8 (brackets lane).
   - Sums equal to powers of `p` preserve `v_p` (Theorem V).
   - Kohl's `3 | q` criterion (Kuratowski lane).
4. **Topological.** A minor appears: `Q_n` is non-planar from 25; in Althöfer's `U_N`, `K_{3,3}` appears at 52, `K_5` at 68 and the Petersen graph at 104 (Kuratowski lane, CITED).

The new instance below combines 1 and 2. Its mechanism is leaf counting (type 1), but where the leaves sit is decided by the positions of the powers of 2 among the powers of 3, which is type-2 data.

### 3.2 Hitting times (FINITE-EXACT, T3.F)

**Proposition (PROVED + finite).** `delta(Q_n) > (sqrt2 - 1) sqrt n - 2`, since `sqrt(x+n) - sqrt x` decreases in `x`. So `delta(Q_n) >= k` for all `n >= h_k`, with

```text
h_1 = 7, h_2 = 31, h_3 = 71, h_4 = 97, h_5 = 161, h_6 = 241, h_7 = 287, h_8 = 391
```

(the finite part is checked up to the point where the bound takes over).

| family (`n <= N`) | last `n` with a vertex of degree `<= 1` | last `n` without Hamiltonian cycle | last local path obstruction | last `n` without chain |
|---|---|---|---|---|
| squares (120) | 30 | 31 | 18 | 24 |
| triangular (120) | 11 | 14 | 8 | 8 |
| pentagonal (160) | 56 | 56 | 27 | 44 |
| perfect powers (120) | 16 | 16 | 14 | 14 |
| `s^2 + 1` (120) | 23 | 31 | 19 | 19 |
| `s^2 - 1` (120) | 38 | 38 | 19 | 30 |
| `s^2 - 5` (120) | 36 | 36 | 22 | 30 |
| cubes (T1.G: cycles `n <= 600`, chains `n < 400`) | 472 | 472 | 384 | 384 |

The cycle threshold equals the minimum-degree-2 threshold in five of the eight families, and lags it by 1 to 8 in the others (squares 31/30, triangular 14/11, `s^2 + 1`: 31/23). This is the square-type (archimedean) behaviour familiar from the random graph process, where the two hitting times coincide almost surely (ANALOGY; the random-graph theorem is an UNCITED-RECOLLECTION of Bollobás and Komlós–Szemerédi). Nothing like it can happen in the Collatz parity graph, whose degrees are constant.

### 3.3 New instance: the Collatz-alphabet sum graph `C_n`

`C_n = G_S(n)` with `S = {2^p} u {3^a}`: arrange `1..n` so that adjacent sums are powers of 2 or of 3. This is the square-sum problem with the squares replaced by the Collatz clock numbers.

**Theorem C1 (PROVED).** For `n >= 2` the largest power of two `2^p <= n` has degree `<= 1` in `C_n`. Hence `C_n` has no Hamiltonian cycle for any `n >= 3`, and `2^p` is an end of every Hamiltonian path.

*Proof.* The targets in `(2^p, 2^p + n]` lie in `(2^p, 3 * 2^p)`. The only power of two there is `2^(p+1) = 2 * 2^p`, which is excluded as the self-pair, and an interval of ratio `< 3` contains at most one power of 3. ∎ (T3.E1: `n <= 20000`; every witness chain below ends at `2^p`.)

Call `n` *admissible* if `C_n` has no isolated vertex and at most two vertices of degree `<= 1`. This is necessary for a Hamiltonian path.

**Theorem C2 (Window Theorem, PROVED).** For `a >= 2` put `B = 3^(a-1)`. The admissible `n` in `[B+1, 3B+1]` form exactly one interval `W_a`:
- if `(B, 3B)` contains two powers of two `P < 2P`, then `W_a = [max(3B - P, 2P - B), 3B - 1]`;
- if it contains exactly one, then `W_a = [2B, 3B]`.

The admissible `n <= 4` are `2, 3`. So the admissible set is `[2,3] u W_2 u W_3 u ...`, with `W_2 = [5,8]`, `W_3 = [18,27]`, `W_4 = [49,80]`, `W_5 = [162,243]`, `W_6 = [473,728]`, `W_7 = [1319,2186]`, `W_8 = [4374,6561]`, `W_9 = [11491,19682]`, and so on (T3.E3: equal to the direct computation for all `n <= 3^13`).

*Proof* (`a >= 3`; `a = 2` and `n <= 4` by inspection). Let `n in [B+1, 3B+1]` and `x in [n]`, and write `Q(x)` for the unique power of two in `(x, 2x]`. There are four kinds of `x`.

- **(alpha)** `x` is not a power of 2 and some power of 3 lies in `(x, 2x]`. Then `Q(x)` and that power of 3 are two targets in `(x, x+n]`, other than `2x`, so `deg >= 2`.
- **(beta)** `x` is a power of 2 with `3x <= n`. Then `4x` and the power of 3 in `(x, 3x]` give `deg >= 2`.
- **(gamma)** `x` is not a power of 2 and no power of 3 lies in `(x, 2x]`. Then `3^(b-1) <= x < 3^b/2` for the least power `3^b` above `x`, the first two targets above `x` are `Q(x)` and `min(2Q(x), 3^b)`, and `x` is a leaf iff `n < g(x) := min(2Q(x), 3^b) - x <= 2 * 3^(b-1)`. With `x <= n`, this needs `3^(b-1) <= n < 2 * 3^(b-1)`.
- **(delta)** `x` is a power of 2 with `3x > n`. Its only possible partners come from powers of 3 in `(x, x+n]`, since `4x > x + n`.

Now go through the values of `n`.

1. **`2B <= n < 3B`.** No gamma-leaves (no `b` fits). The powers of 3 in `(n/3, 2n]` are `B` and `3B`, so a delta vertex has `deg = [x < B] + [x >= 3B - n] >= 1`: it is never isolated, because `x >= B` implies `x >= 3B - n`. The leaves are therefore among the at most two powers of 2 in `(n/3, n]`. **Every such `n` is admissible.**
2. **`n = 3B`.** The leaves are `3B` itself (its only target is `Q(3B)`) and every power of two in `(B, 3B)` (each has only the target `3B`). **Admissible iff there is exactly one such power.**
3. **`n = 3B + 1`.** The vertices `3B`, `3B+1` (a gamma-leaf; `3^a + 1` is `2` or `4` mod 8, so it is never a power of 2 for `a >= 2`) and a power of two in `(B, 3B)` are all leaves. **Not admissible.**
4. **`B+1 <= n < 2B`.** Gamma-leaves lie in `[B, 3B/2)`. A delta vertex has `deg = [x < B] + [x >= 3B - n] <= 1`, because `3B - n > B`, and it is isolated iff `B <= x < 3B - n`. Two cases.
   - **Case one power `P in (3B/2, 2B)`.** Here `g(x) = 3B - x` on `[B, 3B/2)` and `P/2 in (3B/4, B)` is always a delta-leaf. For `n <= 2B - 2` the three leaves `B, B+1, P/2` exclude `n` (`B + 1 = 3^(a-1) + 1` is not a power of 2 since `a >= 3`). For `n = 2B - 1`, the three leaves `B, P/2, P` exclude it. So `W_a = [2B, 3B]`.
   - **Case two powers `P in (B, 3B/2)` and `2P in (2B, 3B)`.** Here `g(x) = 2P - x` on `[B, P)` and `g(x) = 3B - x` on `(P, 3B/2)`. Let `l = max(3B - P, 2P - B)`.
     - For `l <= n < 2B` there are no gamma-leaves, `P` is a leaf that is not isolated, `P/2` is a leaf iff `n < 3P/2`, and `P/4` is out of range. So `n` is admissible.
     - For `B+1 <= n < l` there are three sub-cases. If `P <= n < 3B - P`, then `P` is isolated. If `n < P`, the leaves `B, B+1, P/2` exclude `n`. If `max(P, 3B-P) <= n < 2P - B`, the leaves `B, P, P/2` exclude it.

   So `W_a = [l, 3B - 1]`. ∎

**Reading.** The number of powers of two in `(3^(a-1), 3^a)` is `floor(a log_2 3) - floor((a-1) log_2 3)`. So which case occurs at level `a` is the Beatty (Sturmian) word of `log_2 3`: `2,1,2,1,2,2,1,2,...` for `a = 2,3,...`. The window left ends are differences of two targets (T3.E3b):

```text
5 = 3^2-2^2 = 2^3-3,  18 = 2*3^2,  49 = 3^4-2^5,  162 = 2*3^4,  473 = 3^6-2^8,  1319 = 2^11-3^6,
4374 = 2*3^7,  11491 = 3^9-2^13,  39366 = 2*3^9,  111611 = 3^11-2^16,  347141 = 2^19-3^11,  1062882 = 2*3^12
```

So the admissibility thresholds of `C_n` are the Collatz gaps `|2^q - 3^b|` adjacent to each level (ANALOGY: the same numbers are the denominators `2^p - 3^a` of cycle equations, THM-4471 and THM-4484; no map from chains to orbits is claimed).

**Gersonides zigzags (PROVED, Lemma Z2 below).** The consecutive targets of `S` are exactly `(1,2), (2,3), (3,4), (8,9)` (Gersonides; T3.D4). They give chains of `C_n` at `n = 1, 2, 3, 7, 8`, and together with `n = 5, 6` these are all Hamiltonian `n <= 8` (T3.E6).

**Hamiltonian paths (FINITE-EXACT, T3.E4).** For `n <= 2200`, `C_n` has a Hamiltonian path exactly for

```text
n in 1-3, 5-8, 18-26, 49-63, 65-66, 68-80, 179-194, 224-243, 473-575, 665-728,
     1319-1418, 1536-1616, 1620-1663, 1703-1713, 1792-1802, 1920-2114, 2120-2186.
```

There are 945 local obstructions and 497 non-local "no"s. Every non-local "no" is decided by the exact forcing solver. In the runner, CP-SAT confirms all of those with `n <= 1000` and a deterministic sample (every 8th) above. A development pass with CP-SAT on all 497 agreed everywhere (1812 s, peak RSS 370 MB; reproduce with `NCP = NHAM` in `t3e`). The windows are only partly Hamiltonian, and the finer sub-windows (for example `W_5 = [162, 243]`, Hamiltonian on `179-194` and `224-243`) come from contradictions in the forcing propagation started at the forced ends (the powers of 2). They have no closed description here (OPEN). Two examples:
- The only non-Hamiltonian admissible `n` in `W_4 = [49,80]` are `64` and `67`.
- `n = 27 = 3^3` is admissible but not Hamiltonian, so the right end `3^a` of a one-power window can fail.

### 3.4 Lacunary forests, valuation splitting, consecutive targets

**Theorem L (PROVED).** If `|S cap (M, 2M)| <= 1` for every `M`, then every `G_S(n)` is a forest.

*Proof.* Take the largest vertex `M` of a cycle. Its two cycle neighbours `y != z` are `< M`, so `M + y` and `M + z` are two distinct targets in `(M, 2M)`, a contradiction. ∎ (T3.B1: powers of 2, 3, 5 and five random lacunary sets, `n <= 1000`.)

**Theorem V (PROVED).** If `S = {p^j}`, every edge `{x, y}` of `G_S(n)` has `v_p(x) = v_p(y)`, because `x + y = p^j` with `x, y < p^j`. So `1` and `p` lie in different components for `n >= p`, and there is no chain for `n >= 2` (`p = 2`) or `n >= 3` (`p = 3`); the chain `1,2` exists at `n = 2` for `p = 3`. For `p = 2` the components are exactly the valuation classes, each a tree: the odd class is rooted at 1 through `x -> 2^(ceil log_2(x+1)) - x`, and the other classes are scaled copies of it. ∎ (T3.C1-C3.)

**Lemma Z2 (PROVED).** If `m, m+1 in S`, then `m, 1, m-1, 2, m-2, ...` is a chain of `G_S(m)` with sums alternating `m+1, m`, and dropping the first term gives a chain of `G_S(m-1)`. ∎ (T3.D1, `m <= 600`.)

For perfect powers, `(8,9)` is the only consecutive pair up to `10^7` (Mihailescu for all sizes, CITED). The first chains are `n = 7` (2 chains) and `n = 8` (1 chain, `4,5,3,6,2,7,1,8` with sums `9,8,9,8,...`), and they are exactly these Catalan zigzags (T3.D3).

### 3.5 The requested hybrids (FINITE-EXACT, "no threshold content")

- **Square sums plus Collatz edges.** Add the edges `{x, T(x)}` (`T(x) = x/2` or `(3x+1)/2`) to `Q_n`. Hamiltonian paths then exist for all `n <= 150` except `7` and `12` (T3.G1). The Collatz edges fill exactly the low-degree holes that made 15 the square-sum threshold, and the hybrid encodes neither problem.
- **`Z/2^k` with `x -> 3x+1` and square-residue sums.** Hamiltonian cycles exist for every `3 <= k <= 9` and not for `k = 2` (T3.G2). The square residues mod `2^k` are dense (about `2^k/6` partners per vertex), so no threshold beyond tiny `k` can arise.

Both are recorded as "no content": the brief suggested them, and they fail for a reason (density swamps arithmetic).

## 4. Failures, corrections and hostile checks

- **Two wrong guesses of this lane, both caught by the code before any claim was printed.**
  - I first expected every sign-strategy parity graph to be 2-in/2-out regular. That is false: a flip `s` and its partner `s*` share targets, so in-degrees are `1, 2, 3` (T3.A2; THM-4481's merge pairs). Only constant signs give the de Bruijn-type regular graph.
  - I first expected the Catalan-type end pair to be a general "first possible" law. It is a law only inside the zigzag family, and it is REFUTED as a general law by `c = 1, 2, -1, ...` and by triangular, pentagonal, cube and perfect-power data (T1.D5, T1.G).
- **Solver hygiene.**
  - Every "no Hamiltonian path/cycle" beyond a hand obstruction is decided by the exact forcing solver (edge forcing plus branching, `lib.HCForcing`) and confirmed by CP-SAT. For the square sums a third method, the exhaustive C DFS with sound pruning (`ham.c`, rules P1-P5), also agrees.
  - The forcing solver agrees with the DFS on 264 graphs (6 families, `n <= 45`, paths and cycles; T0.1). The uniqueness test agrees with exhaustive counts on 213 graphs (T0.2).
  - The square-sum counts reproduce A071983/A071984.
  - Witnesses are always re-verified in Python.
- **An abandoned approach.** An exhaustive C count of chains at the first threshold was too slow for some `c` (millions of chains, e.g. `A116980(35) = 69,328,860` for triangular sums). It was replaced by the exact uniqueness test.
- **Tool incidents, no effect on results.**
  - A `pkill -f` pattern matched its own shell and killed a validation run. It was rerun.
  - A watcher loop used a self-matching `pgrep` and was killed immediately.
  - Neither touched another lane's processes.
- **The records are not optima.** Every fence statement uses the table only as a set of lower bounds. The "trivial" values `4, 7, 12` are grids, but their optimality is not proved here (only `A(3) = sqrt(3)/4` follows from F1).
- **Not done.**
  - The Hales-type asymptotic bound for fences.
  - A proof of uniqueness for the simultaneous Pell system.
  - A closed description of the Hamiltonian sub-windows of `C_n`.

## 5. Reproduction

```text
cd <worktree>
python3 -u 04-computation/experiments/procgen_smallgraph_20260926_run.py \
    > 05-knowledge/results/procgen_smallgraph_20260926.out
```

- **Dependencies:** python3 with numpy, scipy, ortools (CP-SAT, 2 workers, 400 MB cap); a C compiler. The C helper is compiled into `scratch/procgen_smallgraph/bin/` on first use.
- **Output:** every line `[OK] ...` is a `check(...)` that raises on failure.
- **Wall time and memory:** 119 checks, all passing. Wall time 779 s on the 8-core macOS host, run alongside other lanes. Section times: T0 48 s, T1.D 152 s, T1.G 104 s, T2.H 38 s, T3.E 415 s, everything else under 12 s each. Peak RSS 368 MB as self-reported (python, CP-SAT in-process; `/usr/bin/time -l`: 386 MB); the largest C child uses 2 MB. This lane never ran more than two heavy processes: the runner, and the development CP-SAT pass described in 3.3.
- **Script sha256 (LF bytes):**

```text
8e1a1744f279119079e20e1c49e5d8f8df018425c32b14f00835e2fa9199037d  procgen_smallgraph_20260926_lib.py
b7a5cc5c58a0dd87be1f13638bae8c631a82a17a4306b901526e10727f6e3d1e  procgen_smallgraph_20260926_ham.c
79c7a043bb448836e780d7b85935b62d1d984a585667e60f731098648fb7a228  procgen_smallgraph_20260926_t1.py
7a465fb8720fd99db10a5b353c6b671c3ec8f1118cb5c9573b7e9755be09afa0  procgen_smallgraph_20260926_t2.py
c81ffd6ca20402ab15075ef6584ceef53198b0f99520fdd45d8151d05adfac80  procgen_smallgraph_20260926_t3.py
e6e06d815d3407ab104996d42f157faac89537bd5ac39a980c1965b170025f23  procgen_smallgraph_20260926_fence6.py
62b76fc2bbb3beb10c98745f96c536797760ed72eb2cc0c4788463932caec6d5  procgen_smallgraph_20260926_run.py
```

- **Output sha256 (raw bytes of the committed `.out`):** `7ef1341de08d565b5609d3dda75a21a2e3f7ee7f476191e344a7247b00350357`. The section-time lines and the `RUN resources` line vary between runs; everything else is deterministic (fixed seeds, exact solvers).
- **Web sources:** Friedman's page and its images `3..17, 24, 50.gif` were fetched once into `scratch/procgen_smallgraph/web/` (not committed), with UA `Mozilla/5.0 (research; math-repo)`. OEIS entries A304120, A399783 and A116980 were read through the JSON API, with the same UA.
