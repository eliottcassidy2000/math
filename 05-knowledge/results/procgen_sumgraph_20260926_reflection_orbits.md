# Sum graphs as reflection orbits: two and three targets exactly, Pythagorean zigzags, the Collatz-alphabet sub-windows as forcing conflicts (with families proved at every level), square sums at 18-24, and the rotation bridge

**Status block.**

- **PROVED.**
  - *Reflections (Lemma R).* `G_S(n)` is the edge-disjoint union of the reflection matchings `M_t` (`r_t(x) = t - x`), with
    `|M_t| = floor((n - |t-n-1|)/2)`; composites of two reflections are translations.
  - *Theorem A2 (two targets).* `M_s u M_t` is always a union of paths, and it is a Hamiltonian path of `[n]` iff the targets are `{n, n+1}`,
    `{n+1, n+2}`, or an odd pair `{n-1, n+1}`, `{n, n+2}`, `{n+1, n+3}`. This recovers THM-4505's Lemma Z2 (gap 1, Gersonides/Catalan) and adds
    the odd gap-2 case. `C_n` has a two-target chain only for `n in {3, 7, 8}` (`n >= 3`).
  - *Theorem A3 (three targets).* `M_a u M_b u M_c` is a Hamiltonian path iff it is tight (maximum degree 2 and `n - 1` edges, listed
    explicitly in Corollary A3') and `gcd(b-a, c-b) = 1` (or 2 with odd targets). Its alternate vertices follow the rotation by `c - b` on
    `Z/(c - a)`.
  - *Consequences of A3.*
    - THM-4505's zigzag (`2k+1, 4k, 6k+1` at `4k-1`) is the rotation by `2k+1` on `Z/4k`.
    - The whole first square-sum window `{15, 16, 17}` is the single tight triple `(9, 16, 25)`.
    - **Pythagorean zigzags.** Every primitive triple `s^2 + t^2 = u^2` (`s < t`) gives three-square chains of `Q_{t^2-1}` and `Q_{t^2}`,
      with ends `s^2` and `t^2/2` for even `t`. "8 at one end and 9 at the other" is the triple `(3, 4, 5)`.
  - *For `C_n`.*
    - Lemma T: the top layer is a forced two-target zigzag with translation `delta = |3^a - 2P|`.
    - Lemma K: the core of a propagation conflict contains the free end.
    - The Choke Lemma.
  - **Propositions W1-W3.** These are three families of non-Hamiltonian ranges, valid at every level `a` with `rho_a = 2^p/3^(a-1)` in
    `(4/3, 12/7)`, `(4/3, 8/5)`, `(4/3, 3/2)` respectively.
    - Their endpoints are linear in `2^p` and `3^(a-1)`.
    - They account exactly for four failing runs of `W_5` and `W_7` (`162-178`, `1419-1535`, `1714-1791`, `1803-1919`) and for the part
      `207-210` of a fifth.
    - They prove new ones, e.g. `n in [40960, 42664]` (`a = 10`) and `n in [334833, 393215] u [419830, 465904]` (`a = 12`).
    - They occur at levels of density `log_2(9/7) = 0.363`, `log_2(6/5)` and `log_2(9/8)`.
  - *Hand certificates* for every failing run of `W_5` and `W_6`, and for 4 of the 6 failing runs of `W_7`.
- **FINITE-EXACT (own solver, a second code path).**
  - The Hamiltonian set of `C_n` for `n <= 2200` equals THM-4505's list (945 local, 758 paths, 497 non-local NONE).
  - `W_8 = [4374, 6561]`: `C_n` is Hamiltonian for all 2188 `n` in `W_8`, the first fully Hamiltonian window since `W_2`.
  - Propagation completeness: every NONE in these ranges is a unit-propagation contradiction; no search is ever needed to refute.
  - Paths forced by propagation alone, hence unique, occur exactly at `n in {3, 5, 8, 65, 243}` among two-leaf `n <= 2200`.
  - Square sums 18-24 typed:
    - 18 is local;
    - 19 and 20 are first-order conflicts (a choke at 3, and an over-saturated 5);
    - 21 and 22 are a lock at the defect `2`;
    - 23 is where the lock opens (`23 = 25 - 2`);
    - 24 is a global forcing contradiction.
- **EMPIRICAL.** The 17 inner sub-window endpoints are `{2,3}`-unit sums `T - w(s)`, `|w| <= 2`. The short forms are strongly enriched:
  - pure connections `T - s`: 6 observed, about 1.2 expected from the per-window base rates;
  - `|w| <= 1`: 14 observed, about 6.4 expected.
- **CONJECTURE** (section 3.6):
  - (B1) propagation completeness;
  - (B2) the conflict mechanism;
  - (B3) `rho`-scaling, the precise form of a "Beatty rule";
  - (B4) the last run reaches the right end.
- **CITED.** OEIS A090461 (Gerbicz: square-sum chains for all `n >= 25`); Weyl equidistribution of `a log_2 3 mod 1`; THM-4505 (all
  facts about `C_n` windows and `Q_n` degrees used here); THM-4484 (Gersonides, cycle denominators).
- **ANALOGY.**
  - The interval-exchange / linear-involution dictionary. The Keane-type "no connection => minimal" is an UNCITED-RECOLLECTION, used in no
    proof.
  - The top translations `delta = |2^q - 3^a|` equal the Collatz cycle denominators as numbers; there is no map from paths to orbits.
- **OPEN.**
  - A proof of propagation completeness.
  - General families for the deeper (second- and third-order) conflicts.
  - Any positive (Hamiltonicity) statement for all levels.
  - The continuous limit in `rho`.
  - Whether `C_{3^a - 1}` / `C_{3^a}` is always Hamiltonian.

Collatz itself is not addressed. No novelty is claimed for the elementary facts in sections 1-2. Session `collatz-procgen-20260922`, lane
`sumgraph` (2026-09-26).
- Scripts: `04-computation/experiments/procgen_sumgraph_20260926_{solver,theory,run}.py`.
- Output: [procgen_sumgraph_20260926.out](procgen_sumgraph_20260926.out).

## 0. The four tasks in one paragraph each

**A (two and three targets).**
- *Two targets.* The orchestrator's picture is right, and it can be made exact. Two reflection matchings never close a cycle, because two
  reflections compose to a translation. So their union is a linear forest, and it is a single path only for the targets `{n, n+1}` and
  `{n+1, n+2}` (Lemma Z2, the Catalan/Gersonides zigzag) or the odd gap-2 pairs.
- *Three targets.* The union is a single path iff it is tight and the rotation it induces, by `c - b` on `Z/(c - a)`, is transitive.
  THM-4505's zigzag is the rotation by `2k + 1` mod `4k`. The square-sum chain at 15 is the rotation by 9 mod 16.
- *Pythagorean zigzags.* Since `9 + 16 = 25`, the same statement produces three-square chains from every primitive Pythagorean triple
  (Corollary A3''').

**B (`C_n` sub-windows).**
- *The top layer.* It is always the forced zigzag of the two largest targets `2P` and `3^a`, a translation by the gap `|3^a - 2P|`
  (Lemma T).
- *Why a run starts and stops.* The sub-windows start and stop when a *choke* is born or dies: a defect (a power of 2 fixed by `r_{2d}`, or
  a power of 3 at the end of a reflection domain) all of whose partners but one are taken by rigid vertices, mostly top vertices of the
  zigzag.
  - A choke is fatal when both ends are forced, or when two chokes need different free ends.
  - It is born when the last saturating top vertex `T - y` enters, and dies when the defect gains a top partner.
  - So every endpoint is a `{2,3}`-unit sum such as `1419 = 3^7 - 2^10 + 2^8`.
- *Every level.* For the defects `P/4` and `P/8` this is proved at every level (W1-W3), with the level condition a window for `rho_a`.
  This is the exact form of the "fractal recursion across levels", and it gives new proved non-Hamiltonian ranges at `a = 10, 12`.
- *What stays conjectural.* The full pattern is decided by unit propagation alone, computed to 2200 and on `W_8`. The deeper conflicts are
  not yet in closed form.

**C (square sums).**
- 18 fails locally (three leaves).
- 19 and 20 fail by small forcing conflicts (a choke at 3, and an over-saturated 5).
- 21 and 22 fail by a lock at the defect 2, which opens at `23 = 25 - 2`.
- 24 is the only global failure.
- From 25 on the new target 49 adds a matching.

**D (continuous bridge).**
- *Rotations.* For three targets the bridge is a theorem: the discrete path is a rational rotation orbit, and Hamiltonicity is the
  coprime-rational shadow of minimality (irrational rotation number).
- *Top layer.* For `C_n` the top layer gives one Rauzy-type induction step.
- *Prediction.* The continuous picture predicts that status changes only at *connections*, where a singular orbit hits the moving domain
  end `T - n`. The data match: every observed endpoint is `T - w(s)`, and the pure connections are strongly enriched.
- *Collatz.* The gaps `|2^q - 3^a|` enter as translation lengths (proved). Their identity with Collatz cycle denominators is an analogy.

## 1. Setting

For a target set `S` and `n >= 1`, `G_S(n)` has vertices `[n] = {1..n}` and an edge `{x, y}` (`x != y`) iff `x + y in S`.
For a target `t` let `r_t(x) = t - x`. On `[n]` it is defined on the interval `D_t = [max(1, t-n), min(n, t-1)]`, and
`M_t = {{x, t-x} : x in D_t, 2x != t}` is a matching. `G_S(n)` is the edge-disjoint union of the `M_t`, `t in S`.

**Lemma R (PROVED).**
- (R1) `|M_t| = floor((n - |t - n - 1|)/2)` for `3 <= t <= 2n-1` (and `M_t` is empty otherwise).
- (R2) `r_{t2} o r_{t1}` is the translation `x -> x + (t2 - t1)`. So along any walk `x_0, x_1, ...` of `G_S(n)` with sums `t_i = x_{i-1} + x_i`,
  `x_{i+1} = x_{i-1} + (t_{i+1} - t_i)`: the even-indexed and the odd-indexed vertices move by translations by target differences.
- (R3) The group generated by `r_t, t in T` acts on `Z` with orbits `(x + gZ) u (t_0 - x + gZ)`, `g = gcd(t - t_0 : t in T)`. Every connected component of
  `U_{t in T} M_t` lies in one such orbit.
- (R4) The reversal `x -> n + 1 - x` maps `M_t` onto `M_{2n+2-t}`.

*Proof.* (R1) `x` runs over `[max(1, t-n), ceil(t/2) - 1]`; for `t <= n+1` this has `floor((t-1)/2)` elements, for `t >= n+1` it has `n - floor(t/2)`, and both equal
`floor((n - |t-n-1|)/2)`. (R2) `t2 - (t1 - x)`. (R3) even words are translations by elements of `gZ`, odd words are `x -> t_0 - x + gZ`. (R4) `(n+1-x) + (n+1-y) = 2n+2-(x+y)`. ∎

## 2. Task A: unions of two and three reflection matchings

### 2.1 Two targets

**Theorem A2 (PROVED).** Let `n >= 3`, let `s < t` be targets in `[3, 2n-1]`, `d = t - s`, and `G = M_s u M_t` on `[n]`.
- (a) `G` has no cycle, so it is a disjoint union of paths, with `n - |M_s| - |M_t|` components.
- (b) `G` is a Hamiltonian path of `[n]` iff
  - `d = 1` and `s in {n, n+1}` (targets `{n, n+1}` or `{n+1, n+2}`), or
  - `d = 2`, `s` odd, and `s in {n-1, n, n+1}` (targets `{n-1, n+1}`, `{n, n+2}` or `{n+1, n+3}`).

*Proof.* (a) Each vertex has at most one edge in each matching. A cycle would alternate `s`- and `t`-edges, and by (R2) its vertices two steps apart satisfy
`x_{2i+2} = x_{2i} + d` (orienting the cycle so that it starts with an `s`-edge), so `x_{2k} = x_0 + kd != x_0`. Hence `G` is a forest of paths and has
`n - |E|` components.

(b) By (a), `G` is a Hamiltonian path iff `|M_s| + |M_t| = n - 1`. Connectivity needs `d <= 2` (for `d >= 3` the vertices `1, 2, 3` lie in three classes
mod `d`, and by (R3) a component meets at most two), and for `d = 2` it needs `s` odd (else `x` and `t - x` have the same parity). Now use (R1).
- `d = 1`. If `t <= n+1`: `|M_s| + |M_t| = floor((s-1)/2) + floor(s/2) = s - 1`, which is `n - 1` iff `s = n`. If `s >= n+1`: the sum is
  `(n - floor(s/2)) + (n - floor((s+1)/2)) = 2n - s`, which is `n - 1` iff `s = n+1`.
- `d = 2`, `s = 2m+1`. If `t <= n+1`: the sum is `m + (m+1) = s`, so `s = n - 1`. If `s >= n+1`: `(n-m) + (n-m-1) = 2n - s`, so `s = n+1`. Otherwise
  `s = n` and the sum is `floor((n-1)/2) + (n - floor((n+2)/2)) = m + m = n - 1`. ∎

(Runner A.1: the brute-force list of two-target Hamiltonian unions for `n <= 150` equals the prediction.)

**Corollary A2' (Lemma Z2 and the Gersonides chains, recovered).** The case `d = 1` is exactly Lemma Z2 of THM-4505 (targets `m, m+1` give a chain of `[m]` and
of `[m-1]`) and says there is nothing else with two consecutive targets. For `C_n` (targets `{2^p} u {3^a}`) the consecutive pairs are `(3,4)` and `(8,9)` (Gersonides),
and no two odd targets `>= 3` differ by 2, so for `n >= 3` a two-target Hamiltonian path of `C_n` exists exactly for `n in {3, 7, 8}`. For the squares (gaps `>= 5`
between squares `>= 4`) it never exists: every square-sum chain with `n >= 3` uses at least three targets.

### 2.2 Three targets: tight unions are rotation orbits

**Theorem A3 (PROVED).** Let `n >= 3`, targets `a < b < c` in `[3, 2n-1]`, `G = M_a u M_b u M_c` on `[n]`, and `g = gcd(b-a, c-b)`. Then `G` itself is a
Hamiltonian path of `[n]` iff
- (i) every `x` in `D_a cap D_b cap D_c = [max(1, c-n), min(n, a-1)]` is a midpoint, `2x in {a, b, c}` (equivalently `G` has maximum degree `<= 2`);
- (ii) `|M_a| + |M_b| + |M_c| = n - 1`, i.e. `sum_{t in {a,b,c}} floor((n - |t-n-1|)/2) = n - 1`;
- (iii) `g = 1`, or `g = 2` and `a` is odd.

Moreover, given (ii), condition (i) holds iff one of the following holds:
- `c - a >= n` (*generic case*);
- `c - a = n - 1` and `b = 2a - 2` (*one midpoint*: the triple point `x_0 = a - 1 = c - n = b/2`);
- `n in {6, 7}` and `(a,b,c) in {(4,6,n+2), (n, 2n-4, 2n-2)}`.

*Proof.* The domains are intervals `[l_t, r_t]` with `l_t, r_t` non-decreasing in `t`, so `D_a cap D_c = [l_c, r_a]` is contained in `D_b`; a vertex has
degree 3 iff it lies in all three domains and is not a midpoint. This gives (i) <=> `Delta(G) <= 2`.

(=>) A Hamiltonian path has maximum degree 2 and `n - 1` edges, so (i), (ii) hold. By (R3) every component lies in `(x + gZ) u (a - x + gZ)`. For `g >= 3` the
vertices `1, 2, 3` are in three classes; for `g = 2` and `a` even, `a - x = x` mod 2. Either way `G` is disconnected, so (iii) holds.

(<=) By (i), `G` is a disjoint union of paths and cycles, and by (ii) it has `n - (n-1) + #cycles = 1 + #cycles` components. So it suffices to show that `G` has
no cycle.

*Generic case `c - a >= n`.* Then `D_a cap D_c` is empty. A vertex of degree 2 therefore lies in `D_b` and in exactly one of `D_a, D_c`, and it is not `b/2`,
so its two edges are one `b`-edge and one `a`- or `c`-edge. Along a cycle `x_0, x_1, ..., x_{2k-1}` the edges alternate `b, e_1, b, e_2, ...` with
`e_i in {a, c}`, and by (R2) `x_{2i+2} = x_{2i} + (e_i - b)`. Closing up: `k_c (c-b) = k_a (b-a)`, where `k_a, k_c` count the `a`- and `c`-edges. Both are
positive, and with `(b-a)/g`, `(c-b)/g` coprime we get `k_a >= (c-b)/g`, `k_c >= (b-a)/g`, so the cycle has length `2k >= 2(c-a)/g >= 2n/g`. For `g = 1`
this exceeds `n`. For `g = 2` the cycle would be Hamiltonian, with `n` edges, while `G` has `n - 1` edges. So there is no cycle.

*One-midpoint case.* Here `D_a cap D_c = {x_0}` and `x_0 = b/2`, so `x_0` carries an `a`-edge and a `c`-edge, and every other vertex of degree 2 carries one
`b`-edge. Since `b = 2a - 2` is even, (iii) forces `g = 1` (for `g = 2`, `a` would be even). A cycle avoiding `x_0` has length `>= 2(c-a) = 2n - 2 > n`
as above. A cycle through `x_0` has edge word `c, (b, e_1), ..., (b, e_j), b, a` (length `2j + 3`), and following it with (R2) from `x_0` gives
`2x_0 = a + c - b + sum_i (e_i - b)`, i.e. `(k_c + 1)(c - b) = (k_a + 1)(b - a)`. With `g = 1`: `k_a + 1 >= c - b` and `k_c + 1 >= b - a`, so the length is
`>= 2(c - a) - 1 = 2n - 3 > n` for `n >= 4`, and for `n = 3` it would be a Hamiltonian cycle with more than `n - 1` edges.

*The remaining (i)-cases.* If `D_a cap D_b cap D_c` has two or three points, each is a midpoint. Solving `{2x, 2x+2} subset {a,b,c}` (resp. `{2x, 2x+2, 2x+4}`)
with the interval `[c-n, a-1]` gives four families:
- `(4, 6, n+2)`;
- its mirror `(n, 2n-4, 2n-2)` under (R4);
- `(4, 5, 6)` at `n = 4`;
- `(6, 8, 10)` at `n = 7`.

The last two fail (ii). For the first two, (ii) reads `3 + floor((n-1)/2) = n - 1`, i.e. `n in {6, 7}`. These four instances are checked
directly (runner A.2). ∎

**The rotation behind it (PROVED, generic case).** On the vertices carrying a `b`-edge, the "return to the `b`-matching" map is
`F(x) = r_e(r_b(x)) = x + (c - b)` or `x - (b - a)`, according to whether `r_b(x)` lies in `D_c` or `D_a`. Both are `x + (c - b)` modulo `c - a`, because
`-(b - a) = (c - b) - (c - a)`. So the alternate vertices of every component follow the rotation `rho: x -> x + (c-b) mod (c-a)`, cut open at the degree-1
vertices, and (iii) says that this rotation is transitive (`g = 1`), or has two cosets exchanged by the odd reflection `r_b` (`g = 2`, `b` odd).

**Corollary A3' (the explicit tight triples, PROVED).** In the generic case (`c - a >= n`), condition (ii) holds iff
- `n` odd: `|b - n - 1| <= 1`, and `c - a = n`, or `c - a = n + 1` with `a = n (mod 2)`;
- `n` even: `b = n + 1` and one of `c - a = n + 1`; `c - a = n` with `a` even; `c - a = n + 2` with `a` odd. Or `1 <= |b - n - 1| <= 2`, `c - a = n` and
  `a` odd.

*Proof.* Put `u_t = t - n - 1`. In the generic case `u_a <= -2 < 2 <= u_c` and `|u_a| + |u_c| = c - a =: m >= n`, and (R1) gives
`floor((n-|u_a|)/2) + floor((n-|u_c|)/2) <= floor((2n-m)/2) <= floor(n/2)` and `floor((n-|u_b|)/2) <= floor(n/2)`. Reaching `n - 1` leaves slack at most one
(n even) or none (n odd), and the parity of `n - |u_a|`, `n - |u_c|` decides each floor. ∎ (Runner A.3 compares this table with (ii) for all `n <= 150`.)

**Corollary A3'' (THM-4505's zigzag and the first square-sum window).**
- The Zigzag Lemma of THM-4505 (`2k+1, 4k, 6k+1` at `n = 4k - 1`) is the tight triple with `c - a = b = n + 1`, `n` odd, `a` odd, and
  `g = gcd(2k-1, 2k+1) = 1`. Its alternate vertices follow the rotation by `2k + 1` on `Z/4k`. For `k = 4` (squares `9, 16, 25`) this is the drift
  `+7, -9` of the note's Proposition 2A.
- For three consecutive squares `((j-1)^2, j^2, (j+1)^2)` (`j >= 2`), `M u M u M` is a Hamiltonian path of `[n]` iff `j = 4` and `n in {15, 16, 17}`:
  `n = 15` (generic, `c - a = n + 1`), `n = 16` (generic, `c - a = n`) and `n = 17` (the one-midpoint case: `b = 16 = 2a - 2`, `x_0 = 8`). Here
  `g = gcd(2j-1, 2j+1) = 1` always. Tightness needs `|j^2 - n - 1| <= 2` and `n <= c - a + 1 = 4j + 1`, which is impossible for `j >= 5`; `j = 2, 3`
  are excluded by the runner (A.4, all `j <= 40`, `n <= 2000`). Since `Q_n` for `n <= 17` has no target besides `4, 9, 16, 25`, **the whole first window
  `{15, 16, 17}` of the square-sum problem is the single tight triple `(9, 16, 25)`**, and it closes at 18 because the triple intersection
  `D_9 cap D_16 cap D_25 = [7, 8]` then contains `7`, which is not a midpoint (so `7` has degree 3).

**Corollary A3''' (Pythagorean zigzags, PROVED).** Let `s^2 + t^2 = u^2` be a primitive Pythagorean triple with `s < t`, and put
`(a, b, c) = (s^2, t^2, u^2)`. Then `M_a u M_b u M_c` is a Hamiltonian path of `[n]` iff
- `n in {t^2 - 1, t^2}`, or
- `n = t^2 + 1` and `t^2 = 2s^2 - 2`, which among Pythagorean triples happens only for `(3, 4, 5)` (runner, `t <= 10^4`).

So every such `Q_{t^2-1}` and `Q_{t^2}` contains a Hamiltonian path that uses only three squares. Its alternate vertices run through one orbit
of the rotation by `s^2` on `Z/t^2`. At `n = t^2 - 1` its ends are:
- `s^2` and `t^2/2` if `t` is even;
- `s^2/2` and `s^2` if `t` is odd.

*Proof.* `c - a = t^2 = b`.
- (iii): `g = gcd(t^2 - s^2, u^2 - t^2) = gcd(t^2 - s^2, s^2) = gcd(t^2, s^2) = 1`.
- (i): the generic case needs `n <= t^2`. The one-midpoint case (`n = t^2 + 1`) needs `b = 2a - 2`.
- (ii), by Corollary A3': at `n = t^2 - 1` we have `b = n + 1 = c - a`, so it holds for `n` even, and for `n` odd because then `s` is odd.
  At `n = t^2` we have `|b - n - 1| = 1` and `c - a = n`, with `a` odd when `n` is even. At `n = t^2 - 2`, `c - a = n + 2` would need
  `b = n + 1`, which fails.
- The ends: `x = s^2` lies only in `D_b` (because `D_c = [s^2 + 1, n]`); the other end is the midpoint `t^2/2` of `r_b`, resp. `s^2/2` of
  `r_a`. ∎

(Runner A.4b checks all 18 primitive triples with `t <= 100`.) This is the cleanest reading of THM-4505's "8 at one end and 9 at the other":
- `9 + 16 = 25` is the triple `(3, 4, 5)`;
- the chain is the rotation by `9` on `Z/16`, i.e. rotation number `(3/4)^2`;
- the ends are `s^2 = 9` and `t^2/2 = 8`.

The next triples give three-square chains of `Q_143` (ends `25, 72`), `Q_224` (ends `32, 64`), `Q_575` (ends `49, 288`), and so on. What is
special about `n = 15` is only that it is the *first* such `n`; for `n <= 14` there is no square-sum chain at all (THM-4505).
The Pell reading of THM-4505 asks for the zigzag *shape* `(2k+1, 4k, 6k+1)`, which always satisfies `a + b = c`. Its Pythagorean instances are
the solutions of the simultaneous Pell system, and there is only one (CITED there).

## 3. Task B: the fine structure of `C_n`

`C_n = G_S(n)` with `S = {2^p} u {3^a}`. Write `P = 2^p <= n < 2P` and `Q = 3^b <= n < 3Q`. In a window `W_a` (Window Theorem C2 of
THM-4505) we have `Q = 3^(a-1)` and `3Q = 3^a` (except at `n = 3^a` itself). THM-4505 (C1) shows that `P` is a leaf, hence an end of every
Hamiltonian path.

### 3.1 The top layer is a forced two-target zigzag

**Lemma T (PROVED).** Let `n >= 3` and `M = max(P, Q)`.
- Every target in `(M, 2n-1]` is `2P` or `3Q`.
- Hence every *top vertex* `x in (M, n]` has exactly the neighbours `2P - x` and (if `1 <= 3Q - x <= n`) `3Q - x`.
- In every Hamiltonian path, every top vertex of degree 2 that is not an end uses both of its edges.
- So the top layer lies on the forced zigzag `M_{2P} u M_{3Q}`, a two-target union (Theorem A2(a)). Its alternate vertices move by the
  translation `delta = |3Q - 2P|`, and the top vertices fall into forced chains along the residue classes mod `delta`.

*Proof.* `(M, 2n-1]` lies in `(P, 4P)`, which contains the single power of two `2P`, and in `(Q, 6Q)`, which contains the single power of
three `3Q`. For `x > M` the neighbours are `t - x` with `t in (x, x+n] subset (M, 2n]`, `t != 2x`. Here `2x != 2P` because `x != P`, and
`2x != 3Q` by parity. Finally `2P - x` lies in `[2P - n, P) subset [1, n]`. ∎ (Runner B.1: the target statement for `n <= 20000`, the
neighbour lists for `n <= 1500`.)

In `W_a` the top translations are `delta = |3^a - 2P|`:

| window | `delta` |
|---|---|
| `W_3` | 5 |
| `W_4` | 17 (`n < 64`), 47 |
| `W_5` | 13 |
| `W_6` | 217 (`n < 512`), 295 |
| `W_7` | 139 (`n < 2048`), 1909 |
| `W_8` | 1631 |

These are the gaps `|2^q - 3^a|`.

### 3.2 Conflicts, cores, and where the free end can be

A vertex is *rigid* if it has degree 2 and is not a prescribed end: both its edges are in every Hamiltonian path, unless it is the free end.
The unit rules are the ones in the solver's docstring:
- (R1) a vertex whose remaining edges number its need uses all of them;
- (R2) a vertex with as many path edges as its need loses its other edges;
- (R4) no forced segment closes;
- contradictions: too few edges (LOW), too many (OVER), a forced cycle, or a forced end-to-end segment that misses vertices.

**Lemma K (PROVED).** Fix the forced ends. Give need 2 to every other vertex and run the rules. Suppose a contradiction is derived, and let
`K` (the *core*) be the set of vertices whose need is used in its derivation. Then every Hamiltonian path has its free end in `K`.

*Proof.* The derivation is a proof that no Hamiltonian path realises the need vector. Each step is sound: if a path `H` realises the needs,
then every forced edge is in `H`, every deleted edge is not, and the contradiction cannot occur. Each step uses only the needs of vertices in
`K`. For ends `{P, e}` with `e` outside `K`, the needs agree with the base run on `K`, so the same proof excludes them. ∎

**Choke Lemma (PROVED, first order).** Say a vertex `v` is *saturated away from* `d` if `v` has two rigid neighbours other than `d` (one, if
`v` is a prescribed end). Then, unless the free end is one of those rigid vertices, the path uses them at `v` and never uses the edge `vd`.
- **Choke.** If `d` is not a prescribed end and at most one neighbour of `d` is not saturated away from `d`, then the free end lies in
  `R(d) = {d} u {the saturating rigid vertices}`. If no neighbour of `d` is left and at least one end is forced (as in `C_n`), there is no
  Hamiltonian path at all: freeing one neighbour would make `d` a third end.
- **Over-saturation.** If `v` has `need(v) + 1` rigid neighbours, the free end is one of them.

Consequently: with two forced ends, any choke or over-saturation refutes Hamiltonicity. With one forced end, conflicts whose resolving sets
are disjoint refute it.

*Proof.* A rigid vertex that is not the free end uses both its edges, and a vertex uses at most `need` path edges. ∎

Runner B.5:
- on every admissible `n` checked, first-order refutations never occur at a Hamiltonian `n`;
- of the 497 non-local NONE verdicts up to 2200, 435 already follow from first-order conflicts, and 62 need deeper propagation.

The *defects* where chokes happen are the singular points of the reflection family:
- the powers of two `d = 2^j`, the half-targets, which are fixed by `r_{2d}`;
- the powers of three `3^j`, which lie just past the end `3^j - 1` of the domain of `r_{3^j}`.

Both are typically local minima of the degree `|S cap (x, x+n]| - [2x in S]` (Lemma 1 of THM-4505).

(A choke can also sit at an ordinary vertex. Example: for `615 <= n <= 647` the propagation also finds conflicts centred at `81 = 3^4` and at the ordinary vertices `82, 83, ...`, inside the longer
failing run `[576, 664]` that the defect `64` already explains. Defects are simply where chokes live longest.)

### 3.3 Exact results (FINITE-EXACT, second code path)

**Theorem B (FINITE-EXACT).**
- (a) For `n <= 2200` the Hamiltonian-path set of `C_n` is exactly THM-4505's list
  `1-3, 5-8, 18-26, 49-63, 65-66, 68-80, 179-194, 224-243, 473-575, 665-728, 1319-1418, 1536-1616, 1620-1663, 1703-1713, 1792-1802, 1920-2114, 2120-2186`
  (945 local obstructions, 758 paths, 497 non-local NONE). This is re-derived by a second code path (own solver, section 6), not by
  THM-4505's.
- (b) On `W_8 = [4374, 6561]`: `C_n` has a Hamiltonian path for **every** one of the 2188 admissible `n`, and every path is verified. `W_8` is the first window since `W_2 = [5, 8]` with no failure at all.
- (c) **Propagation completeness.** Every non-local NONE in (a) and (b) is reached by unit propagation alone: the forced leaves, plus the free
  end ranging over the base core of Lemma K, with zero branching nodes. So in this range every failure of Hamiltonicity has a finite forcing
  certificate, a core in the sense of Lemma K.
- (d) Among two-leaf Hamiltonian `n <= 2200`, propagation alone forces the whole path, which is then unique, exactly for
  `n in {3, 5, 8, 65, 243}`. In particular `C_{243} = C_{3^5}` has exactly one Hamiltonian path.

### 3.4 Why each sub-window starts and stops (PROVED certificates + FINITE-EXACT)

Notation: `r_t(x) = t - x`. "Saturated by `w, w'`" means that `w, w'` are rigid neighbours. Every hypothesis below is re-checked by the runner
for every `n` of the stated range (B.4).

**`W_5 = [162, 243]`, `P = 128`, top targets `256, 243`, `delta = 13`.** Hamiltonian exactly on `179-194, 224-243`.

| run | status | certificate |
|---|---|---|
| `162-178` | NONE | Ends forced `{64, 128}`: `64 = P/2` has the single partner `17` until `179`. **Choke at the defect 32:** `N(32) = {49, 96}`, and `96 = r_128(32)` is saturated by the top vertices `147 = r_243(96)`, `160 = r_256(96)`. So `32` would be an end. |
| `179-194` | PATH | `179 = r_243(64) = 3^5 - 2^6` gives `64` a partner. Only `128` stays forced, and the free end resolves the choke: survivors `{32, 147, 160}`. |
| `195-206` | NONE | **Second-order choke at the defect 16**, born at `195 = r_243 r_64(16)`. `95` is saturated by `148, 161`; so `33` (`N = {31,48,95}`) is rigidly joined to `31, 48`; `48` (`N = {16,33,80,195}`) is saturated by `33` and the new top vertex `195`; `65` is saturated by `178, 191`. Hence `16` (`N = {11,48,65,112}`) is joined to `11` and `112`, and `112`, also saturated by `131, 144`, is over-full. The resolving vertices are disjoint from `{32, 147, 160}`. |
| `207-210` | NONE | `49 = r_81(32)` is now saturated too (`194`, `207 = r_256(49)`): `32` has no usable partner. |
| `211-223` | NONE | `32` gains the rigid top partner `211 = r_243(32)`. Two first-order chokes with disjoint resolving sets: `R(32) = {32,147,160,194,207}` and `R(16) = {16,131,144,178,191,195,208}`. |
| `224-243` | PATH | `224 = r_256(32) = 2^8 - 2^5`: `32` gains a second top partner, so the choke at `32` dissolves and the free end resolves `R(16)`. At `243 = 3^5` a second leaf appears and propagation forces the unique path. |

**`W_6 = [473, 728]`, `P = 256` then `512`, `delta = 217` then `295`.** Hamiltonian exactly on `473-575, 665-728`.

| run | status | certificate |
|---|---|---|
| `576-664` | NONE | Ends forced `{256, 512}`. **Choke at the defect 64:** `N(64) = {17, 179, 192, 448}`. `179 = r_243(64)` is saturated by `333` (rigid) and `550` (top); `192 = r_256(64)` by `320, 537`; `448 = r_512(64)` by `281` and the top vertex `576 = r_1024(448)`. Born at `576 = 2^10 - 2^9 + 2^6`, dies at `665 = r_729(64) = 3^6 - 2^6`. |

**`W_7 = [1319, 2186]`, `P = 1024` then `2048`, `delta = 139` then `1909`.** Hamiltonian exactly on
`1319-1418, 1536-1616, 1620-1663, 1703-1713, 1792-1802, 1920-2114, 2120-2186`.

| run | status | certificate |
|---|---|---|
| `1419-1535` | NONE | Ends forced `{512, 1024}`. **Choke at 256:** `N(256) = {473, 768}`, and `768 = r_1024(256)` is saturated by `1280 = r_2048(768)`, `1419 = r_2187(768)`. Born `1419 = 3^7 - 2^10 + 2^8`. |
| `1536-1616` | PATH | `1536 = r_2048(512) = 2^11 - 2^9` frees the leaf `512`; the free end resolves `R(256) = {256, 1280, 1419}`. |
| `1617-1619` | NONE | A third-order choke centred at `81 = 3^4` (saturation tree through `162, 175, 431, 350, 567, 457, 272, ...`); certificate by propagation only. Born `1617 = r_2048 r_512(81)`, dies `1620 = r_2187 r_729 r_243(81)`, when `567` gains a third partner. |
| `1664-1702` | NONE | **Second-order choke at 128**, disjoint from `R(256)`. `679` is saturated by `1369, 1508`, so `345` (`N = {167,384,679}`) is rigid; `384 = r_512(128)` is saturated by `345` and `1664 = r_2048(384)`; `601`, `896` are saturated by top pairs. Born `1664 = 2^11 - 2^9 + 2^7`, dies `1703 = r_2048(345)`. |
| `1714-1791` | NONE | `473 = r_729(256)` is saturated by `1575`, `1714 = r_2187(473)`, and `768` by `1280, 1419`: `256` has no usable partner unless the free end frees one, and then `256` must itself be an end (three ends). Born `1714 = 3^7 - 3^6 + 2^8`, dies `1792 = r_2048(256) = 2^11 - 2^8`. |
| `1803-1919` | NONE | First-order chokes at `128` (`384` saturated by `1664` and `1803 = r_2187(384)`, `601` by `1447, 1586`, `896` by `1152, 1291`) and at `256` (`1792` rigid, `473`, `768` saturated). Resolving sets `{128,1152,1291,1447,1586,1664,1803}` and `{256,1280,1419,1575,1714}` are disjoint. Born `1803 = 3^7 - 2^9 + 2^7`, dies `1920 = r_2048(128) = 2^11 - 2^7`. |
| `2115-2119` | NONE | `n >= 2048`, ends forced `{1024, 2048}`. A conflict with a 32-vertex core around `32` and `72 = r_128 r_64(8)`; certificate by propagation only. Born `2115 = 3^7 - 72`, dies `2120 = 3^7 - 67 = 2^11 + 72`, symmetric about `(2^11 + 3^7)/2`. |

**`W_8 = [4374, 6561]`, `P = 4096`, top targets `8192, 6561`, `delta = 1631`, `rho = 1.873`.** Hamiltonian at **every** one of its 2188
admissible `n` (FINITE-EXACT, own solver; every path verified). Since `rho > 12/7`, W1-W3 do not apply.
- There is still a choke for `5025 <= n <= 5214`: the base propagation finds one at the defect `512 = P/8` (runner B.3c).
  - Its partners `1536`, `3584`, `1675` are saturated by the rigid pairs `(2560, 5025)`, `(4608, 2977)`, `(2421, 4886)`.
  - Three of these six (`2560, 2977, 2421`) are rigid vertices below the top layer; they have degree 2 because `8192 - x > n`.
  - It is a single choke, so the free end resolves it. In the development scan the paths found in this range end at `512` (189 of 190
    values) or at the saturator `2560` (once).
  - It dissolves at `5215 = 8192 - 2977`, when `2977` gains a third partner.
- So the window is Hamiltonian throughout, even though the defect mechanism is active.

**What the certificates share.**
- Every failing run is the lifetime of one or two conflicts centred at a defect (`2^j` or `3^j`). Each conflict is a choked vertex whose
  partners are eaten by rigid vertices, most of them top vertices of the forced zigzag.
- A conflict is born when the last saturating top vertex `T - y` enters (`n = T - y`, with `y` in the saturation tree). It dies when the
  defect, or a rigid vertex of its tree, gains a top partner.
- A second forced leaf (`P/2` before its first top partner) turns a single conflict into a failure. Two conflicts with disjoint resolving
  sets always fail.
- So every endpoint is an alternating sum of a top target and at most three further powers of 2 and 3, a `{2,3}`-unit sum.

### 3.5 The same conflicts at every level: three families with endpoints linear in `2^p` and `3^(a-1)` (PROVED for all `a`)

Let `a >= 3`, `B = 3^(a-1)`, let `P` be the power of two in `(B, 2B)`, and put `rho = P/B in (1, 2)`. Then `rho = 2^(1 - {(a-1) log_2 3})`.
Every vertex named below is an integer, since `P >= 16`.

**Proposition W1 (two forced ends and a choke at `P/4`).** If
`max(5P/4, 3B - 3P/4) <= n < min(3P/2, 3B - P/2)`,
then `C_n` has no Hamiltonian path. The range is nonempty iff `4/3 < rho < 12/7`.

*Proof.* In the range, `P <= n < 2P` and `B <= n < 3B`, so Lemma T applies with the top targets `2P, 3B`.
- **`P` is a leaf.** Its partners are `t - P` for targets `t in (P, P + n]` other than `2P`. Only `t = 3B` qualifies (`3P > n`), so `P` has
  the single partner `3B - P`.
- **`P/2` is a leaf.** Its only partner is `B - P/2`: `3P/2 > n` and `3B - P/2 > n`, and `P` is its double. So the ends are forced to be
  `{P/2, P}`.
- **`N(P/4) = {3P/4, B - P/4}`.** The candidates `7P/4` and `3B - P/4` exceed `n`, and `B/3 < P/4` because `rho > 4/3`.
- **`3P/4` is saturated.** Its partners `5P/4 = r_{2P}(3P/4)` and `3B - 3P/4 = r_{3B}(3P/4)` are top vertices (the second because
  `rho < 12/7`). Both are `<= n` and have degree exactly 2 (Lemma T), and neither is an end. So the path uses both of them at `3P/4` and never
  the edge `{P/4, 3P/4}`.
- **Conclusion.** `P/4` is not an end but has only one usable partner. ∎

The range is nonempty iff `5P/4 < 3B - P/2` (`rho < 12/7`) and `3B - 3P/4 < 3P/2` (`rho > 4/3`).

**Proposition W2 (the double choke at `P/4`).** If
`max(5P/4, 3B - 3P/4, 2P - B + P/4, 2B + P/4) <= n < min(7P/4, 3B - P/4)`,
then there is no Hamiltonian path. The range is nonempty iff `4/3 < rho < 8/5`.

*Proof.* As in W1, `P` is a leaf, `N(P/4) = {3P/4, B - P/4}`, and `3P/4` is saturated.
- **`B - P/4` is saturated.** It is saturated by the rigid top vertices `2P - B + P/4 = r_{2P}(B - P/4)` and `2B + P/4 = r_{3B}(B - P/4)`.
- **Conclusion.** `P/4` has a usable partner only if the free end is one of the four saturating vertices. Even then it has exactly one, so it
  must itself be an end. That gives three ends. ∎

**Proposition W3 (chokes at `P/4` and `P/8` with disjoint resolving sets).** If
`max(7P/4, 3B - 3P/8, 2B + P/4, 2P - B + P/4) <= n < min(15P/8, 3B - P/4)`,
then there is no Hamiltonian path. These four lower bounds dominate the other saturators' entry times, and `3B - P/8` exceeds the upper bound.
The range is nonempty iff `4/3 < rho < 3/2`.

*Proof.*
- **The choke at `P/4`.** `P/4` now has the rigid top partner `7P/4`, and its other two partners are saturated as in W2. So the free end lies
  in `R(P/4) = {P/4, 5P/4, 3B - 3P/4, 2P - B + P/4, 2B + P/4}`.
- **The choke at `P/8`.** Its partners are `3P/8, 7P/8, B/3 - P/8, B - P/8`, because `15P/8 > n` and `3B - P/8 > n`. The first, second and
  fourth are saturated by the rigid top pairs `(13P/8, 3B - 3P/8)`, `(9P/8, 3B - 7P/8)` and `(2P - B + P/8, 2B + P/8)`. All six are top
  vertices of degree 2 for `4/3 < rho < 3/2`. So the free end lies in
  `R(P/8) = {P/8, 13P/8, 3B - 3P/8, 9P/8, 3B - 7P/8, 2P - B + P/8, 2B + P/8}`.
- **Conclusion.** `R(P/4)` and `R(P/8)` are disjoint. An equality between one of their linear forms and another forces `rho` either outside
  `(4/3, 3/2)` or to one of `24/17, 32/23, 16/11`, and `rho = 2^p/3^(a-1)` is never one of these (runner B.4b2). ∎

All conditions used are linear inequalities in `(n, P, B)`, and every neighbourhood condition is monotone in `n`. The runner (B.4b) checks,
for every `a <= 200` with exact integers:
- each range is nonempty exactly under the stated `rho` criterion;
- the hypotheses hold at both ends of the range, hence on all of it, and fail just outside;
- for `a <= 10`, the hypotheses hold at every `n`;
- for `a <= 8`, the ranges agree with the solver.

**Instances.**

| `a` | `rho` | W1 | W2 | W3 |
|---|---|---|---|---|
| 5 | 1.580 | `[160, 178]` | `[207, 210]` | - |
| 7 | 1.405 | `[1419, 1535]` | `[1714, 1791]` | `[1803, 1919]` |
| 10 | 1.665 | `[40960, 42664]` | - | - |
| 12 | 1.480 | `[334833, 393215]` | `[419830, 458751]` | `[458752, 465904]` |

- For `a = 5` and `a = 7` the ranges are runs or sub-runs of section 3.4 with the same certificates: `162-178`, `207-210` (`W_5`) and
  `1419-1535`, `1714-1791`, `1803-1919` (`W_7`).
- `a = 10` and `a = 12` are new, and PROVED. The solver confirms the `a = 10` endpoints by propagation (B.4d). There the core is
  `{P/4, 3P/4, 3B - 3P/4, 5P/4} = {8192, 24576, 34473, 40960}`.
- The levels are `a` with `{(a-1) log_2 3}` in a fixed interval. By Weyl equidistribution (classical, CITED) the three families occur at
  levels of density `log_2(9/7) = 0.363`, `log_2(6/5) = 0.263` and `log_2(9/8) = 0.170`. The observed counts for `a <= 200` are `71, 52, 33`
  of 198 levels.
- **This is the "fractal recursion across levels" made exact for the first defects.** Whether a level has the failing run is a condition on
  `rho_a`, that is on the Beatty/Sturmian data of `log_2 3` at level `a`. Its endpoints, divided by `3^(a-1)`, are affine in `rho_a`: for
  example W1 is `[max(5rho/4, 3 - 3rho/4), min(3rho/2, 3 - rho/2))`.

### 3.6 Conjecture, and what stays open for `C_n`

**Conjecture B.**
- **(B1) Propagation completeness.** For every admissible `n`, `C_n` has no Hamiltonian path iff unit propagation refutes every candidate
  free end of Lemma K. FINITE-EXACT for all admissible `n <= 2200`, and consistent with all of `W_8`, which has no failure at all.
- **(B2) Conflict mechanism.** Every failing run inside a window is a union of lifetimes of chokes centred at defects, each with a
  saturation tree of bounded depth. Consequently every inner endpoint is a `{2,3}`-unit sum `n* = T - w(s)`, with `T in {2P, 3^a}` a top
  target, `s` a power of 2 or 3, and `w` a product of at most two reflections.
  - FINITE-EXACT for all 17 inner endpoints with `n > 150` up to `W_8` (runner B.6).
  - The significant part is the short forms: 6 pure connections `n* = T - s` where about 1.2 are expected, and 14 with `|w| <= 1` where
    about 6.4 are expected.
  - W1-W3 prove the depth-one part for every level.
- **(B3) `rho`-scaling.** The Hamiltonian set of `W_a`, divided by `3^(a-1)`, depends only on `rho_a = P/3^(a-1)`. It is constant between
  critical values of `rho` at each depth, with endpoints affine in `rho_a`.
  - PROVED for the conflicts of W1-W3; the full statement is OPEN.
  - This is the precise form of "a finite rule on the Beatty word of `log_2 3`". The Beatty letter at level `a` only records whether
    `rho_a > 3/2`. `W_6` and `W_7` share the letter (two powers of two) but have 2 and 7 Hamiltonian runs, so the letter alone does not decide
    the pattern, while `rho_a` does in every proved case.
- **(B4) Right end.** For `a >= 4`, the last Hamiltonian run of `W_a` ends at the right end of `W_a`. FINITE-EXACT for `a = 4..8`; `a = 3`
  (`n = 27`) is the one exception observed.

**OPEN.**
- A proof of (B1).
- General families for the deeper conflicts (`W_5`: 195-206; `W_7`: 1617-1619, 1664-1702, 2115-2119).
- Any positive statement for all `a`, i.e. Hamiltonicity between the conflicts, which is only computed.
- The continuous limit of the scaled pattern as a function of `rho`.
- Whether the right end of `W_a`, `C_{3^a - 1}` (two-power windows) or `C_{3^a}` (one-power windows), is always Hamiltonian. It is for
  `a = 4..8` and fails at `a = 3` (`C_27`); `C_{3^5}` has a unique, fully forced path.

## 4. Task C: square sums at 18-24 through the same lens

The first window `{15, 16, 17}` is the single tight three-square union `(9, 16, 25)` (Corollary A3''); its chains are rotation orbits
(`x -> x + 9 mod 16`, runner D.1). What follows, typed:

| `n` | status | mechanism (reflection language) | type |
|---|---|---|---|
| 18 | LOCAL | three leaves `16, 17, 18` (THM-4505, Theorem 1) | local (degree) |
| 19 | NONE | ends forced `{16, 18}`; **choke at 3**: `N(3) = {1, 6, 13}`, and `1` is saturated by the rigid `8` (`N = {1,17}`) and `15` (`N = {1,10}`), `6` by the rigid `10` (`N = {6,15}`) and `19` (`N = {6,17}`); so `3` would have to be an end | first-order forcing |
| 20 | NONE | one leaf `18`; **over-saturation at 5** (rigid neighbours `4, 11, 20`) puts the free end in `{4, 11, 20}`, the **choke at 3** puts it in `{3, 8, 10, 15, 19}`: disjoint | first-order forcing |
| 21, 22 | NONE | **lock at the defect 2** (`2 = 4/2`, the fixed point of `r_4`): the leaf `18` (`N = {7}`) and the rigid `9` (`N = {7,16}`) saturate `7 = r_9(2)`, so `2` (`N = {7, 14}`) keeps one partner. Base core `{2, 7, 9, 18}`, so the free end is `2, 7` or `9`, and each is refuted by propagation (cores of 3 to 8 vertices) | second-order forcing |
| 23 | PATH | the lock opens: `23 = r_25(2)` gives the defect `2` a third partner | - |
| 24 | NONE | the base core has 18 vertices (18 itself is not in it); all 18 candidate ends are refuted by propagation, with cores of 3 to 19 vertices | global forcing |
| `>= 25` | PATH | the target `49` enters (a new matching `M_49`); Gerbicz's theorem (A090461, CITED) | - |

All rows are FINITE-EXACT (runner C.1-C.4, own solver; the status list for `n <= 40` reproduces A090461). The certificates for 19 and 20 are
checkable by hand from the neighbour lists above. Answer to "local vs global": 18 is local, 19-22 are small forcing conflicts around the
low vertices `2` (a defect: the half-target `4/2`), `3` and `5`, the same choke mechanism as in `C_n` (section 3). 24 is the only genuinely
global failure, and no failure needs search.

## 5. Task D: the continuous bridge

**5.1 The piecewise isometry (definition).** On `J_n = [1/2, n + 1/2]` let `sigma_t(x) = t - x`, defined on `J_n cap (t - J_n)`. It maps the unit
cell around an integer `x` onto the cell around `t - x`, so `G_S(n)` is the graph of the partial involutions `sigma_t` on cells. A Hamiltonian
path chooses at each cell two of its available reflections so that the chosen partial involutions have a single orbit. By (R2) every return
map is a composite of an even number of reflections, i.e. a piecewise translation: an interval exchange of a subinterval. The orchestrator's
picture is therefore correct as a definition. It becomes a theorem where the return maps can be computed:

**5.2 Three targets: rotations (PROVED, section 2.2).** For a tight triple `a < b < c` (`c - a >= n`), the return map to the `b`-matching is
the rotation by `c - b` on `Z/(c - a)`, cut open at the two ends of the path. Its continuous analogue is the rotation by
`alpha = (c - b)/(c - a)` of `R/Z`. Scaling the tightness conditions (i), (ii) by `n` gives `c - a = n + O(1)`, `b = n + O(1)`: the tight triples
form a one-parameter family, parametrised by `alpha in (0, 1)`. So:
- continuous prediction: a single dense orbit iff `alpha` is irrational (a full-measure, dense set of parameters);
- discrete truth (Theorem A3): a single orbit, i.e. a Hamiltonian path, iff `gcd(c - b, c - a) = 1` (or 2 with `b` odd), the coprime-rational
  shadow of minimality.

The square-sum window `{15, 16, 17}` is `alpha = 9/16` three times (at `n = 17` with the cut at the midpoint `x_0 = 8`); THM-4505's zigzag
family is `alpha = (2k+1)/4k`.

**5.3 Two top targets and the first renormalization step (Lemma T PROVED; the Rauzy reading is ANALOGY).** Lemma T makes the top layer of `C_n` a forced two-target zigzag
`M_{2P} u M_{3Q}`, a translation by `delta = |3Q - 2P|`. Contracting each forced top chain to one super-edge between its two exit points below
`M` leaves a problem on `[1, M]` in which the super-edges pair `y` with `y +- k delta` (`k` = number of top visits). This is one step of an
induction of Rauzy type. The next pair of large targets acts in the same way one level down, which is where the "fractal recursion across
levels" lives.

**5.4 What the continuous picture predicts, and the test.** For interval exchanges the orbit structure can change, as parameters move, only when
the orbit of one singular point runs into another singular point (a connection). Here the singular points are the half-targets `t/2` (fixed
points; for `C_n` the powers of 2) and the domain ends (the targets themselves, and the moving ends `t - n`). The discrete counterpart is that
the propagation outcome can change between `n - 1` and `n` only through the two new edges of the new top vertex `n`, namely
`{n, 2P - n}` and `{n, 3^a - n}` (Lemma T; this part is trivial). Conjecturally (B2), it changes only when `y = T - n` sits in the saturation
tree of a singular vertex.
- Test on `C_n`: each of the 17 inner endpoints in `W_5`-`W_7` (`W_8` has none) is `T - w(s)` with `s` a power of 2 or 3 and `|w| <= 2`
  (runner B.6).
- Six of them are pure connections `n = T - s`, `|w| = 0`: `179, 224, 665, 1536, 1792, 1920`. So are the window starts
  `473 = 3^6 - 2^8` and `1319 = 2^11 - 3^6` of the Window Theorem.
- The per-window base rates of `|w| = 0` among all `n` are 18%, 9%, 3.7% and 1.5% in `W_5` to `W_8` (B.6b). Over the 17 endpoints that
  predicts about 1.2 pure connections; 6 are observed.
- For `|w| <= 1` the prediction is about 6.4 and 14 are observed.
- `|w| <= 2` is too common to carry weight (base rate 76% in `W_7`).

Typing:
- The IET dictionary is ANALOGY. The Keane-type statement "no connection => minimal" is an UNCITED-RECOLLECTION and is used in no proof.
- The rotation of 5.2 is PROVED. The connection pattern is FINITE-EXACT / EMPIRICAL (Conjecture B2, section 3.6).

**5.5 Collatz typing.**
- The translation lengths `delta = |2^q - 3^a|` of the top zigzags (`5, 17, 47, 13, 217, 295, 139, 1909, 1631, ...`) are the Collatz cycle
  denominators of THM-4484 as numbers.
- The sub-window endpoints are `{2,3}`-unit sums of at most four terms, for example `1419 = 3^7 - 2^10 + 2^8` and `576 = 2^10 - 2^9 + 2^6`.
- What is proved is internal to the sum graph: `delta` is the translation `r_{3^a} o r_{2P}` of the forced zigzag.
- There is no map from Hamiltonian paths of `C_n` to Collatz orbits or cycles: ANALOGY.

## 6. Methods, independence, reproduction

**Second code path.** `procgen_sumgraph_20260926_solver.py` was written from scratch for this lane and does not call THM-4505's solvers
(`HCForcing`, `ham.c`, CP-SAT). The two share only the definition of the graph. Its method:
1. Local test: an isolated vertex, or three vertices of degree `<= 1`.
2. Ends: two leaves are both ends. With one leaf, the free end ranges over the base core of Lemma K; if the base run finds no contradiction,
   it ranges over all vertices.
3. Unit propagation (R1, R2, R4 with LOW/OVER/CYC/SHORT) for each candidate end pair.
4. Exhaustive DPLL on survivors: force or delete one edge at a most constrained vertex, with a connectivity cut every 16th node. The
   schedule is:
   - three deterministic branching rules with small budgets, in adaptive order;
   - eight randomized restarts with fixed seeds 0-7;
   - one exhaustive run.

A NONE is final only when every candidate is refuted by propagation or by an exhaustive search. In this range no search was ever needed for
a NONE (Theorem B(c)). Every PATH is re-verified edge by edge.

**Cross-checks.**
- The lane's development scan (`scratch/procgen_sumgraph/scan_2200.jsonl`) and the runner both reproduce THM-4505's list to 2200 exactly,
  including the counts `945/758/497`.
- The smallgraph library's node-limited forcing solver (imported unmodified, development only) was run on `W_8` for `4374 <= n <= 4946`. It
  found paths at 477 values, and my solver agrees at every one of them. It left 96 values UNKNOWN at 30000 nodes; my solver decides all of
  them (PATH).

**Runner.** `python3 -u 04-computation/experiments/procgen_sumgraph_20260926_run.py --full > 05-knowledge/results/procgen_sumgraph_20260926.out`
(from the worktree root).
- It re-checks every finite claim of this note as `[OK]` lines; any failure raises.
- Without `--full`, `W_8` is sampled (every 16th `n`).
- It is single-process pure Python and writes only to stdout.
- The committed `--full` run passed 47 checks in 1390 s wall (1329 s of it on `W_8`), with peak RSS 67 MB.
- The section-time and `RUN` lines vary between runs. The verdicts are deterministic; the paths found are not printed.
- sha256 (raw bytes):

```text
ce6479c12c104703eca0b150716089d4ef0073bfd90be0142878ca9ec4fbd041  procgen_sumgraph_20260926_solver.py
85ead1882e362da316e01ebffdf22f7cd03a9b668fa4d1dc130c96d987c76532  procgen_sumgraph_20260926_theory.py
db5cd230ef9af5df202fac07573f2ef987d16234f4c63e5b0bf612461d05f8f2  procgen_sumgraph_20260926_run.py
4e1465d2beaf6bc974b451dc79732030fff6e3be632d8cff72ffcce0049b624b  procgen_sumgraph_20260926.out
```

**Files.**
- `04-computation/experiments/procgen_sumgraph_20260926_solver.py`: solver, Lemma K cores.
- `04-computation/experiments/procgen_sumgraph_20260926_theory.py`: matchings, Theorems A2/A3 checks, Lemma T, Choke Lemma predictor,
  boundary words.
- `04-computation/experiments/procgen_sumgraph_20260926_run.py`: the runner.
- `05-knowledge/results/procgen_sumgraph_20260926.out`: the runner output.
- Scratch (not for commit): `scratch/procgen_sumgraph/`, holding the development scans, traces and drafts.

## 7. Failures and corrections (own)

- **The deterministic DPLL heuristics.**
  - The first one aborted at `n = 1418` (a Hamiltonian `n`) after 200000 nodes.
  - Later all three aborted at `n = 4976` and `n = 5037`. At `5037` the free end must lie in the core of the choke at the defect
    `512 = P/8`, and the path eventually found ends at `512` itself.
  - Fixed by a portfolio with randomized restarts: a random edge choice finds the `5037` path in about 500 nodes.
  - A development check with an own CP-SAT model (2 workers, 500 MB cap, 45 s) found a path at `5037` with the same ends `(4096, 512)`.
  - No verdict ever depended on these heuristics, since every PATH is verified and every NONE came from propagation.
- **The first over-saturation rule.** The first-order predictor wrongly put the saturated vertex itself into its resolving set. If that
  vertex is the free end, its need drops and the saturation gets worse. Corrected before any claim was drawn from it.
- **"Every conflict is centred at a target".** REFUTED by the data: the chokes at `82, 83, ...` in `W_6` sit at ordinary vertices. The
  claim is weakened to "defects are where conflicts live longest" (EMPIRICAL).
- **A pure-reflection word test.** Asking only for short words from powers of 2 or 3 is too weak: most `n` have one with `|w| <= 2`. The
  significant part is the short forms, `|w| <= 1` (see 5.4).
- **Resource slip.** A sample run of the runner was started while the `W_8` scan was running: two CPU-bound processes for about two minutes.
  It was killed at once. Every other heavy run of this lane was alone. Peak RSS stayed below 100 MB for every run except the single
  CP-SAT development check (469 MB, under the 700 MB cap).
