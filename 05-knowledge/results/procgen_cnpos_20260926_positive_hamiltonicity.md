# Positive Hamiltonicity for the Collatz alphabet C_n: the top zigzag returns as a reflection modulo its translation |3^a - 2^k|, an exact reduction at every n of every window, constructions verified at the clean n of every level a <= 15, and why infinitely many levels need the residual problem at every scale

**Status block.**

- **PROVED.**
  - *Theorem F (first return; the construction lemma).* Let `T1 < T2 < 2 T1` be a power of 2 and a power of 3 (in either order),
    `m = T2 - T1`, and `n = T1 - 1 - j` with `j >= 0` and `2n >= T2 - 1`. The reflection matchings `alpha = M_T1` and
    `beta = M_T2` on `[1, n]` return to the residual `R = [1, T2 - 1 - n]` as the reflection `phi(x) = T1 - x (mod m)` of the
    window `W = [j + 1, j + m]`. For every set `mu` of pairs inside `R` (one pair at each vertex of `W`, two at each vertex of
    `F = [1, j]`, one fewer at one end `e`), `alpha u beta u mu` is a Hamiltonian path of `[1, n]` iff `Phi u mu` is a
    Hamiltonian path of `R`, where `Phi` is the set of `phi`-pairs. So every solution of the *residual problem*
    `RP_j(m, T1 mod m)` by C-edges is a Hamiltonian path of `C_n`.
  - *Corollary F1 (every window).* For every `a >= 2` and every `n` in `W_a` with `n != 3^a`, the two top targets of `C_n`
    satisfy the hypotheses. The modulus `m` is the translation `delta = |3^a - 2^k|` of THM-4510's Lemma T.
  - *Proposition F2 (exactness).* In the lower part of every two-power window (`T1 = 2P_n`) the reduction is exact: `C_n`
    has a Hamiltonian path with free end in `R` iff `RP_j` has a solution by C-edges. When `T1 = 3^a` (one-power windows,
    upper parts of two-power windows), the alpha-edges inside `I = [3^a - P_n, P_n]` are not forced. The exact form
    (*Theorem F'*, contraction of the forced chains) lets each interior vertex of `I` choose between `alpha` and `M_P`, and
    Theorem F is its all-alpha choice.
  - *Lemma Z (zone step).* The same first-return argument for an arbitrary forced involution and an arbitrary zone matching.
  - *Proposition U (unit rotations).* Modulo `m` we have `3^a = 2^q`, so the residual targets `2^(q-1), 2^(q-2), 3^(a-1),
    3^(a-2)` compose with `phi` to rotations by `2^(q-1), 3 * 2^(q-2), 2 * 3^(a-1), 8 * 3^(a-2)`. These are units mod `m`, so a
    first zone step with one of them never closes a cycle. The factors `1, 3, 2, 8` are units because of the Gersonides
    relations `2 - 1 = 4 - 3 = 3 - 2 = 9 - 8 = 1`.
  - *Proposition R (rotations).* The residual problem at `n = T1 - 1` can be solved by a single rotation (Theorem A3's shape)
    only if the two linear pieces of one circular reflection of `Z/m` are C-edges: two targets `c < c + m`, or `m + 1` or
    `m + 2` a target. FINITE-EXACT for `a <= 300`: this happens only at `a = 2, 3, 5`, the Pillai coincidences
    `4 - 3 = 9 - 8`, `9 - 4 = 32 - 27`, `16 - 3 = 256 - 243`. `C_7`, `C_26`, `C_242` have Hamiltonian paths that are these
    rotations.
  - *Proposition S (2-exchange rigidity).* Exchanging two edges `{x,y}, {z,w}` of a Hamiltonian path (or of `mu`) for two
    edges on the same four vertices needs `tau1 + tau2 = tau3 + tau4` among targets (PROVED). With targets `< 2^200` the only
    relations are `3 + 9 = 4 + 8`, `3 + 32 = 8 + 27`, `3 + 256 = 16 + 243`, `4 + 32 = 9 + 27` (FINITE-EXACT; completeness
    for all exponents is not claimed and not used). Each uses target 3 or 4, that is, the edge
    `{1, 2}` or `{1, 3}`.
- **FINITE-EXACT** (own code; every path re-verified edge by edge from the definition of `C_n`).
  - *Clean n.* `C_(T1 - 1)` is Hamiltonian at every level `2 <= a <= 15`, where `T1 = min(3^a, 2^q)` and
    `2^q / 3^a in (2/3, 4/3)`: `n = 7, 26, 63, 242, 511, 2047, 6560, 16383, 59048, 131071, 524287, 1594322, 4194303,
    14348906`.
  - *Right ends (Conjecture B4).* The construction also gives `C_(3^a - 1)` at the two-power levels
    `a = 2, 6, 7, 9, 11, 12, 14` and `C_(3^a)` at the one-power levels `a = 8, 10, 13`. At `a = 15` the heuristic did not
    settle `C_(3^15)`. New: `n = 19682, 59049, 177146, 531440, 1594323, 4782968` (`a = 9..14`). With THM-4510
    (`a <= 8`), the right end of `W_a` is Hamiltonian at every level `4 <= a <= 14`.
  - *Whole windows `W_3..W_8`.* Every one of the 2937 Hamiltonian `n` of THM-4510's list (`n != 3^a`, all of `W_8`) gets
    a verified path:
    - 2394 by the all-alpha construction (Theorem F);
    - 540 by the exact form (Theorem F');
    - 3 (`n = 5116, 5117, 5248`, at and just above the choke range of the defect 512) only by the direct search. There the
      exact-form search ran out of budget; the paths end at 512 and 2560, which Theorem F' covers.
  - *`W_9`.* Verified at 42 values (`16343..16383` and `19682`).
- **EMPIRICAL.** The zone-step renormalisation solves the residual problem at every level tried. It applies Lemma Z from the
  top, pairs the power-of-2 fixed point of an even target through its power-of-3 partner, and finishes with an exact search
  below 2500 (or 8000) vertices. Each level takes seconds, with peak RSS at most 316 MB.
- **CONDITIONAL.** If the residual problem `RP(m_a, K_a)` at `n = T1 - 1` has a solution by C-edges for every `a >= 2`
  (Conjecture RP, OPEN), then `C_(T1 - 1)` is Hamiltonian at every level (Theorem F). At two-power levels Conjecture RP is
  equivalent to the existence of a Hamiltonian path of `C_(2^q - 1)` with free end in `[1, m_a]`.
- **ANALOGY.** The modulus `m = |3^a - 2^k|` is the Collatz cycle denominator of THM-4484. The residual problem lives on
  `Z/m`, but there is no map from Hamiltonian paths to Collatz orbits.
- **OPEN.**
  - Conjecture RP, and with it a positive theorem for infinitely many levels. The question is reduced here to one residual
    instance per level and is not proved. Section 5 lists what blocks the direct routes.
  - A proof that the zone-step renormalisation always succeeds.
  - Whole windows beyond `W_8`.

Collatz itself is not addressed. Session `collatz-procgen-20260922`, lane `cnpos` (2026-09-26).
- Scripts: `04-computation/experiments/procgen_cnpos_20260926_{core,run}.py`.
- Output: [procgen_cnpos_20260926.out](procgen_cnpos_20260926.out).

## 0. The four tasks in one paragraph each

**A (structure of the actual paths).**
- The solver's paths are search output with no visible pattern. At `n = 6560`, for example, the interval `I = [2465, 4095]`
  splits into 157 runs of `M_3B`-edges and 156 runs of `M_P`-edges.
- So the structure was read off the forced part instead. Every Hamiltonian path contains the top zigzag, and the top zigzag,
  followed out of the residual, is a reflection modulo its translation (Theorem F).
- At the small levels the whole path is one rotation. The paths of `C_26` and `C_242` found by the solver use only `27, 32`
  and `243, 256` at the top and `4, 9` and `3, 16` in the residual: the rotations `x -> x + 3` of `Z/5` and `x -> x + 6` of
  `Z/13`. This works because `9 - 4 = 32 - 27` and `16 - 3 = 256 - 243` (Proposition R).

**B (construction lemma and theorem).**
- Theorem F is the construction lemma: any solution of the residual problem on `[1, T2 - 1 - n]` lifts through the forced top
  zigzag. It applies at every `n` of every window, the residual modulus being the zigzag translation `|3^a - 2^k|`.
- At `n = T1 - 1` the residual is smallest: it is `RP(m, K)`, a circular reflection of `Z/m` plus a C-matching. There the
  construction succeeds at every level `a <= 15`.
- A theorem for infinitely many levels would follow from Conjecture RP (CONDITIONAL). We could not prove it. The residual
  problem has the same two-choice zigzag at its own top (Lemma T for `C_m`), so a proof must control the renormalisation at
  every scale.

**C (verification).**
- Every constructed path is checked edge by edge.
- The construction is compared with the Hamiltonian list of THM-4510 on `W_3..W_8` (all of `W_8`).
  - Every Hamiltonian `n` gets a verified path: 2394 from the all-alpha construction, 540 from the exact form, and 3 from the
    direct search, where the exact form's search ran out of budget.
  - The all-alpha construction finds exactly the Hamiltonian `n` in `W_5`, `W_6` and in every two-power lower part.
- `W_9` is sampled at 42 values.

**D (what fails).**
- The all-alpha construction loses paths where `T1 = 3^a`, because there the edges `x <-> P - x` of the interval `I` are
  optional and some Hamiltonian paths need them (`C_80`, `C_243`).
- Rotation constructions exist only at `a = 2, 3, 5` (Proposition R).
- Two-edge exchanges are impossible away from vertex 1 (Proposition S), so local cycle-merging proofs are unavailable.
- `RP(m, K)` is not solvable for every `(m, K)`: `RP(11, 1)`, `RP(59, 22)` and `RP(383, 1)` have no solution. A proof must
  therefore use the arithmetic of `(m_a, K_a)`.
- The Collatz reading of `m` is ANALOGY.

## 1. Setting

`C_n` has vertices `1..n` and an edge `{x, y}` (`x != y`) iff `x + y` is a *target*, a power of 2 or of 3 that is at least 3.
Throughout:
- `M_t` is the reflection matching `{x, t - x}` on `[1, n]` (THM-4510, Lemma R);
- `P_n`, `Q_n` are the largest powers of 2 and 3 not above `n`;
- the top vertices are `(max(P_n, Q_n), n]`.

By Lemma T (THM-4510) a top vertex has exactly the two neighbours `2P_n - x` and `3Q_n - x`. `P_n` is a leaf (THM-4505, C1).
Write `T1 < T2` for the two top targets `{2P_n, 3Q_n}`. With `B = 3^(a-1)` and `P` the power of 2 in `(B, 2B)`, the windows of
THM-4505 (C2) split into three regimes:

| part of `W_a` | `T1` | `T2` | `m = T2 - T1` | reduction |
|---|---|---|---|---|
| one-power window `[2B, 3B - 1]` (`P > 3B/2`) | `3B = 3^a` | `2P` | `2P - 3^a` | Theorem F (all-alpha); exact form F' |
| two-power, lower part `[L_a, 2P - 1]` (`P < 3B/2`) | `2P` | `3^a` | `3^a - 2P` | Theorem F, exact |
| two-power, upper part `[2P, 3B - 1]` | `3^a` | `4P` | `4P - 3^a` | Theorem F (all-alpha); exact form F' |

In every row `m` is the translation `delta = |3^a - 2^k|` of the forced top zigzag (THM-4510, table in section 3.1). The
residual `R = [1, T2 - 1 - n]` has `m + j` vertices, where `j = T1 - 1 - n` is the distance to the clean value `n = T1 - 1`.

## 2. Theorem F: the top zigzag returns as a reflection modulo its translation

**Theorem F (PROVED).** Let `T1 < T2 < 2 T1` be targets, one a power of 2 and the other a power of 3, and put `m = T2 - T1`
(odd). Let `j >= 0` and `n = T1 - 1 - j` with `2n >= T2 - 1`, and put `r = T2 - 1 - n = m + j`. Define
- `alpha = M_T1` and `beta = M_T2` on `[1, n]`, the reflection matchings of the two targets;
- `F = [1, j]`, `W = [j + 1, r]`, `R = F u W`, and the zone `Z = [r + 1, n]`;
- `phi : W -> W`, `phi(x) = T1 - x (mod m)`, an involution with the single fixed point `u0` (`2 u0 = T1 mod m`);
- `f` = the half of the even target (a power of 2);
- `Phi = {{x, phi(x)} : x in W, x != u0}`.

Let `mu` be a set of pairs inside `R` such that every vertex of `W` lies in exactly one pair and every vertex of `F` in exactly
two, except one vertex `e in R` that lies in one fewer. Then:
- (a) `G = alpha u beta u mu` has maximum degree 2, `n - 1` edges, and its vertices of degree 1 are `f` and `e`;
- (b) `G` is a Hamiltonian path of `[1, n]` iff `Phi u mu` is a Hamiltonian path of `R` from `u0` to `e`;
- (c) if every pair of `mu` sums to a target, `G` is a Hamiltonian path of `C_n`.

*Proof.* (a) The domain of `M_t` on `[1, n]` is `[t - n, n]` minus `t/2`, so:
- `alpha` covers `[j + 1, n]` except `T1/2`, and `beta` covers `Z = [r + 1, n]` except `T2/2`;
- exactly one of `T1`, `T2` is even, and its half is `f`;
- every vertex of `Z` other than `f` has one `alpha`-edge and one `beta`-edge;
- a vertex of `W` has only its `alpha`-edge (the pair `T1 - x` lies in `[n + 1 - m, n]`);
- a vertex of `F` has neither.

Adding `mu` gives degree 2 everywhere except at `f` and `e`. (If `f` lies in `W`, then `f = u0`, since `T1 - f = f`.) The degree
sum is `2n - 2`. Hence `G` is a Hamiltonian path iff it is connected, iff it has no cycle.

(b) *The first return.* Follow `G` from `b in W` along its `alpha`-edge. By Lemma R2 the walk is
`b, T1 - b, b + m, T1 - b - m, b + 2m, ...`:
- the `beta`-steps land at `b + km`, which lie in `Z` (as `beta`-images of zone vertices);
- the walk stops at the first `k` with `T1 - b - km <= r`.

That vertex lies in `W`. For `k = 0` it is at least `n + 1 - m >= j + 1`. For `k >= 1` it exceeds `r - m = j`, because
`T1 - b - (k - 1)m > r`. It is congruent to `T1 - b`, so it is `phi(b)`, unless the walk meets `f` first.

The walk meets `f` only if `b = u0`. Indeed `f` is an `alpha`-departure `b + km` (if `f = T1/2`) or a `beta`-departure
`T1 - b - km` (if `f = T2/2`). In both cases `b = f (mod m)`, and `2f = T1 (mod m)` gives `f = u0 (mod m)`. Conversely the walk
from `u0` reaches `f` before leaving `Z`. If `f = T1/2`, the `beta`-departures before it are at least `T1 - f + m = f + m > r`.
If `f = T2/2`, `f` is one of the decreasing terms `T1 - u0 - km`, all of which are `> r` up to `f`, since
`u0 <= (T1 - m)/2 = T1 - f`.

The walks cover `Z`:
- every zone vertex has degree 2 in `alpha u beta` (1 for `f`);
- `alpha u beta` has no cycle, because `beta o alpha` is the translation by `m != 0` (Theorem A2(a) of THM-4510);
- so every zone vertex lies on a maximal path whose ends have `alpha u beta`-degree at most 1, i.e. lie in `W` or equal `f`.

Contracting each walk to one edge turns `G` into `Phi u mu` plus the edge `u0 - f`. So `G` is one path iff `Phi u mu` is one
path from `u0` to `e`.

(c) `alpha`, `beta` and `mu` consist of C-edges. ∎

**Corollary F1 (every window, PROVED).** For every `a >= 2` and every `n in W_a`, `n != 3^a`, the top targets
`T1 < T2` of `C_n` satisfy the hypotheses of Theorem F with `j = T1 - 1 - n`. Indeed:
- `T1 > n`;
- `T2 <= 2n` (with equality only at `n = P_n`, where `M_(2n)` is empty);
- `T2 / T1 < 2`, being below `4/3`, `3/2` and `2` in the three rows of the table of section 1.

So `C_n` is Hamiltonian as soon as its residual problem `RP_j` on `[1, T2 - 1 - n]` has a solution by C-edges.

**Proposition F2 (exactness, PROVED).** If `T1 = 2P_n` (the lower part of a two-power window), `C_n` has a Hamiltonian path
with free end in `R` iff `RP_j` has a solution by C-edges.

*Proof.* In a window `P_n > Q_n`, so the top vertices are `(P_n, n]`.
- Every `alpha`-pair `{x, 2P_n - x}` and every `beta`-pair `{x, 3Q_n - x}` has an endpoint above `P_n`. Top vertices have
  degree 2 (Lemma T) and are not ends, since `f = P_n` is a leaf and `e` lies in `R` below `P_n`. So a Hamiltonian path `H`
  contains `alpha u beta`.
- Then every zone vertex is saturated, so the remaining edges of `H` join vertices of `R`, with the degrees of Theorem F. By
  Theorem F(b) they solve `RP_j`. ∎

When `T1 = 3^a`, the pairs `{x, 3^a - x}` with both ends in `I = [3^a - P_n, P_n]` have no top endpoint and are optional.

**Theorem F' (the exact form, PROVED).** Take the *forced* edges: both edges of every top vertex, and the edge of the leaf
`P_n`. Contract every forced chain to one virtual edge between its ends, and call the reduced graph `H_n`. Then `C_n` has a
Hamiltonian path whose free end is not saturated by forced edges iff `H_n` has a Hamiltonian path using every virtual edge,
starting at the far end of the chain of `P_n`.

In the one-power case with `n >= 3P - 3^a`:
- the forced chains run from each `b in W` up by `m` to the unique `c(b) = b (mod m)` in `[3^a - P, P - 1]`;
- each interior vertex `x` of `I` keeps exactly two optional edges, `alpha` (to `3^a - x`) and `M_P` (to `P - x in [1, m - 1]`).

Choosing `alpha` everywhere gives back Theorem F.

**Remark (rho-scaling of the residual data, PROVED).** At `n = T1 - 1` the residual is `RP(m, K)` with `K = T1 mod m`. So
`K/m = frac(1/(rho' - 1))` with `rho' = T2/T1`, and the targets near `m`, divided by `m`, are functions of `rho'` alone (for
example `2^(q-2)/m = 1/(4(rho' - 1))` when `T1 = 2^q`). Here `rho'` is determined by the fractional part of `a log_2 3`. The
continuum data of the residual problem thus follow THM-4510's rho-scaling (Conjecture B3). Solvability, however, also depends
on the integers (section 5.4).

## 3. The residual problem, the zone step, and the renormalisation

**The residual problem `RP_j(m, K)`.** Take `[1, m + j]` with the free part `F = [1, j]`, the window `W = [j + 1, j + m]`, and
the reflection `phi(x) = K - x (mod m)` on `W`. The task is to find C-edges inside `[1, m + j]` giving one edge to each vertex
of `W` and two to each vertex of `F` (one fewer at the end `e`), such that `Phi u mu` is a Hamiltonian path from `u0`. For
`j = 0` it is a circular problem on `Z/m`: the reflection `phi` together with a matching of `C_m`.

**Lemma Z (zone step, PROVED).** Let `psi` be an involution of a finite set `V` with one fixed point `s`, and let `nu` be a
fixed-point-free involution of a subset `Z`. The maximal `psi u nu`-walks through `Z` define an involution `psi'` of
`V' = V \ Z`, with one fixed point `s'`: the vertex whose walk reaches `s`, or `s` itself. Then:
- if some vertex of `Z` lies on no such walk, `psi u nu u mu'` is never a Hamiltonian path;
- otherwise `psi u nu u mu'` is a Hamiltonian path of `V` iff `psi' u mu'` is one of `V'`.

*Proof.* Theorem F(b), with `phi, alpha, beta` replaced by `psi, nu`. ∎ The runner (A.3) confirms Lemma Z on 300 random
instances, enumerating every completing matching.

**Proposition U (unit rotations, PROVED).** Let `n = T1 - 1` and `{T1, T2} = {3^a, 2^q}`, so that `3^a = 2^q (mod m)`. For
`tau in {2^(q-1), 2^(q-2), 3^(a-1), 3^(a-2)}`, the composite `phi o M_tau` is the rotation by
`T1 - tau = 2^(q-1)`, `3 * 2^(q-2)`, `2 * 3^(a-1)`, `8 * 3^(a-2)` respectively. These are units mod `m`, since `m` is prime
to 6. So the first zone step on `RP(m, K)` itself, with such a `tau` and its full zone, never closes a cycle: the walks are the
orbits of one transitive rotation of `Z/m`, and each of them meets the residual.

The four factors `1, 3, 2, 8` are `{2,3}`-units because `2 - 1 = 1`, `4 - 3 = 1`, `3 - 2 = 1` and `9 - 8 = 1`, the solutions
of `3^x - 2^y = +-1` (Gersonides; elementary, section 6). For `j >= 3` neither `2^j - 1` nor `3^j - 1` is a `{2,3}`-unit, so
among the targets `2^(q-j)` and `3^(a-j)` only these four give a rotation that is a unit modulo every such `m`. The targets
lie in `(m, 2m)`, i.e. are available at the top of `C_m`, for:
- `rho' = T2/T1` in `(18/17, 9/8) u (9/8, 3/2)` when `T1 = 2^q`;
- `rho'` in `(19/18, 10/9) u (8/7, 4/3)` when `T1 = 3^a`.

**The renormalisation (EMPIRICAL).** `core.gp_solve` repeats zone steps from the top:
- take a target `t` of the largest vertex and its full zone;
- if `t` is even, pair its fixed point `t/2` through a power-of-3 partner, or make it the end;
- if the step leaves a vertex without partners, re-pair it along a short chain;
- stop with an exact search once at most 2500 vertices remain (8000 if that fails).

It is a heuristic; every result is verified. In the runs it is fast (section 4).

**Proposition R (rotations, PROVED + FINITE-EXACT).** Suppose `RP(m, K)` has a solution that is one circular reflection,
`mu(x) = c - x (mod m)`. Then `Phi u mu` is a single orbit of the rotation `x -> x + K - c` (Theorem A3 of THM-4510).

On `[1, m]` such a `mu` consists of the linear pieces `x + y = c` (nonempty iff `c >= 3`) and `x + y = c + m` (nonempty iff
`c <= m - 1`). So it is made of C-edges only if
- `c` and `c + m` are both targets, or
- `c in {1, 2}` and `m + c` is a target.

For `m = m_a` and `a <= 300` this happens only at `a = 2` (`m = 1`), `a = 3` (`9 - 4 = 5`) and `a = 5` (`16 - 3 = 13`),
runner C.1. These are the Pillai coincidences `9 - 4 = 32 - 27` and `16 - 3 = 256 - 243`. The Hamiltonian paths of `C_26`
and `C_242` found by the solver are exactly these rotations: `x -> x + 3` on `Z/5` and `x -> x + 6` on `Z/13`.

## 4. Results (FINITE-EXACT; every path verified edge by edge)

### 4.1 The clean n of every level (runner B.1)

`n = T1 - 1`, `T1 = min(3^a, 2^q)`, residual `RP(m, T1 mod m)` on `[1, m]`:

| `a` | window | `n` | `T1` | `T2` | `m` | `rho' = T2/T1` | unit targets in `(m, 2m)` (Prop. U) | residual solved by | time |
|---|---|---|---|---|---|---|---|---|---|
| 2 | two-power | 7 | 8 | 9 | 1 | 1.1250 | none (`m = 1`) | trivial (the Gersonides chain) | 0.0s |
| 3 | one-power | 26 | 27 | 32 | 5 | 1.1852 | 8, 9 | zone steps + exact search | 0.0s |
| 4 | two-power | 63 | 64 | 81 | 17 | 1.2656 | 32, 27 | zone steps + exact search | 0.0s |
| 5 | one-power | 242 | 243 | 256 | 13 | 1.0535 | none | zone steps + exact search | 0.0s |
| 6 | two-power | 511 | 512 | 729 | 217 | 1.4238 | 256, 243 | zone steps + exact search | 0.0s |
| 7 | two-power | 2047 | 2048 | 2187 | 139 | 1.0679 | 243 | zone steps + exact search | 0.0s |
| 8 | one-power | 6560 | 6561 | 8192 | 1631 | 1.2486 | 2048, 2187 | zone steps + exact search | 0.2s |
| 9 | two-power | 16383 | 16384 | 19683 | 3299 | 1.2014 | 4096, 6561 | zone steps + exact search | 0.1s |
| 10 | one-power | 59048 | 59049 | 65536 | 6487 | 1.1099 | 6561 | zone steps + exact search | 0.0s |
| 11 | two-power | 131071 | 131072 | 177147 | 46075 | 1.3515 | 65536, 59049 | zone steps + exact search | 0.3s |
| 12 | two-power | 524287 | 524288 | 531441 | 7153 | 1.0136 | none | zone steps + exact search below 8000 | 4.0s |
| 13 | one-power | 1594322 | 1594323 | 2097152 | 502829 | 1.3154 | 524288, 531441 | zone steps + exact search | 1.5s |
| 14 | two-power | 4194303 | 4194304 | 4782969 | 588665 | 1.1403 | 1048576 | zone steps + exact search | 4.7s |
| 15 | one-power | 14348906 | 14348907 | 16777216 | 2428309 | 1.1692 | 4194304, 4782969 | zone steps + exact search below 8000 | 26s |

Times are for the residual problem; lifting and verifying the path of `C_n` takes longer at the top levels. At two-power
levels (Proposition F2) the residual solutions are exactly the Hamiltonian paths of `C_n` with free end in `[1, m]`.
`a = 16` (`m = 9492289`) is left out, since its arrays would exceed the 500 MB cap.

A unit target of Proposition U sits at the top of `C_m` at every level of the table except `a = 2` (`m = 1`) and
`a = 5, 12`.
- The last two are the levels with `rho'` below the thresholds `19/18` and `18/17` of Proposition U, i.e. with `2^q / 3^a`
  closest to 1.
- They are the continued-fraction convergents `8/5` and `19/12` of `log_2 3` (`2^8 = 3^5 + 13`, `2^19 = 3^12 - 7153`).
- At `a = 12` the heuristic needed the larger exact-search base.

### 4.2 The right ends of the windows (Conjecture B4)

| level | right end | construction | status |
|---|---|---|---|
| two-power `a = 2, 6, 7` | `8, 728, 2186` | `T1 = 3^a`, `T2 = 2^(q+1)`, residual solved | Hamiltonian (agrees with THM-4510) |
| two-power `a = 4` | `80` | all-alpha residual has no solution (exact) | Hamiltonian by THM-4510 (a construction loss, section 5.1) |
| two-power `a = 9, 11, 12, 14` | `19682, 177146, 531440, 4782968` | residual solved | **Hamiltonian (new)** |
| one-power `a = 3` | `27` | no solution with the end at `m` | `C_27` is not Hamiltonian (THM-4505) |
| one-power `a = 5` | `243` | no all-alpha solution with the end at `m` | Hamiltonian by THM-4510 (unique path, uses `M_P`) |
| one-power `a = 8` | `6561` | residual solved with the end at `m`, then the edge `{m, 3^a}` appended | Hamiltonian (agrees with THM-4510) |
| one-power `a = 10, 13` | `59049, 1594323` | the same | **Hamiltonian (new)** |
| one-power `a = 15` | `14348907` | heuristic did not finish (residual too large for exact search) | undecided |

So the right end of `W_a` is Hamiltonian at every level `4 <= a <= 14` (`a <= 8` by THM-4510, `9 <= a <= 14` here).
Conjecture B4 is thereby confirmed through `a = 14`.

### 4.3 Whole windows `W_3..W_8` (runner B.2, all of `W_8`)

| window | regime | `n` checked (`n != 3^a`) | Hamiltonian (THM-4510) | all-alpha construction (Theorem F) | exact form (Theorem F') on the rest | direct search only |
|---|---|---|---|---|---|---|
| `W_3 = [18, 27]` | one-power | 9 | 9 | 6 | 3 (`20, 21, 22`) | 0 |
| `W_4 = [49, 80]` | two-power | 32 | 30 | 18 | 12 (upper part `68..80`, not `78`) | 0 |
| `W_5 = [162, 243]` | one-power | 81 | 35 | 35 | 0 | 0 |
| `W_6 = [473, 728]` | two-power | 256 | 167 | 167 | 0 | 0 |
| `W_7 = [1319, 2186]` | two-power | 868 | 509 | 441 | 68 (upper part, from 2055) | 0 |
| `W_8 = [4374, 6561]` | one-power | 2187 | 2187 | 1727 | 457 | 3 (`5116, 5117, 5248`) |
| total | | 3433 | 2937 | 2394 | 540 | 3 |

- **No false positives.** Every constructed path lies at an `n` of THM-4510's Hamiltonian list, as it must, being verified.
- **Exact form.** Wherever the all-alpha construction fails, the exact form (Theorem F': forced chains contracted, exact
  search on `H_n`) finds a verified path, except at `n = 5116, 5117, 5248`.
  - There the search on `H_n` exhausted its budget twice (2e5 and 2e6 nodes).
  - The direct search of THM-4510 (Lemma K cores, both ends prescribed) finds paths ending at 512 and at the rigid vertex
    2560. This matches THM-4510 section 3.4 for the choke range `5025..5214` of the defect `512 = P/8`; 5248 lies just above
    that range.
  - Neither end is saturated by forced edges (neither is a top vertex or a forced neighbour of one), so Theorem F' covers
    these paths.
  - These are search failures, not structural ones.
  - So on `W_3..W_8` the reduction reproduces the whole Hamiltonian list positively, with three values supplied by the
    direct search. The exact search engine is shared with THM-4510; the reduction and the verification are new.
- **All-alpha construction.**
  - It finds exactly the Hamiltonian `n` of `W_5`, `W_6` and of every two-power lower part (Proposition F2).
  - Its losses lie in the regimes with `T1 = 3^a`: `20..22` (`W_3`), `68..80` without `78` (the upper part of `W_4`), 68 values
    in the upper part of `W_7` (from 2055 on), and 460 values of `W_8`.

### 4.4 `W_9` (runner B.3)

- The construction is verified at `n = 16343..16383`: the exact regime, `T1 = 2^14`, `m = 3299`, `j = 0..40`.
- It is also verified at the right end `n = 19682` (`T1 = 3^9`, `T2 = 2^15`, `m = 13085`).
- The whole window was not attempted. At its left end the residual has 8191 vertices, 4892 of them free, and the exact
  search did not finish in 144 s (development run).

## 5. What fails, and why a theorem for infinitely many levels is out of reach here

**5.1 The all-alpha choice (PROVED mechanism + FINITE-EXACT losses).**
- Where `T1 = 3^a`, the reduction of Theorem F fixes every `alpha`-pair inside `I = [3^a - P, P]`. These pairs are optional,
  and the vertices of `I` could instead take `x <-> P - x` (Theorem F').
- Some `n` have Hamiltonian paths only of the second kind. At `n = 80` and `n = 243` the all-alpha residual has no solution,
  by exact search. The unique path of `C_243` uses the `M_P`-edges `10 - 118` and `125 - 3` (runner B.1d).
- In the exact regime (`T1 = 2^k`), all alpha-edges are forced and nothing is lost.

**5.2 Rotation constructions reach only `a = 2, 3, 5` (Proposition R).**
- The cleanest positive statement would be a rotation, as in THM-4510's Theorem A3: one circular reflection in the residual,
  so that the whole path is one orbit of a rotation of `Z/m`.
- This needs the two linear pieces of that reflection to have target sums, i.e. a second representation of `m_a` as a
  difference of targets (or `m_a + 1`, `m_a + 2` a target). For `a <= 300` that happens only at `a = 2, 3, 5` (runner C.1).
  (Finiteness for all `a` would be a Pillai-type S-unit statement. It is not claimed and is used nowhere.)
- So `C_7` (the Gersonides chain), `C_26` and `C_242` are rotation orbits, and no other level `a <= 300` has one at
  `n = T1 - 1`.
- Richer constructions with a fixed number of pieces are not excluded by this. What they meet is 5.3 to 5.5.

**5.3 No two-edge exchanges (Proposition S).**
- A proof that merges the cycles of `Phi u mu` into one path by local switches needs exchanges of two `mu`-edges. These
  require `tau1 + tau2 = tau3 + tau4` among targets.
- The only such relations below `2^200` are `3 + 9 = 4 + 8`, `3 + 32 = 8 + 27`, `3 + 256 = 16 + 243` and
  `4 + 32 = 9 + 27` (runner C.2), and they involve the edges `{1, 2}`, `{1, 3}` only.
- So Hamiltonian paths of `C_n` are rigid under 2-exchanges away from vertex 1. Only rotations at the free end and exchanges
  of three or more edges remain. This rules out merging cycles by two-edge switches, the usual first step of absorption
  arguments.

**5.4 The residual problem is not solvable for every `(m, K)` (FINITE-EXACT examples + EMPIRICAL scan).**
- `RP(11, 1)`, `RP(11, 5)`, `RP(59, 22)` and `RP(383, 1)` have no solution. Propagation or exhaustive search refutes them
  (runner C.3).
- A development scan (scratch, not in the runner) of all odd `m <= 401` and all `K` prime to `m` found, for most odd `m`
  between 65 and 401, between one and nine values of `K` without solution. Many of them have `K` or `K + m` a target, where
  `phi` shares edges with `C_m`.
- The level instances avoid them at every `a <= 15`. For example `m_6 = 217` has unsolvable `K` but not `K_6 = 78`.
- So a proof must use the arithmetic of `(m_a, K_a)`, not only the size or the reflection shape.

**5.5 The renormalisation has to go down through every scale.**
- The residual problem inherits the whole difficulty. The vertices of `[1, m]` above `max(P_m, Q_m)` have exactly two
  C-options (Lemma T for `C_m`), so the residual has its own two-choice zigzag, one scale down.
- The zone steps peel these layers one at a time. Each step's shape is set by new arithmetic: which power of 2 or 3 is the
  top target of the current residual, and whether `phi` shares its edges. At the `a = 12` right end, for example,
  `M_(3^12)` coincides with `phi` on most of `[1, m]`.
- Proposition U controls only the first step.
- The heuristic needed the larger exact-search base (8000 vertices) at `a = 12` and `a = 15`. It could not settle the
  one-power right end at `a = 15`, where the end is forced at `m`.
- Away from `n = T1 - 1` the residual carries `j` free vertices and is as hard as `C_n` itself. The left end of `W_9` is an
  example (section 4.4).

**5.6 Where this leaves the positive question.**
- Theorem F reduces "is `C_n` Hamiltonian" to a residual problem on `Z/m` plus a free part, at every `n` of every window, and
  exactly in the two-power lower parts.
- A positive theorem for infinitely many levels would follow from Conjecture RP (one residual instance per level).
- By 5.2 to 5.5 it does not come from a single rotation or from two-edge repairs, and the residual reproduces the difficulty
  one scale down. It would need a renormalisation argument that controls the residual at every scale. That is OPEN.
- The evidence for Conjecture RP: every level `a <= 15` (section 4.1), plus the right ends through `a = 14`.

## 6. Collatz typing

- **(ANALOGY)** The modulus `m = |3^a - 2^k|` of the residual problem is the translation of the top zigzag (THM-4510) and,
  as a number, the Collatz cycle denominator (THM-4484). The residual problem is a circular reflection on `Z/m`. No map
  sends Hamiltonian paths of `C_n` to Collatz orbits or cycles.
- **(PROVED)** The unit rotations of Proposition U come from `2 - 1 = 1`, `4 - 1 = 3`, `3 - 1 = 2`, `9 - 1 = 8`.
  - For `tau = 2^(q-j)` the rotation `T1 - tau = 2^(q-j) (2^j - 1)` is a unit modulo every `m` prime to 6 exactly when
    `2^j - 1` is a power of 3. That forces `j <= 2`: for `j >= 3`, `3^i = -1 (mod 8)` is impossible.
  - Likewise `3^j - 1` is a power of 2 only for `j <= 2`. Odd `j` gives `2 (mod 4)`, and even `j` gives two powers of 2
    differing by 2.
  - These are Gersonides' solutions of `3^x - 2^y = +-1`, the pairs `(2,3), (3,4), (8,9)` of THM-4484.
  - This is internal to the sum graph.
- **(ANALOGY)** The rotations of `Z/m` by these units are the discrete analogue of THM-4510's continuous bridge. Nothing in
  this note maps them to the Collatz map.

## 7. Methods, independence, reproduction

**Code.**
- `04-computation/experiments/procgen_cnpos_20260926_core.py`:
  - the reduction data (`top_targets`, `reduced_data`, `phi_of`);
  - the lift (`lift`) and the verifier (`verify_path`, which uses only the definition of `C_n`);
  - the exact residual search (`solve_residual_dpll`);
  - the exact form (`forced_structure`, `solve_exact_reduced`);
  - the zone-step renormalisation (`rp_gp`, `zone_step`, `gp_solve`, `check_gp`).
- The exact searches call the sumgraph lane's propagation/DPLL engine (`procgen_sumgraph_20260926_solver.py`), imported
  read-only. Its answers are never trusted: every path is lifted by this lane's code and verified by `verify_path`.
- `04-computation/experiments/procgen_cnpos_20260926_run.py`: the runner, sections A (structure), B (constructions),
  C (obstructions).

**Checks of the proved statements.**
- A.1: the parameters of Theorem F at all 3437 `n` of `W_2..W_8`.
- A.2: the first-return map, by explicitly walking `alpha u beta` from every `b in W`, at 1419 values of `n`.
- A.3: Lemma Z on 300 random instances, by enumerating every completing matching.
- C.1: Proposition R's necessary condition for `a <= 300`.
- C.2: Proposition S's relations for targets below `2^200`.

**Runner.**
- Command: `nice python3 -u 04-computation/experiments/procgen_cnpos_20260926_run.py --full > 05-knowledge/results/procgen_cnpos_20260926.out`
  (from the worktree root).
- It is single-process and writes to stdout only.
- The committed `--full` run passed 13 checks in 4784 s wall, most of it on the whole of `W_8`, with peak RSS 316 MB (at
  `a = 15`).
- The section-time and `RUN` lines vary between runs, and so can the set of values that need the direct-search fallback in
  B.2. The exact-search engine orders its branching rules adaptively across calls (`_ORDER` in the sumgraph solver), so
  whether a search finishes within a fixed budget depends on the call history: a development run needed the fallback at two
  values, the committed run at three. The verified paths and the NONE verdicts do not depend on it.
- Without `--full`, `W_8` is sampled (every 8th `n`).

**sha256 (raw bytes).**

```text
837c00a5ddf65217902a720fd9ad218cf8ba226507a4527068921b3390913968  procgen_cnpos_20260926_core.py
1b2ba537ea8ed4bb443b14fbed6b19296ad9bed3ee189dd8a47cd51df8205f79  procgen_cnpos_20260926_run.py
99bdb281ee273547530d6e6e2234318312547117408aec9ba263b0c7e93b901e  procgen_cnpos_20260926.out
```

## 8. Failures and corrections (own)

- **The first renormalisation was too greedy.**
  - It left the fixed point `t/2` of an even target isolated, and failed from `a = 12` on.
  - Fixed by pairing `t/2` through its power-of-3 partner and adding a short repair chain.
  - The isolation test then missed high leftovers from earlier steps, which made `a = 15` fail. Fixed by testing all
    vertices above the new residual interval.
- **The plain exact search on the residual problem does not scale.** It had not finished at `a = 11` after 12 minutes (the run
  was stopped), and it aborted at the left end of `W_9` after 144 s. Hence the zone steps.
- **The all-alpha construction was at first read as exact everywhere.** The window comparison showed the losses (`C_80`,
  `C_243`, parts of `W_3`, `W_4`, `W_7`, `W_8`). Theorem F' and the exactness statement (Proposition F2) are the corrected
  form.
- **Runner bugs caught before the final run.**
  - A.2 first walked every walk twice and reported a false revisit.
  - My hand list of two-edge relations missed `3 + 9 = 4 + 8` and `4 + 32 = 9 + 27`.
  - C.1 first ignored circular reflections with a single linear piece (the full reflection `M_(m+1)`).
  - The first `--full` run failed its B.2 check. The exact form's search ran out of budget at two values of `W_8`, although
    both are Hamiltonian. The runner now retries with a larger budget and then falls back to the direct search, reporting
    each such value; the rerun passed.
- **Resource slips.**
  - A dict-based development version of the renormalisation reached 1.2 GB peak RSS for about 26 s on the `a = 15` instance,
    above the 500 MB cap. All later code uses arrays: the final runner peaks at 316 MB.
  - Twice a test of a few seconds ran next to a background run: two CPU-bound processes, briefly.
