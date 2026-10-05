# The two-sheet receipt calculus (2026-10-05): reversed `3x-1` paths are the multipliers, a `3x-1`-rooted multiplier costs exactly one fusion relation, the trunk move `m_j = (2^j+1)/3` pays every `-1 mod 2^k` cell, a 16-bit plain-first prefix code with `3x-1`-rooted multipliers covers every odd residue class but `-1 mod 2^16` (the Applegate--Lagarias multipliers `5, 7, 13, 23` are `3x-1`-cyclic and not needed), universal two-sheet receipts cost `0.75` fusion relations per integer below `2^16`, and the first trunk multiplier is the sheet commutator `U_+(3x) = U_-(U_+(x))` whose shadow orbit merges with the orbit of `x` through the quarter-child relation with probability exactly `1/3` (`1/2` in every deep cell, decided at step `k-3`)

**Session:** opus, `pascal-lyapunov-20261005` (third task; worktree `codex/session-pascal-lyapunov-20261005`), 2026-10-05.

**2026-10-05 targeted scope/timing correction:** the sibling automaton is a
sufficient common-future detector, not a necessary one. At x=31 its deep-cell
test fails, yet the actual orbit of31 reaches its shadow35 and then53.
In a successful deep cell the quarter-child relation is reached at
simultaneous time k-3 and the two states become equal one step later.
The probability1/2 uses conditional Haar measure; under the computable
atomic source measure of the [prefix-mass follow-up](collatz_effective_prefix_mass_20261005.md),
section7, it is11/12 for even k and1/12 for odd k. The local algebra survives;
an unguarded claim that structural non-detection prevents every repair does
not. Other statements and finite censuses have not all been independently audited.
Owner's directive: "build the two-sheet receipt calculus and test the trunk multiplier move" (the direction recorded at
the end of the [four-questions note](collatz_receipts_four_questions_20261005.md)).
**Inherits (cited):** the [universal weak receipts note](collatz_universal_weak_receipts_20261005.md) (receipts,
boundary, fusion relation `R`, rule G1); the four-questions note (Theorems 1, 6; the semigroup proof dissected);
Applegate--Lagarias, [The 3x+1 semigroup](https://arxiv.org/abs/math/0411140), Lemmas 2.1-2.3; the
[dependency kernel](collatz_recursive_dependency_kernel_20261004.md) (the quarter-child relation `F_v(n) = 4h + 1`,
`U(F_v(n)) = U(h)`); [THM-4517](../../01-canon/theorems/THM-4517-no-root-uniform-collatz-basin-density.md) (the
three `3x-1` basins); [THM-4512](../../01-canon/theorems/THM-4512-coefficient-descent-classes-one-member.md)
(exact and coarse cylinders).

**Status: PROVED (Propositions A, B, Theorems C, D, E below; elementary); FINITE-EXACT (one script, section 8:
`1,156` move checks, codes to `16` bits, `32,767` universal receipts, `250,000` sheet pairs with `16,660` structural
merges confirmed against the orbits, `10,000` deep-cell members); author-audited only, audit OWED; Collatz OPEN.
The calculus changes the FORM of the universal obligations (one checkable fusion relation per multiplier event),
not their nature; the one dynamical content found is the sheet commutator and its sibling automaton.**

> **Targeted correction, 2026-10-05:** the sibling automaton detects one
> sufficient common-future repair. Its BREAK outcome does not exclude other
> actual-edge repairs or later orbit intersections. The former section5
> elimination iff is refuted by x=7 below; TheoremE's merge predicate is
> scoped to the maintained-relation automaton. The exact counterexample is
> replayed in [the selector audit](vitali_selector_certificate_20261005.md).
> This is not an independent audit of the other probabilities or censuses.

---

## 0. What is new, in one screen

1. **The calculus (section 1).** Edges of both sheets `U_+(u) = oddpart(3u+1)`, `U_-(u) = oddpart(3u-1)`, values
   `u/U_s(u)`; `+` edges with nonnegative stock, `-` edges only reversed. A reversed `-` path from `m` to `1` has value
   `1/m`: **the multipliers of the calculus are exactly the `3x-1`-rooted integers**, and for such `m` the packet
   `(+ path from m x) - P_-(m)` has value `x/x'` and transport defect exactly `R(m, x) = [m x] + [1] - [m] - [x]`
   (Proposition A). The trunk move is `m = m_j = (2^j+1)/3`, a single `-` edge, with `U_+(m_j x) = oddpart(x +
   (x+1)/2^j)` on `x = -1 mod 2^j` (`756/756`; the landing point is itself odd when `j < v_2(x+1)`, `720/720`).
2. **Which multipliers exist.** `3x-1`-rooted odd `m <= 100`: `3, 11, 15, 29, 39, 43, 53, 57, 59, 65, 69, 71, 77, 79,
   85, 87, 95, 97` (the others fall into the cycles `{5,7}` and `{17,...,91}`). Of the paper's wild multipliers,
   `11, 29, 43` are rooted (`11 = m_5`, `43 = m_7`) and `5, 7, 13, 23, 25, 35` are cyclic: the two-sheet calculus
   cannot realize them, and does not need them.
3. **Prefix codes (section 3).** Plain Collatz descent (THM-4512's exact and coarse cylinders) pays `93.55%` of the
   odd residues modulo `2^16`; with one `3x-1`-rooted multiplier on the remaining classes every class but `-1 mod
   2^16` is paid (`6.45%` by multipliers: `3` on `5.35%`, `11` on `0.78%`, `15` on `0.24%`, six others below
   `0.05%`); the trunk numbers `{3, 11, 43, 171, 683}` alone leave four deep cells, the paper's set leaves eight,
   all of the form `-1 mod 2^k`, `k >= 12`. At twelve bits the code pays `88.96%` plainly and `10.99%` by rooted
   multipliers against the paper's `79.69% / 20.26%`, with the same single uncovered class `-1 mod 4096`.
4. **Universal two-sheet receipts (section 4).** Built stage by stage from the code (deep cells by the trunk move),
   every odd `x <= 65536` gets a receipt whose defect is a sum of fusion relations `R(m_i, x_i)`, `0.747` of them
   per integer on average (`46.9%` need none, `36.2%` one, maximum `6`); the greedy code in the style of the
   paper's tables costs `2.68`.
5. **The trunk move as dynamics (sections 5-6).** For `m_j` with `j >= 5` the orbit of `m_j x` never meets the orbit
   of `x` before `x` descends (`0` of `500` in every cell tried, `k = 8..20`). For `m_3 = 3` it does, and the reason
   is an identity: **for `v_2(3x+1) = 1`, `U_+(3x) = U_-(U_+(x))`** (Theorem C; `250,000/250,000`, and never for
   `x = 1 mod 4`). The move's orbit is the `3x-1` shadow of `x`'s own orbit; the two orbits then run in a
   relation `B = 2^c S +- 1` whose exact automaton (Theorem D) merges them, through the quarter-child relation
   `B = 4S + 1`, with probability `1/3` for a sheet pair and `1/2` from the state `(+1, 1)`: measured `0.3333` on
   the `250,000` odd `x = 3 mod 4` below `10^6` with every entry class at its dyadic weight, and exactly `1/2` in
   every deep cell `x = t 2^k - 1`, where the merge is decided at step `k - 3` by one bit of `t` (Theorem E). The
   merge point is the first point after `x`'s climb, so no descent is gained; the fusion relation `R(3, x)` is
   eliminable for those `x` only because their own orbit already carries the shadow's descent.
6. **The mixed trunk (section 6).** `x = (4^e - 1)/(2^j + 1)` with `j | e` has `m_j x` on the `+` trunk, so the
   trunk move reaches `1` in one step; `e = j` gives `x = 2^j - 1`, the representatives of the deepest cells
   (`7, 31, 127, 511, 2047` with `m = 3, 11, 43, 171, 683`): the paper's `x11` at `31` and `x43` at `127` are these.

---

## 1. The calculus

**Definition.** Labels `[u]`, `u` odd positive. Edges `(u, s)`, `s in {+, -}`, from `[u]` to `[U_s(u)]`, value
`u/U_s(u)`. A two-sheet receipt is a finite integer-valued `c` on edges with `c(u, +) >= 0` and `c(u, -) <= 0`;
`V(c) = prod value^c`, `d c = sum c ([u] - [U_s u])`. For a ratio `x/x'` the transport defect is
`Delta(c) = d c - ([x] - [x'])`. A nonnegative `+`-only receipt with `Delta = 0` is an actual Collatz path (G2).

**Proposition A (one multiplier, one relation).** Let `m` be `3x-1`-rooted with `-` path `P_-(m) = (m -> ... -> 1)`,
and let `Q` be any `+` path from `m x` to `x'`. Then `c = Q - P_-(m)` is a two-sheet receipt with `V(c) = x/x'` and
`Delta(c) = R(m, x)`. *Proof.* `V(Q) = m x/x'`, `V(P_-(m)) = m/1`; `d Q = [m x] - [x']`, `d P_-(m) = [m] - [1]`;
subtract and compare with `[x] - [x']`. □ (`400/400` random checks, `m <= 400`, `x < 10^6`, depths `<= 11`.)

**Proposition B (the trunk move).** For odd `j`, `m_j = (2^j+1)/3` is an odd integer with `3 m_j - 1 = 2^j`, so
`P_-(m_j)` is the single edge `m_j -> 1`; for `x = -1 mod 2^j`, `3 m_j x + 1 = 2^j (x + (x+1)/2^j)`, hence
`U_+(m_j x) = oddpart(x + (x+1)/2^j)`, and the landing point `x + (x+1)/2^j` is odd iff `j < v_2(x+1)`. The move's
receipt `{m_j x -> x'} - {m_j -> 1}` has value `x/x'` and defect `R(m_j, x)`. □ This is the paper's Lemma 2.2 in
the accelerated map; its only obligation is one fusion relation whose multiplier is a trunk number of the other
sheet.

**Example.** `x = 7`, `j = 3`, `m_3 = 3`: `21 -> 1` on `+`, `3 -> 1` on `-`: value `7`, boundary `[21] - [3]`, defect
against `[7] - [1]` equal to `[21] - [3] - [7] + [1] = R(3, 7)`.

## 2. Which multipliers the calculus has

A reversed `-` path from `m` has value `1/m` only if it reaches `1`; the `-` map has the fixed point `1` and the
cycles `{5, 7}` and `{17, 25, 37, 55, 41, 61, 91}` (THM-4517: basin densities about `0.327 / 0.325 / 0.348`), so
roughly a third of the odd integers are `3x-1`-rooted: below `100`, `3, 11, 15, 29, 39, 43, 53, 57, 59, 65, 69,
71, 77, 79, 85, 87, 95, 97`. A `-` path into a cycle gives a rational multiplier `c/m` (`c` on the cycle), which
certifies `c x`, not `x`. The paper's `H = {5, 7, 11, 13, 23, 29, 43}` (and the products `25, 35` in its tables)
splits as rooted `{11, 29, 43}` and cyclic `{5 -> 7, 7 -> 5, 13 -> 19 -> 7, 23 -> 17, 25 -> 37, 35 -> 13 -> ...}`.

## 3. Prefix codes with `3x-1`-rooted multipliers (FINITE-EXACT; Test B)

A class `s mod 2^K` is paid by `(m, pos, steps)` if multiplying by `m` at the odd index `pos` and following the
valuation word determined by the `K` known bits gives, after `steps` odd steps, an affine image `P x + Q` with
`P < 1` and `P s + Q < s` (the smallest member; then every member descends, the paper's worst-case criterion);
when the next valuation is not determined but is at least the remaining precision, the coarse bound
`(3y+1)/2^(K-A)` is used (THM-4512's coarse cylinder; it pays e.g. `21845 mod 2^16`, where every member has
valuation at least `16`). Codes are built by refining classes until they are paid (depth `<= 20` odd steps,
multiplier at index `<= 3`, at most one multiplier), either greedily (a multiplier tried before refining, the
style of the paper's tables) or plain-first (refine to `K_max = 16`, multipliers only on the leaves).

| code (`K_max = 16`) | plain | multiplier | uncovered |
|---|---|---|---|
| plain descent only | `93.55%` | -- | `6.45%` (`2,114` classes) |
| plain-first, `3x-1`-rooted `m <= 100` | `93.55%` | `6.45%` (`3`: `5.35`, `11`: `0.78`, `15`: `0.24`, `29`: `0.04`, `39`: `0.02`, `53, 69, 43`: `< 0.01`) | **`-1 mod 2^16` only** |
| plain-first, trunk numbers `{3, 11, 43, 171, 683}` | `93.55%` | `6.44%` (`3`: `5.35`, `11`: `0.78`, `43`: `0.22`, `171`: `0.06`, `683`: `0.03`) | `4` cells `-1 mod 2^k`, `k = 13..16` |
| plain-first, the paper's `{5, 7, 11, 13, 23, 29, 43, 25, 35}` | `93.55%` | `6.43%` (`5`: `5.55`, `7`: `0.55`, `11`: `0.22`, `29, 13, 23, 35, 43`) | `8` cells `-1 mod 2^k`, `k = 12..16` |
| greedy, `3x-1`-rooted (`17` classes) | `68.75%` | `31.25%` (`3`: `25.0`, `11`: `4.7`, `15`: `1.2`, ...) | `-1 mod 2^16` only |
| plain-first, `3x-1`-rooted, `K_max = 12` | `88.96%` | `10.99%` | `-1 mod 4096` only (the paper: `79.69 / 20.26`, same uncovered class) |

**Reading.** The classes that plain descent cannot pay in `16` bits are paid by one `3x-1`-rooted multiplier,
almost always by `3 = m_3`, and the only class left is the deepest `-1` cell, which the trunk family pays by
construction (Proposition B with `j <= k - 5`). The paper's cyclic multipliers `5, 7, 13, 23` are replaceable
everywhere; where the paper multiplies by `11` at `31 mod 64` or by `43` at `127 mod 256` it is already making the
trunk move (section 6). The paper's table also inserts multipliers at classes that a deeper plain code pays, which
is why its multiplier density (`20.26%`) is twice the plain-first value at the same depth (`10.99%`).

## 4. Universal two-sheet receipts and their defect rank (FINITE-EXACT; Test C)

For odd `x`, look up its deepest class in the code, apply the certificate to the integer (multiplying where the
certificate says), reach `x' < x`, repeat; a deep cell `-1 mod 2^k`, `k >= 11`, not in the code takes the trunk move
with `j = 1 mod 6`, `k - 10 <= j <= k - 5` (the paper's choice) and continues from `x + (x+1)/2^j`. The receipt
is the sum of the stage packets; its defect is `sum_i R(m_i, x_i)` over the multiplier events.

| code | mean relations per `x <= 65536` | distribution `0, 1, 2, 3, 4, 5, 6` | maximum |
|---|---|---|---|
| plain-first, `3x-1`-rooted | **`0.747`** | `15383, 11849, 4256, 1033, 214, 20, 12` | `6` (`x = 11983`) |
| plain-first, trunk only | `0.748` | `15383, 11873, 4218, 1039, 214, 26, 14` | `6` |
| plain-first, the paper's set | `0.733` | `15383, 11897, 4415, 992, 78, 2` | `5` |
| greedy, `3x-1`-rooted | `2.680` | `1599, 5690, 8681, 7912, 5212, 2404, 888, ...` | `10` |

Every `x` below `2^16` received a receipt; the trunk fallback fired once (`x = 65535`). The `46.9%` with no
relation are the integers whose stage classes are all plainly paid, so their universal receipt is their own orbit;
the count measures the symbolic blindness of a `16`-bit code, not a property of the integers (all are rooted).

## 5. The trunk move as dynamics: the sheet commutator and the sibling automaton (PROVED + FINITE-EXACT)

**Test D (coarse).** In the cells `x = t 2^k - 1` (`t` odd, `500` values each, `k = 8..20`), the `+` orbit of `m_j x`
meets the `+` orbit of `x` above `1` in about `88%` of cases (the hubs), but BEFORE `x` descends below itself in
`0` cases for every `j >= 5` (`m = 11, 43, 171, 683, 2731`) and in `45-57%` of cases for `j = 3` (`m = 3`).

**Theorem C (sheet commutator).** For odd `x` with `v_2(3x+1) = 1`, `U_+(3x) = U_-(U_+(x))`. *Proof.* `3x + 1 =
2 U_+(x)`, so `3x = 2 U_+(x) - 1` and `3 (3x) + 1 = 6 U_+(x) - 2 = 2 (3 U_+(x) - 1)`. □ For `x = 1 mod 4` the
identity never holds (`0/249,999`). So the first trunk move is `x -> U_+(x) -> U_-( )`: one `+` step followed by
one `-` step of `x` itself, and `y := U_-(U_+(x))` is the `3x-1` shadow of `x`'s second point.

**Theorem D (sibling automaton).** For odd `n` the sheet pair `(U_+(n), U_-(n))` satisfies `B = 2^c S + eps` with
`{B, S}` the two images, `eps = +1` and `c = v_2(3n-1) - 1` when `v_2(3n+1) = 1`, `eps = -1` and `c = v_2(3n+1) - 1`
otherwise. Under one `+` step on both sides: `(+1, 1) -> (+1, b)` with `b = v_2(3S+1)`; `(+1, 2)`: `B = 4S + 1` and
`U_+(B) = U_+(S)` (MERGE, the quarter-child relation); `(+1, 3) -> (-1, b+1)`; `(+1, >= 4)`: the `+-1` form is
lost (BREAK); `(-1, 1)`: `U_+(B) = U_-(S)`, a fresh sheet pair of `S`; `(-1, >= 2)`: BREAK. *Proof.* Direct
expansion of `3B + 1 = 3 2^c S + 3 eps + 1` case by case. □ Starting from `x = 3 mod 4` with the pair of
`n = U_+(x)` (the `x` side is `U_+^2(x)`, the shadow side `y`), the merge probability under uniform valuations is
`1/2` from `(+1, 1)` and `1/3` from a fresh pair. Measured on the `250,000` odd `x = 3 mod 4` below `10^6`: entry
weights `(+1,1): 1/4, (+1,2): 1/8, (+1,3): 1/16, (+1,>=4): 1/16, (-1,1): 1/4, (-1,2): 1/8, ...` exactly, merges
`0.1250 + 0.1250 + 0.0833 = 0.3333`; every automaton merge below `2 10^5` is a true coincidence of the two orbits
(`16,660/16,662`, the two misses at a `400`-step cap), while the orbits meet somewhere above `1` in `93%` of
cases (hub coincidences).

**Theorem E (deep cells, corrected timing and scope).** For
`x = t 2^k - 1`, t odd and k>=4, put `y = 9t 2^(k-3)-1` and
`X = U_+^2(x) = 2y+1`. After k-4 simultaneous valuation-one steps,
the next valuation of the smaller state is
`b = 1 + v_2(3^(k-1)t-1) >= 2`. This automaton reaches the quarter-child
state iff b=2, equivalently `t = 3^k mod4`. That is half the exact cell
under conditional Haar measure. It reaches the quarter-child state at
simultaneous time k-3 and actual equality at k-2. If b=3 it next enters
a negative-sign state with exponent>=2 and stops certifying; if b>=4 it
stops directly. Neither outcome excludes a later or asynchronous common
future. *Proof.* Expand the two affine recurrences and use that an exact
cell `-1 mod2^r` has r-1 initial valuation-one steps, followed by a valuation
at least2. Apply `U_+(4s+1)=U_+(s)` for the final equality. See the independent
exact audit in the prefix-mass follow-up. □ Measured: detector success `0.5000` and the quarter-child detection time exactly
`k - 3` for `k = 8, 12, 16, 20, 24` (`2,000` values each). The merge point is `x`'s first point after its own
`k`-step climb, which may already be below `x` (`x = 255`: climb to `4373 = 4 * 1093 + 1`, then `205`): the move
reveals that half of each deep cell shares its post-climb orbit with its shadow, and supplies no additional descent beyond that already present in the common-future paths.

**Consequence for the receipts, scope-corrected.** An automaton merge supplies
a particular actual common-future repair. BREAK only loses its tracked
relation; it does not characterize all possible actual-edge repairs. For
example, x=7 gives U_+(x)=11 and sheet pair (17,1), so 17=2^4*1+1 immediately
BREAKs. Nevertheless the checked actual ROOT words Q_21=(6), Q_3=(1,4),
Q_7=(1,1,2,3,4) satisfy

    boundary(Q_21-Q_3-Q_7) = [21]+[1]-[3]-[7] = R(3,7).

Replace the minus-sheet path of3 by its checked plus-sheet path Q_3 and
subtract this signed boundary repair; the resulting receipt is the
nonnegative actual path Q_7. This uses those explicit root certificates,
not an unproved automatic repair rule.

The original restricted census remains: among the `21,676` `x3` events of the universal receipts below `2^16` only `794`
(`3.7%`) are at a stage source `x_i = 3 mod 4` whose chain merges, because the code inserts `x3` at classes without
plain descent, where the merge, when it happens, lands at the end of a climb. The sibling repair is real and
rare.

## 6. The mixed trunk

`x = (4^e - 1)/(2^j + 1)` is an integer iff `j | e` (`2^j = -1 mod 2^j + 1`), and then `m_j x = (4^e - 1)/3` is on
the `+` trunk: the trunk move lands at `1` in one step, with defect `R(m_j, x)`. The first members: `j = 3`: `7,
455, 29127`; `j = 5`: `31, 31775, 32537631`; `j = 7`: `127, 2080895`; `j = 9`: `511`; `j = 11`: `2047`. With
`e = j` these are `x = 2^j - 1`, the smallest members of the deepest cells, whose own orbits need `5, 39, ...`
odd steps (`x = 32537631`: `104`). The paper's table entries `31 -> x11 = 341`, `127 -> x43 = 5461`, `63 -> x11 =
693`, `255 -> x43 = 10965` are trunk moves.

## 7. What the calculus changes, typed

* **Form of the obligations.** A universal receipt's defect is `sum R(m_i, x_i)` with every `m_i` carrying a
  checked `-` route to `1`; no wild packet, no induction on wild numbers (the paper's hypothesis (3)) is needed.
  This is a cleaner normal form of universal weak coverage, PROVED by Propositions A-B and the code of section 3
  up to `16` bits (and by the paper's Lemma 2.3 argument for the deep cells).
* **Nature of the obligations.** Unchanged: eliminating `R(m, x)` means an actual route from `x` to a rooted vertex
  (four-questions note, Theorem 2), i.e. the plain descent of the class that needed the multiplier -- the cells
  of the partition-cover obstruction.
* **What the move knows about the dynamics.** Only `m_3 = 3` commutes the sheets (Theorem C); its shadow orbit is
  the `3x-1` image of `x`'s second point and merges with `x`'s orbit through the quarter-child relation with the
  exact probabilities of Theorems D-E. This is the repo's `F_v(n) = 4h + 1` common future written on two sheets;
  it explains the `45-57%` of Test D and the paper's choice of `x3`-free tables (they used `5` where `3` works).
  For `j >= 5` the move is purely multiplicative.
* **Not a Collatz step.** Nothing above descends an integer that its own orbit does not; the calculus reorganizes
  the bookkeeping of universal weak coverage around checkable multipliers and exposes one two-sheet identity.

## 8. Reproduction

```bash
python 04-computation/experiments/collatz_two_sheet_receipts_20261005.py 16 65536     # 16 s
```
Parts A (move and Proposition A), B (codes), C (universal receipts), D (coarse meeting), F (commutator and
automaton, ground truth), E (mixed trunk). Output `.out` beside the script.

## 9. Verdicts

| claim | status |
|---|---|
| a `3x-1`-rooted multiplier costs one fusion relation; the trunk move is a single `-` edge (Propositions A, B) | PROVED |
| `3x-1`-rooted multipliers pay every odd class modulo `2^16` but `-1 mod 2^16`; `5, 7, 13, 23` are cyclic and unnecessary | FINITE-EXACT |
| universal two-sheet receipts: `0.747` relations per integer below `2^16` (plain-first), `2.68` (greedy) | FINITE-EXACT |
| `U_+(3x) = U_-(U_+(x))` for `v_2(3x+1) = 1` (Theorem C); the sibling automaton (Theorem D); deep cells merge with probability `1/2` at step `k-3` (Theorem E) | PROVED + FINITE-EXACT |
| the move for `j >= 5` never meets `x`'s orbit before descent; for `j = 3` the merge is the quarter-child relation | FINITE-EXACT + PROVED |
| the mixed trunk `(4^e-1)/(2^j+1)`, `j \| e`, reaches `1` in one move; `2^j - 1` are its first members | PROVED |
| Collatz | OPEN |
