# The signed descent against real tournaments: what adding vertices can reverse, the ladder at the floor, the 3-adic basin function, and the dominance-pruned slow-factor tree

Session `mac-mini-2026-11-15-descent`, 2026-11-15. **The Chan–Pak conjecture is OPEN; nothing here proves it.**
Not independently audited (owner directive of 2026-11-13: ideas, not audits). Every statement carries its status.

`H(T)` is the number of Hamiltonian paths of a tournament `T`; a *count* is an odd number that is `H` of some
tournament; `T'` *contains* `T` means as an induced subtournament. The **half-count** of `T` is
`n(T) = (H(T) - 1)/2`.

## 0. Results at a glance

| # | statement | status |
|---|---|---|
| 1 | In the half-count chart the increment floor reads `n(T') = n(T)` or `n(T') >= (3 n(T) + 1)/2` for every `T'` containing `T`: the least nontrivial growth is the shortcut Collatz odd step. | KNOWN theorem (Session214 forced shuffle, HYP-9469 (1)); the reading is new |
| 2 | No realiser of a rising target `(4d-3)/3` or `(4d+1)/3` of the signed descent contains a realiser of `d`. Adding vertices can only ever reverse the halving branch. | PROVED (two-line corollary of 1) |
| 3 | The halving branch is not always reversible by adding vertices either: for the 54 counts `h = 1 (mod 4)`, `h <= 223`, the value `(h+1)/2` is the count of a subtournament of some realiser of `h` for 38; not for 13 and 41 (value 7, 21) nor for 57, 69, 77, 85, 93, 105, 109, 125, 133, 141, 153, 169, 173, 201. | FINITE-EXACT, complete |
| 4 | Every count `5 <= h <= 223` has a realiser `T` and a vertex `x` with `|H(T) - 2 H(T - x)| <= 5`: defect 1 for 79 counts, 3 for 27, 5 for 41 and 125. For `225 <= h <= 3357`, realisers of order at most 10 give defect at most 7. | FINITE-EXACT (complete to 223) / NUMERICAL beyond; HYP-9812 (a) |
| 5 | Middle vertex on a bridge (section 3): `H(T) = 2ab + Δ - E S`, `H(T + x) = 3ab + 2Δ`, so `2H(T+x) - 3H(T) = Δ + 3 E S`. | PROVED |
| 6 | Ladder at the floor: for `m = φ_a φ_b` (`φ_0 = 1`, `φ_a = 2^a + 1`) there are tournaments `T` inside `T'` with counts `2m - 1` and `3m`; in half-counts `2k -> 3k + 1`, one step above the floor. For `a = 0` or `b = 0` the chain starts at a realiser of `m`: `k -> 2k -> 3k + 1`, the Collatz odd step as two one-vertex extensions. | PROVED (all `a, b`); the three-term chain also holds for `min(a,b) = 1`, `max <= 3` (FINITE-EXACT) and fails for (2,2), (2,3), (3,3) |
| 7 | Over all tournaments of order at most 9 and all their one-vertex extensions, `2H(T+x) - 3H(T) <= 3` (nontrivial) occurs exactly for `H(T) = 3` (value 1) and `H(T) = 2m - 1`, `m` in {1, 3, 5, 9, 15, 17, 25, 27, 33, 45, 51, 65} = all products `φ_a φ_b <= 65`. | FINITE-EXACT; HYP-9812 (b) |
| 8 | The backward tree of the signed descent: dead rays are the half-counts divisible by 3; roots carrying `j` consecutive rising steps are exactly `2^j` classes mod `3^j` (density `(2/3)^j`); the Haar average of level sizes is exactly `(5/3)^k`; the identity `1/2 + (2/3)(3/4) = 1` makes the tree critical (dimension 1). **21 is the record root** (9 rising steps: 21, 27, 35, 47, 63, 83, 111, 147, 195, 259; next record 489), and 7 is a dead ray. | PROVED (classes, average, identity); FINITE-EXACT (records to `n < 2·10^6`) |
| 9 | Basin densities: 0.481 of the odd numbers below `2·10^7` descend through 21 (0.368 through 27, 0.114 through 41). `F(t) = t ρ(t)` has mean 6.129 (renewal value 6.117) and its means over `n mod 9` are 3.94, 6.6, 7.87 against the exact values `(9/14, 15/14, 9/7)·6.117 = 3.93, 6.55, 7.87` forced by the functional equation. An orbit meets a dyadic window 1, 2 or 3 times (frequencies 0.648, 0.293, 0.059). | PROVED (bound 3, class-mean algebra, additivity) / NUMERICAL (densities) |
| 10 | The dominance-pruned slow-factor tree was built in three versions. With inclusion dominance of column spaces, propagated, the undominated nodes at rate 1.72 number 4, 20, 24, 84, 344, 837, 2616, 4971, 11647 at depths 1–9: growth about 2.3 per level. **It does not close.** | NUMERICAL (exploratory), negative |
| 11 | Fence and square-sum readings (section 7). | ANALOGY |

What is **not** claimed: any new count, any closing of the rank certificate, any proof that the near-doubling
parents persist, the existence of basin densities as limits.

## 1. Setting, sources and credit

**The signed descent** is from `six_tutte_collatz_20261116/COUNTS.md` (another session): on odd `h`,

    D(h) = (h+1)/2   (h = 1 mod 4),      (3h+3)/4   (h = 3 mod 8),      (3h-1)/4   (h = 7 mod 8).

It is strictly decreasing with fixed points 1 and 3. In the half-count chart `n = (h-1)/2` it is `n -> n/2`
(`n` even) and `n -> (3n + ε)/4` with `ε = ±1`, `ε = n (mod 4)` (`n` odd): the 3n±1 map with the sign that allows
two halvings. Its inverse branches are `d -> 2d - 1` (all `d`), `d -> (4d-3)/3` (`d = 3 mod 6`) and
`d -> (4d+1)/3` (`d = 5 mod 6`); only 13 lands on 7 and only 27, 41 land on 21. That packet also records the
half-count chart, the ray `6·2^i` above 3 and the fact that the ordinal sum with a directed triangle is
`n -> 3n + 1`.

**Concurrent work on the same directive**, `inverse_fences_landing_20261116` (another session), realises each of the
three inverse branches on explicit infinite one-apex and two-apex families by fixed-size surgery, proves that an
opaque one-port template cannot be a universal gadget, and shows that no realiser of 105 contains a realiser of 53.
This note is complementary: it determines what *adding vertices* can do at all (sections 2–4), describes the
arithmetic tree itself (section 5) and reports the tree asked for in the directive (section 6).

**The increment floor** is HYP-9469 (1), (6), proved in all orders by Session214's forced-shuffle theorem
(`session214_chan_pak_combined_proof_20261107/forced-shuffle/PROOF.md`): a one-vertex extension either keeps `H` or
raises it to at least `(3H + 1)/2`, and the same holds for every induced extension of any size. The sole-exit case
and its equality cases are THM-4876 (7). (I re-derived the floor in this session before finding it; no claim.)

**The bridge of two one-apex tournaments** is `chanpak_target_ribbon_audit_20261116/TWO_APEX.md`
(`H = hj + PQ - 1`); section 3 uses the same count for arbitrary blocks.

## 2. The floor is the Collatz shortcut step, and what it forbids

**Reading (1).** `H(T') >= (3H(T) + 1)/2` is `n(T') >= (3n(T) + 1)/2`. So in half-counts every induced extension
either stays put or moves at least as far as the shortcut step `n -> (3n+1)/2`. For even `n` that value is not an
integer and the least conceivable child is `3n/2 + 1`. The least children found over all tournaments of order at
most 9 (`out/lift_gap_order9.txt`, line `M`): `n = 1 -> 2` (the only case on the floor), `2 -> 4`, `4 -> 7`,
`5 -> 9`, `6 -> 11`, `7 -> 12`, `8 -> 13`, `9 -> 15`, `11 -> 18`, `12 -> 21`, `13 -> 22`, `14 -> 22`, `15 -> 25`,
`16 -> 25`, `17 -> 28`.

**Corollary (2), PROVED.** Let `T'` contain `T`. Then `H(T')/H(T)` is 1 or at least `3/2`. Since
`d < (4d-3)/3 < (4d+1)/3 < (3d+1)/2` for every `d > 3`, no realiser of a rising target contains a realiser of its
descent value. *Whatever reverses the two rising branches must delete or rewire something*, as the explicit
families of `inverse_fences_landing_20261116` do. Census agreement (`out/lift_gap_order9.txt`): the branch
`(4d-3)/3` is lifted by one vertex for 1 of 512 values of `d` (the fixed point `d = 3`), the branch `(4d+1)/3` for
0 of 506.

The same floor sharpens `triangular_backbone_chanpak_20261116/FACTORED_GAPS.md`: its exact formula for the first
induced gap of an ordinal sum is bounded below there by `2H log_3 M / M`; with the floor every component gap is at
least `(h_i + 1)/2`, so the first positive induced gap is at least `(H + 1)/2` with no hypothesis on the components.

**The halving branch (3), FINITE-EXACT.** If `T'` contains `T` with `H(T') = 2H(T) - 1`, at most one step of any
chain of one-vertex extensions from `T` to `T'` can be strict (two strict steps give at least `(9d+5)/4 > 2d - 1`),
so `H(T)` is the count of a one-vertex deletion of a realiser of `H(T')`. One-vertex parents are therefore the
whole story. A realiser of `h <= 223` has strong components of order at most 10 (a strong tournament of order 11 has
`H >= 225`), and a product parent `q·d'` cannot equal `(h ± 1)/2`. The census of every tournament of order at most
10 and every vertex (`scripts/lift_par.c`, closed under ordinal products by `scripts/parents_closed.py`) is thus
complete for `h <= 223`:

* `(h+1)/2` is a parent of `h` for 38 of the 54 counts `h = 1 (mod 4)`; the 16 failures are listed in row 3.
  For 14 of them the descent value is a count, and still no realiser of `h` contains a realiser of it.
* neither `(h-1)/2` nor `(h+1)/2` is a parent for 29 counts: 13, 15, 27, 41, 43, 51, 57, 63, 69, 77, 85, 87, 93,
  105, 109, 125, 133, 141, 147, 151, 153, 167, 169, 171, 173, 175, 187, 201, 223.

## 3. The ladder at the floor

**Port quantities.** For a tournament `A` and a vertex `u`: `E_A(u)`, `S_A(u)` count the Hamiltonian paths ending,
starting at `u`; `U_A(u)` counts ordered pairs `(P1, P2)` of disjoint directed paths covering `A` with `P2` starting
at `u` and `P1` possibly empty (no condition between them); `U*_A(u)` the pairs with `P1` ending at `u` and `P2`
possibly empty. Cutting a Hamiltonian path just before (after) `u` is an injection into these pairs, so

    D_A(u) = U_A(u) - H(A) >= 0        counts the pairs with P1 nonempty and u -> last(P1),
    D*_A(u) = U*_A(u) - H(A) >= 0      counts the pairs with P2 nonempty and first(P2) -> u,

and the pairs with `P2 = (u)` give `D_A(u) >= H(A - u) - E_A(u)`.

**Bridge.** Let `T = A ⟲ B` be the ordinal sum `A -> B` with the single arc `u -> w` reversed (`u` in `A`, `w` in
`B`), and write `a = H(A)`, `b = H(B)`. A Hamiltonian path that avoids the arc `w -> u` is a path of `A` followed
by a path of `B` whose junction is not `(u, w)`. One that uses it reads `P1, (Q1 w), (u P2), Q2`, and every other
junction is a forward arc. Hence

    H(A ⟲ B) = ab - E_A(u) S_B(w) + U_A(u) U*_B(w).                                            (G2)

With `B` a single vertex this is the tail formula `H(A + q) = 2a + D_A(u) - E_A(u)` for a vertex `q` whose only
out-neighbour is `u`.

**Middle vertex (5), PROVED.** Add `x` with `A -> x -> B`. Without the arc `w -> u` the vertex `x` must sit between
the two blocks: `ab` paths. With it, `x` has exactly two legal places, after `P1` or after `u P2`. So

    H(A ⟲ B + x) = ab + 2 U_A(u) U*_B(w).                                                      (G4)

Put `Δ = D_A(u) b + D*_B(w) a + D_A(u) D*_B(w) >= 0`. Then

    H(T) = 2ab + Δ - E_A(u) S_B(w),      H(T + x) = 3ab + 2Δ,      2H(T + x) - 3H(T) = Δ + 3 E_A(u) S_B(w).

`scripts/port_calculus_check.py` verifies (G2), (G4), the tail formula and the glue formula below on 400 random
pairs of blocks by brute force (1,600 identities, no failure).

**Ladder (6).** `Δ = 0` and `E_A(u) S_B(w) = 1` give the pair `(2ab - 1, 3ab)`. By the inequality above,
`D_A(u) = 0` and `E_A(u) = 1` force `H(A - u) = 1`, so `A` is one-apex with `A(w) + 1 = h(w)`, which holds exactly
for the fan `1 0^r` (TWO_APEX (5)) or a single vertex, up to a transitive block placed before it, which changes
no count; dually `B` is the dual fan `1^s 0` or a single vertex. So this mechanism gives exactly `m = φ_a φ_b`:

    T(a,b) = fan ⟲ dual fan,   H = 2m - 1,   n = m - 1 = 2k;        T(a,b) + x,   H = 3m,   n = 3k + 1      (m = 2k + 1).

`scripts/near_floor_family_check.py` builds these tournaments literally for `a, b <= 4`. For `a = 0` the bridge is
the dual fan plus one vertex, so a realiser of `m` sits inside: the chain of half-counts is `k -> 2k -> 3k + 1`, and
the triangle product `n -> 3n + 1` is two one-vertex extensions, the second one step above the floor. The same
three-term chain exists for `min(a,b) = 1`, `max(a,b) <= 3` and does not exist for `(a,b) = (2,2), (2,3), (3,3)`
(`out/ladder_parents.txt`).

**Census (7).** `out/lift_gap_order9.txt` (line `G`) and `out/gap3_wit_analysis.txt`: among all one-vertex
extensions of all tournaments of order at most 9, the nontrivial ones with `2H' - 3H <= 3` are `3 -> 5`
and `2m - 1 -> 3m` for the twelve `m` of row 7. The witnesses are the fans `1 0^a 1 0` (33 -> 51, 65 -> 99,
129 -> 195), the bridged pairs above (17 -> 27 at order 7, 29 -> 45 at order 8, 49 -> 75 and 53 -> 81 at order 9)
and the triangle multiple of `3 -> 5` (9 -> 15). The first product of three such numbers that is not a product of
two is `m = 75`; no pair `149 -> 225` exists up to order 10.

**Glue formula (context, PROVED by the same splitting).** For `A`, `B` identified at `u = w` with
`A - u -> B - w`:

    H = ab + (H(A-u) - E_A(u))(H(B-w) - S_B(w)) + D_A(u)(b - H(B-w)) + D*_B(w)(a - H(A-u)) + D_A(u) D*_B(w),

every term nonnegative; so gluing is multiplicative exactly at ports with `D_A(u) = D*_B(w) = 0`. This is the
split surplus of `chanpak_second_gap_20261115` in closed form; no priority claim.

## 4. Parents and near-doubling

`out/parents_closed_223.txt` lists the complete parent set `P(h) = { H(T - x) < h : H(T) = h }` for the 109 counts
`h <= 223` (all odd numbers except 7 and 21). The number of parents has minimum 1 (`h = 3`), median 28, maximum 67
(`h = 219`); the fewest are 5, 11, 13 (two each), 9 and 47 (three). The largest parent is at most `0.6615 h`
(`129 -> 195`), in line with the floor. `P(13) = {3, 5}`, `P(27) = {1, 3, 9, 15, 17}`,
`P(41) = {3, 5, 9, 11, 13, 15, 17, 23}`.

**Near-doubling (4).** `min |h - 2d|` over `d` in `P(h)` is 1 for 79 of the 108 counts `h >= 5`, 3 for 27 and 5 for
41 and 125. In half-counts: `n = 2n' + j` with `j` in {0, 1} (defect ∓1), {-1, 2} (∓3) or {-2, 3} (∓5), a binary
expansion with a few extra digits, each digit a real vertex. For a tail vertex at a port the digit is read off the
port: `H(A + q) = 2a + D_A(u) - E_A(u)`. Beyond 223 the census (`out/lift_near_order9.txt`, realisers of order at
most 10 only, so upper bounds) gives defect at most 7 for `225 <= h <= 3357` (1,567 counts: 1,351 with defect 1,
194 with 3, 18 with 5, and 399, 1285, 2005, 2569 with 7). Every odd number up to 4409 other than 7 and 21 has an
order-10 realiser with such a parent; above 3600, at the edge of that census, the defects found reach 13
(`out/near_doubling_summary.txt`).

The strongest targets by this measure are 41 and 125 (defect 5), then the 27 counts of defect 3; by scarcity of
parents 11, 13 and 47.

## 5. The tree of the signed descent

In half-counts the children of `n` are `2n` and, when `3` does not divide `n`, `r(n) = (4n - ε)/3` with `ε = ±1`,
`ε = n (mod 3)`.

* **Dead rays (PROVED).** If `3 | n` the only child is `2n`, again divisible by 3. So the basin of 7 (`n = 3`) is
  `{3·2^k + 1}`, and one third of all odd numbers lie on such rays.
* **Rising chains (PROVED).** `r` maps the classes mod `3^j` inside a unit class mod 3 bijectively onto the classes
  mod `3^(j-1)` (`n = ε + 3i -> ε + 4i`). Hence the `n` that carry at least `j` consecutive rising steps form
  exactly `2^j` classes mod `3^j`: density `(2/3)^j`.
* **Level growth (PROVED on average).** Let `b_k(n)` be the number of units mod 3 at level `k` above `n`. Then
  `b_(k+1)(n) = b_k(2n) + [3 ∤ r(n)] b_k(r(n))` and the Haar average over the 3-adic units is exactly `(5/3)^k`.
  From 21 and from 27 the ratio of consecutive level sizes is 1.6667 at level 30 (FINITE-EXACT).
* **Criticality (PROVED identity).** Weighting a child by the reciprocal of its size ratio gives
  `1/2 + (2/3)(3/4) = 1`: the over-branching `1/2 + 3/4` is cancelled exactly by the rate 1/3 at which rising
  children die. With any other death rate `p` the basins would have dimension `α ≠ 1`,
  `2^(-α) + (1 - p)(3/4)^α = 1`.
* **Records (FINITE-EXACT, `n < 2·10^6`).** The least `n` carrying more rising steps than every smaller one:
  2 (1 step), 4 (3), **10 (9)**, 244 (15), 2660 (17), 4648 (18), 5212 (24), 34628 (27), 58130 (28), 502750 (30),
  737596 (40). The count 21 is `n = 10`: its chain is 21, 27, 35, 47, 63, 83, 111, 147, 195, 259, and every
  depth record of `D` from 35 to 259 lies on it. The two exceptional counts sit at the two extremes of the tree.
* **Landing multiplicity (PROVED / NUMERICAL).** `h/2 < D(h)` and `D(h) - 3 <= (3/4)(h - 3)`, so an orbit that
  starts above `2Y` meets the dyadic window `[Y, 2Y)`, `Y >= 12`, in 1, 2 or 3 points; a landing point has one dipper (dead ray) or two.
  Sampled frequencies 0.648, 0.293, 0.059, mean 1.412 against `1/(1/2 + (1/2) log_2 (4/3)) = 1.413`. This is the
  contrast with THM-4506, where the worst case grows like `0.631 (k - D)`: for a monotone descent the unbounded
  quantity is the basin, not the multiplicity.
* **Basin densities (NUMERICAL, `X = 2·10^7`, `out/dbasin_2e7.txt`).** `basin(t) = {t} ∪ basin(2t-1) ∪ basin(rise t)`
  exactly. Fractions of odd numbers: 0.4811 through 21, 0.3676 through 27, 0.1136 through 41. For units `t`,
  `F(t) = t ρ(t)` ranges over [3.03, 10.16] for `t <= 2000` with mean 6.129; the renewal value is
  `3 / (ln 2 · (1/2 + (1/2) log_2 (4/3))) = 6.117`. Asymptotically `F(n) = F(2n)/2 + (3/4) F(r(n))`, and with
  the symmetry `n -> -n` this forces the means over `n mod 9` to be `9/14` (classes 2, 7, whose rising child is
  dead), `15/14` (1, 8) and `9/7` (4, 5) times the mean: 3.93, 6.55, 7.87; measured 3.94, 6.63 and 6.53, 7.87.
  So `ρ` is a function of the 3-adic expansion of the half-count, and 21, 27 are heavy (`F` about 10).

## 6. The dominance-pruned slow-factor tree

Setting of THM-5252: span-3 ribbon words, the 119-state transfer, its minimal 61-dimensional representation. A
factor `u` is *slow* at rate `ρ` if its operator has rank at least 2 and `H(T[u']) < ρ^|u'|` for every prefix; a
finite tree of slow factors at some rate above `5^(1/3)` is the literal form of HYP-9742. A node is
(column space of the operator of `u`, raw state vector of `T[u]`); a node dominates another at the same depth if
its column space contains the other's and its state is coordinatewise smaller, since the raw matrices are
nonnegative. Three versions at rate 1.72 (column spaces modulo a prime, exploratory):

| depth | slow factors, duplicates removed (`domanalyse.py`) | essential-suffix classes (`domtree.py`) | inclusion dominance, propagated (`domtree2.py`) |
|---:|---:|---:|---:|
| 5 | 21,372 | 528 | 344 |
| 6 | 82,514 | 1,590 | 837 |
| 7 | 329,834 | 6,110 | 2,616 |
| 8 | | 16,088 | 4,971 |
| 9 | | 43,825 | 11,647 |
| 10 | | 157,910 | |

Inclusion dominance removes all but about one node in 100 of the raw level, but the surviving level still grows
by about 2.3 per letter: the column spaces keep shrinking along a word (ranks 13, 17, 48, 57 for single letters,
3 and 4 for most nodes at depth 9) and smaller spaces are incomparable. That growth cannot be outrun. The least
`H` of a rank-two factor is 1 up to 8 vertices, 3 at 9 to 11 and 5 at 12 (THM-5252), so even if the true rate is
`ρ* = 1.770` with the constant fixed by `g_12 = 5`, a tree at rate 1.72 closes only near depth
`ln(ρ*^12 / 5) / ln(ρ*/1.72)`, about 180 letters. **The certificate is not built, and a tree enumerated level by
level cannot build it**: the tree has to be folded into a finite graph, which needs a comparison between nodes
of different column spaces (a potential on the pair, not an order on the state).

## 7. Fence and square-sum readings (ANALOGY)

* Friedman's fence problem bounds each field's area *above* by 1 and asks for the largest total; its lattice budget
  is `2k + ⌈2√k⌉` fences for `k` unit cells (proved for axis-parallel fences in
  `inverse_fences_landing_20261116/FENCE_AND_SQUARE.md`), a two-dimensional boundary term `2√k`. Tournaments bound
  each new vertex's contribution *below* (the floor) and the boundary term of a concatenation is `O(1)`
  (THM-5252): a one-dimensional isoperimetry. In the table of the page (read 2026-11-15) the best known area
  equals the lattice value at `n = 4, 5, 7, 12` and exceeds it by less than 0.05 at the rectangle budgets 17, 24,
  31, 40.
* The square-sum graph on 1..15 has exactly two vertices with a single neighbour, 8 and 9, which are therefore
  the ends of its Hamiltonian path (THM-4510). The ladder of section 3 needs the same object: a block
  with exactly one Hamiltonian path ending at its port. Such blocks are the fans, and the pairs on the ladder are
  exactly the pairs of fans.

## 8. Open

1. HYP-9812 (a): every count has a parent within a bounded defect of half of it. (b): `2H' - 3H = 3` only for
   `H = 2 φ_a φ_b - 1`. Test of (b): `149 -> 225` at order 11.
2. Existence of the basin densities and of `F` as a function on the 3-adic units (a renewal theorem for a tree
   whose branching is read off 3-adic digits).
3. A dominance relation between different column spaces (section 6).
4. Chan–Pak.

## Files

`scripts/`: `domtree.py`, `domanalyse.py`, `domtree2.py` (need `lit_base.pkl` from
`chanpak_waves_20261113/scripts/w39_lit_setup.py`); `lift_census.c`, `lift_gap.c`, `lift_par.c`, `lift_near.c`,
`gapwit.c` (input: `dyadic_cones_two_block_20261111/out/classes_order8.txt`; usage `prog 8 1 DMAX < classes`);
`gapwit_analyse.py`, `gapsearch.c`, `parents_closed.py`, `port_calculus_check.py`, `near_floor_family_check.py`,
`ladder_parents.py`, `dtree.py`, `dbasin.py`. `out/`: the outputs quoted above.
