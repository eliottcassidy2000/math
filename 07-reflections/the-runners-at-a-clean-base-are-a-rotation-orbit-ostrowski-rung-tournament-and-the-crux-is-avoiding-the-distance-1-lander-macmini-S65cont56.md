# The runners at a clean base are a rotation orbit: Ostrowski rung = three-gap orbit = rotational tournament, and the crux is avoiding the distance-1 lander

*mac-mini-2026-07-09-S65 (cont.56). Synthesis of the owner's base-change reframing with
opus-S250's dictionary, the mac-mini S38 Ostrowski ladder, the cont.54–55 Farey/crux
picture, and the opus-2026-07-01-S14 heptagon-tournament bridge. Everything below is
computationally verified (`04-computation/lrc14_rotation_orbit_tournament_macmini_S65cont56.py`,
`lrc14_crux_mechanism_macmini_S65cont56.py`).*

---

## The owner's reframing, made precise

The owner offered: *a runner's position is an infinite base-`b` expansion `0.d1 d2 …`; each
runner changes by a fixed amount, and there is some base where that amount is "the unit," so
some base makes each perspective cleanest — the problem is a change of bases.*

opus-S250 supplied the dictionary that makes this exact:

> **base = witness modulus `q`. digit of runner `i` at multiplier `a` = residue `v_i · a mod q`.
> a "clean base" for runner `i` = one where its residue is a unit (nonzero, far from 0).
> `M(v) = max_q best_q(v)`, `best_q(v) = max_a min_i ‖v_i a / q‖`.
> LRC(14) ⟺ some base cleans **all** runners simultaneously to margin ≥ 1/14.**

This reflection adds the missing geometric object and closes the loop with the tournament
and Ostrowski threads: **at the optimizing base, the runners are the orbit of a single
rotation, and `M` is that orbit's closest approach to the origin.**

---

## One picture: the runners at their clean base are a rotation orbit

Fix a family `v = {v_1, …, v_13}`. At base `q`, multiplier `a`, runner `i` sits at residue
`r_i = v_i a mod q` on the circle `Z/q`. The margin at that base is

> `best_q = (1/q) · min_i dist(r_i, 0)` — the **closest approach** of the residue set to 0.

`M(v) = max over bases of this closest approach`. So LRC is a **covering-radius / closest-
approach problem for rotation orbits**: find the base whose orbit stays farthest from 0.

Why "rotation orbit"? Because `r_i = a · v_i` is `a·(·)` applied to the speed set — a single
rotation of the speed lattice into `Z/q`. When the speeds are consecutive (`{1,…,13}`) the
residues are a genuine arithmetic orbit `{a, 2a, …, 13a}`, i.e. the first 13 iterates of the
rotation by `a/q`. This is exactly the classical setting of the **three-gap (Steinhaus)
theorem**, which is why the extremal families are three-gap-regular
(cont.44 — that regularity is now *explained*, not just observed).

---

## The Ostrowski ladder IS the sequence of clean bases

mac-mini S38 found the covering-min rungs `M_k = k / (13k + 1)`. Read through the rotation
lens they are transparent:

> **Rung `k`: clean base `q = 13k + 1`, rotation angle `a/q = k/q = M_k`.**
> The extremal family's residues are the orbit `{k · v mod (13k+1)}`.

Verified endpoints:

| rung `k` | base `q=13k+1` | family | residues (rotation orbit) | three-gap lengths | `M_k` |
|---|---|---|---|---|---|
| 1 | 14 | AP `{1..13}` | `{1,2,…,13}` | `{1, 2}` (mult 12,1) | `1/14` |
| 14 | 183 | deep well `{1..12,182}` | `{14,28,…,168,169}` | `{1, 14, 28}` (mult 1,11,1) | `14/183` |

Both orbits have **≤ 3 distinct gap lengths** — the three-gap theorem holds on the nose. At
rung 14 the gaps `{1, 14, 28}` show the mechanism visibly: eleven gaps of 14 (the consecutive
multiples `14·1 … 14·12`), one gap of 1 (`168 → 169`, where the far element 182 lands), and
one double gap `28 = 14+14` (the wrap `169 → 14`, where runner 13 is **missing**).

The whole ladder is realized by the **core-plus-multiple-of-13** families: the census of
`{1,…,12, f}` (`lrc14_landerdodge_floor_macmini_S65cont56.py`) gives
`M({1,…,12, 13k}) = k/(13k+1)` **exactly** — rung `k` is `f = 13k`. So `1/14` (`f=13`, AP),
`2/27` (`f=26`), `3/40` (`f=39`), `4/53`, … ascending to `1/13`. Two facts fall out:

- **The AP is the unique minimizer among core families.** `M` over `{1,…,12,f}` is minimized
  at `f = 13` (the AP), value `1/14`; every other core family is strictly above. This is the
  local, computational face of "the AP is the LRC(14) extremal."
- **Off-ladder is easy.** If `f ≢ 0 mod 13`, then base 13 alone gives `M ≥ 1/13` (the core
  `{1,…,12} mod 13` has closest approach 1, and `f` doesn't land on 0). Only the exact
  multiples of 13 stay on the delicate ladder.

The covering condition is `14 | k` (S38). Ladder rungs `k = 2` (`2/27`) and `k = 3` (`3/40`)
sit *below* `14/183` in value but are **non-covering** (killed by a `t = 1/q` sieve), so they
do not bound the covering class. (Distinct from the cont.54 Farey *mediants* `3/41 = 1/14 ⊕
2/27` etc., which interleave the ladder rungs on the Stern–Brocot tree but are yet other
families.) The smallest **covering** rung is `k = 14`, value `14/183` — this is why "covering
forces the jump from rung 1 to rung 14," and why the crux is `14/183`.

---

## The archimedean ↔ finite bridge opus-S250 flagged is the continued fraction

opus-S250 named the deep obstruction: base-change is a **finite-place** move (choose a
modulus), but the 1/14 loneliness band is an **archimedean** fact. A pure `p`-adic argument
can't see the band; a pure real argument can't see the moduli. LRC couples the two places.

The Ostrowski ladder is exactly that coupling made into an identity:

> **Each finite base `q = 13k+1` delivers an archimedean margin `k/q`, and
> `sup_k k/(13k+1) = 1/13`.**

The rungs `1/14, 2/27, 3/41, …, 14/183, …` are the convergent-margins of a single continued
fraction ascending to the archimedean limit `1/13`. The finite places (the bases `13k+1`) and
the archimedean place (the limit `1/13`, the band `1/14`) are reconciled **by the continued
fraction itself**. This is the concrete form of "the problem is a change of bases": the
sequence of clean bases is the Ostrowski/CF expansion of the loneliness threshold.

---

## The crux mechanism: to beat 1/14, omit the distance-1 lander at your best base

The rotation view exposes the single mechanical reason the deep well beats the AP. At base
`q`, multiplier `a`, the residue that lands at **distance 1** from 0 is `v ≡ ±a^{-1} mod q`.
A family's margin exceeds `1/14` only if it places **no** runner within `q/14` of 0 — in
particular it must **not contain the distance-1 lander**.

Verified:

- **AP `{1..13}` at its best base `q=14, a=1`:** distance-1 landers are `v ≡ {1, 13} mod 14`.
  **Both are present** (runners 1 and 13). They pin the closest approach at exactly 1 ⟹
  `M = 1/14`. The AP is tight *because it contains its own distance-1 landers.*

- **Deep well `{1..12, 182}` at its best base `q=183, a=14`:** distance-1 landers are
  `v ≡ {13, 170} mod 183`. **Neither is present.** The nearest runners, `v=1` and `v=182`
  (`≡ -1`), both land at distance **14** ⟹ `M = 14/183`. The deep well is divisor-rich
  enough to be covering, yet at base 183 its distance-1 lander is precisely the one speed
  (13) it omits; the far element 182 replaces it and lands 14× farther out.

So the crux is one sentence:

> **Covering forces small speeds in; beating 1/14 forces you to a base whose distance-1
> lander is a speed you don't have. The smallest base where a covering set can arrange this
> is `q = 13·14+1 = 183`, where the lander is `13` (absent) and `1, 182` sit at distance 14 —
> giving the floor `14/183`.**

The compressed sub-case falls out of the same principle: a divisor-*closed* family too tight
to reach base 183 is confined to base ≈ 13, where dodging the lander yields the larger margin
`1/13` (extremal near `{1..12,14}`, scan-confirmed). Compressed → `1/13`, non-compressed →
`14/183`, both `> 1/14`, one mechanism.

---

## The tournament connection: the clean-base orbit carries a rotational tournament, and at `q=14` it is the heptagon

A rotation orbit on `Z/q` carries a canonical **rotational (circulant) tournament**: `i → j`
iff `(r_i − r_j) mod q ∈ (0, q/2)`. This is the tournament side of the base-change picture —
the residues are not just points, they are a tournament, and the *loneliness margin is the
tournament's closest edge to the "tie" antipode.*

At the AP's clean base `q = 14`, reduce by the apex prime (`14 = 2·7`, CRT): the residues
`{1,…,13} mod 7` are **all of `Z/7`**, and the induced rotational tournament is the **full
heptagon `R_7` = "beat the next 3"** — self-complementary, `Aut = C_7`, with exactly **14
cyclic triangles = `|D_7|`**. Computationally recovered here (13 vertices → 7 classes → `R_7`,
14 triangles), it is *precisely* the opus-2026-07-01-S14 heptagon tournament, but now
*derived* as the clean-base rotation orbit rather than posited. The heptagon's `D_7` of order
14 is the LRC modulus; its self-complementarity is the parity symmetry `t ↔ 1−t` of the
loneliness band.

At the deep well's clean base `q = 183 = 3·61`, the 13 residues form a rotational tournament
with 90 cyclic triangles — a genuine three-gap orbit tournament, no longer the symmetric
heptagon but the same construction one rung deeper. The tournament thread and the runner
thread are the **same object viewed at two Ostrowski rungs**: the self-complementary heptagon
at the tight non-covering rung `k=1`, and the three-gap orbit tournament at the covering crux
rung `k=14`.

---

## What this buys the endgame

1. **A single geometric object** — the closest approach of a rotation orbit — subsumes the
   base-change reframing, the Ostrowski ladder, the three-gap regularity, and the tournament
   bridge. Every endgame sub-thread is a statement about one orbit at one base.
2. **The crux is mechanical**, not mysterious: `14/183` is "the smallest base at which a
   covering set can omit its own distance-1 lander." This is a *finite* certificate obligation
   (check the lander at each candidate base) for an *archimedean* fact (the `1/14` band) —
   exactly the coupling that makes LRC(14) resist single-place proofs.
3. **The compressed/non-compressed reduction is one principle**: dodge the distance-1 lander;
   how far you must travel to a lander-free base sets the floor (`1/13` if confined near base
   13, `14/183` if you can spread to 183).

The next lever is to turn "omit the distance-1 lander" into a clean lower bound over *all*
covering families — i.e. prove that any covering set, at every base `q < 183`, contains a
runner within `q/14` of 0 (forcing `M ≤ 1/14` to be escaped only at `q ≥ 183`). That is the
rotation-orbit form of the non-compressed peel, and it is the same content as the
far-element floor (cont.55) and klein's compressed floor, now stated as a lander-exclusion
count on rotation orbits.

---

*Cross-links: [[the-M-spectrum-is-a-farey-tree-and-the-crux-is-the-top-ostrowski-rung]],
[[the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends]],
[[the-three-gap-regularity-is-why-the-extremals-are-extremal]],
[[the-heptagon-tournament-6-units-plus-antipode-are-roots-of-z7-minus1]],
[[the-coverage-clearing-duality-one-dichotomy-governs-both-endgame-sides]]. opus-S250
supplied the base=modulus dictionary; this reflection adds the rotation orbit, the
three-gap/tournament identification, and the distance-1-lander mechanism.*
