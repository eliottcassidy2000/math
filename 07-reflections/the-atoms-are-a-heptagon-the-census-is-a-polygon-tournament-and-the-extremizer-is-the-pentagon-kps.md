# The atoms are a heptagon, the census is a polygon-tournament, and the extremizer is the pentagon

*kind-pasteur-2026-07-01-S9. The owner's bridge: take the 6 unit atoms `(Z/14)*` at `a/14`, add the vertex `7/14`, build a tournament on 7 vertices (`k ↔ (2k+1)/14`), apply the tiling model, and bring in dihedral groups. Following it lands on a clean structural correspondence between the LRC census and the project's tournament/tiling pillar — and it recasts the finish as an optimization over regular polygons, minimized at the pentagon. Verified computationally.*

## The heptagon (the owner's construction, made precise)

The seven odd residues `{1,3,5,7,9,11,13}/14` are **equally spaced by `1/7`** — a **regular heptagon**, `vertex k ↔ (2k+1)/14`. Vertex 3 `= 7/14 = 1/2` is the owner's "added vertex," the one non-unit (the runner that always half-laps). Its symmetry group is the **dihedral group `D₇` of order `14`** — the LRC modulus is `|D₇|`.

Two group actions live on these seven points, and their tension is the whole additive-vs-multiplicative theme of the project:

- **Additive** — the heptagon rotation `k ↦ k+1` (`+1/7` on the circle) plus reflections = the geometric `D₇`.
- **Multiplicative** — `(Z/14)* = ⟨×3⟩ = C₆` (orbit `1→3→9→13→11→5`), which **fixes the center 7** and permutes the six units. This is the Frobenius-type action, *not* a heptagon symmetry.
- **Inversion** `ι: t ↦ −t` (`a ↦ 14−a`): `1↔13, 3↔11, 5↔9`, fixing `7` — a reflection of the heptagon. Crucially, on the circle→tournament map `ι` **reverses every arc = tournament complementation `T ↦ Tᵒᵖ`** — the merged-metagraph `ℤ₂` of CLAUDE.md, and the Borsuk–Ulam `ι`-odd antipode of the lonely set. (The arithmetic inversion `a ↦ a⁻¹ mod 14`: `1,13` fixed, `3↔5, 9↔11` — the collar-binding involution of HYP-3793, a *different* dihedral element.)

## The tournament on 7 vertices, through the tiling model

Orienting each edge by the forward arc of the heptagon gives the **rotational tournament `C₇(1,2,3)`** (`i` beats the next three); the quadratic residues give the **Paley tournament `C₇(1,2,4)`**. Both are regular (score `3` everywhere). Via the tiling model (fix the base Hamiltonian path `0–1–…–6`; the 15 tiles are the gaps `d≥2`, forward iff `d ∈` the connection set):

| tournament | forward tiles | `H` (Ham paths) | `c₃` | `|Aut|` | self-converse |
|---|---|---|---|---|---|
| rotational `C₇(1,2,3)` | gaps `{2,3}` | **175** | 14 | 7 | yes |
| Paley `C₇(1,2,4)` | gaps `{2,4}` | **189** (the max) | 14 | 21 (`C₇⋊C₃`) | yes |

`H=189` is the project's canonical `n=7` maximum (Paley = the `H`-maximizer). **Both tournaments are self-converse** — i.e. `ι`-invariant — so the dihedral reflection of the atom-heptagon is exactly tournament complementation. The atoms' antipodal symmetry and the tournament's self-complementarity are the *same* `ℤ₂`.

## The general bridge: census class `(Z/N)*` = regular polygon = `D_p` = tournament

The heptagon is one instance of a rule that makes the whole census geometric. A full-orbit census class `(Z/N)*` with `N=2p` is the **regular `p`-gon** of odd residues, with dihedral symmetry `D_p` of order `N`, carrying a `p`-vertex rotational tournament:

| `N=2p` | polygon | `|D_p|` | min `meas(L_C)` | rotational `H(p)` | status |
|---|---|---|---|---|---|
| 6 (`p=3`) | triangle | 6 | `0.04413` | 3 | loose |
| **10 (`p=5`)** | **pentagon** | 10 | **`0.03226`** | 15 | **loose — the extremizer** |
| 14 (`p=7`) | heptagon | 14 | `0` | 175 | **tight boundary** (`M=1/14`, the AP `{1..13}`) |

The mechanism is sharp: a `(Z/N)*` core has covering-min `M=1/N`, so it is **loose iff `N < 14`**. Thus the loose full-orbit classes are exactly `N = 3..13`, the heptagon `N=14` is the *tight boundary*, and `p≥7` is covering. Among the clean polygon classes only the **triangle (`p=3`)** and the **pentagon (`p=5`)** survive as loose, and the pentagon is smaller. Across *all* eleven loose full-orbit classes `N=3..13` (from S8), the minimum is `N=10` — the pentagon — at `313/9702`, and every class clears `1/36`. **The census infimum is the pentagon because it is the largest regular polygon still loose.**

## The finish, made concrete

The tournament/dihedral frame turns the finish into a finite, structured statement:

> `inf meas(L_C)` over 11-cores `=` min over **(a)** the eleven full-orbit polygon classes `N=3..13` — attained at the **pentagon `N=10`**, `313/9702` — and **(b)** the two-clash *sporadic* family (`N=19,17,…`, non-polygonal, 2-atom), which competes at `0.03238` but stays above. Both clear `1/36`.

This is the tournament image of THM-523's tight-locus dichotomy `{AP dilate, Goddyn–Wong sporadic}`: the **full-orbit = circulant/Paley-like regular polygon** (the "AP"), the **two-clash = sporadic non-polygon** (the "GW"). The `{AP, GW}` split at 13 speeds and the `{polygon, two-clash}` split at 11 speeds are the same `ℤ₂`-symmetric, apex-prime-7 dichotomy — one seen in the lonely measure, one in the self-converse tournament.

So the finish is now: **(1)** the loose full-orbit classes are the finite list `N=3..13` (`M=1/N>1/14`), min at the pentagon — *verified*; **(2)** the pentagon core `{1..13}\{6,10}` has `meas = 313/9702 ≥ 1/36` — *exact*; **(3)** the two-clash sporadic family stays above `1/36` — the remaining Diophantine piece (the 11-speed Goddyn–Wong bound). The infinite search has become a two-family check with the polygon side pinned down.

## Honest status

**Verified:** the heptagon/`D₇` structure, the two tournaments and their tiling/`H`/self-converse data, the `N=2p ↔ p`-gon correspondence, the loose-classes-are-`N=3..13` mechanism, and the pentagon as the census minimum (`≥1/36`). **The bridge is a structural correspondence and a finite reframe, not a proof** — the two-clash sporadic family's uniform `≥1/36` bound is still open (the 11-speed analogue of the Goddyn–Wong Diophantine estimate), and I have not found a closed form expressing `meas` as a tournament invariant (`H`, independence polynomial) — the correspondence is of *symmetry and family structure*, not (yet) a numerical identity. But the finish is now concrete: a finite list of polygon classes (pinned) plus one sporadic family (open), unified by the self-converse/`ι` = complementation symmetry.

— Related: HYP-3793 (atoms `=(Z/N)*`, Ramanujan moments), opus/mac-mini HYP-3795/3799 (the loop-map / group-action dictionary — this adds the tournament/tiling pillar and the dihedral `D_p`), THM-523 (`{AP, Goddyn–Wong}` — the dichotomy this mirrors), THM-255/027 (Paley `H`-maximizer, BIBD), HYP-2576 (the difference-winding map *is* the circular tournament), THM-024 (SC tournaments have an involutive anti-automorphism = the `ι` here), the merged-metagraph `G_n/ℤ₂` (complementation), OPEN-Q-108. Scripts: `04-computation/lrc14_heptagon_tournament_tiling_kps.py` (+ `.out`). Extends HYP-3793; no new HYP (converging frontier).
