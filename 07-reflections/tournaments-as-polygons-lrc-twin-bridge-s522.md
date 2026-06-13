---
source: oracle-2026-06-01-S522
status: reflection (seeding) + one verified computational anchor
tags: [LRC, twin-primes, goldbach, roots-of-unity, regular-polygon, rotational-tournament, paley, A000568, necklace, mod-6]
---

# Tournaments as polygons on the unit circle: the root-of-unity substrate shared by LRC and twin-Goldbach

**Prompt (user):** spend this session thinking how the twin-Goldbach necklace
structure (S521) relates to LRC; think of tournaments as regular polyhedra on the
unit circle.

This is the right frame. Both problems live on a circle, both are an **additive
orbit forced to meet a target that is structured by roots of unity**, and the
"regular polygon" is the common extremal object. One verified computation pins
the dictionary; the rest is the seed.

## 1. The dictionary (tournament = inscribed polygon)

Put `m` runner positions on the unit circle: `x_i = e^{2πi θ_i}`. Orient the
**half-turn tournament** `i -> j  iff  frac(θ_j - θ_i) ∈ (0, 1/2)` — "j is in the
forward semicircle of i." This is exactly the S511/S512 runner clock. Then:

- An arbitrary configuration = an inscribed `m`-gon with a half-turn orientation.
- **Equally spaced points `θ_i = i/m` = the regular `m`-gon = the `m`-th roots of
  unity.**

**Verified (`tournaments_as_polygons_lrc_s522.py`, m=3..7):**
- The regular `m`-gon half-turn tournament **is the rotational (circulant)
  tournament `R_m`** with connection set `{1,…,(m-1)/2}` — checked equal for
  m=3,5,7.
- For **even `m`** the antipodal pairs sit at `frac = 1/2` *exactly* → a TIE, no
  arc. The regular even-gon is **degenerate: a WALL** of the arrangement, not a
  tournament. (m=4 square: 2 tie-pairs; m=6 hexagon: 3 tie-pairs.)
- So among regular polygons, **odd-gons are clean rotational tournaments,
  even-gons are walls.**

## 2. The regular polygon is a reachable LRC *source* — and it is geometric, not Paley

Cross-checking the regular polygon against the S520 LRC reachable source-menu
(HYP-1987):

```
m=7 (n=8): regular heptagon R_7 has H=175, score (3,3,3,3,3,3,3)  -> IN the menu.
           Paley P_7 has H=189, same regular score              -> NOT in the menu.
```

So **the literal regular polygon (evenly spaced runners) is a reachable lonely
configuration**, and the menu distinguishes it from the *arithmetic* regular
tournament (Paley). The LRC clock reaches the geometric regular `m`-gon, not the
quadratic-residue one. This matters: it says the "multiplicatively-shaped target"
of HYP-1987 contains the **rotational** roots-of-unity tournament as its most
symmetric source, while the **Paley** roots-of-unity tournament — the one number
theory loves (THM-369) — is excluded. The arc geometry, not the field structure,
decides reachability.

(Why it should be a source: the classical LRC-extremal speeds `{1,…,n-1}` reach,
at the tight time, equally spaced positions — a regular polygon with the smallest
possible largest-gap `1/n`. The regular polygon is LRC's *barely-lonely* extreme;
finding it inside the reachable menu confirms the menu's extremal vertex is the
roots of unity.)

## 3. The shared template

Strip both problems to one sentence each:

- **twin-Goldbach (S521):** the additive sumset `K+K` must cover every large `m`,
  where `K = {k : 6k±1 prime}` lives on the **mod-6 wheel** — the 6th roots of
  unity, sieved to the two primitive beads `±1`.
- **LRC (S512/S520):** the additive clock-walk `{v_i t mod 1}` must reach a
  **source polygon** in the arc-confined menu inside `A000568(n-1)`, whose
  extremal vertex is the regular `m`-gon — the `m`-th roots of unity, oriented to
  `R_m`.

> **Common form:** an *additive* orbit (a sumset / a linear flow) must intersect
> an *arithmetically thinned, rotation-symmetric* target on a circle of roots of
> unity. The hard core is that the additive object cannot dodge the
> multiplicative/root-of-unity target forever.

Twin-Goldbach is the **finite, solved-shaped** instance: the circle is `Z/6`, the
target `K` is explicit, and the failure set is a finite 11-hole comb ending at
`m=701`. LRC is the **infinite, open** instance: the circle is the real track,
the target is the regular-polygon menu, and "no counterexample" is exactly "the
clock-walk cannot avoid a source polygon."

## 4. The hexagon coincidence (seed, not claim)

The twin-Goldbach wheel is the **hexagon** (`m=6`, mod 6). On the LRC side `m=6`
is `n=7` — **the exact point where the S520 source-menu collapses** (menu =
`A000568(n-1)/2` for n≤6, then 6 « 28 at n=7) — and the regular hexagon is a
**degenerate wall**. Three facts colliding at `m=6`:

1. the hexagon is the twin-prime modulus (`6 = 2·3`, first wheel with two live
   beads `±1`);
2. the regular hexagon is an LRC wall (even-gon, antipodal ties);
3. the LRC reachable menu collapses at `n=7`.

I do **not** claim these are causally one thing. But "the modulus where twin
primes organize" = "the regular polygon that degenerates" = "the `n` where the
LRC menu stops being half of everything" is the kind of triple coincidence this
project treats as a reflection, not noise. The disciplined version of the
question: *does the even/odd (wall/clean) parity of the regular `m`-gon control
the `A000568(n-1)/2 → collapse` transition of the LRC source menu?* (S520 already
shows the collapse is at `m=6`; this gives a candidate mechanism — antipodal
degeneracy killing the complement-pairing that produced the clean `/2`.)

## 5. Seeds to chase (→ HYP-1995)

- **R_m reachability:** is the rotational regular polygon `R_m` a reachable LRC
  source for *every* odd `m`, always the unique regular class in the menu, always
  `≠ Paley_m`? Test by extending S520's menu to m=9,11 and checking `R_m`/`P_m`
  membership. (Anchor verified at m=7.)
- **Complement = antipodal:** the S520 `/2` menu for small `n` smells like the
  complement involution (`T ↔ T^op`), which on the circle is the **antipodal map**
  `θ ↦ θ + 1/2`. Even-gons are fixed by it (→ walls); the `/2` count = orbits of
  antipodal pairing. Does the collapse at `m=6` = the antipodal map acquiring
  many fixed walls?
- **An LRC "necklace reduction":** twin-Goldbach reduced to binary Goldbach on the
  necklace `K` by quotienting the ±2 thickening. Is there an analogous quotient
  that reduces the LRC arc-walk to a covering statement on the rotational
  tournaments `{R_m}` alone (the polygon skeleton), with the non-regular menu
  classes as the "thickening"?
- **Roots of unity as the literal bridge:** both targets are sieved root-of-unity
  sets. Twin primes = `e^{2πik/6}` surviving the `2,3` sieve; LRC sources =
  `e^{2πik/m}` surviving the half-turn arc constraint. Is "survives the arc" a
  sieve in the same Hardy–Littlewood sense that "is twin-prime" is?

## Anchor
`04-computation/tournaments_as_polygons_lrc_s522.py` (+ `.out`): regular `m`-gon =
`R_m`; even-gon = wall; `R_7` (H=175) reachable LRC source, `P_7` (H=189) not.
Builds on HYP-1987 (S512/S520) and HYP-1994 (S521).
