# The change of base has three faces — arithmetic, tournament, Ostrowski — and they are one object

*klein-2026-07-12-S269. The owner reframed LRC as a change of bases and asked to push it into the
project's tournament machinery and its relation to Ostrowski. opus-S250 developed the arithmetic
face (base = modulus, digit = residue). This adds the two the owner named — the **tournament** face
and the **Ostrowski** face — and the finding is that all three are the same object, sampled at three
resolutions, with the apex prime 7 (via `14 = 2·7`) the seam that runs through every one of them.*

---

## Face 1 — arithmetic (opus-S250, recap)

Sample at `t = p/q`; runner *i* sits at digit `(v_i·p mod q)` on a *q*-cell dial; lonely = all
digits in the middle band, `M = ⌈q/14⌉/q` at the best base. The AP `{1..13}` is clean at base `n=14`
(`t=1/14`, exactly `1/14`); a divisor-complete family has a multiple of 14, so digit 0 is occupied
at base 14 — **base 14 is blocked** — and the extremal deep well `{1..12,182}` is forced up to base
`q = 183 = Φ₆(14)`, giving `14/183`. opus's caveat stands: cleaning the *motion* (choice of base)
distorts the *danger band* (an archimedean fact), which is why "just a change of base" reframes but
doesn't close. What the two new faces add is that the *same* object is a tournament and a continued
fraction, and both make the obstruction — the apex 7 — visible in a new way.

## Face 2 — tournament: at the clean base the winding tournament is a circulant (and the even base degenerates)

At the optimum `t* = a/q` the runners land on the *q*-grid (residues mod *q*), and the project's
**winding tournament** — the THM-373 runner-phase clock, `u→v` iff `frac((s_u−s_v)t) ∈ (0,1/2)`,
observer = speed 0 — becomes a **circulant** on those residues (mac-mini-S57), its invariant the
Rédei count `H = #directed Hamiltonian paths = I(Ω,2)`. Computed exactly:

| family | base *q* | digit gaps | winding tournament | cyclic △ `c₃` |
|---|---|---|---|---|
| AP `{1..13}` | **14 (even)** | {1} (one-gap) | `C₁₄({1..6})` **degenerate: 7 tied antipodal arcs** (each strict outdeg 6) | 154 |
| deep well `{1..12,182}` | 183 (odd) | {1, 14} (two-gap) | non-regular circulant, outdeg split 6/7 | 112 |

The key subtlety the exact computation surfaces: at the **even** base `14 = 2·7`, the apex difference
`7` places the 7 diameter pairs `(0,7),(1,8),…,(6,13)` at distance *exactly* `1/2` — the order-2
antipodal symmetry `x ↦ x+1/2`. Those 7 arcs are **tied**, so this is **not a single tournament**; it
is 7 undecided arcs, `2⁷` honest resolutions. The genuinely clean, maximally-regular object lives one
level **down, at the odd prime base 13**: the AP `{1..12}` at `t=1/13` gives `C₁₃({1..6}) = R₁₃`, the
rotational tournament — vertex-transitive, constant score 6, and the **H-maximizer** at n=13
(`H = 3{,}711{,}175`, essentially unique; *not* literally Paley — the half-turn/Dirichlet circulant
beats Paley for n ≥ 13, kps-S13 / opus). So "regularity is extremal" is exact at the **prime** base;
the **composite** base 14 is where it degenerates into the 7 apex ties.

**Why this is the obstruction — the honest lever.** A tournament has odd `|Aut|` (2 arc-states ⇒ no
order-2 automorphism), so it *cannot* carry the antipodal symmetry `x↦x+1/2`. The 7 tied
diameter-arcs must therefore resolve, and **every one of the `2⁷` resolutions gives `M ≥ 1/14`, never
below** — "you cannot beat 1/14" = "you cannot realize the order-2-symmetric configuration." This is
the load-bearing statement. (It is *tempting* to say "the winding tournament forbids `H = 7`" — Ω=K₃
is impossible, THM-029/200, and `14 = 2·7`, the apex prime — but that avoidance is **vacuous**: every
real tournament avoids `H=7`. The real content is the odd-`|Aut|` / antipodal tie-resolution, not a
mystical `H=7` exclusion. HYP-3099.) The zero-divisor `7` — `2·7 ≡ 0 (mod 14)`, doubling kernel
`{0,7}` — is the same `7` on both sides: it ties the tournament and it blocks the clean base.

**What H does and doesn't certify.** `H` is a *selector* of the extremal iso-class (the tight optimum
is the regular circulant), not a loneliness scalar: over all `t` the achievable H-spectrum does **not**
separate tight from loose (order forgets the distance `M`, mac-mini-S57), and finitely many
iso-classes cannot encode `M` over unbounded speeds (HYP-2926). The one exact reading is **THM-374:
`H=1 ⟺ empty semicircle ⟺ max gap > 1/2`** — the *coarse* `n=2` threshold. And it points the *right*
way for the frame: anchored `1/n` loneliness lives at the **high-H regular pole** (the tight AP row at
`t=1/14` has `H ≈ 2.4·10⁷` with every vertex lonely), while `H=1` is the bunched pole — so the tight
case being the most-regular circulant is consistent, it just isn't measured by a scalar.

## Face 3 — Ostrowski: the clean base is a continued-fraction convergent

The clean bases are not arbitrary moduli — they are the **CF convergent denominators of the optimal
rotation**. This is the 2-term Ostrowski ladder (mac-mini-S38, HYP-4078/3738):

$$M_k = [0;\,n-1,\,k] = \frac{k}{k(n-1)+1},\qquad \text{AP} = \text{rung }k{=}1,\quad \text{deep well} = \text{rung }k{=}n.$$

Verified continued fractions:

| object | `t*` | CF | convergent denominators (= clean bases) |
|---|---|---|---|
| AP (tight) | `1/14` | `[0;14]` | 1, **14** |
| compressed floor (`k→∞` limit) | `1/13` | `[0;13]` | 1, **13** |
| **deep well** (covering-min) | `14/183` | **`[0;13,14] = [0;n−1,n]`** | 1, 13, **183** |

The extremal time's continued fraction **spells out `n−1` and `n`**. Its convergent denominators
`1, 13, 183` are exactly the base tower: base 13 is the compressed floor `1/13`, base `183 = Φ₆(14)`
is the deep well's clean base. The digits I computed at base 183 — `{0,14,28,…,168,169}` — are the
`{kα}` progression `{14k mod 183}` (step 14 = the regular gap `n/Φ₆`) *plus one killer point* 169,
which sits one unit above the top core point 168. That solitary **unit gap is exactly the `+1` in
`Φ₆ = (n−1)n + 1`** — the cyclotomic value is the arithmetic shadow of the killer's unit gap
(HYP-3738, proved n=5..9). By the three-gap (Steinhaus) theorem a `{kα}` set has ≤3 gaps; the AP is
the degenerate **one-gap** (perfectly regular) case, the deep well the **two-gap** case `{1/183,
14/183}` — the `{1, n, 2n}/Φ₆` Ostrowski place-values with the observer at the center of the deep
`2n` hole. LRC lives in the 2-term Ostrowski world the way Zeckendorf/Fibonacci lives in the golden
one `[1;1,1,…]` — same numeration framework, a different continued fraction (ratio 14, not φ).

## The three are one object

The dictionary closes on itself:

> **base *q*  =  the resolution at which the winding tournament is a circulant  =  the
> continued-fraction convergent denominator of the optimal rotation.**

- The AP is: the clean base `n` / the regular rotational tournament `C_n({1..(n−1)/2})` / the
  one-gap `{kα}` at `[0;n]`.
- The deep well is: the forced base `Φ₆(n)` / the non-regular circulant / the two-gap `{kα}` at
  `[0;n−1,n]`.

And the open crux is **one statement in all three languages**: *`covering ⟹ M ≥ 14/183`*
(arithmetic) `⟺` *the extremal optimum tournament is a circulant* (tournament, = the achievable-class
census mac-mini-S57) `⟺` *the extremal residues are a `{kα}` progression, so `g(n) ≤ 3`*
(Ostrowski/three-gap, HYP-2913). The three faces don't multiply the difficulty — they're three
handles on the single fact that the extremal configuration is a rotation orbit. The seam in all
three is the apex prime 7: the zero-divisor blocking base 14, the 7 apex ties / H=7-forbidden, the
even `n` whose CF must route through `2·7`. That is exactly why `n = 13` (prime, a field, a clean
`{kα}`) was tractable and `n = 14` is the first open case — the same apex, seen three ways.

## What it buys, honestly

It doesn't close LRC(14) — opus-S250's archimedean/finite coupling is the genuine content, and
mac-mini-S57's caveat stands (the winding tournament *over all t* forgets the metric; only the
*at-optimum* circulant carries `M`). But it converts the crux into a claim about *rotation orbits*
with three independent certificates to attack — and it says the obstruction is not three separate
walls but one prime, `7`, wearing three masks. The owner's instinct — "the problem is a change of
bases" — lands on the deepest available frame: the right base is the Ostrowski base of the extremal
rotation, at which the runners become a Steinhaus orbit and their comparison graph becomes the
regular circulant.

*Files: `04-computation/lrc14_base_tournament_ostrowski_klein_S269.py` (+out). HYP-6210. Extends
opus-S250 (arithmetic face, HYP-6200) with the tournament + Ostrowski faces; consumes mac-mini-S57
(winding tournament = at-optimum circulant), mac-mini-S38/HYP-4078/3738 (Ostrowski ladder + base
tower), kps-S13 (regularity=extremal / AP=Paley), THM-374 (H as loneliness meter), HYP-2913
(three-gap census), klein-S267/S268 (covering-min 14/183, the blocked witness). Connects
[[LRC14-as-change-of-base-archimedean-vs-finite-places-opus-S250]],
[[the-covering-min-is-an-ostrowski-ladder-and-the-ap-and-deep-well-are-its-ends]],
[[lrc-as-tournament-analysis-the-winding-realization-and-where-the-metric-lives]],
[[the-covering-min-is-the-most-blocked-family-not-the-small-q-one-klein-S268]]. HOUSEKEEPING:
HYP-6200 double-claimed (klein-S268 first-push keeps 6200; opus-S250 → please renumber to 6205);
klein-S269 takes HYP-6210.*
