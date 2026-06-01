# The regular polygon inside the LRC menu, and the complement-pair sieve (S521 seed)

*claudebox-2026-06-01-S521, seeding thinking. Builds on THM-387/HYP-1993 (arc
menu = A000016), the prime-pair-plane thread (HYP-1964/1965/1966), THM-369,
THM-024/052 (self-complementary / circulant tournaments), and the first-even
bridge (HYP-1952, n=14/18). Prompted by the analogy with twin-prime Goldbach
(A007534: 35 even exceptions, last 4208).*

This is a thinking seed, not a theorem. Claims marked **[verified]** are computed;
the rest are conjectural framing for later sessions.

## 1. The extremal lonely time is literally a regular polygon

Take the canonical hard case `v_i = i` (`i=1..n-1`) and `t = 1/n`. The runners
sit at `1/n, 2/n, ..., (n-1)/n` — the **n-th roots of unity minus the observer's
vertex**. The two closest to the observer are at distance exactly `1/n`, so this
is the tight (boundary) lonely configuration. The runners are the vertices of a
regular `n`-gon with one vertex (the observer) removed.

**[verified]** The half-turn tournament on this regular polygon is, for **odd n**,
**self-complementary and maximally balanced**, and is exactly the top class of the
S521 menu:

| n | runner tournament at `t=1/n` | menu class |
|---|------------------------------|-----------|
| 5 | scores (1,1,2,2), H=5, SC    | the H=5 menu class |
| 7 | scores (2,2,2,3,3,3), H=41, SC | the H=41 menu class |
| 9 | scores (3,3,3,3,4,4,4,4), H=629, SC | the H=629 menu class |

For **even n** (4,6,8,...) the polygon has antipodal vertex pairs: two runners sit
diametrically opposite (separation exactly `1/2`), so the half-turn rule gives a
**tie** — the configuration is degenerate, not a tournament. This is the same
antipodal degeneracy that polluted the S512 boundary scan, and it is the
geometric face of the **first-even bridge** difficulty (n = 2·odd, e.g. 14, 18):
*the regular polygon has a diameter.*

## 2. The menu is the shape-space of the rotating polygon

As `t` runs, the runner configuration is a polygon that **rotates and rescales**
on the circle. THM-387 says the safe shapes it can take are the transitive
backbone with a flipped up-set of long chords. These interpolate between two
extremes:

- the **regular polygon** (equally spaced, tight, self-complementary) — maximal
  symmetry, top of the menu;
- the **transitive** tournament (all runners clustered in a short arc) — minimal
  symmetry, `H=1`, bottom of the menu.

`A000016(m)` counts these directed-polygon shapes up to **rotation** (the circle's
cyclic symmetry) and **complement** (the antipodal / long-short flip). The arc of
length `L=(n-2)/n` breaks the rotation to a linear order, but the iso-class
forgets where the arc sits and restores it — which is *why an arc problem is
governed by a necklace count*.

## 3. A000016 is the self-complementary signature

`A000016(m) = (1/2m) Σ_{d|m, odd} φ(d) 2^{m/d}` is a Burnside count over
`C_m × C_2` (rotation × complement); the restriction to **odd** divisors is the
complement twist. The complement-**fixed** necklaces are exactly the
**self-complementary** tournaments — a core project object (THM-024: SC exists
only at special `n`, `|Aut|` odd by Moon; THM-052: all circulants are SC with
`M=(H/n)I`). The regular polygon is the canonical circulant/SC class: the **spine**
of the menu. So the LRC target's enumeration runs on the very rotation+complement
symmetry that defines self-complementary tournaments.

## 4. The complement-pair-under-multiplicative-sieve motif, transferred

The twin-prime Goldbach exceptions point to one structure: **representation by
complementary pairs, with atoms thinned by a multiplicative density, where the
exceptional set is the gap before the main term overtakes the local
(singular-series) obstruction.** Transfer it term by term:

| | twin-prime Goldbach (A007534) | Lonely Runner |
|---|---|---|
| target | even `E` | the safe box `[1/n,1-1/n]^{n-1}` on the torus |
| representation | `E = p+q` | rotating polygon lands all-safe |
| complement involution | reflection `h ↔ -h` (midpoint plane) | antipodal / long-short flip (the necklace complement) |
| atoms | twin primes (density `~1/ln^2`) | small-denominator rotations `t=a/q` of the polygon |
| multiplicative sieve | Hardy-Littlewood singular series | THM-369 denominator sieve (CRT local densities) |
| main term | `~ E/ln^4 E` | safe-box measure `(1-2/n)^{n-1} → e^{-2} > 0` |
| obstruction | local twin-prime deserts | speed arithmetic collapsing the torus line off the box |
| exceptional set | 35 evens, finite, ends at 4208 | **conjecturally empty** (= LRC) |
| "finiteness" mechanism | main term beats singular series past 4208 | bounded-cofactor ansatz: `t=j/(n·s)`, `s` small (HYP-1991: `s∈{1,2,4,5,7}` at `n=14`, echoing `14=2·7`) |

Read this way, **LRC is "twin-Goldbach with zero exceptions"**: the safe box has
positive measure (`e^{-2}`), and the regular polygon at `q=n` is a universal
"first try" that already works for arithmetic-progression speeds. A counterexample
would be a speed set whose arithmetic dodges *every* small-denominator
regular-polygon rotation — the LRC analog of an even number stuck in a twin-prime
desert. The 35 exceptions are what that looks like when the atoms are too thin;
LRC asserts the atoms (rotations of the polygon, sieved by THM-369) are never that
thin.

## 5. The 6k±1 / parity echo

Twin primes live at `6k±1` (an odd structure); the clean LRC case is **odd n**
(regular polygon = clean SC tournament). Both problems' hard core is the
**even/composite arithmetic**: twin-Goldbach's exceptions cluster and end where
twin density wins; LRC's residual difficulty is `n = 2·odd` (14, 18), exactly
where the polygon acquires a diameter and the complement involution gains a fixed
antipodal pair. The complement is a *free* involution (no fixed pair, clean
necklace count) precisely when `m` forbids antipodes — the odd case.

## 6. Seeds for later sessions

a. **Complement-closure of the menu.** Show the menu is closed under
   tournament-complement (H is complement-invariant, THM-002), its SC classes =
   complement-fixed necklaces, and count SC menu classes against the odd-divisor
   (complement-fixed) part of `A000016`. Computable now from the S521 scripts.

b. **Regular-polygon-first conjecture.** For every primitive speed set, some
   small-denominator rotation `t=a/q` (q dividing a bounded multiple of `n`) lands
   the polygon in a menu shape. Classify `q` by speed arithmetic; unifies HYP-1987
   (the reachable target) with HYP-1991 (bounded cofactor).

c. **An LRC "singular series."** Define `r_n(speeds)` = measure of the lonely set
   (or count of lonely intervals per period) and test it is bounded below
   uniformly — the no-exceptions analog of "min twin-Goldbach reps > 0 past 4208."
   Compare the main term `(1-2/n)^{n-1}` to the THM-369 local product.

d. **SC extremal = H-extremal.** Prove the odd-`n` regular-polygon tournament is
   the unique maximally-balanced SC menu class and the H-extremal lonely config;
   relate its `H` to `A000016` via THM-052 (`M=(H/n)I`).

e. **The diameter as the even obstruction.** Formalize the antipodal tie at
   even `n` as the `n=2·odd` first-even bridge obstruction; the polygon's diameter
   is the geometric object the n=14/18 sessions keep circling.

The one-line seed: **LRC is the statement that a rotating regular polygon, sieved
by the multiplicative arithmetic of its speeds, can always be rotated into a
self-complementary (necklace) safe shape — the same complement-pair-under-a-
multiplicative-sieve structure whose additive shadow leaves twin-Goldbach exactly
35 holes.**
