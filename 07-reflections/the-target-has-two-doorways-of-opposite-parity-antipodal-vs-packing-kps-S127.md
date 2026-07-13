# The target box has two doorways of opposite parity: antipodal (tight) vs packing (covering)

*kind-pasteur-2026-07-11-S127 cont.56. Owner: "get understanding of the shape related to the target and come up
with creative geometric/topological arguments toward proofs." Reading the fleet's geometric picture
(opus-S250 archimedean-vs-finite-places, klein-S269's three faces with the order-2 antipodal tournament lever,
opus-S252 Chebyshev-equioscillation), I worked the shapes and found a parity duality the fleet had not
connected — and it says something sharp about which tool can reach the open problem.*

---

## The target, as a shape

LRC(14) is: **the closed geodesic in direction `v = (v₁,…,v₁₃)` on the 13-torus `Tⁱ³` hits the lonely box**

> `L = [1/14, 13/14]¹³` — a 13-cube of side `6/7`, centered at `(½,…,½)`.

`M(v)` is the largest `μ` such that the geodesic reaches the shrunk box `[μ, 1−μ]¹³` — the `L^∞` inradius of the
geodesic relative to the coordinate hyperplanes `xᵢ = 0`. The two-bucket dispatch (THM-366) is really a
statement about **which doorway (rational base `q`) the geodesic uses to enter `L`**:

- **non-DC / tight bucket** (misses a multiple of some `d ≤ 14`): enters at a **small** base `q = d`. The
  hardest one — the AP `{1..13}`, which misses only `d=14` — enters at `q = 14`. *Proved:* `M ≥ 1/14`.
- **DC / covering bucket** (has a multiple of every `d ≤ 14`): every small base is blocked (a runner sits at
  `0`), so it enters at a **large** base. The hardest one — the deep well `{1..12, 182}` — enters at
  `q = Φ₆(14) = 183`. *Open:* `M ≥ 14/183` (klein's ILP certifies `speeds ≤ 182`).

## The duality: the two doorways have opposite parity

Line up the two doorways and their symmetry:

| bucket | hardest family | doorway base `q` | parity | grid symmetry available |
|---|---|---|---|---|
| non-DC (tight) | AP `{1..13}` | `N = 14` | **even** | `x ↦ x+½`: **7 antipodal pairs** `(0,7),…,(6,13)` |
| DC (covering) | deep well `{1..12,182}` | `Φ₆(N) = 183` | **odd** | **none** — 0 antipodal pairs |

The right column is forced, not incidental:

- **`Φ₆(N) = N²−N+1` is always odd** (`N(N−1)` is even), for *every* `N`. So the covering doorway *never*
  admits the order-2 antipodal symmetry `x ↦ x+½` — there is no pair of grid points at distance exactly `½`.
  Verified `N = 3..20` (all odd), and the deep-well residues at `N = 8,10,12,14` carry **0** antipodal pairs.
- **The even doorway `N=14` carries exactly `N/2 = 7`** antipodal pairs `(j, j+7)` — precisely klein-S269's
  tied diameter arcs. At the odd tight case (`N` odd, e.g. `13`) there are `0`, and the AP is instead the
  rotational tournament `R_N`.

So the same word "regular extremal circulant" splits by parity: **even doorway → antipodally-symmetric (the 7
tied arcs, an order-2 degeneracy); odd doorway → rotational `R_q` regularity.** An individual DC family can
clear at some even base, but then it is *loose* (`M > 14/183`, not binding); the **binding case — the
covering-minimum itself — sits at the odd base `Φ₆(N)`**, and that is the case a proof must handle.

## The shape of each doorway

**Even doorway (tight, base 14).** The 14 grid points `{0,…,13}/14` fold into 7 antipodal pairs. The AP config
*wants* to be `x ↦ x+½`-symmetric, but the winding tournament (THM-373) has **odd `|Aut|`** — two arc-states,
no order-2 automorphism — so it cannot realize the symmetry; the 7 ties must resolve, and **every one of the
`2⁷` resolutions gives `M ≥ 1/14`** (klein-S269's load-bearing lever). This is a genuine Borsuk-Ulam-flavored
obstruction: the barrier is an *unrealizable symmetry*.

**Odd doorway (covering, base 183).** No antipode exists to obstruct. The 13 residues `14·{1,…,12,−1} mod 183 =
{14, 28, …, 168, 169}` are a **one-sided band-packing**: a scaled AP `{14i}` filling the band `[μ, q−μ] =
[14, 169]`, with only the two AP-endpoints `(14, 169)` related (by reflection across `0`, *not* antipodally).
The floor `14/183` is a *packing* fact — the tightest the AP-plus-outlier fits the band — driven by
`13μ ≤ q` (`13·14 = 182 ≤ 183`, tight), not by any symmetry.

## Toward proofs — the payoff, and the honest limit

The parity duality is not decoration; it says **where each tool can and cannot reach.**

> The order-2 antipodal / Borsuk-Ulam lever is intrinsically an **even-base** tool. The AP it obstructs is in
> the **non-DC bucket — which THM-366 already closes elementarily.** The genuinely **open** part of LRC(14) is
> the DC bucket, and its **binding case — the covering-minimum — enters at the odd base `Φ₆(N)`**, where the
> antipodal symmetry *structurally does not exist*. So the beautiful antipodal lever, however sharp, sits on
> the **already-closed side** of the two-bucket divide and cannot by itself reach the open problem.

Concretely this **redirects effort**: do not try to push the antipodal/tournament-`|Aut|` argument onto DC
families — there is no `x↦x+½` at odd base for it to bite on. The DC bucket must be closed by the **odd-base
tools**: three-gap band-packing (the residues fill `[μ,q−μ]`), rotational-regularity rigidity at `R_q`, and the
`13μ ≤ q` packing inequality — the shapes that actually live at the odd doorway. This is the same lesson as
opus-S250's archimedean/finite-place coupling, read through parity: the tight archimedean scale `1/14` couples
to an **even** finite place for the tight bucket and to an **odd** one for the covering bucket, and only the
odd-place packing story is still open.

What this does *not* do: it does not prove `DC ⟹ M ≥ 14/183` (still klein's finite ILP for `q ≤ 183`). It
sharpens the *map* — it explains why two different geometric obstructions (antipodal vs packing) are needed,
pins each to a parity, and rules the symmetry tool out of the open bucket on structural grounds.

*Files: lrc14_two_doorways_geometry_kps_S127.py (+.out). Builds on klein-S269 (antipodal tournament lever, the
7 tied arcs), opus-S250 (archimedean/finite-place coupling), the ladder/lcm reflections (Φ₆(N)=(N−1)N+1 the
denominator, (N−1)N=lcm the outlier). HYP-6220.*
