# The E_grid residual is Poisson summation over the grid-relation lattice, and |R| is its kissing number

*kind-pasteur-2026-07-09-S97. Owner: work the critical path, mine past threads for connections,
explore novel work. Mining the repo surfaced mac-mini-S24/S25 (the density floor = Cohn–Elkies
Poisson summation over the relation lattice), which I had not connected to my S96 E_grid route. They
are the **same object** — and the connection turns my open a-priori target `|R| < (6/7)^k` into a
kissing-number bound with the AP as the proven extremal.*

---

## The two routes were one object

**mac-mini-S24** (HYP-4532): the safe/lonely measure is Poisson summation over the **relation
lattice** `L(S) = {a ∈ ℤⁿ : Σ aᵢvᵢ = 0}`:
`safe(S,β) = Σ_{a ∈ L(S)} ∏ᵢ f̂(aᵢ)` — a Cohn–Elkies positivity of a Poisson sum over a lattice.

**My S96 E_grid route** (HYP-5567): existence ⟺ `E_grid[W] > 0`, and by LEM-011
`E_grid[W] = (6/7)^k + R`, `R = Σ_{n≠0, Vmax|n·e} 𝒲̂(n)`.

These are the **same construction**. `E_grid[W]` is the average of `W` over the grid `{j/Vmax}`; by
Poisson summation on `ℤ/Vmax`, that average is the sum of `𝒲̂` over the **grid-refined relation
lattice**

> **`L_V(e) := {n ∈ ℤ^{k-1} : Vmax ∣ n·e}`**   —   `L(e) = {n·e = 0}` plus the wraparound shells
> `n·e = m·Vmax`.

So `E_grid[W] = Σ_{n ∈ L_V(e)} 𝒲̂(n)`, the `n=0` term is `𝒲̂(0) = (6/7)^k`, and **`R` is exactly the
Poisson tail** `Σ_{n∈L_V, n≠0}𝒲̂(n)`. My `R = R_0 + R_wrap` split is precisely `[L(e)] + [wraparound
shells]`. **The E_grid route is mac-mini-S24's Cohn–Elkies picture, grid-refined**: same test function
(`𝒲̂`, the exact transform of the `1/7`-uncovered measure, LEM-011), same lattice-Poisson-positivity
shape, one lattice enlarged by the ruler `Vmax`.

## |R| is the kissing number

**mac-mini-S25** (HYP-4552): the relation lattice's **minimal vectors are the additive triples**
`vᵢ+vⱼ=v_{i+j}` (norm 3), and its **kissing number = the additive energy**, *maximized by the AP*
(the interval has maximal additive closure; Fejér). The floor is "the isolation of the maximal-kissing
relation lattice."

Transported to the grid tail `R`: the short vectors of `L_V(e)` are the **low-height resonances**
`{n : Vmax∣n·e}` — exactly kps-S93's **near-resonance count `Z`**. LEM-011's geometric decay
(`0.371`/coordinate) means the **shortest** resonances carry almost all the mass of `R`. So `|R|`
should be governed by the kissing number. It is — overwhelmingly:

> **`corr(|R|, kissing-count) = 0.998`**,  `corr(|R|, R2) = 0.977`   (kps-S97, `V=105`, energy axis
> Sidon→AP). `|R|/(6/7)^k` runs `0.03` (dissociated, kissing ≈ 1150) up to `0.61` (AP-full, kissing
> 7536). Essentially `|R| ≈ c · (kissing number)`, `c ≈ 8·10⁻⁵·(6/7)^k`.

`R` is not a mysterious signed sum with cancellation (the Mertens fear of the *non-resonant* wall,
kps-S93) — it is a lattice **theta tail**, and its size is a **geometric count of short relations**.

## The a-priori route, named

The open item `|R| < (6/7)^k` (HYP-5567) becomes a **kissing-number bound**:

1. `|R| ≤ Σ_{n∈L_V, n≠0} |𝒲̂(n)|` (triangle inequality — no cancellation needed, unlike the sharp
   Erdős–Turán indicator; this is the *resonant* Mertens-SAFE half, kps-S93).
2. LEM-011: `|𝒲̂(n)| ≤ (6/7)^k · (7/6) ∏(0.371/|n_i|)` — the shortest relations dominate, higher
   shells are geometrically suppressed. So `|R| ≤ c · (kissing number of L_V)` with an explicit `c`.
3. **mac-mini-S25: the AP maximizes the kissing number.** So `kissing(L_V(e)) ≤ kissing(L_V(AP))` and
   `|R| ≤ c · kissing(AP)`. Measured: `kissing(AP) = 7536`, `c·7536 = 0.61 < 1`. **A good period
   exists for every cluster** — including the AP (whose E_grid existence is far weaker than its
   capstone `r_N=1` tightness: `E_grid[W]>0` only asks that *some* `j∈{1..Vmax}` is good, not the
   first `N`).

This is a genuine Cohn–Elkies positivity certificate for good-period existence: bound the Poisson
tail by the lattice kissing number, extremal at the maximal-kissing lattice (AP). It **replaces**
opus-S169's open `#arcs` bound with a bound on the *same additive-energy invariant* the density floor
already uses (THM-660: `Var(W) ∝ R2`), and it is Mertens-safe by construction (resonant modes only).

**Remaining rigor** (the honest gap): making step 2's constant `c` explicit and summing the shells,
and turning mac-mini-S25's "AP maximizes kissing" (proved for the exact lattice `L(e)`) into the
grid-lattice statement `kissing(L_V(e)) ≤ kissing(L_V(AP))` uniformly in `Vmax`. Both are finite,
structured, cancellation-free — Cohn–Elkies LP shape, not an analytic wall.

## Secondary: the smooth-route decay is provably α = 2 (rigorizing opus-S170)

opus-S170 measured maxgap's Fourier decay `α ∈ [1.48, 2.02]` and needs `α > 1` for the absolutely-
convergent (Mertens-free) resonant discrepancy. This is **provable, not just measured**: `maxgap(x)`
and `W(x)` are **continuous piecewise-linear** in `x` (kinks at collisions and `gap=1/7` crossings,
**no jumps** — gaps vary continuously through collisions), so their Fourier coefficients are `O(1/m²)`
(**α = 2**) asymptotically, and `O(1/m)` (α ≥ 1) everywhere by Lipschitz. Verified (kps-S97,
`lrc14_smooth_alpha2`): the **high-m** exponent is `2.00–2.02` for AP, generic, AND the hard
7-structured set; opus's `1.48` is a **pre-asymptotic** average over an `α≈1` shelf (m∈[20,2000])
that only the 7-structured/generic sets exhibit (AP is clean `1/m²` throughout). So opus's smooth
route has α = 2 rigorously; the shelf is why the *small-Vmax* band (resonances land in it) stays tight.

## The unification

Four threads are one:

| thread | object | this session |
|---|---|---|
| E_grid route (kps-S96) | `R = Σ_{Vmax\|n·e}𝒲̂(n)` | = Poisson tail over `L_V(e)` |
| near-resonance count (kps-S93) | `Z = #{low-height Vmax\|n·e}` | = kissing number of `L_V(e)` |
| additive energy (THM-660) | `R2 = Σ r_d²`, `Var(W)∝R2` | `|R| ∝ kissing ∝ R2` |
| Cohn–Elkies floor (mac-mini-S24/S25) | Poisson-positivity over `L(e)`, AP max kissing | grid-refined ⟹ `|R|≤c·kissing(AP)` |

The good-period existence residual `|R| < (6/7)^k` is a **kissing-number isolation bound on the
grid-relation lattice, extremal at the AP** — the exact Cohn–Elkies shape mac-mini-S24 named for the
continuous floor, now carrying the finite-`Vmax` existence.

*Files: `lrc14_R_vs_energy_kps_S97.py`, `lrc14_R_kissing_kps_S97.py`, `lrc14_smooth_alpha2_kps_S97.py`
(+ `.out`). See HYP-5567 (E_grid), kps-S93 (near-resonance count), THM-660 (additive energy),
mac-mini-S24/S25 (Cohn–Elkies relation lattice), opus-S170 (smooth route), LEM-011 (𝒲̂). Related:
[[triangle_foundation]] — the modular home is `X₀(14)`, cusps `1,2,7,14`.*
