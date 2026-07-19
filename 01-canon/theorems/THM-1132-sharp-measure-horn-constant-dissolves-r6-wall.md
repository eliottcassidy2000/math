---
id: THM-1132
title: The closed-component r=6 measure horn has the optimal threshold 1/(7L); sampled R_sharp is 0.8011 and the uniform tail remains open
status: PROVED local closed-component horn and exact equal-step AP G(sigma) formula; broad finite search gives max sampled R_sharp = 0.8011, but no uniform all-core/all-killer theorem follows
source: death-star-2026-07-18-S58
depends_on:
  - THM-1102   # source of the conservative min(N/(6mu),1/(3L)) coordinate
related:
  - THM-1061   # where 1/(3L) was introduced ("a killer leaves gaps 6/(7k), so 1/(3L)=7k/18")
  - THM-724    # deep-well covering-min; the sharp horn certifies the deep well directly
  - THM-1121   # exact bounded atlas, independent of the local improvement here
  - THM-1134   # exact multiplier cone and all-scale step-two closure
  - THM-1145   # refutes extrapolation from a fixed-width sequential max-T bank
  - THM-1135   # harmonic all-scale finite reduction
  - MISTAKE-168
script:
  - 04-computation/r6_sharp_horn_deathstar_S58.py
  - 04-computation/r6_sharp_horn_robust_deathstar_S58.py
  - 04-computation/r6_sharp_horn_search_deathstar_S58.py
  - 04-computation/r6_Gsigma_exact_bands_deathstar_S58.py
result:
  - 05-knowledge/results/r6_sharp_horn_search_deathstar_S58.out
  - 05-knowledge/results/r6_Gsigma_exact_bands_deathstar_S58.out
---

# THM-1132 — the sharp closed-component horn and the equal-step AP functional

## Setup (kind-pasteur's covering-killer stratum)

A covering 13-family in the r-ladder splits as **core P** (a size-`13−r` subset of `{1,…,12}`)
plus **r killers** `k_1<⋯<k_r` with `k_i > 13·max(P)`. For `r=6`: `P` is a 7-subset (792 of
them), killers `≥ 92`. To prove the family is `1/14`-lonely it suffices to exhibit **one**
real time `t` with `‖v t‖ ≥ 1/14` for all 13 speeds.

## The measure horn, sharpened

Remove the 5 smaller killers exactly:
`E = S(P) ∖ ⋃_{i<6} danger(k_i)`, where `S(P) = {t : ‖p t‖ ≥ 1/14 ∀ p∈P}` and
`danger(k) = {t : ‖k t‖ < 1/14}`. Let `L` = length of the largest component of `E`.

> **Sharp certification (PROVED).** If `L >= 1/(7 k_6)` then the family is `1/14`-lonely.

*Proof.* Each `danger(k_i)` is open and `S(P)` is closed, so `E` is closed. Its connected
components are therefore closed arcs (or points). The set `danger(k_6)` is a disjoint union of
open arcs of width exactly `w=1/(7k_6)`, separated by positive safe gaps. Let `J` be a largest
component of `E`. If every point of `J` were `k_6`-dangerous, connectedness would put `J` inside
one open danger arc. A closed arc contained in an open arc of length `w` has length strictly less
than `w`; this contradicts `|J|=L>=w`. Thus some `t* in J` is safe for `k_6`, and by construction
it is also safe for `P,k_1,...,k_5`. ∎

Equivalently, `k_6 >= 1/(7L)` suffices. This threshold is optimal from component length alone:
for every `L<1/(7k)` a closed interval of length `L` can lie strictly inside one open
`k`-danger arc. At equality it cannot; the danger-arc endpoints are safe. The earlier assertion
that equality could be covered had the open/closed topology backwards.

## What the factor `7/3` does—and does not—rescale

THM-1102 used

`T_old = min(A,C)`, where `A=N/(6mu)` and `C=1/(3L)`.

With `k_5` the largest removed killer, distinguish

`R_old=T_old/k_5`, `R_comp=C/k_5`, and `R_sharp=(1/(7L))/k_5`.

Then the exact identity is

> `R_sharp = (3/7) R_comp`.

It is **not** generally `(3/7)R_old`, because the counting branch `A` can be the minimum.
Only on configurations where `T_old=C` do the two ratios coincide. On THM-1102's displayed
window row, core `[1,2,4,7,9,11,12]` and killers `(158,160,162,164,166)`, the component branch
is active, so that row's `1.85794` rescales to about `0.796`. This does not rescale the maximum
over a different search domain and is not a uniform r=6 theorem.

## Verification (this session)

The floating-point search `r6_sharp_horn_search_deathstar_S58.py` directly evaluated
`R_sharp`, rather than deriving it by rescaling `T_old`. Its reported sampled maxima are:

| config family | max R_sharp |
|---|---|
| (a) consecutive small killers, **all 792 cores**, offsets 0–23, steps 1–5 | 0.8011 |
| (b) consecutive, moderate/large scale (300–3000) | 0.8011 |
| (c) **clustered large killers** (≥2 killers ≥333, tight) — *codex's flagged worry* | 0.8011 |
| (d) random spread, widths up to 1500 | 0.8011 |

The displayed sampled row is recomputed exactly by `r6_sharp_horn_deathstar_S58.py`:
core `[1,2,4,7,9,11,12]`, removed killers `[171,173,175,177,179]`,
`L=72/72275`, threshold `1/(7L)=72275/504`, and
`R_sharp=10325/12888=0.80113...`. Exact arithmetic certifies this row; it does not prove that
the row is the global maximizer. Its old component threshold is `1/(3L)=72275/216=334.60...`.

**Robustness** (`r6_sharp_horn_robust_deathstar_S58.py`): the sharp horn certifies the deep
well `{1..12,182}` (`L=1/168 > 1/(7·182)`), its tower `{1..12,364}`, and the worst `|core|=1`
body `{1..11,13,84}` — all genuinely lonely. It correctly does **not** fire on the small-`k`
tight families `{1..13}`, `{1..11,13,24}` (which are non-covering, outside the far-killer
stratum). These are regression examples; soundness comes from the topological proof above.

## What is closed, and what remains

- **Proved here:** the local component-only certification, including its equality case.
- **Proved elsewhere:** THM-1121's exact weighted atlas closes the old bounded branch
  `92<=k_i<333`. THM-1135 gives a different uniform finite mixed-scale reduction; its
  remaining box is not closed. THM-1145 refutes the separate scalar max-`T`
  extrapolation in THM-1102, without refuting this component method.
- **Proved for the former worst shape:** THM-1144 gives a search-free proof on the old
  worst core, and THM-1134 strengthens it to all 792 cores and every legal scale for the
  consecutive step-two family. THM-1134 also supplies a multiplier-chart cone and an
  exact separated-ratio gate.
- **Not proved:** `R_sharp<=1` for every core and every five-killer configuration. The
  global `0.8011` number is a sampled maximum, not a 25% uniform margin. The honest
  residual consists of non-step-two shapes outside THM-1134's cone and ratio gate.

## Toward the uniform bound: the equal-step G(σ) functional

For the special equal-step family `{b,b+d,...,b+4d}`, the phases at time `t` are
`bt+m sigma`, where `sigma=dt`. Thus the five phase offsets form an arithmetic progression.
The previously studied step-two family is the case `d=2`. This reduction does **not** apply
to an arbitrary five-killer configuration. In the AP slice,
"all 5 killers safe at `t`" is "`φ = bt` avoids 5 danger-arcs of width `1/7` centered at
`{0,-σ,-2σ,-3σ,-4σ}`". Reflection gives the equivalent positive-centre set
`{0,σ,2σ,3σ,4σ}`. As `t` ranges over a core-safe arc, `φ = bt` sweeps `b·(width) ≫ 1`
turns only in a favorable scale regime. Define the exact frozen-offset functional

  **`G(σ) = ` largest gap in `[0,1) ∖ ⋃_{m=0}^{4}(mσ−1/14, mσ+1/14)`.**

Let `H(sigma)` be the largest circular gap between the five centres. Removing open arcs of
radius `1/14` trims `1/14` from both ends of a largest centre gap. Since five centres force
`H>=1/5>1/7`, a complementary gap remains and

> `G(sigma)=H(sigma)-1/7`.

Taking `sigma` modulo one and reflecting lets us put `u=min(sigma,1-sigma)` in `[0,1/2]`.
Sorting the five centres gives

```
H = max(u,1-4u)                 for 0 <= u <= 1/4,
H = u                           for 1/4 <= u <= 1/3,
H = max(3u-1,1-2u)             for 1/3 <= u <= 1/2.
```

Indeed, the respective cyclic orders have gap multisets generated by
`{u,1-4u}`, `{u,4u-1,1-3u}`, and `{3u-1,1-2u}`; the displayed interval inequalities identify
their maxima. Solving `G>1/7`, equivalently `H>2/7`, proves the exact strict-good set

```
[0,5/28) union (2/7,5/14) union (3/7,4/7)
  union (9/14,5/7) union (23/28,1].
```

Its circle measure is exactly `9/14`. Equality `G=1/7` holds at the eight finite band
endpoints, and the closed-component horn accepts equality; an eventual landing/drift argument
may still need interior margin. The global minimum is `2/35` at
`sigma in {1/5,2/5,3/5,4/5}`. Also `G(1/3)=4/21` and `G(1/2)=5/14`.

**Worst-shape telemetry.** For the displayed exact row, `L*b=12312/72275>1/7`.
`r6_worst_shape_finite_check_deathstar_S58.py` is a floating-point scan on the single core
`[1,2,4,7,9,11,12]`: step two uses `157<=b<=4000`, while steps one, three, and four use
`157<=b<=2000`. It reports ratios below one on precisely those finite ranges. It is neither
an exact finite certificate nor an all-core or all-offset theorem.

The tempting relation `L*b approximately G(dt)` freezes `sigma` while `t` moves and was
only a heuristic before fixed landing rectangles were supplied. THM-1144 resolves that
drift on the old worst core: an exact finite head handles `157<=b<=399`, while the
rational witness interval `(71/154,13/28)` handles `b>=400`. THM-1134 strengthens this
to all cores: an exact 792-core rectangle atlas closes `b>=164`, and an independent
12,771-row endpoint bank closes every legal `b<=164`. Thus the step-two family no longer
has a landing/drift gap.

Arbitrary non-AP killer offsets remain outside this one-variable model. For them the
multiplier chart is the more flexible object: THM-1134 closes its cone
`B>=17 max(A,80)` and its exact separated-ratio `Q5` gate, but the complementary shapes
still form part of the uniform `r=6` frontier.

## Methodological note (for the fleet)

The local horn really does improve `1/(3L)` to `1/(7L)`, but a local improvement is not a
uniform tail theorem. Keep three quantifiers separate: component length alone, a sampled bank,
and all cores/all killer tuples. For the AP functional, the useful vertices are the five centres
or the five cyclic gaps. A tournament orientation discards the metric gap lengths that determine
`G`, so no clean tournament quotient is being claimed here.
