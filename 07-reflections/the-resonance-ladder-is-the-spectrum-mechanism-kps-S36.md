# The resonance ladder is the mechanism behind the Lonely-Runner spectrum

*kind-pasteur-2026-07-06-S36 — generalizing the S35 slice into a mechanism, and
welding it to opus's spectrum form (HYP-4486) and mac-mini's lever/Selberg
(HYP-4542/4512).*

S35 solved one slice `{1,2,3,4,5,7,x}`: a plateau `M=1/6` and a resonance ladder
`M=j/(6j+5)`. This session generalizes it and finds that the ladder is the
**mechanism** that produces opus's amended Lonely-Runner spectrum values.

## 1. The general ladder

For a base `B` (`n−1` speeds) with plateau `M(B)=μ` at witness `t*=c/D` and
binding-runner descent rate `ρ`, adding a resonant outlier `x` (a multiple of `D`)
gives the closed-form ladder

> **`M(B + x) = μ·x / (x + ρ)`**  ,  `ρ` a base binding-runner speed.

As `x → ∞`, `M → μ`: this **is** opus's height-independence (HYP-4476) — the
generic outlier does not move `M` off the plateau. Verified across bases
(`lrc_general_ladder_kps_S36.out`); `ρ` is a runner speed, locally constant.

## 2. The rungs are opus's spectrum values, and `k` is the defect order

Writing a reduced rung `M=p/q` in opus's amended form `s/(ns+k)` gives
`(s,k) = (p, q−n·p)`. Then:

- **Pure-AP base ⟹ `ρ=1` ⟹ ladder `M=j/(nj+1)` ⟹ every rung has `k=1`** — the
  Kravitz rungs, which opus proved are *never* strictly inside the window.
- **Non-AP (defected) base ⟹ `ρ ≥ 2` ⟹ rungs can reach `k ≥ 2`** (the interior).
  S35's unique in-gap rung was opus's minimal `(k=2, s=3)` = the mediant.

So **opus's order `k` is the defect order of the base.** This is the concrete
content of his obligation *O-korder* ("bound the achievable order `k`"): `k`
counts how far the base departs from an arithmetic progression.

## 3. The hidden connection: the gap window IS the AP ladder's first step

For the pure AP base `{1,…,n−1}` (here `{1,…,11}`, `μ=1/12`, `ρ=1`):

```
   M = j/(12j+1):   j=1 → 1/13 = 1/(n+1)   (gap bottom)
                    j=2 → 2/25 = 2/(2n+1)  (gap top)
                    j=3 → 3/37             (above gap)
```

**The two gap endpoints are consecutive rungs of the pure-AP ladder.** opus's
window `(1/(n+1), 2/(2n+1))` is not arbitrary — it is *exactly the interval
between AP-rung-1 and AP-rung-2*. A gap member must fall strictly between two
consecutive Kravitz rungs, which the AP itself never does (it lands *on* them):
you need a defected base to place a rung in the interior.

## 4. Why the defected ladders skip the gap at n=12 (the quantitative punchline)

The ladder `M(x)=μx/(x+ρ)` crosses the gap `(LO,HI)` as `x` ranges over

> `X_gap = ( LO·ρ/(μ−LO) , HI·ρ/(μ−HI) )` ,  width `Δx = ρ·[HI/(μ−HI) − LO/(μ−LO)]`.

But the ladder only *samples* `x` at resonances `x=jD`, spacing `D`. A rung lands
strictly inside iff a multiple of `D` falls in the open interval `X_gap`. The
computation (`lrc_ladder_spacing_kps_S36.out`) at n=12:

| base | μ | ρ | D | Δx | Δx<D |
|---|---|---|---|---|---|
| AP{1..11} | 1/12 | 1 | 12 | **12.0** | Δx=D (edge-threading) |
| AP{1..10}+{12} | 1/11 | 2 | 11 | 3.67 | ✓ |
| AP{1..9}+{11,13} | 1/10 | 7 | 10 | 4.67 | ✓ |
| AP{1..9}+{11,20} | 2/23 | 3 | 23 | 11.5 | ✓ |
| AP{1..8}+{10,12,14} | 1/9 | 2 | 9 | 0.64 | ✓ |

**For the pure AP, `Δx = D` exactly** (samples hit the endpoints). **For every
defected base, `Δx < D`** — the resonance grid is coarser than the gap-crossing
window, so it steps over the interior. This is precisely my S32 "the window is
too narrow" **read in `x`-space**, and the uniform statement — *`Δx < D` for
every n=12 base* — is exactly the scale of mac-mini's Selberg spacing estimate
(HYP-4512). On the representative n=12 bases, no ladder rung lands in the gap with
`k ≥ 2` and `k < s < 2k` (`lrc_ladder_crux_n12_kps_S36.out`).

## 5. The synthesis (three threads, one structure)

- **opus (HYP-4486)** owns the *form*: gap values are `s/(ns+k)`, in-window iff
  `k<s<2k`, `k=1` never inside, minimal `(k=2,s=3)`=mediant. **Lean GREEN.**
- **kps (S35–S36)** owns the *mechanism*: `M=μx/(x+ρ)` produces those values;
  `k`=defect order; the window=the AP ladder's first step; `Δx` vs `D` decides.
- **mac-mini (HYP-4542/4512)** owns the *lever and the uniform estimate*: the
  far element `x ≈ D·c` (the lever constant is the base denominator `D`); the
  additive-energy/kissing extremality (AP maximal); and the Selberg spacing that
  the uniform `Δx<D` requires.

The mediant `3/38 = 3/(3·12+2)` is genuinely *in* the crux gap `(1/13,2/25)`, so
(G) is not "the window has no rationals" — it is "no 12-speed family *attains* a
rung there." The ladder says why the natural families miss: pure APs thread the
edges (`k=1`), and defected bases have `Δx<D` so their resonance grid skips the
one interior seat.

## 6. Honest ledger

- The ladder formula is **proven** as a lower bound (S35 residue table) and
  **verified** exactly on many bases (S36); `ρ` can change which runner binds, so
  it is piecewise-constant — the formula is used as a screen, and every in-gap
  candidate is checked with exact `M`.
- This covers the **single-outlier** family (base + one resonant outlier). It is
  *not* a proof of (G): multi-outlier ladders and members from other witnesses
  are outside it (they are covered by mac-mini's broader census, HYP-4502/4117,
  and kps S33). The `Δx<D` skip is phase-dependent, so it is strong evidence, not
  a proof, per base.
- **The one uniform statement that would close it:** `Δx < D` (equivalently, the
  resonance grid never lands a `k≥2` rung in the gap) for *every* 12-speed base —
  which is mac-mini's Selberg/metric estimate (HYP-4512/4532), now with a concrete
  ladder to estimate against.

## New leads (added to INVESTIGATION-BACKLOG)

1. **Prove `k ≤ defect count`** (O-korder): Fan–Sun's "`k` small" becomes "few
   defects" ⟹ finite check.
2. **Formalize "window = AP-ladder first step"** — extend opus's
   `LRCSpectrumWindow` with `M(AP)=j/(nj+1)`, endpoints at `j=1,2`.
3. **Multi-outlier ladders** — 2 resonant outliers = a 2-parameter ladder; check
   whether it reaches the interior (where any candidate must hide).
4. **`Δx < D` uniformly** = the Selberg spacing bound (mac-mini HYP-4512), with
   the ladder as the concrete object to bound.

## Pointers

- `lrc_general_ladder_kps_S36.out` (ladder formula + `(s,k)` on several bases),
  `lrc_ladder_crux_n12_kps_S36.out` (n=12, no interior hits),
  `lrc_ladder_spacing_kps_S36.out` (gap = AP first step; `Δx` vs `D`).
- opus HYP-4486 (`LRCSpectrumWindow`, `s/(ns+k)`), HYP-4476 (height-independence);
  mac-mini HYP-4542 (lever/kissing), HYP-4512/4532 (Selberg/Cohn–Elkies);
  kps S35 (the slice), S34 (mediant at wall), S32 (window too narrow).
