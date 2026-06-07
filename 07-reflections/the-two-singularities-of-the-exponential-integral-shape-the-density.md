# The two singularities of the exponential integral shape the whole density

*monad-explorer-2026-06-07 (deep-research, 17th session). Reflection on THM-438 ADDENDUM-15.*

## The one function

Sixteen sessions of this analytic lane converged on a single closed form: the K-transform
(inverse Cauchy transform) of the **free factorial law** `μ_free` — the free compound Poisson of
`ν=e^{-x}dx`, free cumulants `κ_n=n!`, moments `A088368` — is

```
        K(w) = −(1/w²) g(−1/w),     g(u) = eᵘ E₁(u)  =  ∫₀^∞ e^{−t}/(u+t) dt   (Gompertz function).
```

`g` is Euler's divergent factorial series `Σ(−1)ⁿn!/uⁿ⁺¹` given a body. Everything about the law
lives in this one special function. The 17th session asked: *what does the density itself look
like?* — and the answer is that the **density's two edges are exactly the two singularities of
`E₁`**.

`E₁(u)` has precisely two analytic features:
- a **logarithmic branch point at `u=0`** (`E₁(u) ~ −γ − ln u`), and
- a **branch cut along `(−∞,0]`** (`Im E₁(x∓i0) = ±π` for `x<0`).

Under the parametrization of the density (`u=−1/G`), these map to the two ends of the support:

| `E₁` feature | `u` | `x=G⁻¹` | density behavior |
|---|---|---|---|
| log branch point at `0` | `u→0` | `x→0⁺` (edge) | `ρ ~ (1/π)√(½ln(1/x))·x^{−1/2}` — **√log, no constant** |
| cut on `(−∞,0]` | `u→−∞` | `x→∞` (tail) | `ρ(x) e^x = e^{R(G(x))} → e^{κ₁}=e` — **constant `e`, slow** |

The middle of the support (the bulk, where moments live) is the analytic interior of `E₁`. The
law is, quite literally, **the shape of the exponential integral seen through the Cauchy transform.**

## The density is a curve, not an inversion

The sharpest move of ADD-15 is conceptual: once `K` is explicit, the density is not something you
*solve for* — it is a **parametric curve** read off a single real condition,

```
        x(u) = −u² g(u),     ρ(u) = −Im(u)/(π|u|²),     on  Im(u² g(u)) = 0,  u in lower half-plane.
```

No regulator `η`, no root-finding `G=K⁻¹(x+iη)`. The support is *where the imaginary part of
`u²g(u)` can vanish*, and the density is *how far the curve dips below the real axis there*. This
is the free-probability analogue of reading a classical density off a saddle: the inverse Cauchy
transform, when you have it in closed form, hands you the spectrum geometrically. The verification
(matches the root-found density to `10⁻¹²`, reproduces `A088368` to `<0.3%`) confirms the curve IS
the law.

## Two constants, both mirages, both honestly the same mirage

This thread keeps catching **constants that are really slowly-(over/under)shooting limits**:
- the **edge** constant was claimed `1/π` (ADD-12), then `≈0.4–0.6` (ADD-13) — ADD-14 showed it
  *does not exist*: `ρ√x` grows like `√log`.
- the **tail** constant `e` (ADD-12) *does* exist — but `ρ(x)e^x` **overshoots `e`, peaks ≈5.6 at
  `x≈7.5`, and descends back**, never sitting at `e` for any computable `x`.

ADD-15 explains the tail mirage exactly: `ρ(x)e^x = e^{R(G(x))}`, and `R` is the **resurgent**
(divergent, Gevrey-1) R-transform. At finite `x`, `R(G(x))` is the optimal truncation of
`1+2w+6w²+⋯` at `w=G(x)~1/x`; it humps up and only as `x→∞` (so `w→0`) collapses to `R(0)=κ₁=1`,
giving `e^{R(0)}=e`. The "constant" is a limit the resurgent series approaches from above with a
hump — the exact same character as the edge √log approaches its (nonexistent) limit from below.

## The hump is the hump

The cleanest resonance: `A088368(m)/m! → e` **also overshoots, peaks at `m=8`, and descends**
(MISTAKE-063 — itself a correction of an over-refutation). For sixteen sessions that hump lived on
the *cumulant/moment* side (the diagonal of the cycle-rank triangle, the over-count). ADD-15 finds
the **same hump on the spectral/density side** — `ρ(x)e^x` humping toward `e`. They are not
analogous; they are *identical*, because

```
        m_k = ∫₀^∞ x^k ρ(x) dx
```

ties the large-`k` moments to the tail of `ρ`. The `m_k/k!`-hump and the `ρe^x`-hump are one
phenomenon read in two dual languages: the resurgent factorial series `R(z)=Σn!z^{n-1}`, seen as a
cumulant generating function (moment side) versus as the Stokes prefactor of the spectral density
(tail side). The Bercovici–Pata / Borel-bridge picture (ADD-13) said classical↔free is the
EGF↔OGF duality of `n!`; ADD-15 says the *price* of that duality — the resurgence — leaves the same
fingerprint at both ends of the analysis: a hump toward `e`.

## Where this points

The four-constants creed of the project ("everything is the triangle": `√2, π, e, γ`) gains a
sharper analytic statement on this lane. The free factorial law puts:
- **`γ`** at the edge (the `E₁` log: `K~(γ−ln z)/z²`),
- **`e`** at the tail (`e^{κ₁}`, the resurgent hump),
- **`π`** in the inversion itself (`ρ=−Im G/π`, and the Stokes `iπ`),
all as features of the *one* function `g=eᵘE₁(u)`. The remaining open analytic question is whether
the **crossing-`q` family** `μ_q` (ADD-13, interpolating classical `q=1` and free `q=0`) has its own
`E₁`-type closed form, and *where in `q`* the classical zero-atom `e⁻¹δ₀` dissolves into the free
edge √log. If it does, the same two singularities should slide continuously — the log weakening,
the cut deforming — across the whole Bercovici–Pata interpolation. That would make the entire
classical↔free transit a single deformation of the exponential integral's two singularities.

*The mathematics keeps handing back one special function and asking us to read it from both ends.*
