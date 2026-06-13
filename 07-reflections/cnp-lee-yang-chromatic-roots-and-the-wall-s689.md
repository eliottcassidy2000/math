---
source: claude-2026-06-06-S689
status: NO elimination of {5,6,7} (open) — but the Lee-Yang lens the user asked for: VERIFIED Lee-Yang circle of the HN gadget + the outward march + an honest analytic-wall meta-statement
tags: [hadwiger-nelson, chromatic-plane, lee-yang, fisher-zeros, chromatic-roots, potts, deletion-contraction, wheel, moser-spindle, zero-free-region, barvinok, sokal, integrality-gap, honest-negative]
---

# χ(ℝ²) through the Lee-Yang lens: the HN gadget's zeros are a circle, and where the analytic methods die

**Prompt (user):** a long session on the chromatic number of the plane; attempt any
novel statement, even small, like eliminating one of {5,6,7}; *"still think
Lee-Yang."*

**Honest headline: I did not eliminate any of 5, 6, 7** (each is a major open
result). What the Lee-Yang lens *did* give: a clean reformulation (with an honest
note on what's tautological in it), a **verified** Lee-Yang circle for the basic HN
gadget, a verified "outward march" of the zero cloud as χ rises, and a meta-statement
on *why the analytic Lee-Yang method stalls at the same wall* as the measure method
(HYP-2278). This is the chromatic-roots companion to HYP-2278 (which took the
measure / Loeschian / field-tower angles).

## The object: chromatic roots ARE Lee-Yang zeros

The chromatic polynomial is the zero-temperature antiferromagnetic Potts partition
function,
```
P(G,q) = Σ_{S⊆E} (−1)^{|S|} q^{c(S)} = Z_Potts(G; q, v = −1),
```
so its complex roots are the chromatic / **Lee-Yang–Fisher zeros** of that model.
Computing `P` by this Whitney-rank sum is the same deletion–contraction the repo
already lives in (tie-induction = DC = Tutte/Potts, S682/HYP-2263). The chromatic
number is where this partition function's real-axis behaviour turns on.

## The reformulation — and what is tautological in it

By de Bruijn–Erdős, `χ(ℝ²) = sup_G χ(G)` over finite unit-distance graphs. In
Lee-Yang language one is tempted to say "`χ(ℝ²) = 1 + (largest integer that is a
chromatic root of a UD graph)`." **This is true but tautological**: for *every*
graph, `P(G,k)=0` for all integers `0 ≤ k < χ(G)` (you simply can't `k`-colour),
so the integers `0,…,χ−1` are automatically roots. The integer layer just restates
dBE.

The **substantive** Lee-Yang content is the *complex* zero-free region:
> **Eliminate 7  ⟺  `q = 6` lies in a zero-free region of the entire UD-graph family**
> (`P(G,6) > 0` for all UD `G`) — a Barvinok/Sokal-type analytic statement.
> **Eliminate 5  ⟺  some UD graph's zero locus pinches the real axis at `q = 5`.**

So the lens turns CNP into: *does the chromatic-zero locus of the unit-distance
family reach the real integers 5 and 6?*

## Verified: the Lee-Yang circle of the HN gadget

The minimal Eisenstein gadget is the wheel `W₆` = unit hexagon + center (`χ=3`).
Wheels have
```
P(W_n,q) = q·P(C_n,q−1) = q[(q−2)^n + (−1)^n(q−2)] = q(q−2)[(q−2)^{n−1} + (−1)^{n+1}],
```
so the "rim" roots solve `(q−2)^{n−1} = ±1`, i.e. **`|q−2| = 1`**. The chromatic /
Lee-Yang zeros of every wheel lie **exactly on the circle of radius 1 centred at
`q = 2`** — an antiferromagnetic Lee-Yang circle. Verified for `W₆`: the 5 rim roots
sit on `|q−2|=1` to `1e-6` (`q = 1`, and `2 + e^{±iπ/5}, 2 + e^{±3iπ/5}`). Because
the locus is this small circle, every real `q ≥ 3` is *outside* it — which is
exactly why the HN gadget 3-colours, and the real-axis crossing sits at `q = 2`.

**The zeros are shifted roots of unity.** `(q−2)^{n−1} = ±1` means `q − 2` is a root
of unity, so `W₆`'s Lee-Yang zeros are `q = 2 + ζ` with `ζ⁵ = −1`. The **same
roots-of-unity locus** the FTA keystone `z⁷ − z = z(z⁶−1)` put in the *geometric*
plane (HYP-2276 — the hexagon as 6th roots of unity) reappears here, shifted by
`+2`, as the gadget's Lee-Yang zeros in the *q*-plane. Cyclotomic zeros on both sides
of the dictionary: the gadget's geometry (`z⁶−1`) and its colouring partition
function (`(q−2)⁵+1`). The triangular patch, by contrast, has its zeros fill *inside*
the disk — circularity is special to the hub-symmetric wheel.

## Verified: the outward march

As gadgets force higher `χ`, the zero cloud expands and its real-axis crossing
moves right (exact chromatic polynomials, Whitney-rank):

| gadget | χ | real-axis crossing | zero-cloud radius `max|q−2|` |
|---|---|---|---|
| `K₃` (triangle) | 3 | 2 | 1 |
| `W₆` (hexagon+center) | 3 | 2 | 1 (on the circle) |
| triangular patch `R=1` (9 pts) | 3 | 2 | 1 (cloud fills *inside*) |
| **Moser spindle** | **4** | **3** | **√2** |

The χ=3 family keeps its zeros on/inside `|q−2|=1`; forcing `χ=4` (the Moser
spindle, which needs the non-cyclotomic `√−11` rotation — HYP-2276/2277) pushes the
zeros out to `|q−2| = √2` and the real crossing to `q = 3`. So
> **`χ(ℝ²) − 1` = the rightmost point where the UD family's Lee-Yang zero locus
> crosses the real axis.** de Grey's `χ=5` graph is the statement that this crossing
> reaches `q = 4`; `{5,6,7}` is the question of whether it reaches `5` and `6`.

## The wall: why the analytic Lee-Yang method stalls (the honest negative)

Sokal's theorem confines chromatic roots to `|q| ≤ 7.96·Δ` for max degree `Δ`, and
Barvinok-style interpolation gives zero-free disks around integers for bounded-degree
graphs. **These cannot cap `χ(ℝ²)` at 6**, because unit-distance graphs have
**unbounded degree** — a single point can have arbitrarily many unit-neighbours
(any number of points on its unit circle) — so no uniform `Δ`, no uniform zero-free
region. This is *precisely* the wall HYP-2278 hit from the other side: the
fractional/measure bound gives `χ_f ≤ 4.36 < 5`, so it cannot even *reach* 5.

> **Both analytic frontiers stall at the same place.** The measure/spectral bound
> pushes *up from below* and dies at `4.36` (can't reach 5); the Lee-Yang/Barvinok
> zero-free-region bound pushes *down from above* and dies at unbounded degree
> (can't cap at 6). The `{5,6,7}` gap is exactly the integrality gap *between* the
> two analytic frontiers — the LRC Vitali wall (THM-406 M2). Narrowing it is
> irreducibly combinatorial: you must exhibit a finite gadget (a zero at `q=5`) or
> a finite colouring (zero-freeness at `q=6`), not an inequality.

## What would be a genuine next step (concrete, Lee-Yang-flavoured)

1. **A bounded-degree reduction.** Is `χ(ℝ²) = sup{ χ(G) : G unit-distance, Δ(G) ≤ D }`
   for some fixed `D`? Moser (`Δ=4`) and de Grey (bounded) suggest yes. If `χ(ℝ²)`
   is attained on degree-≤`D` UD graphs, Sokal/Barvinok *do* apply — and the question
   becomes whether the UD-specific zero-free region around `q=6` survives. (Sokal's
   `7.96D` is far too weak, but UD geometry might force a much smaller region.)
2. **The UD-specific root-free interval.** Jackson: no chromatic roots in `(1,32/27]`
   for any graph. Is there a *larger* real root-free interval just left of an integer
   for unit-distance graphs (a geometric Sokal bound using the Szemerédi–Trotter
   `O(n^{4/3})` unit-distance edge bound)? A root-free `(6−ε,6]` for all UD graphs
   would eliminate 7.
3. **The zero-cloud radius law.** The march `1 → √2` (χ: 3→4) — is the cloud radius a
   monotone function of `χ` with a clean form? If `max|q−2|` for a `χ=k` UD graph were
   bounded by `r(k)` with `r(k)+2 < k+1` impossible, that would bound `χ` — the
   quantitative version of the outward march.

## Where it sits

Companion to **HYP-2278** (measure/Loeschian/field-tower, the wall from below) and
the **FTA/Lee-Yang** thread (HYP-2275, the LRC face; HYP-2276/2277, the HN
field-tower — the `√−11` that makes the Moser spindle's zeros march out). The honest
sum: the Lee-Yang lens reframes CNP cleanly, *verifies* the HN gadget's circular
zero-locus and the outward march, and pins *why* analytic methods can't finish — but
the elimination itself remains open and combinatorial.
