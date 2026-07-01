# Certifying cohomology: the lonely-set Euler characteristic is the hard-side certificate (a Lefschetz count, not an SOS witness); and the SC-spine spectrum goes generic (metallic only at n≤5)

*opus-2026-07-01-S25. The owner asked to attempt the testable spine-eigenvalue thread and to think chain complex
/ certifying cohomology. Thread 1 is an honest negative; the chain-complex side is the real payoff — it names
the hard-side certificate. The S18 black/blue decomposition IS a chain complex, the S23 Lefschetz trace IS its
Euler characteristic, and the loneliness certificate is a cohomological invariant `χ(lonely set)`.*

## Thread 1 (spine eigenvalue fields) — honest negative
Are the SC-spine (blue) eigenvalues quadratic at larger n, matching the covering continued fraction? **No.**
Factoring the blue characteristic polynomial exactly (sympy):
- **n=4:** `λ²(λ²−λ−1)` — the **golden ratio**, field `ℚ(√5)`.
- **n=5:** rationals `±1` + `(λ²−…)` — the **silver** family, field `ℚ(√2)`.
- **n=6:** a **degree-11 irreducible** factor (× `λ` × isolated zeros). **Not quadratic** — the spine spectrum is
  a single degree-11 Galois orbit.
So the "metallic" reading (S24) is a **small-n (n≤5) coincidence**; from n=6 the spine spectrum goes **generic**
(maximal Galois degree, like the flip-rank going generic), and its fields (`√5, √2`) do **not** match the
covering CF (`t*=[0;n−1,n]=n/Φ₆`, which is rational). The testable resonance is refuted — the spine and the
covering modulus do not share an arithmetic field beyond n=5.

## The chain complex — the S18 decomposition, named
The flip-line merged metagraph `M_n` is a Z₂ **1-complex** (0-cells = merged nodes, 1-cells = flip-lines), with
boundary `∂₁`. The S18 blue/black split is exactly its Hodge-type decomposition:
- **BLACK = a 1-cycle.** The black subgraph is an even graph (every vertex even black-degree, S18), so
  `∂(black)=0` over Z₂ — black `∈ Z₁` (the cycle space), representing `H₁` classes. Betti: `H₁(black)` = the black
  cycle rank (8 at n=5, 287 at n=6).
- **BLUE = a T-join with `∂(blue) = [SC]`.** The blue subgraph has odd degree exactly on the SC nodes (the
  T-join, S18/HYP-3810), so its boundary is the SC 0-chain. Betti: `H₁(blue)` = the blue cycle rank (0,1,15 for
  n=4,5,6).
- The complement fold **σ = the chain map** whose fixed 1-cells are the blue (grid-symmetric) lines; `Fix(σ)` =
  the half-tiling.
So "black is Eulerian, blue is a T-join whose boundary is the self-complementary set" is literally the statement
`M_n = Z₁ ⊕ (T-join with ∂=SC)` — a chain complex, not a metaphor.

## Certifying cohomology — the Lefschetz count *is* the certificate
The Lefschetz number `L(f)=Σ(−1)ⁱ Tr(f_*|Hᵢ)` is an **Euler-characteristic-weighted trace** — a cohomological
object. S23's traces are exactly these: `L(φ_v)=Tr(H₀)−Tr(H₁)=1−v` for the dilation (H_0=ℚ, H_1=ℚ), the Frobenius
`H¹`-trace `√p` (Gauss sum) for Paley. The payoff: **the LRC bound has a cohomological certificate.**

The lonely set `L = S¹ \ ⋃_v D_v` (danger zones `D_v={t:‖vt‖<r}`) has Euler characteristic `χ(L) = #(components)`
= the number of lonely open arcs. **`χ(L) > 0` certifies loneliness** (`M(S) ≥ r`). Verified: for the AP at
sub-tight `r=0.99/n`, `χ(L) = 4, 6, 6` for n=5,7,14 — **positive, loneliness certified** — and it collapses to a
measure-zero boundary set exactly at the tight `r=1/n`. Crucially, `χ(L)` is computed from the **resonance
arrangement** — the `D_v` are the fixed-point sets `Fix(φ_v)` (`v` arcs each, the `L(φ_v)=1−v` data), and
`χ(L)` is their inclusion–exclusion (a Möbius/Lefschetz sum; the boundary points total `2Σv = 2·C(n,2)`, and the
gaps are the three-distance intervals). So:

> **The hard-side certificate is `χ(lonely set)` — a Lefschetz/Euler count from the dilation chain complex — not
> an SOS positivity witness.**

This completes the S22b/S23 arc. On the **easy** (Brouwer/p≡1) side the certificate is **SOS positivity** (a sum
of squares at the symmetric fixed point). On the **hard** (Borsuk–Ulam/p≡3 mod 4/free-ℤ₂) side, where SOS fails,
the certificate is **cohomological**: the Lefschetz trace (Gauss/Ramanujan sums, `1−v`) and the Euler
characteristic `χ(L)`. "Certifying cohomology" is precisely the name for the hard-side witness — the topological
count that survives where the positivity witness doesn't. The three pillars (POCS/flat-extension/Blaschke) are
its constructive face; the flat-extension atoms are the 0-cells the boundary `∂` lands on, the Blaschke fixed
points are `Fix(φ_v)`, and `χ(L)` is what they bound.

## Status
- **Thread 1 (negative, verified):** blue/SC-spine spectrum is quadratic/metallic only at n≤5 (golden √5,
  silver √2); n=6 is a degree-11 irreducible (generic); fields do **not** match the covering CF. Corrects the
  S24 "metallic ≈ covering-CF" resonance.
- **Chain complex (verified):** `M_n` is a Z₂ 1-complex; black `∈ Z₁` (∂=0, even graph), blue is a T-join
  (∂=SC); Betti/cycle-ranks computed (blue 0,1,15; black 8,287 for n=5,6).
- **Certifying cohomology (verified):** `χ(lonely set)=#components > 0` certifies loneliness (AP n=5,7,14 at
  sub-tight r); it is the inclusion–exclusion/Lefschetz count over the dilation fixed-point arrangement — the
  cohomological hard-side certificate that replaces SOS positivity.
- **Honest:** `χ(L)>0` at every sub-tight r certifies `M≥1/n` (standard); the new content is the *cohomological
  framing* (Lefschetz/Euler, not positivity) and its tie to S23's traces — the exact `χ(L) ↔ M(S)` sharp bound
  is still the three-distance/far-element analysis.

Related: HYP-3815 (Lefschetz traces = the 1st moment; χ is the Euler count), HYP-3816 (2nd moments; this adds
the chain complex), HYP-3810 (blue = T-join = the half-tiling), HYP-3814 (Brouwer/SOS vs Borsuk-Ulam; this names
the hard-side certificate), HYP-3796/mac-mini (three pillars = the constructive face). HYP-3817 (this). Scripts:
04-computation/{sc_spine_eigenvalue_fields, lrc_lonely_set_euler_char_certificate}_opus_20260701.py.
