# The Two Staircases — How Tournaments Encode Their Own Recursion

**Filed by:** opus-2026-04-04-S4

## The Core Insight

A tournament on n vertices has C(n,2) = n(n-1)/2 arcs. These arcs live on the staircase Young diagram δ_{n-1}, which has exactly C(n,2) cells.

When we fix a Hamiltonian path (the base path P₀: n→n-1→...→1), we pin n-1 of these arcs, leaving C(n-1,2) = (n-1)(n-2)/2 free binary choices — the TILES. These tiles live on the SMALLER staircase δ_{n-2}.

So a tournament is simultaneously:
- **δ_{n-1}** when all arcs are free (the "full" tournament)
- **δ_{n-2}** when the base path is fixed (the "tiling" model)

The difference between these two staircases is the HYPOTENUSE — the n-1 base-path arcs that sit on the anti-diagonal of δ_{n-1}.

**Fixing the Hamiltonian path = removing the hypotenuse of the staircase.**

## The n → n+2 Recursion

The staircase δ_n (for an (n+2)-tournament with base path) decomposes as:

```
δ_n = δ_{n-2} ⊔ bottom_leg ⊔ top_leg ⊔ apex

     ·
    /|            apex = tile (n+2, 1), connecting endpoints
   / |
  /  | δ_{n-2}   overlap = n-tournament's tiles
 /   |
·----·
↑    ↑
top  bottom      legs = wiring of endpoint vertices
```

Verified computationally (THM-291):
  **H_{n+2}(t_inner, 0_boundary) = H_n(t_inner)**

Setting boundary tiles to 0 makes the endpoints source/sink, and every Hamiltonian path must traverse: source → [inner HP] → sink. The HP count reduces exactly to the inner tournament's count.

## The Two Legs Are Mirror Images

The bottom wiring tiles (x, 1) have skips 2, 3, ..., n.
The top wiring tiles (n+2, y) have skips n, ..., 3, 2 (reversed).

By grid reflection symmetry, the polynomials on these two legs are IDENTICAL:
H(0_inner, t_bottom, 0_top, 0_apex) has the same values as
H(0_inner, 0_bottom, t_top, 0_apex) under the natural skip-preserving map.

The bottom leg controls where vertex 1 (the sink) fits in the hierarchy.
The top leg controls where vertex n+2 (the source) fits.
Grid reflection swaps source↔sink, bottom↔top.

## The Degree Hierarchy

| Component | Tiles | Max degree contribution |
|-----------|-------|------------------------|
| Inner (δ_{n-2}) | C(n-1,2) | 2⌊(n-1)/2⌋ (THM-260) |
| Bottom leg | n-2 | 2 (from 3-cycles through vertex 1) |
| Top leg | n-2 | 2 (from 3-cycles through vertex n+2) |
| Apex | 1 | 1 (linear only for pure boundary) |
| Mixed inner-boundary | — | up to 2⌊(n+1)/2⌋ (from Hamiltonian cycles) |

The PURE BOUNDARY polynomial has degree exactly 4 at every n tested (n=5,6,7).
This is 2 × 2 = (contribution from bottom 3-cycle) + (contribution from top 3-cycle).

The DEGREE GROWTH from n to n+2 comes entirely from MIXED interactions:
Hamiltonian cycles of the (n+2)-tournament that traverse both boundary and interior.

## Why n → n+2 (Not n → n+1)?

The n → n+2 step adds TWO vertices (one at each end of the base path). This is the natural recursion because:

1. **Geometric**: Both legs of the right triangle are added simultaneously, preserving the triangle's symmetry.

2. **Algebraic**: Adding one vertex would break the source/sink structure. With two vertices (one source, one sink), the inner tournament is completely enclosed.

3. **Parity**: The degree cap 2⌊(n-1)/2⌋ changes by 2 when n→n+2 for odd n, and by 0 for n→n+1 (even to odd). The n→n+2 step always increases capacity.

4. **The Cayley-Dickson connection**: The project's "Mode B" recursion (documented in CLAUDE.md) maps n to n-2, removing both legs. This is the REVERSE of the n→n+2 construction. In the Cayley-Dickson tower (R→C→H→O→S), each step doubles the algebra dimension — analogous to the staircase growing by two legs.

## The Boundary Polynomial as a New Object

The pure boundary polynomial at n+2 is:
  P_n(t_bottom, t_top, t_apex) = H_{n+2}(0_inner, t_boundary)

This polynomial lives on 2(n-2)+1 = 2n-3 variables and has degree 4.
It captures how the two endpoint vertices interact with a TRANSITIVE inner tournament.

The structure:
- Linear terms: 2^(skip-1) for each boundary tile (THM-284)
- Same-leg quadratics: -2^max(1,|Δs|-1) (antiferromagnetic within each leg)
- Cross-leg quadratics: complex (the bottom-top interaction, controlled by which 
  inner vertices are "shared" between the bottom and top 3-cycles)
- Degree-3 and degree-4 terms: involve at most 2 tiles per leg

The boundary polynomial is the "wiring problem" — given a fixed inner tournament,
how do the endpoint connections affect the HP count?

## Recursive Construction of the Full Polynomial

H_3(t₁) = 1 + 2t₁     [1 tile, 1 boundary]
H_5(t_{inner}, t_{bnd}) = H_3(t_inner) + Δ₃     [1+5 = 6 tiles]
H_7(t_{inner}, t_{bnd}) = H_5(t_inner) + Δ₅     [6+9 = 15 tiles]
H_9(t_{inner}, t_{bnd}) = H_7(t_inner) + Δ₇     [15+13 = 28 tiles]

At each step:
- The inner polynomial is PRESERVED EXACTLY
- A boundary correction Δ is added
- Δ involves 2n-1 new boundary tiles and their interactions with inner tiles
- ~97% of coefficients at each level are boundary-touching
- The degree increases by 2 (for odd n)

This recursive construction builds the multilinear polynomial H(t) layer by layer,
like adding successive shells to a growing crystal. The inner structure is never
disturbed — it's encapsulated inside the boundary wiring.

## Connection to the Independence Polynomial

From the OCF perspective: adding boundary tiles creates new odd cycles.
New 3-cycles through vertex 1: (1, a, b) where a,b are inner vertices.
New 3-cycles through vertex n+2: (n+2, a, b).
New 5-cycles through both: (1, ..., n+2, ...).
New Hamiltonian cycles: the highest-degree terms.

The independence polynomial I(Ω_{n+2}, 2) = I(Ω_n, 2) × (1 + corrections from new cycles).
The corrections are NOT simply multiplicative (new cycles conflict with inner ones),
but the leading-order structure is: I(Ω_n) gets multiplied by wiring factors.
