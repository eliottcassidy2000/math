# HYP-2310 — Degree = #roots: the n+1 ↔ n ↔ n−1 off-by-one map, the master correspondence of the arc

**Session:** claudebox-2026-06-04-S636. **Prompt (user):** the mapping between n−1 and n things; a degree-n polynomial
has n+1 coefficients (with the constant) mapped to n roots (with multiplicity). **Threads:** HYP-2305 (AG_n/cube root),
HYP-2270 (Vieta sum/product), HYP-2245 (H/indep poly), HYP-2299 (Tutte).

## The master off-by-one map
A degree-`n` polynomial has **`n+1` coefficients** (including the constant — the base / "+1") and, over `ℂ`, exactly
**`n` roots with multiplicity** (FTA). For a monic polynomial the `n` non-leading coefficients are the `n` elementary
symmetric functions of the roots — **Vieta is a bijection** between the `n` monic coefficients and the `n` roots (the
S629 sum/product is `n=2`). The **constant term** (the `n+1`-th coefficient) `= (−1)^n·(product of roots)` is the
ground state / vacuum / base; the `n` roots are the excitations.

## The off-by-one across the arc
| structure | `n` things | the `±1` partner |
|---|---|---|
| polynomial | `n` roots | `n+1` coefficients (the constant = base); derivative → `n−1` critical points (Gauss–Lucas) |
| Hamiltonian path | `n` vertices | `n−1` edges (the path / spanning tree) — `H` counts these |
| LRC | `n` runners | gap `1/(n+1)` (the stationary runner / origin = the +1) |
| perspective key | `n` | `persp(n) = #structures(n+1)` (small n) |
| independence/partition poly (`H`) | degree `n` | `n+1` counts `α₀,…,α_n`, `α₀ = 1` (empty set = base); roots = resonance/Lee–Yang spectrum |
| cyclotomic `Φ_n` | `φ(n)` roots = primitive n-th roots | `Φ₃ = x²+x+1`, 2 roots `ω,ω²` (the cube root) |
| Euler char / Tutte | `V − E + F` alternating | deletion–contraction (S633) |

**The deep principle:** the `+1` (constant term / vacuum / base / apex / stationary runner / `α₀=1` / identity) is the
ground state; the `n` roots/objects are the excitations. The **partition function** (independence / `H` / chromatic /
Tutte) IS this polynomial; the **coefficients↔roots (Vieta) duality is the counts↔spectrum duality** the whole arc
runs on.

## The cube-root instance (formalized + verified)
`X³ − 1` (degree 3) has exactly the **3 roots `1, ω, ω² = e^{±2πi/3}`** — the eigenvalues of a 3-cycle, the AG_n
generators (S635), the `π/3`/`Φ₃` resonance. Sum `= 0` (the `x²` coefficient — the `n` roots sum to the trivial,
Vieta), product `= 1` (the constant — the base). `x³−1 = (x−1)·Φ₃`. **Formalized** `Math/Polynomial/RootCount.lean`:
`card_nthRoots_complex` (`X^n−1` has `n` roots over `ℂ` — FTA, degree = #solutions), `card_cube_roots` (`= 3`).

## Synthesis
The degree-`n` / `n`-roots / `n+1`-coefficients correspondence is the off-by-one map underlying everything: the `+1`
is the base/vacuum (the constant term, the empty set, the lonely origin, the apex, the identity), the `n` is the
excitation count (roots, runners, cycles, non-identity elements), and the polynomial that ties coefficients to roots
is the partition function. FTA = "excitations match the degree"; Vieta = "their symmetric functions are the
coefficients"; the cube-root case is the resonance.

## Formalized (math-lean, sorry-free) — `Math/Polynomial/RootCount.lean`
`card_nthRoots_complex`, `card_cube_roots`.

## Open
- The Vieta bijection `n` monic coefficients ↔ `n` roots (symmetric functions) formalized in general.
- The independence/partition polynomial's constant term `= 1` (the base) and degree `= #objects`, as the off-by-one
  for `H`.
- Gauss–Lucas (the `n → n−1` derivative/critical-point map) and its meaning for the resonance spectrum.
