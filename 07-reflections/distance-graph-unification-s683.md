---
source: claude-2026-06-06-S683
status: synthesis + demonstrated shared spectral bound; the LRC / Hadwiger-Nelson / unit-distance unification
tags: [LRC, hadwiger-nelson, unit-distance, distance-graph, cayley-graph, fourier, delsarte, hoffman, bessel, eisenstein, de-grey, moser-spindle, chromatic, unification]
---

# LRC, Hadwiger–Nelson, unit distance: one object, transferable keys

Prompt: unify these three; see them as the same underlying thing; let insights
of one be keys to the other. They are the **chromatic / independence theory of
distance Cayley graphs**, and the shared key is the **Fourier transform of the
connection set**.

## The one object

A **distance Cayley graph** `Cay(A, S)`: an abelian group `A` (or `ℝ^d`) with a
**connection set** `S` ("the forbidden distances"); vertices adjacent iff their
difference lies in `S`. Two parameters:

- **independence ratio** `m(S)` = max density of an `S`-avoiding set;
- **chromatic number** `χ(Cay) ≥ 1/m(S)`.

**Shared spectral bound (Hoffman/Delsarte):**
```
χ ≥ 1 − λ_max/λ_min,    λ = Fourier transform of the connection set S.
```
The **negativity** of the connection-set Fourier transform is what forces `χ`
up. This single inequality is the spine of all three problems.

## The three problems are instances (demonstrated)

| problem | group `A` | connection `S` | `λ` = Fourier(S) | bound |
|---|---|---|---|---|
| **Hadwiger–Nelson** | `ℝ²` | unit circle | `J₀(2π|ξ|)` (Bessel) | `λ_min≈−0.4028 ⇒ χ_m ≥ 3.48` |
| **LRC** | `ℝ/ℤ` or `ℤ_m` | speeds / arcs | cosine-sum (arc → sinc) | AP ⇒ complete graph, `χ = n` (tight) |
| **unit distance** | `ℤ[ζ₆]`-lattice | unit vectors | (finite) | `u(n)` = max edges (HYP-2170, `n=22`) |

- **Hadwiger–Nelson:** the plane's unit-distance operator has "eigenvalues"
  `J₀(2π|ξ|)`; the spectral bound gives `χ_m ≥ 3.48`. de Grey's `χ ≥ 5` is the
  *same method* sharpened by a finite rigid gadget (below).
- **LRC** *(Barajas–Serra, n=7)*: LRC for a speed set `D` is equivalent to the
  **(circular/fractional) chromatic number** of the distance graph `G(ℤ, D)`
  being large enough. Computed: the AP `{1,…,n−1}` mod `(2n−1)`-type gives the
  **complete** distance graph (`χ = n`) — the extremal/tight case (the AP is the
  unique tight config, S679); non-AP speed sets give sparser graphs with smaller
  `χ`.
- **Unit distance:** the unit-distance graph is a finite Cayley graph on an
  Eisenstein lattice (HYP-2170 identified `n=22` as `Cay(ℤ[ζ₆], U₆)`).

## How the insights transfer (the keys)

**LRC → Hadwiger–Nelson.** My LRC covering-depth Poisson formula
`p₀ = Σ_{c∈L(V)} ∏ κ(c_i)` (HYP-2154, `κ` = the arc's Fourier coefficient) **is**
the Delsarte/LP method that bounds the plane's measurable independence ratio —
with `J₀` in place of `κ`. The **resonance lattice `L(V)`** is exactly the dual
lattice the LP sum runs over. So the machinery I built for LRC (Poisson sum over
the relation lattice) is, verbatim, the linear-programming bound for `χ(ℝ²)`.
The transferable lesson: *the connection-set Fourier transform, summed over the
dual/relation lattice, bounds the independence ratio.*

**Hadwiger–Nelson → LRC.** The Hadwiger–Nelson lower bounds come from **finite
rigid gadgets** — the Moser spindle (`χ ≥ 4`), de Grey's 1581-vertex graph
(`χ ≥ 5`): a finite unit-distance subgraph that *cannot* be properly colored with
few colors. This is precisely LRC's **tight-configuration** method: the AP and
`V*` are finite extremal configurations that *pin* the loneliness bound. So
"build a finite rigid certificate" transfers both ways. And — the key geometric
import — **de Grey's gadget lives on the Eisenstein / triangular lattice
`ℤ[ζ₆]`**, the *same* prime-3 / `π/3` geometry as LRC's resonance shells
(HYP-2170) and the `Cl₂(π/3)` tropical constant (HYP-707). So the *lattice* on
which the hard configurations live is shared.

**The shared 7.** The Hadwiger–Nelson **upper** bound `χ ≤ 7` is the hexagonal
7-coloring — which is the Eisenstein lattice modulo a **norm-7 prime** (`7 ≡ 1
(mod 3)` splits in `ℤ[ζ₃]`, giving an index-7 sublattice that 7-colors the
triangular tiling). This is plausibly the *same* prime-3 root as the forbidden
tournament H-value `7 = |Fano PG(2,2)|` and the metagraph's perfectness break at
`n=7`. (Flagged as suggestive, not proven — the kind of `7`-coincidence worth one
careful arithmetic check via splitting in `ℤ[ζ₃]`.)

## Why they are "the same"

All three ask: **for a distance Cayley graph, how large is the chromatic number
/ how small the independence ratio** — and the answer is controlled by where the
**connection-set Fourier transform goes negative**. LRC lives on `ℝ/ℤ` (1-D,
arc → sinc), Hadwiger–Nelson on `ℝ²` (circle → Bessel), unit distance on a finite
Eisenstein lattice. The *group* changes; the *method* (Delsarte/LP on the dual
lattice) and the *extremal-gadget* technique (tight configs / Moser–de Grey) are
identical. The repo's LRC tooling (resonance lattice, Poisson formula, tight
configs, the prime-3 shells) is a ready-made toolkit for the plane, and vice
versa.

## Concrete payoffs / next

1. **Run the LRC Poisson/Delsarte LP for the plane:** transcribe HYP-2154's
   relation-lattice sum to `ℝ²` with `J₀` and reproduce the known `χ_m` LP bound
   — a literal code reuse, confirming the unification operationally.
2. **An LRC "Moser spindle":** find the *smallest* rigid LRC sub-configuration
   that forces `M < 1/n` would-be (the tight gadget) — the analogue of the
   spindle — and see whether de Grey-style amplification gives sharper LRC bounds.
3. **Pin the prime-3 `7`:** check whether `7` splitting in `ℤ[ζ₃]` (hexagonal
   7-coloring) and the forbidden tournament-`7` share a root, or are independent.
4. **de Grey lattice → LRC@19 apex:** the sieve-apex residual (S679) lives at the
   negation fixed point; test whether placing it on the Eisenstein lattice (the
   de Grey geometry) gives the off-grid witness the apex needs.
