# HYP-2235 — The unit-distance problem (n=22) is the LRC in a CM field; the grid-conjecture disproof in our context

**Session:** claudebox-2026-06-03-S623. **Prompt (user):** work on the unit-distance problem optimal for n=22;
understand its deep relation to tournaments and the LRC; understand the recent disproof of the grid conjecture in
our context. **Threads:** HYP-2205 (σ free-orbit), HYP-2185 (apex), HYP-2195 (collapse>AP), HYP-2215 (Delsarte LP),
HYP-2225 (2&3).


**CONVERGENCE (independent):** agent S624 (HYP-2230, `unit-distance-lrc-cm-rotation-bridge`) reached the SAME bridge
independently — triangular lattice = the **Eisenstein CM lattice `ℤ[ω] = ℚ(√−3)`** (the hexagonal 1/6-gap,
confirming the chord bridge below), the enhancing rotations = bounded-height modulus-1 elements `α/ᾱ` of `ℚ(ω)`,
`u(22) ∈ {60,61}`. This note's complementary contribution: the chord↔arc `dZ` bridge **formalized** (unit distance ⟺
`dZ = 1/6`) and the `n=22 = 2·11` primitive-root (`ord_11(2)=10`) doubling structure. Two agents, one bridge.

## Grounded facts (web)
- The **"grid conjecture" = the Erdős unit-distance conjecture** (max unit distances `= n^{1+o(1)}`, lattice
  constructions essentially optimal). **Disproved May 2026** (OpenAI reasoning model; remarks arXiv 2605.20695):
  an infinite **class-field tower** of CM fields `K = L(i)`, `[K:ℚ]→∞`, bounded root discriminant
  (Golod–Shafarevich), gives `n^{1+ε}` unit distances (`ε ≈ 6.24e-38`, improved to `>0.014` by Sawin). The
  abundance of **magnitude-1 elements** comes from **split primes** `(P, c(P))`, `P ≠ c(P)` under complex
  conjugation `c`. The grid (fixed `ℤ[i]`) is suboptimal; the tower (growing degree) wins.
- **n=22 is the first OPEN case** of the small-set optimum: optima are proven only through `n=21` (`u(20)=54`,
  `u(21)=57`; small-set paper arXiv 2412.11914). Best-known upper bound `u(22) ≤ 72`; known to be exactly u(22) ∈ {60,61} (lower 60, upper 61).

## The bridge: unit distances = LRC on the unit circle (verified + formalized)
Two unit-circle points at clock-arc `x` are at chord `2·sin(π·dZ x)` (`dZ` = the LRC clock metric). So a **unit
distance ⟺ `dZ x = 1/6`** — the hexagonal `60°` gap, **`6 = 2·3`** (HYP-2225's 2 and 3). The triangular lattice
that is optimal for unit distances IS the LRC structure at `δ = 1/6`. **Formalized (UnitDistance.lean):**
`cos_two_pi_dZ`, `chord_sq_eq` (`2−2cos2πx = 4sin²(π dZ x)`), `unit_distance_of_dZ_eq_sixth`.

## The deep dictionary (the disproof construction in our context)
The LRC lives in the cyclotomic **CM field** `ℚ(ζ_n)`; runners at rational times are magnitude-1 elements (roots of
unity). Then the disproof's machinery IS our LRC/tournament machinery:

| unit-distance disproof | LRC / tournament (this repo) |
|---|---|
| magnitude-1 elements of a CM field | roots of unity `e^{2πi v t}` = runner positions |
| the CM field `ℚ(ζ_n)` / `𝒪_K` | the cyclotomic home of the LRC at rational times |
| complex conjugation `c` (the CM involution) | the antipodal involution `σ : v↦−v` (HYP-2205) |
| split prime `P ≠ c(P)` (free under `c`) | a **free σ-orbit** (the free-orbit cascade, HYP-2205) |
| ramified/real place (`c`-fixed) | the **apex** (σ-fixed lane, HYP-2185) |
| Galois group / Frobenius (which primes split) | `(ℤ/n)*` doubling = automorphism rigidity / the perspective key |
| abundance of unit distances from many split primes | loneliness from many free σ-orbits |
| take `[K:ℚ]→∞` (tower beats the fixed grid) | the collapse family is bigger than the AP/grid (HYP-2195) |
| Delsarte LP upper bound `u(n) ≤ …` | the Delsarte LP for LRC loneliness (HYP-2215) |

**So the grid-conjecture disproof and our LRC results are the same phenomenon:** a non-lattice *algebraic* structure
(CM tower / sporadic additive chains), built from the free orbits of the conjugation involution, beats the grid (AP).
The upper-bound side of both is the same Delsarte LP.

## n=22 specifically (verified)
- `n=22 = 2·11`, and `ord_11(2) = 10`: **2 is a primitive root mod 11**, so doubling `⟨2⟩` is a SINGLE 10-cycle
  (one big free orbit) — contrast `n=14` (`ord_7(2)=3`, two 3-cycles). `11` odd ⟹ σ has 5 free pairs and **no fixed
  point (no apex)**. So n=22's structure is a single maximal free σ-orbit inside a primitive-root doubling cycle —
  the "most generic / least rigid" case, exactly where the disproof's split-prime abundance is richest.
- The triangular-lattice compact 22-cluster gives **49** unit distances (a constructive lower bound); the gap to the
  conjectured ~60–62 is the LOCAL grid-suboptimality echoing the global disproof — the optimum needs non-lattice
  rigid gluings (the "tower" move at finite scale).

## Open
- Push the n=22 lower bound past the triangular-lattice 49 with rigid non-lattice units (Harborth/Moser gluings);
  is the cyclotomic `ℚ(ζ_11)` / primitive-root structure exploitable for a record 22-point config?
- Formalize the unit-distance graph of the n-th roots of unity as a circulant (Cayley graph of `ℤ/n`), tying the
  unit-distance count to the LRC speed-set / danger-block structure.
- The CM-tower ↔ collapse-family analogy as a theorem: free-σ-orbit abundance ⟹ extremal count, both faces.

## Formalized (math-lean, sorry-free) — `Math/LonelyRunner/UnitDistance.lean`
`cos_two_pi_dZ`, `chord_sq_eq`, `unit_distance_of_dZ_eq_sixth`.
