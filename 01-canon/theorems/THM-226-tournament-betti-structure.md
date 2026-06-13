---
id: THM-226
name: Tournament Betti Structure Theorem
status: PARTIALLY-PROVED
proved_parts: [β₂=0, β₁≤1]
conjectured_parts: [β₁·β₃=0, β_p≤1, χ∈{0,1}]
verified_computationally: n=6 exhaustive (32768), n=7,8 sampled
proved_by: opus-2026-03-15-S72d
related: [THM-103, THM-108, THM-100, THM-095]
---

# THM-226: Tournament Betti Structure Theorem

## Statement

For any tournament T on n ≥ 3 vertices in GLMY path homology:

1. **β₀(T) = 1** (connected) — TRIVIAL
2. **β₁(T) ≤ 1** — PROVED (THM-103)
3. **β₂(T) = 0** — PROVED (THM-108)
4. **β₁(T) · β₃(T) = 0** (mutual exclusivity) — CONJECTURED
5. **β_p(T) ∈ {0, 1} for all p** — CONJECTURED
6. **χ(T) ∈ {0, 1}** — CONJECTURED

## Consequences

The path homology Betti vector of a tournament has one of these forms:
- **(1, 0, 0, ...)** — contractible (point-like)
- **(1, 1, 0, ...)** — one 1-hole (circle S¹-like)
- **(1, 0, 0, 1, 0, ...)** — one 3-hole (sphere S³-like)
- **(1, 1, 0, 0, 1, 0, ...)** — both 1-hole and 4-hole (χ = 1)
- Other patterns with each β_p ∈ {0,1} and χ ∈ {0,1}

The mutual exclusivity of β₁ and β₃ means: a tournament cannot simultaneously
have a "free directed cycle" (β₁) and a "directed 3-dimensional cavity" (β₃).
These represent fundamentally incompatible topological structures.

## Computational Evidence

### n=5 EXHAUSTIVE (1024 tournaments)

| Betti vector | Count | χ |
|-------------|-------|---|
| (1,0,0,0,0) | 720 | 1 |
| (1,1,0,0,0) | 304 | 0 |

### n=6 EXHAUSTIVE (32768 tournaments)

| Betti vector | Count | χ | Fraction |
|-------------|-------|---|----------|
| (1,0,0,0) | 27648 | 1 | 84.4% |
| (1,0,0,1) | 320 | 0 | 1.0% |
| (1,1,0,0) | 4800 | 0 | 14.6% |

**NOTE**: β₃ first appears at n=6 (corrects previous claim of n=7).
320 tournaments at n=6 have β₃ = 1 — ALL have β₁ = 0.

### n=7 SAMPLED (3000 tournaments)

| Betti vector | Count | χ |
|-------------|-------|---|
| (1,0,0,0,0,0) | 2658 | 1 |
| (1,0,0,1,0,0) | 222 | 0 |
| (1,1,0,0,0,0) | 119 | 0 |
| (1,1,0,0,1,0) | 1 | 1 |

**KEY DISCOVERY**: One tournament has BOTH β₁ = 1 AND β₄ = 1:
- bits = 1251547, score = (2,2,2,3,4,4,4), c₃ = 11, SC = true
- Euler characteristic = 1 (consistent with χ ∈ {0,1})
- This shows strict Betti concentration (≤1 nonzero β_p) is FALSE
- But β₁ · β₃ = 0 still holds (0 violations)

## Key Observations

### Euler Characteristic Pattern
χ(T) = 1 - β₁ + 0 - β₃ + β₄ - β₅ + ...

All observed values: χ ∈ {0, 1}.
- χ = 1: contractible tournaments, or β₁ = β₄ = 1 (cancellation)
- χ = 0: tournaments with exactly one odd-dimensional hole

### Mutual Exclusivity Mechanism — Seesaw (THM-095)

**Algebraic identity** (from β₂ = 0):
  β₁ + β₃ = ker(d₁) + dim(Ω₃) - dim(Ω₂) - im(d₄)

The im(d₂) terms cancel! This sum depends only on Ω dimensions and im(d₄).
**Computationally verified**: this sum is always 0 or 1 (n=5,6 exhaustive).

Since β₁ ∈ {0,1} (THM-103), β₁ + β₃ ≤ 1 implies β₁·β₃ = 0.

**CORRECTION**: The claim "β₃=1 requires ALL 3-cycles dominated" is FALSE.
At n=6, 240 of 320 β₃=1 tournaments have 2 free 3-cycles (score [2,2,2,3,3,3]).
However, β₁ = 0 for these tournaments — free cycles being present does NOT
force β₁ = 1. The mechanism is subtler than domination structure.

### LES Induction Approach (partial)
For β₁(T)=1: pick v not in a specific free cycle σ (exists for n ≥ 7).
σ remains free in T\v → β₁(T\v) = 1 → β₃(T\v) = 0 by induction.
LES gives β₃(T) = dim H₃(T,T\v). Proof complete IF H₃(T,T\v) = 0.
**Status**: H₃(T,T\v) = 0 NOT YET PROVED.

### β₃ Generator Structure (n=7)
H₃ generators have 45-67 terms in A₃, use ALL n vertices,
and have irrational coefficients. They are truly GLOBAL objects —
not localizable to any subset of vertices.

## Proof Status

### Part 4 (β₁ · β₃ = 0): Proof approaches

**Approach A — Seesaw (THM-095, strongest):**
1. β₂ = 0 gives conservation: im(d₂) + im(d₃) = dim(Ω₂)  [PROVED]
2. ker(d₁) = C(n,2) - n + 1  [PROVED, from β₀ = 1]
3. β₁ ≤ 1  [PROVED, THM-103]
4. Algebraic identity: β₁+β₃ = ker(d₁)+Ω₃-Ω₂-im(d₄)  [PROVED]
5. This sum ≤ 1  [VERIFIED n=5,6 exhaustive; NOT YET PROVED for general n]
Status: REDUCES TO STEP 5. The bound im(d₄) ≥ ker(d₁)+Ω₃-Ω₂-1.

**Approach B — LES induction (for n ≥ 7):**
β₁(T)=1 → free cycle σ → pick v∉σ → σ free in T\v → β₁(T\v)=1
→ β₃(T\v)=0 (induction) → β₃(T) = dim H₃(T,T\v).
Need: H₃(T,T\v) = 0. Status: NOT YET PROVED.
Note: fails at n=6 because all vertices are in free cycles.

**Approach C — Cup product/Massey products:**
If β₁ = 1, the cocycle ζ ∈ H¹(T) has ζ² = 0 (since H² = 0).
The Massey triple product ⟨ζ,ζ,ζ⟩ ∈ H³ might relate to β₃.
Status: NOT DEVELOPED.

### Part 5 (β_p ≤ 1): Likely requires Part 4 + higher analog of THM-103.

### Part 6 (χ ∈ {0,1}): Would follow if the full Betti structure is constrained.

## Relation to Existing Work

- THM-103 alternative proof: all directed 3-cycles are homologous in H₁.
  This is the template for proving β_p ≤ 1 at higher dimensions.
- THM-108: β₂ = 0 proof via LES + good vertex existence.
  The same framework could prove β₁ · β₃ = 0 if relative H₃ vanishes.
- The spectral identity C^TC = nP (THM-224) operates at the edge-triple
  level and directly controls β₁. β₃ involves higher-dimensional structure
  not captured by this identity.

## Chain Complex Data (n=6 exhaustive)

### β₁=1 (4800 tournaments)
im(d₂) = 9 always (= ker(d₁) - 1). β₁ = 10 - 9 = 1.
ker(∂₃) = im(∂₄) EXACTLY for all 4800. β₃ = 0.
Free 3-cycles: always {4, 6, or 8} free cycles.

### β₃=1 (320 tournaments)
im(d₂) = 10 always (= ker(d₁)). β₁ = 10 - 10 = 0.
ker(∂₃) = im(∂₄) + 1 EXACTLY. β₃ = 1.
Two score types: [1,1,1,4,4,4] (80, 0 free) and [2,2,2,3,3,3] (240, 2 free).

### Euler characteristic
χ ∈ {0, 1} for all 32768 tournaments.
χ=1: 27648 (contractible). χ=0: 5120 (β₁=1 or β₃=1).

## Files

- `04-computation/betti_concentration_v2.py` — main verification script
- `04-computation/beta3_generator_anatomy.py` — H₃ generator analysis
- `04-computation/beta13_mechanism.py` — chain complex dimension analysis
- `04-computation/seesaw_identity.py` — seesaw identity verification
- `04-computation/domination_fix.py` — corrected domination check
- `05-knowledge/results/betti_concentration_v2_S72d.out`
- `05-knowledge/results/beta3_anatomy_S72d.out`
- `05-knowledge/results/beta13_mechanism_S72d.out`
- `05-knowledge/results/seesaw_identity_S72d.out`
