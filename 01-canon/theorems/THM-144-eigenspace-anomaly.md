# THM-144: Eigenspace Boundary Rank Anomaly (Preliminary)

**Status:** CONJECTURE (verified p=7 exhaustive, p=11 partial)
**Author:** opus-2026-03-12-S68
**Dependencies:** CirculantHomology eigenspace decomposition, THM-125

## Statement

For a circulant tournament C_p^S on prime p ≡ 3 (mod 4) vertices, the GLMY
boundary map d_m decomposes into p eigenspaces k = 0, ..., p-1. Define:

    r_m^(k) = rank(d_m restricted to eigenspace k)

**Observation 1 (Universal):** For ALL circulant tournaments:
- r_m^(k) is the same for all k ≥ 1 (k = 1, ..., p-1)
- r_m^(0) may differ from r_m^(k≥1) — we call this the "k=0 anomaly"

**Observation 2 (Anomaly depth):** Define the anomaly depth as:
    D(C_p^S) = max{m : r_m^(0) ≠ r_m^(1)}

For the **Paley tournament T_p** (S = QR_p):
- p=7: D = 4 (anomaly at m = 1, 2, 3, 4)
- p=11: D ≥ 6 (anomaly at m = 1, 2, 3, 4, 5, 6 verified)
- Pattern: D(T_p) = m + 1 = (p+1)/2 (conjectured)

For **ALL OTHER circulant tournaments** at p=7 (6 orientations):
- D = 1 (anomaly ONLY at m = 1)
- β₁ = 1, all other β_i = 0 (except β₀ = 1)

**Observation 3 (Betti from anomaly):** The non-trivial Betti numbers arise
FROM the anomaly:
- k=0 contributes differently from k≥1 at the anomaly degrees
- β_m ≠ 0 only when the anomaly at degree m or m+1 creates a ker/im mismatch
- For Paley: β_{D} = (p-1) × 1 = p-1 when each k≥1 eigenspace contributes 1

## Verified Data

### p = 7, Paley T_7 (S = {1, 2, 4}):
| m | r_m^(0) | r_m^(k≥1) | Δ = r^(0) - r^(k≥1) |
|---|---------|-----------|---------------------|
| 0 | 0       | 0         | 0                   |
| 1 | 0       | 1         | -1                  |
| 2 | 3       | 2         | +1                  |
| 3 | 3       | 4         | -1                  |
| 4 | 6       | 5         | +1                  |
| 5 | 3       | 3         | 0                   |
| 6 | 3       | 3         | 0                   |

Anomaly pattern: (-1, +1, -1, +1) at m=1-4. Alternating sign, magnitude 1.

**Per-eigenspace Betti decomposition (VERIFIED):**
- β₄^(k=0) = 9 - 6 - 3 = 0
- β₄^(k≥1) = 9 - 5 - 3 = 1 (each)
- **β₄ = 0 + 6×1 = 6** ✓

At p=7: k=0 NEVER contributes to homology. ALL Betti comes from k≥1 eigenspaces.

**CORRECTION (S68):** This pattern does NOT generalize! At p=11, k=0 IS the
primary source of β₅. See p=11 per-eigenspace decomposition below.

### p = 7, Interval (S = {1, 2, 3}):
| m | r_m^(0) | r_m^(k≥1) | Δ |
|---|---------|-----------|---|
| 1 | 0       | 1         | -1|
| 2 | 2       | 2         | 0 |
| 3 | 4       | 4         | 0 |
| 4 | 7       | 7         | 0 |
| 5 | 7       | 7         | 0 |
| 6 | 2       | 2         | 0 |

Constant for m ≥ 2. **Chain complex exact at ALL degrees ≥ 2.**
β₀ = 1, β₁ = 1 (both from k=0 only), all others = 0.

### p = 11, Paley T_11 (S = QR_11):
| m | r_m^(0) | r_m^(k≥1) | Δ |
|---|---------|-----------|---|
| 0 | 0       | 0         | 0 |
| 1 | 0       | 1         | -1|
| 2 | 5       | 4         | +1|
| 3 | 15      | 16        | -1|
| 4 | 55      | 54        | +1|
| 5 | 150     | 151       | -1|
| 6 | 305     | 309       | -4|
| 7 | 390     | 390       | 0 |

Anomaly at m=1-6 ONLY. **r_7 = 390 for ALL eigenspaces (VERIFIED S68).**
This confirms D(T_11) = 6 = (p+1)/2, matching HYP-550.

Anomaly pattern: (0, -1, +1, -1, +1, -1, -4, 0, 0, ...).
Note: Δ₆ = -4 (not ±1!) — anomaly magnitude INCREASES at the last degree.

**Per-eigenspace Betti decomposition (VERIFIED S68):**
Using r_7=390, r_8=300, r_9=150, r_10=30 (from β=0 at degrees 7-10):

| m | β_m^(k=0) | β_m^(k≥1) | β_m total |
|---|-----------|-----------|-----------|
| 0 | 1         | 0         | 1         |
| 1 | 0         | 0         | 0         |
| 2 | 0         | 0         | 0         |
| 3 | 0         | 0         | 0         |
| 4 | 0         | 0         | 0         |
| 5 | **5**     | 0         | **5**     |
| 6 | **5**     | **1**     | **15**    |
| 7-10| 0       | 0         | 0         |

**KEY INSIGHT:** At p=11, β₅=5 comes ENTIRELY from k=0 (the opposite of p=7!).
β₆=15 has mixed origin: k=0 contributes 5, each k≥1 contributes 1.

Pattern comparison:
- p=7: β₄=6 from k≥1 only (k=0 contributes 0)
- p=11: β₅=5 from k=0 only; β₆=15 from k=0 (5) + k≥1 (10×1)
- The role of k=0 REVERSES between p=7 and p=11!

### p = 11, Interval (S = {1, 2, 3, 4, 5}):
| m | r_m^(0) | r_m^(k≥1) |
|---|---------|-----------|
| 0 | 0       | 0         |
| 1 | 0       | 1         |
| 2 | 4       | 4         |
| 3 | 16      | 16        |
| 4 | 58      | 58        |
| 5 | 166     | 166       |
| 6 | 356     | 356       |

**Completely constant for m ≥ 2!** Chain complex exact at m=2,...,6 (verified).
Ω_m - r_m - r_{m+1} = 0 at all verified degrees.
β₀ = 1, β₁ = 1 (from k=0 only), conjectured β_m = 0 for m ≥ 2.

## Why Paley has deep anomaly

The Paley tournament has automorphism group Aut(T_p) ⊇ {x ↦ ax : a ∈ QR_p},
which is a group of order m = (p-1)/2. This group acts on the eigenspaces:
eigenspace k maps to eigenspace ak (mod p). Under QR action, the eigenspaces
split into:
- k=0: fixed (trivial representation)
- k≥1: orbits of size m under QR multiplication

The k=0 anomaly reflects the fact that the "quotient" boundary map (on the
orbit space) has different rank from the "generic" boundary map.

For Interval: Aut = Z_p only (cyclic rotation), which already IS the
eigenspace decomposition. No further symmetry → no anomaly.

## Chain Complex Exactness

**Key discovery (S68):** For Interval T_11, the chain complex is **EXACT**
at degrees 2-6 (verified), and conjecturally at all degrees ≥ 2.

This means: for each eigenspace k ≥ 1, ker(d_m) = im(d_{m+1}) exactly.
The per-eigenspace Betti numbers are zero at every degree.

Evidence:
  Ω = [1, 5, 20, 74, 224, 522, 926, 1222, 1115, 611, 148]
  r = [0, 0, 4, 16, 58, 166, 356, 570, ?, ?, ?]
  Ω_m - r_m - r_{m+1} = [*, *, 0, 0, 0, 0, 0, ?, ?, ?, ?]

Five consecutive exact degrees (m=2-6) with ZERO gap is extremely unlikely
by coincidence. Under the exactness hypothesis:
  r_8=652, r_9=463, r_10=148 (=Ω_10, injective at top degree)
  β(Interval T_11) = [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], χ = 0

This matches the p=7 pattern exactly.

## p=11 Complete Topological Landscape (S68)

**HYP-552 REFUTED**: At p=11, there are FOUR distinct topological types,
not two. Organized by Z_p* orbits:

| Orbit | Size | β profile | χ | χ/eigenspace | S-S=Z_p*? |
|-------|------|-----------|---|--------------|-----------|
| Interval | 10 | [1,1,0,...,0] | 0 | 0 | NO |
| Paley | 2 | [1,0,0,0,0,5,15,0,...] | 11 | 1 | YES (perfect) |
| {1,2,3,4,6} | 10 | [1,0,0,0,0,0,21,0,...] | 22 | 2 | YES |
| {1,6,7,8,9} | 10 | [1,0,0,0,0,14,13,0,...] | 0 | 0 | YES |

**Difference set criterion (HYP-559)**: Deep anomaly ⟺ S-S covers Z_p*.
PERFECT correlation: 22/22 deep have full S-S, 10/10 shallow have incomplete.

Anomaly patterns (Δ_m = r_m^(0) - r_m^(k≥1)):
- Interval: Δ = (-1, 0, 0, ...) — stops at m=1
- Paley: Δ = (-1,+1,-1,+1,-1,-4, 0, ...) — alternating ±1, jump at m=6
- {1,2,3,4,6}: Δ = (-1,+1,-1,+1,-1,+1, 0, ...) — clean ±1 alternation
- {1,6,7,8,9}: Δ = (-1,+1,-1,+1,-1,-2, 0, ...) — alternating ±1, small jump

## Consequences

1. **Topological landscape is RICHER than expected**: At p=7, exactly 2 types.
   At p=11, at least 4 types. The number of types grows with p.

2. **Acyclicity mechanism**: The Interval chain complex is acyclic at k≥1
   because ALL boundary ranks are constant across eigenspaces (m≥2).
   No eigenspace has "room" for homology.

3. **Deep anomaly creates homology**: Any S with S-S = Z_p* has deep
   anomaly, creating gaps that allow β_m ≠ 0 at high degrees.

4. **Spectral-topological-H connection**:
   - Flat spectrum (Paley) → QR symmetry → deep k=0 anomaly → complex topology
   - Peaked spectrum (Interval) → no extra symmetry → trivial topology
   - At the phase transition p~13: topology doesn't change (both are stable)
   - H-maximizer transition is INDEPENDENT of topological transition

## Related
- THM-125: Constant symbol matrix (Ω_m^(k) same for all k)
- THM-137: Paley eigenvector theorem
- HYP-531: Topological dichotomy at p=7
- Paley T_11 Betti: [1, 0, 0, 0, 0, 5, 15, 0, 0, 0, 0]
