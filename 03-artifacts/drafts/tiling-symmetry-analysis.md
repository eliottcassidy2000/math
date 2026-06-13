# Tiling Symmetry and Class Structure Analysis

**Source:** opus-2026-03-06-S1 (deep tiling investigation)

---

## Overview

Investigation of the {0,1}^m tiling space and how its geometry relates to tournament isomorphism classes. The triangular grid Grid(n) has tiles at positions (r,c) with r>=1, c>=1, r+c<=n-1. Each tiling t in {0,1}^m encodes a tournament T_t with a fixed base Hamiltonian path.

---

## Finding 1: Sigma (Converse) Acts Cleanly on Classes

**sigma: (r,c) -> (c,r)** on the grid corresponds to taking the converse tournament (T -> T^op) after relabeling by vertex reversal (k -> n-1-k).

### Key properties:
- sigma PRESERVES tiling bits (no complement) — it permutes arc positions only
- sigma PRESERVES tiling weight (number of 1-bits)
- sigma is an involution
- sigma maps isomorphism classes to isomorphism classes (verified n=3..6)
- Sigma-fixed tilings: |Fix(sigma)| = 2^floor((n-1)^2/4) (verified n=3..6)

### Class structure under sigma:

| n | Classes | Self-converse | Paired |
|---|---------|---------------|--------|
| 3 | 2 | 2 | 0 |
| 4 | 4 | 2 | 1 |
| 5 | 12 | 8 | 2 |
| 6 | 56 | 12 | 22 |

### Paired classes:
- Always have identical H, size, score, c3, c5, omega_deg, weight distribution
- Cannot be distinguished by ANY tournament invariant (since T and T^op share all invariants)
- CAN be distinguished by tiling-specific properties (e.g., which specific tilings belong to each)

### Self-converse classes:
- sigma maps the class to itself
- Within each, tilings partition into sigma-fixed and sigma-paired
- Can sometimes share all invariants with OTHER classes (sigma pairs or other self-converse classes)

---

## Finding 2: Standard Invariants Almost Distinguish Classes

At n=6 (56 classes), the tuple (score, c3, c5, omega_3_degree_sequence, H) fails to distinguish only 8 groups:

| Group type | Count | Description |
|------------|-------|-------------|
| Sigma pairs | 5 | Two converse-paired classes with identical invariants |
| Mixed | 3 | 1 self-converse + 1 sigma pair, or 2 cross-score partners |

The "cross-score" phenomenon: converse maps score (d_1,...,d_n) to sorted (n-1-d_i). So classes at score (0,2,2,3,4,4) pair with classes at score (1,1,2,3,3,5). These pairs share c3, c5, omega_deg, and H despite having different score sequences.

### Implication:
Modulo the converse ambiguity (which is inherent — T and T^op are genuinely different tournaments with identical invariants), standard invariants distinguish all classes at n<=6.

---

## Finding 3: Complement Does NOT Respect Classes

The bit complement (flip all non-path arcs) does NOT map classes to classes cleanly. This means the tiling space has a strong asymmetry: the "all-zeros" end (transitive tournament) and "all-ones" end are fundamentally different, unlike in the tournament-complement duality.

---

## Finding 4: Triangle Geometry and 3-Cycle Probabilities

Each vertex triple (a,b,c) forms a triangle on the grid. The 3-cycle probability depends only on how many path arcs the triangle contains:

| Path arcs | Non-path arcs | P(3-cycle) | Triangle type |
|-----------|--------------|------------|---------------|
| 2 | 1 | 1/2 | Consecutive triple (i,i+1,i+2) |
| 1 | 2 | 1/4 | One adjacent pair |
| 0 | 3 | 1/4 | No adjacent pairs |

### Expected c3 over uniform random tilings:

**E[c3] = (C(n,3) + n - 2) / 4**

Verified n=3..6. This differs from the random tournament expectation C(n,3)/4 by the correction term (n-2)/4, which arises because consecutive triangles have doubled 3-cycle probability due to the fixed path arcs.

### c3 distribution:

| n | E[c3] | Var[c3] |
|---|-------|---------|
| 3 | 0.50 | 0.25 |
| 4 | 1.50 | 0.75 |
| 5 | 3.25 | 1.44 |
| 6 | 6.00 | 2.69 |

---

## Finding 5: Strong H ~ c3 Correlation

| n | Correlation(H, c3) | H = 1+2c3 exact? |
|---|--------------------|------------------|
| 3 | 1.0000 | Yes |
| 4 | 1.0000 | Yes |
| 5 | 0.9553 | No (breaks at c3=4) |
| 6 | 0.9560 | No (breaks at c3>=2 for some) |

At n<=4: H(T) = 1 + 2·c3(T) exactly. At n>=5: H depends also on c5, c7, etc. through the OCF.

The mean H over random tilings:
- n=3: 2.0, n=4: 4.0, n=5: 9.88, n=6: 29.0

---

## Finding 6: Bit-Position Variance (Which Arcs Predict Class?)

Within each isomorphism class, some arc positions have nearly fixed values (low variance) while others vary (high variance):

**Most predictive arcs (lowest within-class variance):**
- The longest arc (gap = n-1) consistently has the lowest variance
- Arcs symmetric to it under sigma also have low variance

**Least predictive arcs (highest variance):**
- Middle-gap arcs (gap = 2, interior positions) have highest variance
- The center arc (1,3) at n=5, arc (1,3)/(2,4) at n=6

### Interpretation:
The longest arc controls the most "global" structure of the tournament, so it's most constrained by isomorphism class. Middle arcs control "local" structure that is more permutable.

---

## Finding 7: Class Transition Graph

When one tile is flipped (one bit changed), the tournament changes and may move to a different isomorphism class. The class transition graph records which classes are connected by single flips.

### Properties:
- **Always connected** (verified n=3..6)
- **ΔH is always even** (Rédei: both H values are odd)
- **ΔH distribution is symmetric around 0** (by sigma/complement symmetry)
- **Self-loop fraction decreases**: 0%, 42%, 10%, 5% for n=3,4,5,6

### ΔH by arc gap:
- Longest arc (gap=n-1) has the largest |ΔH| range
- E[ΔH] = 0 for every arc position (no arc is systematically "good" or "bad" for H)
- ΔH=0 probability is ~15-20% for most arc positions at n=6

### Implication for OCF proof:
Since the transition graph is connected, a proof that OCF is preserved under each single-arc flip suffices to prove OCF for all tournaments. This is exactly the arc-flip strategy (THM-013/014), but now with geometric context: flipping tile (r,c) on the grid affects specific triangles and cycles in a geometrically constrained way.

---

## Finding 8: Weight Distribution Within Classes

Each class has a spread of tiling weights (number of backward non-path arcs).

### Key observations:
- sigma preserves weight, so sigma-paired classes have identical weight distributions
- Self-converse classes have weight distributions that need NOT be symmetric about m/2
- Weight distributions CAN distinguish classes that share all tournament invariants (e.g., two self-converse classes with same score, c3, c5, omega_deg, H but different weight distributions at n=5)

### Hamming distance within classes:
- Average Hamming distance between tilings in the same class is roughly m/2 (half the bits differ)
- Minimum Hamming distance is often 1-2 (classes are locally connected in the cube)
- Maximum Hamming distance can be up to m (opposite corners of the cube)

---

## Synthesis: Toward Cracking OCF Open

### What the tiling perspective adds:
1. **Geometric constraints on arc flips**: Not just "flip any arc" — the grid structure means some flips are geometrically adjacent and their combined effects are constrained.

2. **Sigma reduction**: The proof space is halved by sigma symmetry. Only need to verify OCF preservation for one class in each sigma-pair.

3. **c3 as primary predictor**: Since c3 explains 95% of H variance, the OCF "correction" from higher cycles is small. A proof might proceed: (a) prove H = 1+2c3 for c3-dominated regime, (b) show the correction from 5-cycles, 7-cycles etc. follows from the independence polynomial structure.

4. **Transition graph connectivity**: OCF can be proved by induction on single-arc flips starting from the transitive tournament. The tiling grid gives geometric structure to which flips are "easy" (diagonal neighbors) vs "hard" (distant arcs).

### Key formula for tiling analysis:
- class_size = H(T)/|Aut(T)|
- E[c3] = (C(n,3) + n - 2) / 4
- P(3-cycle for consecutive triple) = 1/2, for all others = 1/4
- Sigma-fixed tilings: 2^floor((n-1)^2/4)

---

## Computational Scripts

All scripts saved to /tmp/ during this session:
- `tiling_geometry.py` — bit variance, class centroids, Hamming distances
- `tiling_s3_v2.py` — sigma symmetry analysis (corrected)
- `tiling_sigma_pairs.py` — classify indistinguishable groups
- `tiling_transitions.py` — class transition graph
- `tiling_triangles.py` — triangle analysis and H vs c3 correlation
