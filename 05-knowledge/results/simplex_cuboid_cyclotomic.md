# Simplex-Cuboid-Cyclotomic Connection to Forbidden H Values

**opus-2026-03-14-S71e**

## Setup

User directive: "k-nacci approaches 2 and weighted k-nacci approaches 3. Think of simplices as (x+1)^n and cuboids as (x+2)^n, and think about packing them inside each other."

## Key Results

### 1. k-nacci and weighted k-nacci limits

- **k-nacci**: a_{n+1} = a_n + ... + a_{n-k+1}. Ratio → root of x^{k+1} - 2x^k + 1 = 0. Limit = **2**.
- **Weighted k-nacci** (weight 2): a_{n+1} = 2*(a_n + ... + a_{n-k+1}). Ratio → root of x^{k+1} - 3x^k + 2 = 0. Limit = **3**.

These are exactly (x+1)|_{x=1} = 2 and (x+2)|_{x=1} = 3.

### 2. OCF as simplex/cuboid evaluation

At n=5 (where α₂=0): I(Ω, x) = 1 + ((H-1)/2)·x

| Point | Name | Value | Limit of ratio |
|-------|------|-------|----------------|
| x=1 | Simplex | (H+1)/2 | - |
| x=2 | OCF | H | OCF/simplex → 2 (k-nacci!) |
| x=3 | Cuboid | (3H-1)/2 | cuboid/simplex → 3 (weighted k-nacci!) |

### 3. Gap formula: (x+2)² - (x+1)² = 2x+3

Setting 2x+3 = q²+q+1 = Φ₃(q):
- x = (q+2)(q-1)/2

| q | Φ₃(q) | x | (x+2)²-(x+1)² |
|---|--------|---|----------------|
| 2 | **7** | 2 | 7 |
| 3 | 13 | 5 | 13 |
| 4 | **21** | 9 | 21 |
| 5 | 31 | 14 | 31 |

**Both forbidden H values are gaps between adjacent squares at x=2 and x=9=3².**

### 4. Cyclotomic polynomial connection

- H=7: The only decompositions (α₁=3, α₂=0) or (α₁=1, α₂=1) correspond to independence polynomial I(Ω,x) = 1+3x or 1+x+x² = **Φ₃(x)**.
- H=21: The key decomposition α₁=2, α₂=4 corresponds to I(Ω,x) = 1+2x+4x² = **Φ₃(2x)**.
- Both Φ₃(x) and Φ₃(2x) are **irreducible** over Q.

### 5. Simplex-cuboid sandwich

Define the "simplex band" m as the integers in [log₄(H), log₃(H)].
- If 3^m ≤ H ≤ 4^m: H is "in band m" (achievable via m independent cycle VSs).
- Gaps between bands: 4^m < H < 3^{m+1}

| Gap | Range | Odd values | Forbidden |
|-----|-------|------------|-----------|
| 1→2 | [5, 8] | {5, 7} | **{7}** |
| 2→3 | [17, 26] | {17,19,21,23,25} | **{21}** |
| 3→4 | [65, 80] | {65,...,79} | ∅ |

The forbidden values live in the gaps, near the midpoints (6.5 and 21.5).

### 6. T = C(q+1, 2) — triangular numbers

For H = Φ₃(q): T = (H-1)/2 = (q²+q)/2 = C(q+1, 2)
- q=2: T=3 = C(3,2), 2 decompositions to block
- q=4: T=10 = C(5,2), 6 decompositions to block
- q=8: T=36 = C(9,2), 19 decompositions to block (TOO MANY → achievable)

### 7. Independence polynomial factorization at n=6

Exhaustive computation shows:

| H | (α₁, α₂) | I(Ω,x) | Factors |
|---|-----------|---------|---------|
| 9 | (2, 1) | x²+2x+1 | **(x+1)²** ← simplex! |
| 17 | (6, 1) | x²+6x+1 | irreducible |
| 29 | (10, 2) | 2x²+10x+1 | irreducible |
| 45 | (14, 4) | 4x²+14x+1 | irreducible |

**Key**: H=9 with (2,1) gives I = (x+1)² — a PERFECT simplex polynomial.
This is the "simplex tournament": 2 disjoint independent cycles, each with d=1.

### 8. Why only {7, 21} are permanently forbidden

1. **Gap location**: Both are in simplex-cuboid gaps (between 3^m and 4^m).
2. **Small gap width**: Gaps 1→2 and 2→3 have few odd values (2 and 5 respectively).
3. **Few decompositions**: T=C(q+1,2) gives only 2 or 6 cases to block.
4. **Tournament structure constraints**: HYP-1142 (internal triple forcing) and the Binary Phase Theorem block all decomposition routes.

For Φ₃(2³)=73: the gap 3→4 has 8 odd values and T=36 has 19 decompositions. Tournament structure cannot block all routes simultaneously.

## The Packing Interpretation

- **Simplex = (1+x)^m**: Tournament with m mutually disjoint cycle VSs, each with d=1.
  - H = 3^m (powers of 3)
  - This is "m simplices packed with no overlap"

- **Cuboid = (1+2x)^m**: Tournament with m disjoint cycle VSs, each with d=2.
  - H = 5^m (powers of 5)

- **Forbidden values**: Require independence polynomial I(Ω,x) = Φ₃(q·x) for some q.
  - Cyclotomic = **unpacked** = cannot decompose into independent components
  - The tournament structure FORCES decomposability of the independence complex

## Conjecture (HYP-1179)

**The permanent moat of odd H values consists of exactly {7, 21} = {Φ₃(2), Φ₃(4)}.**

No other odd value is permanently forbidden. This is because:
1. Only gaps 1→2 and 2→3 are small enough for tournament constraints to block
2. The number of decomposition routes (T//2 + 1) grows quadratically with q
3. Tournament structure constraints (HYP-1142, Binary Phase) have fixed blocking power
