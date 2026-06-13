# Terminology Evolution: From Confusion to Clarity

**opus-2026-03-24-S297**

## The Definitive Naming Convention

All connections between tilings are **waggly lines**. They decompose by Hamming distance d = 1, ..., m into layers. Every old term maps to this framework.

## Term Map

### Current (use these)

| Term | Definition | Layer |
|------|-----------|-------|
| **Waggly line** | ANY connection between two distinct tilings | All (d=1..m) |
| **Wiggly line** | Flip exactly 1 tile (d=1 waggly) | d=1 |
| **Blue/black line** | Flip all m tiles = complement tiling (d=m waggly) | d=m |
| **Blue line** | Blue/black line where the tiling is grid-symmetric | d=m, grid-sym |
| **Black line** | Blue/black line where the tiling is NOT grid-symmetric | d=m, not grid-sym |
| **Silent mutation** | Any waggly line that preserves the iso class (self-loop) | Any d |
| **Expressive mutation** | Any waggly line that changes the iso class (cross-class) | Any d |
| **Wiggly class** | Which tile was flipped (A, B, C, ...) | d=1 only |
| **Neutrality** | Fraction of a tiling's d=1 neighbors in the same class | d=1 |
| **Reach** | Number of distinct classes among d=1 neighbors | d=1 |

### Deprecated (do not use, mapped to current)

| Old term | What it meant | Current equivalent |
|----------|--------------|-------------------|
| **Translucent** (S270-S273) | A single-arc flip preserving the iso class | **Silent mutation** at d=1 |
| **Opaque** (S270-S273) | A single-arc flip changing the iso class | **Expressive mutation** at d=1 |
| **"BB is purely diagonal"** (S295) | WRONG claim from MISTAKE-033 | Blue/black has 18 cross-class edges at n=5 |
| **"BB lives outside Q_m"** (S295) | WRONG — confused complement tiling with T^op | Blue/black IS in Q_m (distance m) |
| **"BB and wiggly are disjoint"** (S275) | True but misleading — they're different layers of same spectrum | d=1 ≠ d=m, both ⊂ waggly |
| **"Blue/black = complement pairing"** (old CLAUDE.md) | Partially correct but imprecise | Blue/black = complement TILING (not T^op) |

### Legacy terms with correct meaning (keep using)

| Term | Meaning | Status |
|------|---------|--------|
| **SC/NS** (self-complementary) | Tournament iso class where T ≅ T^op | Correct, keep |
| **Spine/Ribs/Sea** | SC-SC / SC-NS / NS-NS edge types | Correct, keep |
| **Fiber** | Number of tilings in an iso class = H/\|Aut\| | Correct, keep |
| **Meta-graph** | Quotient of tiling space by isomorphism | Correct, keep |
| **Merged meta-graph** | Meta-graph with complement pairs identified | Correct, keep |

## The Key Insight That Was Nearly Lost

The old "translucent/opaque" terminology from S270-S271 contained a beautiful discovery: the **translucent subgraph** of the transitive class is the **(n-1)-permutohedron**. This is because transitive tournaments correspond to permutations, and neutral single-arc flips correspond to adjacent transpositions. This result is VALID — only the name needed updating. The correct statement:

> The d=1 silent-mutation subgraph within the transitive class (H=1) is the permutohedron Π_{n-1}.

## The Five Mistakes to Never Repeat

1. **MISTAKE-033**: Complement tiling ≠ complement tournament (T^op). The tiling complement preserves base-path arcs.
2. **Old S275 framing**: "Disjoint" was technically true (d=1 ∩ d=m = ∅) but created the false impression that blue/black and wiggly are fundamentally different kinds of thing. They are layers of the same waggly spectrum.
3. **Old S295 claim**: "Blue/black is purely diagonal in the merged meta-graph" — WRONG. It has 18 cross-class edges at n=5.
4. **Sessions S211-S231**: Used "blue/black" to mean SC-type preservation at the class level. This is DIFFERENT from the tiling-level blue/black definition. Use spine/ribs/sea for class-level types.
5. **Confusing "all tiles" with "all arcs"**: There are m = C(n-1,2) tiles and M = C(n,2) arcs. Flipping all tiles ≠ flipping all arcs. The difference is the n-1 base-path arcs.
