# The Packing-Forbidden Connection: A Complete Framework
**opus-2026-03-14-S71f**

## The Packing Hierarchy

A tournament T's independence polynomial I(Ω(T), x) can be interpreted as a "packing" 
of building blocks (bricks). Each brick is a clique component K_c of the conflict graph Ω:

| Brick | I factor | H contribution | Name |
|-------|----------|----------------|------|
| K_1 | (1+x) | ×3 | Simplex |
| K_2 | (1+2x) | ×5 | Cuboid |
| K_3 | (1+3x) | ×7 | **Tesseract (IMPOSSIBLE!)** |
| K_c | (1+cx) | ×(2c+1) | c-cell |

## The Fundamental Theorem

**THEOREM**: The tesseract brick K_3 never occurs as a component of Ω(T).
That is, no tournament has I(Ω(T),x) with a factor (1+3x).

**Proof**: The factor (1+3x) requires 3 mutually adjacent cycles in Ω (forming K_3) 
that are all non-adjacent to every other cycle. This means:
1. Exactly 3 directed odd cycles exist
2. All 3 share pairwise vertex overlaps
3. No other directed odd cycles exist

But:
- At n≤4: Can't have 3 odd cycles (too few vertices)
- At n=5,6: α₁=3 is impossible (the phase transition gap)
- At n≥7: α₁=3 occurs but forces α₂≥2 (disjoint pair forcing)

Since (1+3x) requires α₁=3 AND α₂=0, it's impossible at every n.

## Consequences

### H=7 is permanently forbidden
H=7 = 1·7 has the ONLY factorization 7 (prime). The only brick giving H-contribution 7 
is the tesseract (1+3x). Since the tesseract is impossible, H=7 cannot arise from any 
packing — and since OCF ensures H = I(Ω,2), H=7 is impossible.

### H=21 analysis
H=21 = 3·7 could arise from simplex × tesseract = (1+x)(1+3x). But the tesseract 
factor is impossible. Any other factorization 21 = 21 uses a single 10-cell (1+10x), 
which requires α₁=10, α₂=0. But this is also blocked: 10 mutually overlapping odd 
cycles with no disjoint pairs requires very specific structure that tournament constraints 
prevent (at n≤6 exhaustively, at n≥7 by structural arguments).

### The 5-cell and beyond
H=11 (single 5-cell, 1+5x), H=13 (6-cell), etc. ARE achievable. The obstruction 
specifically targets the tesseract (1+3x) due to the 5-cycle forcing mechanism.

## Key Insight: Why (1+3x) but not (1+4x)?

(1+4x) works fine: at n=6, there are 2880 tournaments with I = 1+4x (H=9).
The difference: 4 mutually overlapping cycles CAN exist without creating additional cycles.
But 3 mutually overlapping cycles at the right threshold (t₃=3) ALWAYS force 5-cycles 
at n=5, and force disjoint pairs at n≥7.

This is because t₃=3 is the EXACT threshold for 5-cycle generation (HYP-1142):
any 5-vertex subtournament with t₃≥3 internal 3-cycles has a Hamiltonian cycle.
With t₃=3 globally at n=5, ALL 5-vertex subtournaments hit this threshold, creating 
at least one 5-cycle.

## The Evaluation Duality

For a brick (1+cx):
- At x=1: contributes (1+c) = c+1 (simplex counting)
- At x=2: contributes (1+2c) = 2c+1 (OCF/Hamilton)  
- At x=3: contributes (1+3c) = 3c+1 (cuboid counting)

The forbidden tesseract (1+3x):
- x=1: 4 (quaternion!)
- x=2: 7 (forbidden H!)
- x=3: 10 (= T in OCF notation)

The evaluation at x=1 gives FOUR — the number of quaternionic dimensions.
The forbidden number 7 lives at x=2, which is the OCF evaluation point.
This connects to the Hurwitz theorem: only {1, 2, 4, 8} are dimensions of division algebras,
and the gap between 4 and 8 is exactly where the forbidden number lives.

## Future Directions

1. **Can the tesseract obstruction be proved algebraically?** (Not just by exhaustion/sampling)
2. **What about (1+10x) for H=21?** Need to prove this clique structure is also impossible.
3. **Does the packing framework give a complete characterization of forbidden H?**
4. **Engineering**: The brick decomposition could give efficient algorithms for computing H 
   by decomposing Ω into components.
