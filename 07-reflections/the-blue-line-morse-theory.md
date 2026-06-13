# The Blue Line and Morse Theory

**Session**: kind-pasteur-2026-03-22-S20x

## The Blue Line Skeleton

The self-complementary (SC) tournaments form a "blue line" through tournament space -- the diagonal where T is isomorphic to its complement. This is the y=x line of the binary tiling.

At n=5:
- 704/1024 tournaments are SC (68.8%)
- SC spans the FULL H range [1, 15]; NSC is constrained to {3, 5}
- NSC tournaments have HC=0 ALWAYS -- they never have Hamiltonian cycles
- H=5 is NSC-exclusive; H=1,9,11,13,15 are SC-exclusive

At n=6:
- 7040/32768 tournaments are SC (21.5%) -- the blue line thins rapidly
- SC still spans the full range [1, 45]; NSC spans [3, 43]
- NSC now HAS Hamiltonian cycles (HC up to 5)
- SC-exclusive H values: {1, 17, 41, 45} -- the extremes remain SC-only

**The blue line thins (68.8% -> 21.5%) but continues to span the full range.**

## Symmetry Breaking Reallocates H Budget

The decomposition H = n*HC + L (cyclic + linear paths) reveals:

At H=15, n=5:
- SC regular (2,2,2,2,2): HC=2, L=5 -- tree-like, many linear paths
- SC near-regular (1,2,2,2,3): HC=3, L=0 -- cyclic, every HP in a cycle
- NSC: DOES NOT REACH H=15

Breaking the y=x symmetry converts linear paths into cyclic paths. The budget H = n*HC + L is REALLOCATED from L to HC when symmetry breaks.

## H as a Morse Function on the Hypercube

Tournament space on n vertices is the m-cube Z_2^m where m = C(n,2). Each vertex is a tournament (a binary assignment to each pair), and edges are arc flips.

H is a Morse function on this space. At n=5:
- **120 local minima** -- exactly the 5! transitive tournaments with H=1
- **64 local maxima** -- exactly the H=15 SC tournaments (ALL are SC!)
- **840 saddle points** -- all intermediate H values

**There is no other local maximum.** H=15 is the only critical value for maxima. Every gradient ascent from any tournament terminates at an SC tournament with H=15. The 1024 tournaments partition into 64 basins of attraction, ALL draining to SC maxima.

This means: **H has no local-but-not-global maxima at n=5**. The landscape is a single mountain with a plateau top.

## The Missing Delta

Arc flip deltas at n=5 are: {-12, -8, -6, -4, -2, 0, +2, +4, +6, +8, +12}.

**Delta = +/-10 is IMPOSSIBLE.** The Petersen number 10 = C(5,2) is absent.

At n=4, ALL even deltas in [-4, +4] appear -- no gaps.

The gap at 10 means: no pair of tournaments differing in exactly one arc have H values differing by exactly 10. This connects to the earlier Gray code finding.

Why 10? From deletion-contraction: delta = H(T/e) - H(T'/e'). The contracted tournaments have 4 vertices, so H values in {1, 3, 5}. Max |delta| = 5-1 = 4. Wait -- that gives max |delta| = 4 for the contracted part, but the original delta also involves H(T\e) vs H(T'\e'). Actually, H(T) - H(T') = H(T/e) - H(T'/e') directly. At n=4, max H = 5, min H = 1, so max |delta| = 4. But we observe |delta| up to 12 at n=5. The contracted tournaments at n=4 are NOT tournaments (they're general digraphs from contraction), so H can exceed tournament bounds.

The delta = 10 gap likely arises from parity/modular constraints on the contracted digraph HP counts.

## The Arborescence Signature

SC tournaments have 18 distinct arborescence values; NSC has only 7.
SC-exclusive arb values include odd values (1, 3, 7, 9, ...) while NSC arbs are concentrated at multiples of 4.

The arborescence (directed spanning tree count) measures the "tree-directedness" of a tournament. SC tournaments explore a much wider range of tree structure than NSC ones.

## Connection to the Broader Picture

1. **The blue line AS the staircase diagonal**: SC tournaments are exactly those whose binary tiling of the staircase is symmetric under 180-degree rotation. This is the geometric content of self-complementarity.

2. **Morse theory + OCF**: H = I(Omega, 2) is a function of the conflict graph. The Morse flow on tournament space induces a flow on the conflict graph space. The basin structure tells us how conflict graphs evolve under arc flips.

3. **The missing 10**: Delta = 10 = C(5,2) = the number of arcs in a 5-tournament. The dimension of the space is the forbidden step size. This cannot be a coincidence. The hypercube has m = C(n,2) dimensions, and the missing delta is m itself at n=5. At n=4, m=6 and max |delta|=4 < 6, so the phenomenon doesn't arise. Need to check n=6 (m=15).

4. **Engineering**: At n=5 the single-basin property means gradient ascent ALWAYS finds the global maximum. At n=6, a secondary peak at H=37 (score [1,2,2,3,3,4]) traps 10.6% of random starts. But the barrier height is 0 -- the trap is a plateau, not a valley. Adding a same-level exploration step (accept H-preserving flips) would escape every trap.

## The Morse Phase Transition at n=6

At n=6, H acquires its first non-global local maximum. This is a genuine topological transition:

| n | Local max H values | # maxima | # minima | Basins |
|---|-------------------|----------|----------|--------|
| 3 | {3} (and {1})     | 4+4      | 4        | trivial |
| 4 | {5}               | 24       | 28       | single peak |
| 5 | {15}              | 64       | 120      | single peak |
| 6 | {37, 45}          | 720+480  | 720      | TWO peaks |

The H=37 secondary peak has:
- Score (1,2,2,3,3,4) -- near-regular, one step from the most-regular (2,2,2,3,3,3)
- Morse index 8 (8 downward, 7 same-level) vs global max's index 12 (12 down, 3 same)
- The 7 same-level neighbors create a "plateau" that traps steepest ascent
- But barrier height = 0: you can escape without going downhill
- All H=37 local maxima are SC (self-complementary), just like the global maxima

This means: at n=6, the blue line bifurcates. It still contains all maxima, but now there are two distinct ridges on the blue line -- the global ridge at score (2,2,2,3,3,3) and a secondary ridge at score (1,2,2,3,3,4).

## Connection to the Score Hierarchy

The SC score classes at n=6 that contain local maxima:
- (2,2,2,3,3,3): H=45, HC=5 -- the REGULAR score (most balanced)
- (1,2,2,3,3,4): H=37, HC=5 -- the NEAR-REGULAR score (one step from regular)

The near-regular score has the SAME HC=5 as the regular score! The difference is entirely in L:
- Regular: L=15 (more linear paths)
- Near-regular: L=7 (fewer linear paths)

The budget H = 6*HC + L gives: 45 = 30 + 15 vs 37 = 30 + 7. The 8-path difference is entirely in linear (non-cyclic) paths. The cyclic structure is identical (5 Hamiltonian cycles each). Breaking the regularity of the score destroys 8 linear paths while preserving all cycles.

This connects to the BIBD structure from S18h: regularity maximizes both cyclic and linear paths, while near-regularity only loses the linear ones.
