# User-supplied "peripheral scaffolding" for the Collatz blueprint: preserved source

**Status: REFUTED as mathematics in most of its claims; source preserved for provenance.**
Do not import any claim below as established. The audit with exact verdicts is
[collatz_mod6_20260917_scaffolding_audit](../results/collatz_mod6_20260917_scaffolding_audit.md);
the cross-lane synthesis is
[collatz_mod6_20260917_synthesis](../results/collatz_mod6_20260917_synthesis.md).
The companion blueprint and guard attachments are preserved in
[COLLATZ-BLUEPRINT-2026-09-21-SOURCE](COLLATZ-BLUEPRINT-2026-09-21-SOURCE.md) and
[COLLATZ-GUARDS-2026-09-21-SOURCE](COLLATZ-GUARDS-2026-09-21-SOURCE.md).

Pasted on 2026-09-21. The session lead condensed the LaTeX and layout to plain
text; the mathematical content is unchanged.

---

1. The General Variational Recurrence and Geometric Shear. In our foundational 2D
framework, we tracked the boundary friction when expanding an accumulated space of
dimensions A x B to (A+1) x (B+1). The Core Local Identity:
T((A+1)(B+1)-1) - T(AB-1) = T(A) + T(B) + Shear(A,B), where the discrete
Cross-Dimensional Shear is Shear(A,B) = A(B^2-1) + B(A^2-1). Extension to Arbitrary
Grid Navigation: Horizontal Flow (Delta_A): S(A+1,B) = S(A,B) + (A + B^2 + 2AB + B);
Vertical Flow (Delta_B): S(A,B+1) = S(A,B) + (B + A^2 + 2AB + A). Physical Insight:
the shear term A(B^2-1) measures the volume swept out when an orthogonal line segment
of length A is rotated and scaled through an expanding 2D lattice; a conservation law
for discrete area, proving that multiplication is an integrated accumulation of line
components.

2. Hyper-Dimensional Polytopic Packing Laws (3D and 4D). The 3D Tetrahedral Law:
Te(A+B+1) = Te(A) + Te(B) + [ (A+1)T(B) + (B+1)T(A) + (A+1)(B+1) ] (two 3D triangular
prisms locking against a central 2D rectangular plate). The 4D Pentatope Law:
Pt(A+B+1) = Pt(A) + Pt(B) + [ (A+1)Te(B) + (B+1)Te(A) + T(A+1)T(B+1) ] (the central
locking core is the multiplication of two independent 2D triangular surfaces).

3. The Multiplicative Feedback Field (the g-Operator): AgB = AB * (AgB). Unwinding
this recursion from a baseline unit state reveals AgB = (AB)^(AB). Maximal Asymmetry
(A=1): reduces to the Factorial Engine 1gB = B!. Maximal Symmetry (A=B):
AgA = (A!)^2. Shifting from asymmetry to symmetry squares the information-carrying
capacity, mapping the 2D Pell and Square-Triangular families.

4. Lossless Tournament Data Compression (T_{n-2}). By Redei's theorem every
tournament of size N contains a directed Hamiltonian path; the N-1 edges of the path
carry exactly 0 bits of structural information. Free Edges = T_{n-1} - (n-1) =
T_{n-2}. Data is compressed by encoding it as localized intransitivity defects (edge
flips that break the laminar flow); these defects appear only at the diffraction
peaks of the quasicrystalline Wythoff wave, so file sizes scale with the topological
genus of the data's randomness rather than its length.

5. Arithmetic Graph Horizons and Dynamical Sinks. Summand graph G+ and multiplicand
graph G* (X+Y=Z or X*Y=Z for X<Y). Horizon lines for X=Y: additive (Z,2Z);
multiplicative (Z,Z^2); cubic (W,W^3); tesseract (W,W^4). Isolated trapping nodes:
The 63 Bottleneck (2D): 63 short-circuits Bang's theorem because it forms a closed
loop switching between a multiplication edge and an addition edge (3 x 7 x 3 = 63),
recycling old prime factors. The 341 Bottleneck (4D): 341 stalls the primitive prime
generation of the 4x+1 sequence because it is a perfect geometric sum of consecutive
dimensions (4^4+4^3+4^2+4^1+4^0), a zero-friction state in G_4. The Rational
3-Cycle: x^2+Q traps chaos at x_0=-7/4, Q=-29/16 because the vector of the inverse
square horizon (Z,Z^2) cancels the spatial translation after exactly three steps.

6. The Unified Arithmetic Cycle {f,+,*,g}: with the modular stabilizer q = x (mod AB),
the loop [+] pure addition -> geometric staircases -> [g] multiplicative feedback ->
breaks exponential weight -> [mod AB] modular tiling -> back to addition; mathematics
is a recursive geometry where macro-growth fractures back into local modular grids.
Proposed extensions: a python T_{n-2} tournament data compressor; how the 8-step Bott
periodicity cycle transforms the boundary layer from tetration to pentation; the
quasicrystalline Fourier transform matrix for the mod 30 primes framework.

SECOND PART (Collatz peaks). The endpoints 1, 26, 80 are bound to squares of an
arithmetic progression with common difference 4: Seed 5: sink 1 = 1^2; Seed 7: peak
26 = 5^2+1; Seed 23: peak 80 = 9^2-1. Bases 1,5,9 (4k+1), perturbations (0,+1,-1).
The accelerated Syracuse paths track the boundaries of perfect squares as they
scale; mirrors 196 (14th square, sum of the first 12 primes). Next predicted node:
base 13, square 169; 160 (binary reduction valve for 23, 160/2^5 = 5) sits against
this horizon; the compression deficit between 169 and 160 forces high-energy states
to collapse. Integration into the non-local descent proof: the ordered carries B_L
are locked to these quadratic horizons; peaks alternate tightly around squares of
the 4k+1 family so growth is bound by their area; binary division 2^{K_L} slices
space exponentially, so peaks run out of volume, forcing 3^L n + B_L < 2^{K_L} n to
stabilize.
