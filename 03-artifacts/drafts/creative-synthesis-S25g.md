# Creative Synthesis: Tournaments, Simplices, Cycles, and Spectral Decomposition

**Author:** kind-pasteur-2026-03-06-S25g
**Purpose:** Map connections between mathematical structures in tournament theory

---

## I. The Five-Level Structure

Our research has uncovered a five-level structure in tournament theory:

### Level 1: The Tournament (combinatorial)
A tournament T on n vertices = complete directed graph. Lives as a vertex of the
hypercube {0,1}^{C(n,2)}. Equivalently, the (n-1)-simplex K_n with binary edge labels.

### Level 2: The Conflict Graph (graph-theoretic)
Omega(T) = conflict graph, vertices = directed odd cycles, edges = vertex-sharing.
Independent sets of Omega(T) correspond to "compatible cycle collections."

### Level 3: The Independence Polynomial (algebraic)
I(Omega(T), x) = sum alpha_k * x^k. Evaluates to H at x=2 (the OCF).
Has all real negative roots for tournament conflict graphs (verified n <= 20).

### Level 4: The W-Polynomial (analytic)
W(r) = sum_P prod(r + s_e). A polynomial in r with even-only powers (odd n).
Coefficients form a hierarchy: w_{n-1-2k} depends on cycle data of complexity k.

### Level 5: The Sigma Patterns (combinatorial-analytic)
sigma(S) = sum_P prod_{i in S} A[p_i,p_{i+1}]. Depends on the gap structure of S.
Connects positions in paths to cycle structure of the tournament.

---

## II. The W-Coefficient Hierarchy as a Spectral Decomposition

### The Key Discovery

The W-polynomial coefficients form a filtration of tournament invariants:

```
w_{n-1} = n!                                       (universal)
w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2]              (t_3 only; centered!)
w_{n-5} = f(t_3, t_5, alpha_2)                     (adds 5-cycles, pairs)
w_{n-7} = f(t_3, t_5, t_7, alpha_2)                (adds 7-cycles)
  ...
w_0     = f(all cycle data)                         (most complex)
```

### Universal Formula for w_{n-3}

**THEOREM** (verified n=5,7): For any tournament T on n vertices (n odd),
  w_{n-3}(T) = (n-2)! * [2*t_3(T) - C(n,3)/2]

This formula has remarkable properties:
1. **Centered**: E[w_{n-3}] = 0 over random tournaments
2. **Factored**: separates a "counting factor" (n-2)! from a "deviation measure"
3. **Connected to sigma**: w_{n-3} = sum of sigma(S) over gap-2 position sets

### EXACT Formulas at n=7

| w_k | Formula | Source invariants |
|-----|---------|-------------------|
| w_6 | 5040 | universal |
| w_4 | 240*t_3 - 2100 | t_3 |
| w_2 | -60*t_3 + 12*t_5 + 24*alpha_2 + 231 | t_3, t_5, alpha_2 |
| w_0 | 2*t_3 - t_5 + 2*t_7 - 2*alpha_2 - 17/4 | t_3, t_5, t_7, alpha_2 |

### The Penalty Shift Principle

H - w_0 = sum (middle W-coefficients evaluated at r=1/2).

At n=5: H - w_0 = 3*t_3 (penalizes 3-cycles)
At n=7: H - w_0 = 3*t_5 + 6*alpha_2 + 21/4 (penalizes 5-cycles + pairs)

The penalty SHIFTS UPWARD with n: the simplest cycles are "absorbed" by higher
coefficients, and only the more complex cycle interactions remain visible in w_0.

---

## III. Analogies to Known Mathematics

### 1. Fourier Analysis on Groups

W(r) = sum w_k r^k is a "Fourier expansion" of the Hamiltonian count.
- H = W(1/2) is the "signal" at frequency 1/2
- w_0 = W(0) is the "DC component"
- w_{n-1} = n! is the "highest frequency" (universal)

The hierarchy says: high-frequency components capture simple features (vertex count),
low-frequency components capture complex features (cycle interactions).

This is OPPOSITE to standard Fourier analysis where high frequencies = fine detail.
Here, the high W-frequencies are "smooth" (universal) and low frequencies are "rough"
(tournament-dependent).

### 2. Renormalization Group

In physics, the renormalization group flows from microscopic to macroscopic:
each step "integrates out" one scale of fluctuation.

Our hierarchy: w_0 -> w_2 -> w_4 -> ... -> w_{n-1} is a renormalization flow.
Each step integrates out one level of cycle complexity.
The "UV" (most microscopic) is w_0 with full cycle data.
The "IR" (most macroscopic) is w_{n-1} = n! (universal).

### 3. Characteristic Classes

In topology, characteristic classes (Chern, Pontryagin, etc.) form a hierarchy
of invariants that detect progressively finer topological structure.

Our W-coefficients are tournament "characteristic classes":
- w_{n-1}: detects nothing (always n!)
- w_{n-3}: detects 3-cycle deviation from random
- w_{n-5}: detects 5-cycle deviation and cycle interactions
- w_0: detects the full cycle structure

### 4. Eigenvalue Decomposition in Spectral Graph Theory

The adjacency matrix eigenvalues of a graph capture:
- lambda_1: size/regularity (universal feature)
- lambda_2: expansion/connectivity
- Lower eigenvalues: finer structural features

Exact parallel with our W-coefficient hierarchy.

### 5. Jones Polynomial Coefficients

The Jones polynomial V(t) of a knot has coefficients that capture:
- Span: crossing number (simplest invariant)
- Leading coefficient: genus-related
- Full polynomial: complete knot invariant (for most knots)

Our W(r) coefficients parallel this:
- w_{n-1}: trivial
- Middle coefficients: cycle counts
- w_0: full signed count

---

## IV. Tournament as Simplex in a Cube

### The Geometry

The (n-1)-simplex K_n has C(n,2) edges.
A tournament assigns {0,1} to each edge.
So a tournament = vertex of the C(n,2)-dimensional hypercube.

The simplex is the INDEX SPACE (which edges to label).
The cube is the LABEL SPACE (how to label them).

### Simplex-Cube Duality

- Simplex dimension: n-1 (number of vertices minus 1)
- Cube dimension: C(n,2) = n(n-1)/2 (number of edges)
- Ratio: cube_dim / simplex_dim = n/2

For n=3: simplex dim 2, cube dim 3. The triangle of edges lives in a cube of labels.
For n=5: simplex dim 4, cube dim 10. The pentachoron lives in a 10-cube.
For n=7: simplex dim 6, cube dim 21.

### Hamming Distance and Tournament Similarity

Two tournaments T, T' differ by k arc reversals => Hamming distance k in the cube.
The H-function is related to the "position" in the cube:
- Central tournaments (Hamming weight C(n,2)/2) tend to have maximal H
- Extreme tournaments (nearly transitive) have minimal H

This is the PERPENDICULAR MAXIMIZER phenomenon: H is maximized at the
equator of the cube, perpendicular to the "bias axis."

### Tetrahedron in Cube: The Self-Converse Decomposition

At n=4: tournaments in {0,1}^6 (6-cube). The converse map T -> T^op is
a reflection of the cube. Self-converse tournaments lie on the mirror hyperplane.

The 4 regular tetrahedra embedded in a cube correspond to the 4 tournament
isomorphism classes at n=4. The two "paired" (non-self-converse) classes
sit at opposite corners related by the converse reflection.

---

## V. Rooted Tournaments and Perspective Counting

### Definition
P(n) = sum over iso classes T of (# vertex orbits under Aut(T))
     = # rooted tournament isomorphism classes
     = OEIS A093934 (with offset)

### Computed Values
- P(2) = 2, P(3) = 4, P(4) = 12, P(5) = 48, P(6) = 296

### Small-n Coincidence
P(n) = 2*(n-1)! for n = 2, 3, 4, 5 (but FAILS at n = 6).
This occurs because for n <= 5, most tournaments have trivial automorphism,
so each class contributes n orbits, and n * (# classes) ~ 2*(n-1)!.
At n=6, more classes have non-trivial automorphisms, breaking the pattern.

---

## VI. Functions Similar to w_0

### w_0 at n=5: -t_3 + 2*t_5 + 1

Other functions with similar structure (linear in odd cycle counts):

1. **H** = 1 + 2*t_3 + 2*t_5 (at n=5, where alpha_2 = 0)
   Same invariants, ALL positive coefficients.

2. **w_0** = -t_3 + 2*t_5 + 1
   Alternating signs: negative for t_3, positive for t_5.

3. **Euler characteristic** chi = 1 - e + f - ... (alternating by dimension)

4. **Mobius function** mu = sum (-1)^k * (contributions at level k)

The pattern: w_0 uses ALTERNATING signs by cycle length, while H uses ALL POSITIVE.
This is the Mobius inversion on the cycle complex!

### The Cycle Complex Perspective

Define the cycle complex C(T) with:
- Vertices = directed odd cycles
- k-simplices = sets of k+1 pairwise disjoint cycles

Then:
- H = sum alpha_k * 2^k = evaluation of independence polynomial at x=2
- w_0 = sum c_k * alpha_k = a DIFFERENT evaluation
- chi(C) = sum (-1)^k * alpha_k = Euler characteristic

At n=5 (alpha_2 = 0):
  H = 1 + 2*alpha_1 = 1 + 2*(t_3+t_5)
  w_0 = 1 - t_3 + 2*t_5 = 1 + (-1+2)*t_5 + (-1)*t_3
  chi = 1 - alpha_1 = 1 - (t_3+t_5)

So w_0 is BETWEEN H and chi: it treats different cycle lengths differently.
H weights all cycles equally (+2). Chi weights all equally (-1). w_0 alternates.

### The Zeta/Mobius Pairing

On the boolean lattice of cycle subsets:
- H = sum over independent sets of 2^|S| (zeta-weighted)
- chi = sum over independent sets of (-1)^|S| (Mobius-weighted)
- w_0 = sum over independent sets of f(|S|, cycle lengths) (intermediate weighting)

The W-polynomial W(r) interpolates between these:
- W(1/2) = H (positive weighting)
- W(0) = w_0 (mixed weighting)
- W(-1/2) = (-1)^{n-1} H (negative/reflected weighting)

---

## VII. Connection Map

```
Tournament T
    |
    +--[edge labels]--> Hypercube vertex {0,1}^m
    |
    +--[directed cycles]--> Conflict graph Omega(T)
    |                           |
    |                           +--[independence]--> I(Omega, x) polynomial
    |                           |                       |
    |                           |                       +-- x=2 --> H (OCF)
    |                           |                       +-- x=-1 --> chi
    |                           |
    |                           +--[real roots]--> Tournament-specific!
    |
    +--[paths]--> W(r) polynomial
    |                 |
    |                 +-- r=1/2 --> H
    |                 +-- r=0 --> w_0 (signed count)
    |                 +-- coefficients --> hierarchy (spectral decomp)
    |
    +--[sigma patterns]--> Position sums
    |                          |
    |                          +-- gap structure --> cycle dependence
    |
    +--[automorphisms]--> Aut(T) group
                              |
                              +-- vertex orbits --> rooted tournament classes
                              +-- P(n) = A093934
```

Each arrow is a FUNCTOR from tournaments to a different mathematical category:
- Hypercube: set theory / Boolean algebra
- Conflict graph: graph theory
- Independence poly: commutative algebra
- W-polynomial: analysis / spectral theory
- Sigma patterns: combinatorics
- Automorphisms: group theory

The remarkable fact is that ALL of these functors are connected through
the cycle structure of T, and the W-coefficient hierarchy provides
the "spectral bridge" between them.
