# The Dimension Axis

**Session:** kind-pasteur-2026-03-21-S18q
**Arising from:** The {p,q} hierarchy (S18p), the cascade g_{p+1}(rho_p) = -1

---

## The Discovery

The k-nacci constants form a ladder on the real line:

1 < phi < tau < rho_4 < rho_5 < ... < 2

Each rho_p is the Euclidean boundary of the {p,q} tessellation theory. This ladder IS a dimension axis. The "tessellation dimension" of a real number x > 1 is its position on this ladder, interpolated linearly between consecutive rho_p.

---

## The Axis

| Constant | Value | D(x) | Meaning |
|----------|-------|------|---------|
| 1 (unit) | 1.000 | 0.00 | Trivial: no structure |
| sqrt(2) | 1.414 | 0.67 | Sub-edge: below all thresholds |
| e^{1/e} | 1.445 | 0.72 | Sub-edge: the "self-root" |
| **phi** | **1.618** | **1.00** | **Edge boundary (Fibonacci)** |
| sqrt(e) | 1.649 | 1.14 | Between edges and triangles |
| sqrt(3) | 1.732 | 1.52 | Hexagonal: between edges and triangles |
| sqrt(pi) | 1.772 | 1.70 | Circle: approaching triangles |
| **tau** | **1.839** | **2.00** | **Triangle boundary (tournaments)** |
| e^{2/pi} | 1.890 | 2.58 | Between triangles and squares |
| **rho_4** | **1.928** | **3.00** | **Square boundary (bipartite)** |
| **rho_5** | **1.966** | **4.00** | **Pentagon boundary (Petersen)** |
| ... | ... | ... | ... |
| **2** | **2.000** | **infinity** | **Universal hyperbolic: ALL dimensions** |
| e | 2.718 | infinity | Super-hyperbolic |
| pi | 3.142 | infinity | Super-hyperbolic |

---

## What This Means

### Integer dimensions are algebraic

The points D = 1, 2, 3, 4, ... are occupied by the k-nacci constants phi, tau, rho_4, rho_5, .... These are all ALGEBRAIC numbers — roots of integer-coefficient polynomials (the k-nacci equations x^p = x^{p-1} + ... + 1). The integer dimensions of tessellation theory are algebraic.

### Fractional dimensions are transcendental

The points between integer dimensions are occupied by transcendental constants: sqrt(e) at D = 1.14, sqrt(3) at D = 1.52, sqrt(pi) at D = 1.70. These constants parameterize "fractal-dimensional" theories — structures whose atoms are between p-gons and (p+1)-gons. A D = 1.5 theory has atoms that are "between edges and triangles" — perhaps something like a fractal lattice or a random graph at a critical percolation threshold.

### The circle lives at D = 1.70

sqrt(pi) = 1.772... has tessellation dimension 1.70. The circle constant parameterizes a structure that is 70% of the way from edges to triangles. This makes geometric sense: a circle is "between" a polygon with 2 sides (degenerate) and a polygon with 3 sides (triangle). The circle is the LIMIT of polygons, and its position at D = 1.70 on the dimension axis reflects its position in the polygon hierarchy.

### The hexagonal lattice lives at D = 1.52

sqrt(3) = 1.732... has dimension 1.52. The hexagonal lattice ({6,3} tiling) involves 120-degree angles and is built from the constant sqrt(3) (the height of an equilateral triangle). Its position at D = 1.52 — about halfway between edges and triangles — reflects its role as the "flat" (Euclidean) version of triangular structure. It is the honeycomb: triangles arranged without curvature.

### The tournament evaluation point is infinite-dimensional

D(2) = infinity. The evaluation point x = 2, where H = I(Omega, 2), has infinite tessellation dimension. This means the tournament's Hamiltonian path count probes ALL cycle lengths simultaneously — it is not specialized to any particular face type but integrates over the entire dimension hierarchy.

This is WHY the SRCP needs all odd cycle lengths (c3, c5, c7, ...) to determine H: each cycle length is one dimension, and the evaluation at x = 2 "sees" all of them equally (because g_p(2) = 1 for all p).

### e and pi are super-hyperbolic

D(e) = D(pi) = infinity. Both transcendentals are above 2, hence hyperbolic for ALL tessellation theories. Evaluating the independence polynomial at x = e or x = pi would give a "super-tournament" invariant that is even more dominated by large independent sets than H = I(Omega, 2).

In the hard-core lattice gas: evaluating at fugacity lambda = e means each independent set element contributes a factor of e instead of 2. This is deeper in the crystal phase — more ordered, more dominated by the maximal packing.

---

## The Cascade Identity on the Dimension Axis

The cascade theorem g_{p+1}(rho_p) = -1 has a beautiful interpretation on the dimension axis:

**Each integer-dimensional constant is exactly one quantum spherical for the dimension above it.**

phi (D=1) is one quantum spherical for the triangle theory (g_3(phi) = -1).
tau (D=2) is one quantum spherical for the square theory (g_4(tau) = -1).
rho_4 (D=3) is one quantum spherical for the pentagon theory (g_5(rho_4) = -1).

The "quantum" is always -1, regardless of dimension. And g_p(2) = +1, regardless of dimension. So the dimension axis has a clean structure:

At each integer dimension d = p-1:
- g_p(rho_p) = 0 (the Euclidean boundary)
- g_{p+1}(rho_p) = -1 (one quantum spherical for the next dimension)
- g_p(2) = +1 (one quantum hyperbolic from this dimension's perspective)

The three values {-1, 0, +1} at each dimension are the complete classification: spherical, flat, hyperbolic.

---

## The Depth of 2

The number 2 is not "just" the smallest prime, or the binary choice, or the tournament fugacity. On the dimension axis, it is the point at infinity — the accumulation point of all k-nacci constants. It is the number that is **one quantum hyperbolic for every dimension simultaneously.**

No algebraic k-nacci constant has this property. Each rho_p is Euclidean for its own dimension and spherical for the next. Only the LIMIT (2) is universally hyperbolic. This is why tournament theory is governed by the modular group PSL(2,Z): the evaluation point 2 sits at the boundary of the hyperbolic plane, which is where the modular group acts.

The interval [phi, 2) is the OPEN interval of "theories that are hyperbolic for some dimensions but not all." The point 2 is the CLOSURE — the unique point that is hyperbolic for all. Tournament theory lives at this closure point. It is the theory that sees everything.

---

*The real number line between 1 and 2 is a dimension axis. The algebraic k-nacci constants sit at integer dimensions like lattice points. The transcendentals sit at fractional dimensions like irrational lattice sites. The golden ratio is dimension 1 (edges). The tribonacci constant is dimension 2 (triangles). The number 2 is dimension infinity (all cycles). Every constant in mathematics — pi, e, sqrt(3), every algebraic and transcendental number — has a tessellation dimension that determines which theories it is spherical for, which it is Euclidean for, and which it is hyperbolic for. The tournament evaluation point x = 2, sitting at dimension infinity, is the unique point that is one quantum hyperbolic for every theory simultaneously. This is why it sees everything. This is why H = I(Omega, 2) is the right formula.*
