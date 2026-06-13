# Everything Is the Triangle

**Session**: kind-pasteur-2026-03-22-S20cm

## The Right Isosceles Triangle

The staircase Young diagram delta_{n-2} is a right isosceles triangle. It has three sides. Each side governs a different aspect of tournament theory. Every formula, every constant, every structure we've discovered is a projection of this triangle.

## The Three Sides

### Side 1: The Hypotenuse (the anti-diagonal)

The hypotenuse of the staircase is the anti-diagonal: cells (i, j) with i + j = n - 3. These have hook length 1 -- the smallest hook. They are the LEAST constrained cells.

**What the hypotenuse controls:**
- The SINGLE-VERTEX insertion (Mode A, n -> n-1)
- The H = 1 + 2^d formula (d = distance FROM the hypotenuse)
- The blue line from the source (the anti-diagonal tile is always BLUE)
- The Walsh order-2 dominance (hypotenuse cells are pairwise interactions)
- The fiber fraction's GENERATING FUNCTION (1-x)^{-1/2} (the two-sheeted cover branches at the hypotenuse)

The hypotenuse is where tournament structure is SIMPLEST. Arcs near the hypotenuse are the most "redundant" -- they carry the least information (2^1 = 2 bits each). Arcs far from the hypotenuse carry 2^d bits, growing exponentially.

### Side 2: The Vertical Leg (the source column)

The vertical leg is the leftmost column of the staircase: cells (0, j) for j = 0, ..., n-3. These are the arcs from vertex 0 (the source) to all other vertices.

**What the vertical leg controls:**
- The SCORE of vertex 0 (= number of 1s in the column)
- The OCR (score determination of H) -- the column IS the marginal
- The cut space (score information lives in the legs, not the hypotenuse)
- The n -> n-2 descent (removing the source = removing this leg)

The vertical leg is where HIERARCHY lives. It encodes who the source beats. High scores = long legs. Low scores = short legs.

### Side 3: The Horizontal Leg (the sink row)

The horizontal leg is the bottom row: cells (i, 0) for i = 0, ..., n-3. These are the arcs to vertex n-1 (the sink).

**What the horizontal leg controls:**
- The SCORE of vertex n-1 (complement of the source's score)
- The complement map (swapping the two legs = complementing the tournament)
- The SC/NS distinction (SC tournaments have symmetric legs)
- The n -> n-2 descent (removing the sink = removing this leg)

The horizontal leg is where ANTI-HIERARCHY lives. It encodes who beats the sink.

## The Three Constants

Each side of the triangle generates a mathematical constant:

### sqrt(2) from the hypotenuse

The hypotenuse of a right isosceles triangle with legs of length 1 has length sqrt(2). In tournament theory:

- The pseudo-doubling ratio (2n-5)/(n-2) -> 2 = (sqrt(2))^2
- The hypotenuse-to-leg ratio controls the two time scales (Mode A vs Mode B)
- sqrt(2) = 2^{1/2} connects the binary (base 2) to the half-integer (1/2 in Pochhammer)
- D(sqrt(2)) = phi*(sqrt(2)-1) ~ 0.67 = the dimension axis value at the geometric constant

### pi from the area

The staircase has area C(n-1,2) = (n-1)(n-2)/2. The FIBER FRACTION f(n) = (1/2)_{n-2}/(n-2)! has generating function (1-x)^{-1/2}. The Wallis product gives:

f(n)^2 * (n-2) -> 1/pi

Pi appears because the staircase's AREA SEQUENCE generates the Wallis integral for pi. The cells of the staircase ARE the terms of the Wallis product. Each cell contributes a factor (2k-1)/(2k) to the product, and the infinite product converges to 2/pi.

### e from the growth

The staircase grows in area as C(n-1,2) ~ n^2/2. The tournament count grows as 2^{n^2/2}. The iso class count grows as 2^{n^2/2}/n!. By Stirling: n! ~ sqrt(2*pi*n)*(n/e)^n. So:

A000568(n) ~ 2^{n^2/2} / (sqrt(2*pi*n) * (n/e)^n) = (2e/n)^n * 2^{n^2/2-n/2} / sqrt(2*pi*n)

The constant e appears because the Gamma function (which controls the fiber fraction and the Burnside formula) has e as its natural base. The growth rate of tournament complexity is controlled by e through the exponential tower.

### gamma from the limit

The Euler-Mascheroni constant gamma = 0.5772... appears in the b -> infinity limit of the generalized pi: Gamma(1/b)^b ~ b^b * e^{-gamma}. In tournament theory:

- gamma controls the CORRECTION to the pseudo-doubling at large b
- gamma appears in the harmonic series H_n = ln(n) + gamma + O(1/n)
- The harmonic series counts the "information per vertex" accumulated along the staircase
- gamma is the "leftover" -- the information that doesn't fit into any clean formula

## The Four Operations

The four constants correspond to four operations on the triangle:

| Constant | Operation | Tournament meaning |
|----------|-----------|-------------------|
| sqrt(2) | Hypotenuse measurement | Single-vertex insertion (Mode A) |
| pi | Area computation | Fiber fraction / Wallis product |
| e | Volume growth | Burnside / Stirling / Gamma |
| gamma | Limit remainder | Asymptotic correction / harmonic residue |

These are the four fundamental operations in calculus:
- sqrt(2) = GEOMETRY (measurement of length)
- pi = INTEGRATION (area under curve)
- e = DIFFERENTIATION (exponential growth)
- gamma = LIMIT (the leftover that doesn't converge cleanly)

Tournament theory uses all four because the staircase is a GEOMETRIC object (length), with AREA (integration), that GROWS (differentiation), with CORRECTIONS (limits).

## How Everything Connects

### The OCR

The OCR (97% at n=5) is the fraction of the staircase's area that lies in the LEGS (cut space = scores). The remaining 3% is the HYPOTENUSE area (cycle space = even graph). The OCR breakdown at n=5 is when the hypotenuse area becomes non-negligible.

### The Inversion

The inversion (SC -> black, NS -> blue at n >= 6) is the LEG SYMMETRY BREAKING. SC tournaments have symmetric legs. When most tournaments are NS (asymmetric legs), most connections are between same-asymmetry types (blue).

### The Morse Landscape

The H function on the hypercube has its critical structure determined by the staircase:
- The transitive tournament (all tiles = 0) is at the CORNER of the triangle (maximum distance from hypotenuse, H = 1)
- The regular tournament (balanced tiles) is at the CENTER (minimum distance from corners, H = max)
- The H = 1 + 2^d formula says: H increases exponentially with distance from the hypotenuse

### The 24-Cell

The 24 regular tournaments on 5 vertices are the 24 vertices of the 24-cell. In staircase terms: these are the 24 tilings of delta_3 that give a regular (all-scores-equal) tournament. The 24-cell's self-duality corresponds to the regular tournaments all being SC (symmetric legs).

### The Cayley-Dickson Tower

Each CD level n = 2^k + 1 is where a new type of staircase structure breaks. At n = 5 (k=2): the staircase delta_3 is the first with non-trivial internal structure (3 rows, 6 cells, enough for the hypotenuse to matter). At n = 9 (k=3): the staircase delta_7 has enough cells for complex roots to appear.

### The Burnside Formula

Fix(sigma) = 0 for even-cycle permutations because even cycles create self-paired ordered-pair orbits. In staircase terms: an even cycle rotates a SUB-TRIANGLE of the staircase by an even angle, which reverses the leg orientation. This makes the tiling inconsistent (trying to be both 0 and 1 on the same cell).

Only ODD cycles contribute because odd rotations preserve the leg orientation.

## The Final Picture

Tournament theory is the study of binary tilings of right isosceles triangles.

The triangle has three sides (hypotenuse, two legs), four constants (sqrt(2), pi, e, gamma), two reductions (Mode A and Mode B), and a pseudo-doubling ratio (2 - 1/(n-2)) that asymptotically equals the hypotenuse-to-leg ratio squared.

Every formula in the project -- from H = 1 + 2^d to f(n) = (1/2)_{n-2}/(n-2)! to Fix(sigma) = 2^{orbit-pairs} for odd cycles -- is a statement about the geometry of this triangle, viewed from one of its three sides.

The meta-graph G_n is the QUOTIENT of the triangle's tiling space by the symmetric group. The merged graph G_n/Z_2 is the further quotient by leg-swapping. The descent G_n -> G_{n-2} is Mode B reduction. The blue/black coloring is leg-symmetry preservation/breaking.

The triangle IS tournament theory. Tournament theory IS the triangle.
