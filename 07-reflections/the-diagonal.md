# The Diagonal

**Session:** kind-pasteur-2026-03-21-S18s
**Arising from:** The dimension axis (S18q/r), the observation that sqrt(2) relates to doubling

---

## D(sqrt(2)) = phi * (sqrt(2) - 1)

The tessellation dimension of the diagonal is the golden ratio times the diagonal excess:

D(sqrt(2)) = phi * (sqrt(2) - 1) = 1.618... * 0.414... = 0.670...

This is exact (within the linear interpolation of the dimension axis on [1, phi]).

### What it says

**sqrt(2) - 1** is the amount by which the diagonal exceeds the side. It is the "cost" of doubling the area. In the unit square, the diagonal is 1 + 0.414... — the side plus 41.4% more.

**phi** is the self-similar ratio. It is the number where the whole is to the larger part as the larger part is to the smaller.

Their product is the DIMENSION OF DOUBLING: 0.670. This is the tessellation dimension of the operation "double the area." It is about two-thirds of an edge — the diagonal is two-thirds of the way from having no structure (D=0) to having edge structure (D=1).

### Why two-thirds matters

2/3 is the cooperation fraction at n=5 — the fraction of the attention matrix that is symmetric (not tournament). The diagonal sits at the cooperation fraction of the edge. Doubling creates a structure that is two-thirds cooperative and one-third competitive.

---

## The Doubling Root Formula

For large n, the n-th root of 2 has tessellation dimension:

D(2^{1/n}) ~ ln(2) * phi / n = ln(2) / (n * (phi - 1))

The constant is **ln(2) * phi = 0.693... * 1.618... = 1.122...**

This means: the tessellation dimension of the n-th doubling root is inversely proportional to n, with proportionality constant ln(2) * phi. The golden ratio mediates between the logarithm (information) and the dimension (geometry).

### The formula

D(2^{1/n}) * n -> ln(2) * phi as n -> infinity

This says: **n dimensions of the n-th root of 2 = one bit of information times the golden ratio.**

The "total dimension" (D * n) of splitting a doubling into n equal parts converges to a universal constant: the bit (ln 2) dressed by the golden ratio (phi). This is the dimensional cost of dividing doubling equally.

---

## What doubling IS on the dimension axis

The number 2 is at D = infinity. It is the universal hyperbolic point. sqrt(2) = 2^{1/2} is at D = 0.670. The cube root 2^{1/3} is at D = 0.421. The fourth root 2^{1/4} is at D = 0.306.

Taking roots of 2 DESCENDS the dimension axis. Each root pulls 2 down from infinity toward zero. The n-th root of 2 sits at D ~ ln(2)*phi/n, approaching zero as n -> infinity.

This is the Cayley-Dickson tower in reverse. The CD tower DOUBLES dimension at each step (1 -> 2 -> 4 -> 8 -> ...), generating primes by adding 1 (2, 3, 5, 9, 17, ...). The n-th root of 2 UNDOES this doubling, descending from the infinite-dimensional evaluation point back toward the trivial.

sqrt(2) undoes ONE doubling: it descends from D=infinity to D=0.670.
The cube root undoes 1.5 doublings: D=0.421.
The fourth root undoes 2 doublings: D=0.306.

Each undone doubling costs about D = 0.33 = 1/3 of an edge dimension. The diagonal is worth two undone doublings: D(sqrt(2)) = 0.670 ~ 2 * 0.335.

---

## Logs of things that double

The user's observation: "logs of things double." ln(2) = 0.693. ln(4) = 1.386 = 2*ln(2). ln(8) = 2.079 = 3*ln(2). The logarithm converts doubling to ADDING.

On the dimension axis:
- ln(2) = 0.693 has D(ln(2)) = 0. The logarithm of 2 is at the ZERO of the dimension axis. It is pure information, no geometry.
- 2 itself has D(2) = infinity. The number 2 is at the top.
- sqrt(2) = e^{ln(2)/2} has D = 0.670. Halfway (in log scale) between 1 and 2.

The logarithm is the map from the number to its information content. On the dimension axis, ln maps:
- D = infinity (the number 2) to D = 0 (the number ln(2) = 0.693)
- D = 1 (the number phi) to D(ln(phi)) = D(0.481) = 0 (sub-structural)

The logarithm sends EVERY point on the dimension axis to the sub-structural region (D < 1). It is a DIMENSIONAL COLLAPSE: the log of any number in [1, 2] is in [0, 0.693], which maps to D in [0, 0].

This means: **the logarithm destroys tessellation dimension.** It converts geometric structure (D >= 1) to pure information (D = 0). The information content of a structure has no geometry of its own.

And the exponential does the reverse: it CREATES tessellation dimension from information. e^{0.693} = 2 (D = infinity). e^{0.481} = phi (D = 1). The exponential inflates the dimension axis.

---

## The square root vs the logarithm

Both sqrt and ln reduce the dimension, but differently:

**sqrt** is the DIMENSIONAL PROJECTOR (from S18q): it halves the exponent, projecting from D = infinity down to finite D. D(sqrt(pi)) = 1.70, D(sqrt(e)) = 1.14. It preserves SOME structure.

**ln** is the DIMENSIONAL ANNIHILATOR: it converts any D to essentially D = 0. ln(2) = 0.693 has D = 0. ln(phi) = 0.481 has D = 0. It destroys ALL structure.

The difference: sqrt(x) = x^{1/2} takes the HALF-POWER. ln(x) takes the... well, ln is not a power function. It is the INVERSE of exponentiation. And on the dimension axis, inversion is annihilation.

This makes sense: the exponential CREATES dimension (e^x goes from D=0 to D=infinity very fast), and its inverse DESTROYS dimension. The square root only partially destroys it because it is a FRACTIONAL power, not an inverse.

---

## The diagonal as the fundamental doubling witness

sqrt(2) witnesses the operation of doubling the area. Its tessellation dimension D = phi*(sqrt(2)-1) = 0.670 says: doubling an area creates a structure with 67% of an edge's worth of geometric content.

This is the most basic geometric operation — constructing the diagonal of a square — and the dimension axis assigns it a precise position: two-thirds of the way to having edge structure.

The ancient Greeks discovered that sqrt(2) is irrational — that the diagonal is incommensurable with the side. On the dimension axis, this incommensurability appears as D(sqrt(2)) being irrational (phi*(sqrt(2)-1) is irrational). The diagonal's tessellation dimension, like its length, does not fit into a rational framework.

And the formula D(sqrt(2)) = phi*(sqrt(2)-1) ties the two most ancient irrationals together: the diagonal (sqrt(2), from doubling the square) and the golden ratio (phi, from dividing a line in extreme and mean ratio). Their meeting point — at tessellation dimension 0.670 — is where doubling and self-similarity first touch.

---

*The diagonal is two-thirds of an edge. Not metrically — in tessellation dimension. The ancient discovery that sqrt(2) is irrational was the discovery that the dimension of doubling is incommensurable with the dimension of proportion. The formula D(sqrt(2)) = phi*(sqrt(2)-1) says it all: the golden ratio (edge boundary) times the diagonal excess (doubling cost) gives the position of the fundamental geometric operation on the axis that runs from nothing (D=0) through edges (D=1) through triangles (D=2) to everything (D=infinity at x=2). The logarithm annihilates this dimension, converting geometry to information. The square root halves it, projecting infinity to finitude. And the diagonal — the simplest non-trivial geometric construction — sits at exactly the position where two-thirds of edge structure has formed, the rest still waiting to crystallize.*
