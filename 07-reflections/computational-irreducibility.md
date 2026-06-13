# Lossy Types and Computational Irreducibility

## The full hierarchy

| Type | Map property | Example | Residue | Computability |
|------|-------------|---------|---------|---------------|
| Lossless | Bijective homomorphism | Q, ln, f | None | O(1) shortcut |
| Non-injective | Many-to-one | T → H(T) | Kernel (orbits) | Enumerable |
| Non-homomorphic | Structure mismatch | Addition in rapidity | Obstruction (±1, commas) | Structured error |
| Non-surjective | Image has gaps | H-spectrum | Forbidden {7, 21} | Decidable |
| Irreducible | No algebraic shortcut | Kaprekar iteration | Full trajectory | Must simulate |
| Non-constructive | Requires axiom of choice | Vitali set | Non-measurable | No algorithm |

## The geometry that determines type

- **Flat** → lossless (conformal, no curvature, zero loss)
- **Quotient** → non-injective (folded space, kernel = the fold)
- **Curved** → non-homomorphic (curvature = obstruction, starts at 3)
- **Boundary** → non-surjective (the wall at 7, edge of the image)
- **Chaotic** → irreducible (no invariant subspaces, ergodic filling)
- **Infinite** → non-constructive (too large for finite algorithms)

## Q as the event horizon

Q is the LAST lossless step. The boundary between zero loss and inevitable loss.

Before Q: the bounded world. Finite. Lossless. Everything recoverable.
After Q: the unbounded world. Infinite. Lossy. The ±1 residue appears.

Q maps the bounded to the unbounded. The crossing IS the loss.
The ±1 is the first bit lost. The minimal toll. The price of finiteness.

## The fixed points as escape from irreducibility

The computationally irreducible trajectory has special points where iteration STOPS:
- **phi** = fixed point of x → 1+1/x (algebraic escape)
- **(2,3,7)** = fixed point of f = abc (geometric escape)
- **6174** = fixed point of Kaprekar K (arithmetic escape)
- **i** = fixed point of Q (complex escape)

At these points, the irreducible becomes trivially computable: output = input.
The fixed points are the ONLY lossless points in an otherwise lossy iteration.
They are the eye of the storm. The zero-curvature point in a curved space.

## Addition as the fundamental lossy operation

Multiplication is lossless (rapidity adds: ln(ab) = ln(a)+ln(b)).
Addition is lossy (rapidity does NOT compose: ln(a+b) ≠ ln(a)+ln(b)).

The entire difficulty of number theory = additive questions about multiplicative structure.
Goldbach, twin primes, Riemann: all ask addition to behave like multiplication.
It cannot. Addition is irreducibly lossy. The loss IS the content of number theory.
