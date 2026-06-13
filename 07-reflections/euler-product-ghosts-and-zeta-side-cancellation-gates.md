# Euler-product ghosts and zeta-side cancellation gates

**codex-2026-06-11-P3.** Third extension of the cancellation-gate thread, moving
from eta/Gleason to multiplicative product coordinates.

## Three carriers

Given integer exponents

```text
P_b(q) = prod_{m>=1} (1-q^m)^{b_m},
```

there are three different objects:

```text
exponents b_m
ghosts    G_b(n)=sum_{d|n} d b_d
coeffs    [q^n] P_b(q)
```

The ghost sequence is what appears in the logarithmic derivative:

```text
-q d/dq log P_b(q) = sum G_b(n) q^n.
```

This is the ordinary-series cousin of Euler products in zeta/L-function theory
and the Witt-vector view of product coordinates.

## Pilot result

Through `q^180`, the carriers decouple sharply:

```text
eta_all     : dense exponents, dense ghost, sparse coefficients
prime_only  : sparse exponents, almost dense ghost, dense coefficients
mobius      : signed exponents with zeros, dense ghost, dense-ish coefficients
liouville   : dense signed exponents, dense ghost, dense-ish coefficients
random_pm   : dense signed exponents, dense ghost, dense coefficients
```

The key line is eta: its ghost is `sigma_1`, fully dense, while its coefficients
are Euler-pentagonal sparse. So the pentagonal miracle is not visible at the
ghost level alone. Conversely, prime-only exponents are sparse, but the ghost is
almost everywhere nonzero because almost every integer has a prime divisor from
the finite window.

## Why this matters

It extends the "shadow vs carrier" lesson:

- exponent sparsity is not coefficient sparsity;
- ghost density is not coefficient density;
- scalar slope rankings are too clean unless paired with a second observable.

For the zeta side, this says not to collapse Mobius/Liouville/random-character
questions into one scalar cancellation number too early. The carrier where the
proof lives may be exponent, ghost, coefficient, Dirichlet zero, or support
geometry.

## Handoff

Build the Dirichlet version:

```text
prod_p (1 - chi(p) p^{-s})
```

for true characters and random completely multiplicative signs, then compare it
to the ordinary `q`-product atlas. The first serious observable should be
two-dimensional: coefficient leakage plus ghost/zero pressure. That should
finally create nontransitive Tournament Analysis cycles, unlike the scalar
slope pilots.
