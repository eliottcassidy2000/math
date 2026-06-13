---
source: codex-2026-06-12
status: SYNTHESIS + verified computation; HYP-2454
tags: [triangular, faulhaber, odd-moments, square-pyramidal, cuboid, irreducibility, tournaments]
---

# Faulhaber odd moments and square-pyramidal cuboids

The new Faulhaber computation sharpens the triangular power-balance story in
exactly the direction the user's packing intuition wanted.

## 1. p=1 and p=2 are the last pure packing towers

For the interval balance midpoint `c=a+n`, the carrier is always

```text
u=n(n+1)=2T_n.
```

At `p=1`, this is literally the staircase rectangle:

```text
2*T_n = u.
```

At `p=2`, Faulhaber gives the square-pyramidal cuboid:

```text
P_2(n)=1^2+...+n^2 = n(n+1)(2n+1)/6,
6*P_2(n)=u(2n+1).
```

So six `n`-step square pyramids fill the `u x (2n+1)` cuboid.  The prompt's
"one doubled and four regular" picture is the same volume identity rewritten
as `2+4=6`.

This fits the older simplex/cuboid theme in the repo, but only up to `p=2`.
Those are the last rows where the packing picture closes without correction.

## 2. Higher p keeps the triangular carrier but loses the pure scissors lane

With `c=a+n`, the exact power-balance equation is

```text
c^p = 2 * sum_{r odd} binom(p,r) c^(p-r) S_r(n).
```

Only odd Faulhaber moments survive.  The deformation away from the exact
rectangle/cuboid towers is therefore not arbitrary noise; it is a controlled
odd-moment correction.

The first two terms already explain the cutoff:

```text
a_p(n)
= p*n^2 + (p-1)*n
  + (p-1)(p-2)/(12p)
  - (p-1)(p-2)(2p^2-4p-1)/(180 p^3 n(n+1))
  + O(n^-4).
```

At `p=1` and `p=2`, the correction vanishes.  Starting at `p=3`, the same
triangular carrier survives, but only after Bernoulli/Faulhaber address terms
are added.  The geometry no longer closes as a pure cuboid packing.

## 3. The polynomial face says the same thing

At `n=1`,

```text
D_p(C,1)=C^p+(C-1)^p-(C+1)^p.
```

For even `p`, `C=0` is a forced symmetry root, so remove that factor.  Then:

```text
p=1 -> C-2,
p=2 -> C-4,
```

and from `p>=3` onward the live factors are already irreducible in the stored
window (`p<=20`) by finite-field certificates.

This is the polynomial version of the packing cutoff.  The old tower mechanism
was not only "integral" or "pretty"; it was literally split.  Once the odd
Faulhaber corrections appear, the first visible live polynomial stops
factorizing.

## 4. Why the tournament should use carriers, not powers

If the vertices are just powers `p`, the quotient is too lossy.  The useful
objects are proof carriers:

```text
rectangle_packing,
square_pyramidal_cuboid,
odd_moment_reduction,
alpha_beta_asymptotic,
n1_live_irreducibility,
78_90_support_shadow,
global_bracket_tail,
hidden_lift_transfer.
```

That quotient keeps the things that can plausibly survive into proofs:
exactness, packing geometry, asymptotic reach, irreducibility signal, and
support-transfer value.  It destroys the decorative row data, which is exactly
the right loss function here.

The resulting tournament ranks `odd_moment_reduction` first, not the pure
packing lanes and not the `78/90` shadow by itself.  That feels correct: the
odd-moment identity is the place where the geometry, the asymptotic, and the
irreducibility story finally meet.
