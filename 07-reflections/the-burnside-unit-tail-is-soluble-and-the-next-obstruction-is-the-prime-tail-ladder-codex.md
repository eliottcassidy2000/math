# The Burnside unit tail is soluble, and the next obstruction is the prime-tail ladder

**Source:** codex-2026-06-15-S12. Prompt: consider the odd-partition split
`lambda = mu U 1^r`, the resulting `B[m,t]` kernel, and push on what that means.

## What the split really does

The quoted formula is not just a rearrangement. It isolates the only completely
universal part of the odd-part Burnside sum.

For A000568, the odd partition `lambda` is a cycle type. A `1`-cycle is special
because it interacts with everything by the same gcd:

```text
gcd(1,L)=1.
```

That makes the `1`-tail soluble. Once the odd partition is written as

```text
lambda = mu U 1^r,
```

all dependence on the tail passes only through:

- its size `r`,
- and the core length `t = ell(mu)`.

The rest of the arithmetic stays in the core.

So the Burnside sum becomes:

```text
hard arithmetic core  *  universal unit-tail background.
```

That is why the outer partition count at `n=100` really collapses from
`444793` odd partitions to exactly `834` active `(m,t)` states.

There is also an arithmetic bonus hiding inside that collapse:

```text
K[m,t] = m! B[m,t] = sum_{mu} (m!/z(mu)) 2^e(mu)
```

is an integer. The factor `m!/z(mu)` is the labelled permutation count of cycle
type `mu`. So the core kernel is not a formal rational residue. It is a genuine
weighted labelled-cycle enumerator.

## The exact support is already clean

There is no hidden sparsity beyond the obvious core constraints. The active
state space is exactly:

```text
(m,t) = (0,0) or [m>=3t and m==t mod 2].
```

Nothing subtler is needed. Every feasible `(m,t)` occurs.

This matters because it says the tail peel is complete: the remaining difficulty
is not support combinatorics, only the arithmetic content of the kernel
coefficients.

## The universal tail function

After the peel, the ordinary generating function factors through one universal
series

```text
F(x) = sum_{r>=0} 2^C(r,2) x^r / r!,
```

and the whole A000568 series becomes

```text
A(x) = sum_{m,t} B[m,t] x^m F(2^t x).
```

So the tail is a fixed analytic background, while the kernel only decides which
rescaled copies `F(2^t x)` appear and with what coefficients.

Even better, the tail satisfies

```text
F'(x)=F(2x).
```

That is a very rigid object. The uncertainty is no longer in the tail at all.

## The next rung is not another partition, but one divisor statistic

The first new progress beyond the quoted formula is that the same mechanism
continues for odd primes.

If we peel `3`-cycles from the core,

```text
mu = nu U 3^s,
```

the interaction is no longer controlled just by length, because

```text
gcd(3,L) = 1 or 3.
```

But it still does not require the full residual partition. It needs only one
extra statistic:

```text
c_3(nu) = number of remaining parts divisible by 3.
```

That is the key refinement:

```text
1-tail  -> needs only t = ell(core),
3-tail  -> needs only t plus c_3,
prime p tail -> needs only t plus c_p.
```

So the obstruction after peeling `1`s is not "full partitions again." It is a
prime-tail ladder of divisor statistics.

At `n=100`, the next rung is still meaningfully compressed. The residual odd
partitions with parts at least `5` number `7551`, while the exact `3`-free
kernel uses only `2049` active `(m,t,c_3)` states. So the ladder does not
collapse to the raw partition DP immediately after `p=1`; there is still a real
quotient at `p=3`.

## Why this is the right next object

This is the partition-side analogue of the FKN / tiling story that just landed
upstream.

- In the tiling cube, the universal easy part is the transitive / level-1 side.
- In the Burnside sum, the universal easy part is the `1`-tail.
- In both cases, the hard part is what remains after that easy background is
  peeled away.

For A000568 the next obstruction is arithmetic, not combinatorial:

```text
which remaining parts are divisible by 3, by 5, by 7, ...
```

That is where the kernel lives.

## The right open question

The true next problem is not "can we peel the `1`-tail?" That is done.

It is:

```text
how far can the odd-prime tail ladder be pushed before the state space becomes
the full divisor-sum DP again?
```

The exact `3`-tail formula shows there is still substantial structure left after
`B[m,t]`. The kernel is not featureless. It is organized by a ladder of prime
divisibility statistics.
