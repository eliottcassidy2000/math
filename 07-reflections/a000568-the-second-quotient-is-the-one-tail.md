# A000568: the second quotient is the one-tail

The cleanest new idea was not a subtle divisor profile. It was admitting that
`1` is different from every other odd part.

Once I write an odd partition as

```text
lambda = mu union 1^r,
```

the `1`-tail stops behaving like generic partition data. Its contribution to the
Burnside exponent is universal:

```text
C(r,2) + r * (#parts(mu)).
```

That means the ones are an endpoint variable, not real core complexity.

So the second quotient is:

- full odd partition `lambda`
- core mass / core length state `(m,t)`
- hidden kernel `B[m,t]`.

This is much better than I expected. It is exact, and by `n=100` it shrinks the
outer sum from `444793` odd partitions to `834` active states.

What it does **not** do is solve the whole problem. It just relocates the
difficulty. The hard object is now the core kernel `B[m,t]`, i.e. odd partitions
with no ones. That is still a lot of structure, but it is the right structure.

Conceptually this feels like endpoint recursion in the strongest sense:
the long tail of ones is not where the entropy lives. It is an address. The real
interaction starts at `3`.
