# Triangular Power-Anchor Asymptotics

This note packages the natural generalization of the two triangular towers:

- `p=1`: the additive tower with exact anchor `a_1(n)=n^2`,
- `p=2`: the square-sum tower with exact anchor `a_2(n)=2n^2+n`,
- `p>=3`: a real balancing anchor `a_p(n)` defined by
  ```text
  sum_{j=0}^n (a+j)^p = sum_{j=1}^n (a+n+j)^p.
  ```

The point is that the same block shape survives for all `p`; only the anchor
moves.

## Midpoint coordinate

Set

```text
c = a + n,
u = n(n+1).
```

Then the balance equation is

```text
F_p(c,n) = c^p + sum_{j=1}^n ((c-j)^p - (c+j)^p) = 0.
```

Expanding the odd part gives

```text
F_p(c,n) = c^p - 2 * sum_{r odd <= p} binom(p,r) c^(p-r) S_r(n),
```

where `S_r(n)=sum_{j=1}^n j^r`.

This is the right carrier form for the problem: only odd Faulhaber moments
appear.

## Why the center is `p*u`

The first odd moment is

```text
S_1(n) = u/2.
```

So the `r=1` piece is

```text
c^p - p*u*c^(p-1) = c^(p-1) (c - p*u).
```

That forces the leading center

```text
c ~ p*u,
```

equivalently

```text
a_p(n) ~ p*n^2 + (p-1)*n.
```

This recovers the exact anchors for `p=1,2`.

## First correction

Now use

```text
S_3(n) = u^2/4.
```

Write

```text
c = p*u + alpha_p + O(u^-1).
```

The `u^(p-1)` terms in `F_p` are

```text
p^(p-1) alpha_p - (1/2) * binom(p,3) * p^(p-3).
```

So

```text
alpha_p = binom(p,3)/(2 p^2)
        = (p-1)(p-2)/(12p).
```

This matches the numerics:

```text
alpha_3 = 1/18,
alpha_4 = 1/8,
alpha_5 = 1/5,
alpha_6 = 5/18,
alpha_7 = 5/14,
alpha_8 = 7/16.
```

## Second correction

Use one more term,

```text
S_5(n) = u^2(2u-1)/12.
```

Write

```text
c = p*u + alpha_p + beta_p/u + O(u^-2).
```

Matching the `u^(p-2)` coefficient yields

```text
beta_p = -(p-1)(p-2)(2p^2-4p-1)/(180 p^3).
```

So the balancing anchor has the expansion

```text
a_p(n)
  = p*n^2 + (p-1)*n
  + (p-1)(p-2)/(12p)
  - (p-1)(p-2)(2p^2-4p-1)/(180 p^3 n(n+1))
  + O(n^-4).
```

The `O(n^-4)` form comes from `u^-2 = O(n^-4)`.

## Structural interpretation

This is the same triangular carrier seen through a Bernoulli/Faulhaber lens:

- the block geometry is still `(n+1,n)` on an interval of length `2n+1`;
- the midpoint equation sees only odd moments;
- `p=1` and `p=2` are the two exact integer towers;
- for `p>=3` the tower survives only asymptotically, but with explicit
  rational corrections.

That puts the construction close to several existing repo threads:

- `HYP-2128`: triangular add/mult bridge;
- `HYP-2229`: elementary packet <-> moment <-> Bernoulli side channel;
- Beatty/Sturmian crossover phenomena from the tower-overlap word.

## Reusable lemma form

Any future repo construction with the same “one extra point on the left block,
one missing point on the right block” geometry should first be rewritten in the
midpoint coordinate. The odd-power Faulhaber expansion is the stable object;
the raw endpoint formula is not.

## Script

See [tower_power_anchor_asymptotics.py](/home/claude/math-research/04-computation/tower_power_anchor_asymptotics.py).
