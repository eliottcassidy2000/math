---
id: THM-561
title: Miss-zeta factorial moments invert to missed-count atoms
status: PROVED
source: codex-2026-06-21
depends_on:
  - THM-406
  - THM-534
related:
  - HYP-2716
  - HYP-2718
  - HYP-2680
---

# THM-561 - Miss-Zeta Factorial Atom Inversion

## Statement

Let `U` be a finite set of size `m`.  Let `M` be a random subset of `U`,
and write `T=|M|`.  For each `R subset U`, define the miss-zeta coordinate

```text
z(R) = Pr(R subset M).
```

For `0 <= j <= m`, set

```text
Z_j = sum_{|R|=j} z(R).
```

Then `Z_j` is the `j`th binomial, or normalized falling-factorial, moment of
`T`:

```text
Z_j = E[binom(T,j)] = E[(T)_j]/j!.
```

Moreover the atom masses `p_t=Pr(T=t)` are recovered by exact binomial
inversion:

```text
p_t = sum_{j=t}^m (-1)^(j-t) binom(j,t) Z_j.
```

Equivalently, the origin atom is the alternating finite difference

```text
p_0 = sum_{j=0}^m (-1)^j Z_j.
```

In falling/rising factorial language this is the same basis twice:
`binom(T,j)=(T)_j/j!` and `(T)_j=(-1)^j(-T)^{overline j}`.

## Proof

For a fixed realization `M=A` with `|A|=t`, exactly `binom(t,j)` subsets
`R` of size `j` satisfy `R subset A`.  Therefore

```text
Z_j
 = sum_{|R|=j} Pr(R subset M)
 = E[# {R subset M : |R|=j}]
 = E[binom(T,j)].
```

Writing `p_t=Pr(T=t)` gives the triangular system

```text
Z_j = sum_{t=j}^m binom(t,j) p_t.
```

Now compute the proposed inverse:

```text
sum_{j=t}^m (-1)^(j-t) binom(j,t) Z_j
 = sum_{s=t}^m p_s
     sum_{j=t}^s (-1)^(j-t) binom(j,t) binom(s,j).
```

Use

```text
binom(j,t) binom(s,j) = binom(s,t) binom(s-t,j-t).
```

The inner sum is therefore

```text
binom(s,t) sum_{r=0}^{s-t} (-1)^r binom(s-t,r)
 = binom(s,t) (1-1)^(s-t).
```

This is `1` when `s=t` and `0` otherwise.  Hence the inverse returns exactly
`p_t`.  The formula for `p_0` is the special case `t=0`.

## LRC14 Corollary

In the seven-sector LRC14 model, take `U={1,...,6}` and let `M(x)` be the set
of missed inner sectors at time `x`.  If two laws are compared, for example
the actual line law and the shared-slow-`x` independent-carrier product law,
write

```text
W_j = Z_j^prod - Z_j^actual,
Q_t = p_t^prod - p_t^actual.
```

Then

```text
Q_t = sum_{j=t}^6 (-1)^(j-t) binom(j,t) W_j,
Q_0 = sum_{j=0}^6 (-1)^j W_j.
```

Thus the all-inner-sector cover discrepancy is exactly the origin atom of the
missed-count law:

```text
Q_0 = ProductCover - p0.
```

This is the same scalar called the `h=6` Krawtchouk/MacWilliams top character
in HYP-2716, because `K_6(j)=(-1)^j`.

## Verification

`04-computation/lrc14_multiblock_miss_zeta_layers_codex_20260621.py`
implements this inversion exactly as `atom_profile_from_factorial`, asserts
`Q_0=ProductCover-p0`, and asserts `sum_t Q_t=0` on the current HYP-2715 row
bank.
