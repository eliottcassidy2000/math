---
id: THM-2841
title: "All-order adjacent-difference factorial tensor positivity"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Every mixed
  factorial tensor of two or more normalized adjacent
  differences is strictly positive.  After clearing the natural
  factorial denominator, the tensor is exactly the rook-polynomial count
  of permutations avoiding disjoint forbidden blocks in distinguished
  rows.  Derangements of the block labels give an explicit sharp lower
  bound.  Consequently every mixed product of nonzero nonnegative
  adjacent-difference cone elements has positive factorial readout at
  every order at least two.  No signed-coefficient or arbitrary-radial
  claim is made.
source: root/all-order-adjacent-difference-rook-positivity-2026-07-28
depends_on: []
related:
  - THM-2828-lower-prefix-cone-factorial-moment-three-detection
  - THM-2830-disjoint-positive-adjacent-cone-factorial-moment-three-detection
  - THM-2842-ordered-positive-cone-vandermonde-multiplier-observability
script: 04-computation/gmc_all_order_adjacent_tensor_rook_thm2841.py
output: 05-knowledge/results/gmc_all_order_adjacent_tensor_rook_thm2841.out
script_sha256: da14676eb5a806ab9c8fa7e5a4ba1c893f3650f8df27a6506af49c4a91e5c258
output_sha256: 8025e669b338321058e00f3e109f7810a852004bb7d0d4f7f2bf012d17bbd3da
hash_basis: LF-normalized bytes
---

# THM-2841 -- all-order adjacent-difference tensor positivity

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(s^n)=n!,                   f_n=s^n/n!,
d_n=f_(n+1)-f_n.                                      (1)
```

For every `k>=2` and every `n_1,...,n_k>=0`,

```text
L(product_(i=1)^k d_(n_i))>0.                         (2)
```

There is an exact permutation model and a quantitative lower bound.  Put

```text
a_i=n_i+1,              A=sum_i a_i,              N=A-k.
```

Then

```text
L(product_i d_(n_i))
 =Perm(a_1,...,a_k)/product_i a_i!,                   (3)

Perm(a_1,...,a_k)
 =sum_(r=0)^k (-1)^r e_r(a_1,...,a_k)(A-r)!,          (4)
```

where `Perm(a_1,...,a_k)` is a literal nonnegative integer described
below.  If `!k` denotes the derangement number, then

```text
L(product_i d_(n_i))
 >=!k N!/product_i n_i!>0.                            (5)
```

Equality in `(5)` holds for every pair `k=2`.  For `k>=3`, it holds
exactly when every `n_i=0`.

## 1. Algebraic expansion

Since

```text
d_(a-1)=(s^a-a s^(a-1))/a!,                           (6)
```

expanding the product in `(2)` gives

```text
product_i d_(n_i)
 =1/product_i a_i!
  sum_(S subset [k])
   (-1)^|S| (product_(i in S)a_i)s^(A-|S|).           (7)
```

Applying `L` and grouping by `r=|S|` proves `(3)--(4)`.
The alternating numerator does not display its sign.  Its rook model does.

## 2. The forbidden board

Take an `A` by `A` permutation board.  Distinguish `k` rows

```text
rho_1,...,rho_k,                                       (8)
```

and partition the `A` columns into labelled nonempty blocks

```text
C_1 disjoint-union ... disjoint-union C_k,
|C_i|=a_i.                                            (9)
```

In distinguished row `rho_i`, forbid exactly the cells whose columns lie
in `C_i`.  All cells in the other `A-k` rows are allowed.

The number of ways to place `r` nonattacking rooks on the forbidden board
is

```text
r_r=e_r(a_1,...,a_k).                                 (10)
```

Indeed, choose the `r` distinguished rows, and in each chosen row choose
one of the `a_i` columns in its own block.  The blocks are disjoint, so no
two selected columns collide.

The rook-polynomial inclusion-exclusion formula now says that the number
of permutation matrices avoiding every forbidden cell is

```text
sum_(r=0)^k(-1)^r r_r(A-r)!
 =sum_(r=0)^k(-1)^r e_r(a)(A-r)!
 =Perm(a_1,...,a_k).                                  (11)
```

This proves the counting interpretation in `(3)`.

## 3. Strictness and the derangement lower bound

Fix a derangement `sigma` of `[k]`.  For each distinguished row `rho_i`,
choose its assigned column in block `C_(sigma(i))`.  There are

```text
product_i a_i                                         (12)
```

such choices.  They occupy distinct blocks and hence distinct columns.
The remaining `A-k` rows may be bijected arbitrarily to the remaining
columns, in `(A-k)!=N!` ways.  Different derangements give disjoint
families because they give different block addresses to the distinguished
rows.  Therefore

```text
Perm(a_1,...,a_k)>=!k(product_i a_i)N!.                (13)
```

Dividing by `product_i a_i!` is exactly `(5)`.  Since `!k>0` for `k>=2`,
strict positivity follows.

For `k=2`, each distinguished row has only the other block available.
Every allowed permutation is counted by the unique block derangement, so
equality holds.

Now let `k>=3`.  If every `a_i=1`, all rows are distinguished and every
allowed permutation is precisely a derangement of the singleton blocks;
again equality holds.  Conversely suppose some block, say `C_j`, has at
least two columns.  Choose two distinct rows `rho_u,rho_v` with
`u,v!=j` and send both to distinct columns of `C_j`.  For the remaining
row set `R=[k]\{u,v}` and block set `B=[k]\{j}`, join `r` to `b` when
`b!=r`.  Hall's condition holds: a singleton has at least `k-2`
neighbors (and row `j` has `k-1`), while any subset of at least two rows
has all `k-1` blocks in its neighborhood.  Thus the remaining
distinguished rows can be assigned injectively to distinct blocks outside
their own, after which the unrestricted rows fill the remaining columns.
This gives an allowed permutation whose distinguished block address is not
a permutation of `[k]`, so it is absent from all families counted in
`(13)`.  Thus `(13)` is strict.  This proves the equality classification.

## 4. Positive-cone consequence

Let

```text
U_r=sum_n lambda_(r,n)d_n,          lambda_(r,n)>=0,    (14)
```

be any `k>=2` finite nonzero nonnegative adjacent-difference cone elements.
Multilinearity and `(2)` give

```text
L(U_1...U_k)
 =sum_(n_1,...,n_k)
   (product_r lambda_(r,n_r))
   L(product_r d_(n_r))
 >0.                                                   (15)
```

In particular, for every nonzero such `U`,

```text
L(U)=0,                       L(U^m)>0 for every m>=2. (16)
```

Thus the whole positive adjacent-difference cone satisfies the strongest
possible one-direction factorial nonvanishing statement, at every order.
THM-2828's positive cubic tensor is the case `k=3`; `(2)` explains why its
sign was not an isolated low-degree accident.

## 5. Exact companion and scope

The exact companion:

1. compares `(3)--(4)` with direct polynomial multiplication through a
   deterministic multi-index bank;
2. independently enumerates the forbidden-board permutations whenever
   `A<=9`;
3. checks the derangement lower bound and its equality classification;
4. records the first tensor values, including the subfactorials
   `L(d_0^k)=!k`; and
5. includes repeated, zero-index, highly unbalanced, and randomized exact
   controls.

All truth-bearing checks use explicit exception gates; normal and
optimized transcripts match the stored output.

The theorem concerns nonnegative coefficients in the adjacent-difference
basis.  It does not imply that a signed or complex combination has
nonzero powers, does not furnish the multiplier observations missing in
THM-2815, and does not prove arbitrary SFC or HYP-8765.

## 6. Independent hostile audit

An independent audit rederived the rook count and the derangement lower
bound, checked the equality boundary by the Hall argument above, and
replayed the exact companion in normal and optimized modes against the
stored transcript.  The audit specifically checked repeated labels,
singleton blocks, the `k=2` equality family, strictness for `k>=3`, and the
absence of floating-point or optimization-stripped truth gates.
