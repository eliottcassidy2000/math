---
id: THM-2199
title: "Effective positive-subspace safe floor and rank lift"
status: >
  PROVED RELATIVE TO LRCUpTo13 through THM-2053. Every rational subspace
  T of Q^13 of dimension at least two that contains a positive integer row
  has 1/14-safe Haar mass at least 182^(-13) on its associated torus.
  Consequently, a zero-Haar positive row in T has an integer relation outside
  the annihilator of T of height at most
  H=78*182^13=18750922831149193194381342621696. Taking T to be the
  orthogonal complement of all relations through height H forces twelve
  independent relations by that same height. Hence every primitive zero-Haar
  row has maximum speed at most 12^6 H^12. The bounds are explicit but
  enormous; the finite primitive locus is not enumerated or emptied, so this
  does not prove LRC(14).
source: codex-2026-07-24-effective-carry-lock-rank-ladder
depends_on:
  - THM-2053-rank-two-parameter-plane-geodesic-terminal
  - THM-2193-uniform-rank-six-safe-torus-floor
  - LRCUpTo13
related:
  - THM-2196-bounded-relation-cone-circuit-atlas-and-carry-lock-rank-ladder
  - THM-780-haar-compactness-safe-measure-floor
  - THM-763-strict-finite-height-for-tight-lrc-instances
---

# THM-2199 -- effective positive-subspace safe floor and rank lift

The qualitative carry-lock argument has a uniform quantitative form. Its
finite object is not a list of rational subspaces: it is one heavy cell of
Haar measure on the phase torus.

Put

```text
J=[1/14,13/14] subset T=R/Z.                            (1)
```

For a rational subspace `T_0 subset Q^13`, define

```text
Lambda_(T_0)=T_0^perp intersection Z^13,
K_(T_0)={x in T^13:a.x=0 for every a in Lambda_(T_0)}. (2)
```

Thus `K_(T_0)` is the connected torus whose rational speed space is `T_0`.

## 1. A uniform positive-subspace safe floor

> **Positive-subspace floor.** Assume `LRCUpTo13`. If `T_0` has dimension
> at least two and contains a positive integer row, then
>
> ```text
> Haar_(K_(T_0))(K_(T_0) intersection J^13)>=182^(-13). (3)
> ```
>
> In fact the set on the left contains a measurable packet on which every
> coordinate clearance is strictly greater than `1/14`.

Choose a positive integer row `v in T_0`, divide by its coordinate gcd so
that it is primitive, and choose an integer `z in T_0` independent of it.
The rational plane

```text
U=span_Q{v,z} subset T_0                              (4)
```

contains a positive row. THM-2053 supplies a point

```text
x_* in K_U,             min_i ||(x_*)_i||>=1/13.     (5)
```

Since `U subset T_0`, its relation lattice contains that of `T_0`, and hence

```text
K_U subset K_(T_0).                                  (6)
```

It remains to thicken one point uniformly inside an arbitrary compact
subtorus. Let `K=K_(T_0)`, with normalized Haar measure `mu_K`, and partition
the ambient `T^13` into the `182^13` half-open coordinate cubes

```text
product_(i=1)^13 [j_i/182,(j_i+1)/182).               (7)
```

Their intersections with `K` partition `K`, so one cube `C` has

```text
mu_K(A)>=182^(-13),             A=K intersection C.  (8)
```

Choose `x_0 in A` and put `B=A-x_0`. Translation invariance gives

```text
B subset K,             mu_K(B)=mu_K(A)>=182^(-13).  (9)
```

The two representatives of each coordinate of `x in A` and `x_0` lie in
the same half-open interval of length `1/182`. Therefore every `b in B`
satisfies

```text
||b_i||<1/182.                                        (10)
```

The circle triangle inequality, (5), and

```text
1/13-1/182=1/14                                      (11)
```

now give, for every `b in B`,

```text
||(x_*+b)_i||>1/14                 for every i.       (12)
```

Since `x_*+B subset K`, equations (9) and (12) prove (3). This is the
subgroup-Haar version of THM-780's heavy-phase-cell argument.

## 2. An explicit bounded relation outside every proper speed subspace

Put

```text
D=182^13=240396446553194784543350546432,
N=39D+1,
H=2N-2=78D=18750922831149193194381342621696.          (13)
```

> **Effective rank lift.** Let `T_0 subset Q^13` satisfy the hypotheses of
> Section 1, and let `v in T_0 intersection Z_(>0)^13` have zero Haar-safe
> mass on its Kronecker line. Then some
>
> ```text
> a in Z^13,       a.v=0,       ||a||_infinity<=H,    (14)
> ```
>
> does not belong to `Lambda_(T_0)`.

Use the normalized squared-Fejer approximant `f_N` to `f=1_J` from
THM-2193, Section 6. Its Fourier degree is `2N-2=H`, and

```text
eta_N=||f_N-f||_1<3/(2N),
26 eta_N<39/N<1/D.                                  (15)
```

Suppose (14) fails. Every Fourier frequency `a` of the tensor

```text
F_N(x)=product_(i=1)^13 f_N(x_i)                     (16)
```

which is resonant on the line `t |-> tv` has
`a.v=0` and `||a||_infinity<=H`, hence belongs to
`Lambda_(T_0)`. Conversely, every `a in Lambda_(T_0)` annihilates `v`
because `v in T_0`. Thus the line average and the `K_(T_0)` Haar average
of (16) agree exactly, frequency by frequency.

The line has zero safe mass. Product telescoping compares its `F_N` average
with its `1_(J^13)` average at cost at most `13 eta_N`. Every coordinate
character is nontrivial on `K_(T_0)`: if `e_i` belonged to
`Lambda_(T_0)`, it would annihilate the positive row `v`, which is
impossible. The same product telescoping on the torus therefore costs at
most another `13 eta_N`. The exact Fourier equality yields

```text
Haar_(K_(T_0))(K_(T_0) intersection J^13)
   <=26 eta_N
   <1/D.                                               (17)
```

This contradicts (3), proving the lift.

## 3. One step forces rank twelve

For a positive integer row `v`, write

```text
W_H(v)=span_Q{a in Z^13:a.v=0, ||a||_infinity<=H}.    (18)
```

Suppose `v` has zero Haar-safe mass and `dim_Q W_H(v)<=11`. Set

```text
T_0=W_H(v)^perp.                                      (19)
```

Then `T_0` is rational, has dimension at least two, and contains the positive
row `v`. Its annihilator is

```text
Lambda_(T_0)=W_H(v) intersection Z^13.               (20)
```

Section 2 produces a relation `a` satisfying (14) but not (20). Yet the
height condition in (14) puts `a` in `W_H(v)`, hence in (20), a
contradiction. Therefore:

> **Explicit rank-twelve harvest.** Every positive thirteen-speed row with
> zero Haar-safe mass satisfies
>
> ```text
> dim_Q W_H(v)>=12,
> H=78*182^13=18750922831149193194381342621696.        (21)
> ```

The relation space of a nonzero thirteen-speed row has dimension twelve, so
equality holds in (21).

Choose twelve independent relations of height at most `H` and place them
in a `12 x 13` matrix. Its signed maximal minors generate the
one-dimensional integer kernel. After dividing their gcd, Hadamard's
inequality bounds the primitive positive kernel generator `v_prim` by

```text
max_i (v_prim)_i
 <=(sqrt(12)H)^12
 =12^6 H^12.                                          (22)
```

This is an explicit finite primitive box and an explicit relation-rank
certificate, but it is vastly larger than the speed ceiling in THM-763.
Neither theorem enumerates or eliminates the box. The contribution here is
that every carry-lock rank lift is simultaneously available at one universal
Fourier height; there is no residual compactness minimum and no uncomputed
sequence `H_8,...,H_12`. QED.
