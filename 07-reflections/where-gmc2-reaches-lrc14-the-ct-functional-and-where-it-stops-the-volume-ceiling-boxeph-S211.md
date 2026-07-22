# The constant-term thread needs an augmentation coordinate

*Scope-repaired audit of boxeph-2026-07-21-S211 / HYP-8840. See
MISTAKE-226 and the corrected exact script/output.*

The original session found real constant-term structure but merged three
different functionals too early. The repair is useful because it identifies
the coordinate that was lost: **augmentation**, or equivalently the stationary
observer at speed zero.

## What survives exactly

For `P(z)=sum_(v in S)z^v`,

```text
E_m(S)=CT[P(z)^m P(z^(-1))^m]
      = sum_x r_m(x)^2.
```

This is the unanchored additive `m`-energy. The single-character first-return
seed is also a constant-term statement. If `S_1+tS_2` is a genuine no-carry
digital sum at order `m`, then its energy factors as

```text
E_m(S_1+tS_2)=E_m(S_1)E_m(S_2).
```

The AP generating polynomial has its standard cyclotomic factorization. And
THM-2047 supplies the exact analytic ceiling: positive safe-set Haar measure is
a strict-exit certificate, while Euler characteristic also detects a finite
tight fiber.

These are four separate facts. None implies LRC(14).

## Why diagonal energy is not the LRC lattice

A monomial contributing to `E_m` has a relation vector `k=alpha-beta` with

```text
k dot v=0,             sum_i k_i=0.
```

The second equation is the hidden affine/augmentation constraint. The LRC
Fourier formula sums over every linear relation `k dot v=0`, with no
augmentation restriction. That is why unanchored energy is translation blind
while LRC is not.

The clean hostile control uses three affine copies of a five-term AP:

```text
A={1,2,3,4,5},   B={2,3,4,5,6},   C={1,3,5,7,9}.
```

They have identical `E_m` at every order, but exact pair-sum enumeration gives

```text
M(A)=1/6,        M(B)=1/4,        M(C)=1/2.
```

Thus even the entire diagonal energy profile cannot be a faithful Wall-A or
LRC invariant. THM-730 counts `a+b=c`; this relation has augmentation one and
is not the `m=2` diagonal-energy theorem claimed in the original synthesis.

## Two equivalent repairs

One repair adjoins the observer:

```text
S -> {0} union S.
```

Every linear relation lifts to an affine relation by setting
`k_0=-sum_i k_i`. The hostile triple above then separates already at anchored
second energy:

```text
E_2({0} union A,B,C)=146,130,106.
```

The more informative repair keeps the mixed bigrading

```text
M_(r,s)(S)=CT[P(z)^r P(z^(-1))^s].
```

It selects augmentation `r-s`; in particular

```text
THM-730's Schur count = CT[P^2 Pbar]=M_(2,1).
```

Every finite LRC relation appears in some bidegree
`(r,s)=(|k^+|,|k^-|)`. The exponential package

```text
F_S(x,y)=CT_z exp(xP(z)+yP(z^(-1)))
        =sum_(r,s>=0) M_(r,s)x^r y^s/(r!s!)
```

is therefore a full relation-lattice transform with the augmentation fugacity
`x/y` retained. Its one-coordinate factors are Bessel-type kernels. LRC safe
volume pairs the same full lattice with sinc coefficients; swapping Bessel and
sinc weights is not a positivity theorem.

## Why the other proposed bridges fail

The AP itself is a carry obstruction:

```text
{0,...,11}={0,1,2}+3{0,1,2,3},
E_2({0,...,11})=1156 != 19*44=836.
```

No-carry tensorization removes precisely the carry relations that make the AP
extremal, so it cannot rank-reduce an arbitrary thirteen-speed core. The
original computation checked AP energy maxima in only nine finite
`(N,k,m)` cases and used a made-up geometric weight as a GMC proxy. Its raw
box-truncated Fourier sum was not the Fejer limit required by THM-2047.

Finally, `M(S)=max_t min_v||vt||` is a one-variable nonsmooth maximum. Calling
it a game/Morse saddle adds no map or theorem.

## Live research thread

The corrected question is now concrete:

> Can the augmentation-graded mixed CT table, or its observer-anchored
> specialization, constrain the small-relation branch supplied by THM-2051
> strongly enough to interact with THM-2047 deletion or THM-2048 fiber tax?

This is a worthwhile wildcard. It retains the exact coordinate the original
bridge discarded and has hostile controls built in. It is not yet a Wall-A
reduction or an implication from GMC(2) to LRC(14).
