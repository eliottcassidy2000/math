---
id: THM-2041
title: "Frobenius stability of exact-period projectors"
status: >
  PROVED in two compatible forms. Every coprime exponent-dilation-invariant
  cyclic twist has its twisted constant term carried to the p-th power, and
  every Fourier-packet idempotent in a good-characteristic finite abelian
  group algebra commutes with injective Frobenius. Ramanujan exact-order
  projectors, phase-function energies, parity packets, and conductor packets
  are corollaries. This supplies the whole-layer preservation half of the
  THM-2022 mechanism, but not LRC's missing seed or pointwise-witness theorem.
  The coprimality hypothesis is sharp.
source: codex-2026-07-21-NC2-transfer
related:
  - THM-2022
  - HYP-8800
  - HYP-2978
  - HYP-2979
  - HYP-3036
  - LEM-032
  - LEM-033
scripts:
  - 04-computation/lrc_frobenius_exact_period_projector_codex_20260721.py
  - 04-computation/frobenius_exact_period_projectors_codex_20260721.py
outputs:
  - 05-knowledge/results/lrc_frobenius_exact_period_projector_codex_20260721.out
  - 05-knowledge/results/frobenius_exact_period_projectors_codex_20260721.out
---

# THM-2041 -- Frobenius stability of exact-period projectors

This theorem isolates the part of the NC2 proof that really does transfer to
the LRC Ramanujan packets. It preserves an entire tied exact-period layer; it
does not choose a single character or assume termwise noncancellation.

## 0. Finite-abelian packet theorem

The cyclic statement is a specialization of a more general fact. Let `G` be
a finite abelian group, let `k` have characteristic `p` with `p` not dividing
`|G|`, and put `A=k[G]`. If `e` is any idempotent in `A`, then for every
`X in A` and `r>=0`,

```text
e X^(p^r) = (eX)^(p^r),                                  (0)
```

and the right side is nonzero exactly when `eX` is nonzero. Indeed, `A` is
commutative, `e^(p^r)=e`, and Maschke plus commutativity makes `A` a finite
reduced algebra, so Frobenius is injective.

After a splitting-field extension, `e` may be the sum of the primitive
Fourier idempotents in any character packet that descends to `k`. Exact
order, parity, and conductor packets are examples. Thus the theorem preserves
the `G=U_q` parity/conductor projections used in LEM-032 and the conductor
projections paired with valuation grades in LEM-033 whenever `p` does not
divide `phi(q)`. The valuation-grade encoding itself is separate data and is
not automatically transported by this theorem.

There is a scope distinction. On the multiplicative group `U_q`, group-ring
Frobenius sends the basis element `[u]` to `[u^p]`, whereas an LRC multiplier
often acts by translation `u -> pu`. Both preserve the relevant character
support when `p` is a unit, but they are not the same physical operation.
An application must state its encoding rather than identify them silently.

## 1. Cyclic stable-twist formulation

Fix a cyclic period `N`, a prime `p` coprime to `N`, and a field `k` of
characteristic `p`. In

```text
A=k[u,u^(-1)]/(u^N-1),
```

let `sigma_p(u^a)=u^(pa)`. This is an automorphism. If an integral twist
`Theta` has coefficients in the prime field after reduction and satisfies

```text
sigma_p(Theta)=Theta,                                     (1a)
```

then Freshman's dream gives `Theta^p=Theta`. Hence, for every `Lambda in A`
and `m>=0`,

```text
Theta Lambda^(pm)=(Theta Lambda^m)^p.                     (1b)
```

Because `p` is invertible modulo `N`, the only exponent class whose p-fold
dilate is zero is zero itself. With `ct_N` denoting the identity coefficient,

```text
ct_N(Theta Lambda^(pm))=ct_N(Theta Lambda^m)^p.            (1c)
```

Iteration propagates every nonzero twisted constant term through levels
`p^r m`. This formulation is useful even when `Theta` is a stable kernel
rather than an idempotent. The primitive-support kernel

```text
sum_(a in (Z/NZ)^*) u^a
```

is stable because multiplication by `p` permutes the units. The Ramanujan
identity `c_N(pn)=c_N(n)` is the character-side version of the same fact.

## 2. Exact-period idempotents

Let `k` be a field of characteristic `p`, let `p` not divide `N`, and write

```text
A = k[C_N] = k[u]/(u^N-1).
```

For every `d|N`, let

```text
c_d(x) = sum_(a mod d, gcd(a,d)=1) exp(2*pi*i*a*x/d)
       = sum_(r|gcd(d,x)) r*mu(d/r)
```

be the integer Ramanujan sum, reduced in `k`, and define

```text
e_d = (1/N) sum_(x mod N) c_d(x) u^x in A.                (1)
```

> **Exact-period decomposition.** The elements `e_d`, `d|N`, are pairwise
> orthogonal idempotents, their sum is one, and `e_d A` is exactly the sum of
> the character spaces whose characters have order `d`.

To see this, extend scalars to a splitting field and evaluate at every
`N`-th root. If the character `u -> zeta_N^j` has exact order `r`, then

```text
e_d(zeta_N^j) = 1 if r=d, and 0 otherwise.               (2)
```

Indeed, expand `c_d(x)` as the sum over primitive `d`-th roots and use finite
geometric-series orthogonality. Because `p` does not divide `N`, `u^N-1` is
separable, so these evaluations detect equality. Formula (2) proves all the
claims.

## 3. Frobenius fixes every whole exact-period layer

Let `Fr_p(F)=F^p` be absolute Frobenius on `A`. Then

```text
e_d^p
 = N^(-p) sum_(x mod N) c_d(x)^p u^(px)
 = N^(-1) sum_(x mod N) c_d(x) u^(px)
 = e_d.                                                   (3)
```

The second equality uses that `N^(-1)` and every reduced integer `c_d(x)` lie
in the prime field. For the last equality, multiplication by `p` permutes
`Z/NZ`, and

```text
c_d(p^(-1)x)=c_d(x)                                      (4)
```

because `p` is a unit modulo `d`; equivalently, it permutes the primitive
residues in the defining sum.

Writing `Pi_d(F)=e_d F`, (3) gives the exact identity

```text
Pi_d(F^p) = e_d F^p = (e_d F)^p = Pi_d(F)^p.              (5)
```

Moreover, `A` is reduced because `u^N-1` is separable. Frobenius on a reduced
ring of characteristic `p` is injective. Consequently, for every `r>=0`,

```text
Pi_d(F) != 0  iff  Pi_d(F^(p^r)) != 0.                    (6)
```

Thus Frobenius preserves not just the set of exact periods, but the complete
possibly-colliding sum in each exact-period component.

## 4. Phase functions and Ramanujan energy

For a function `f:Z/NZ -> k`, encode it by

```text
F_f(u)=sum_x f(x)u^x.
```

The coefficient function of `F_f^p` is

```text
(Fr_p f)(x)=f(p^(-1)x)^p.                                 (7)
```

Hence (5) is the exact-period projector identity for cyclic signals. The
Ramanujan shell energy used in HYP-2979 is

```text
E_d(f)=sum_(x,y mod N) f(x)f(y)c_d(x-y).                  (8)
```

Changing variables by the unit `p` and using (4) gives

```text
E_d(Fr_p f)=E_d(f)^p.                                     (9)
```

Also,

```text
E_d(f)/N = [u^0] e_d F_f(u)F_f(u^(-1)).                  (10)
```

Formula (10) explains both the usefulness and the limitation of the scalar
energy. It is a readout of the exact-period component, but over a finite field
its bilinear pairing can be isotropic. The group-algebra component `e_dF_f`
is therefore the stronger nonvanishing carrier.

## 5. The finite-place form used by seed-and-amplify arguments

Suppose the coefficients of `F` lie in a number field and `Pi_d(F)` is
nonzero in characteristic zero. After clearing denominators, only finitely
many finite places divide `N` or kill every nonzero coordinate of this layer.
At every other place `pfrak|p`, reduction and (5) give

```text
Pi_d(F^(p^r)) = Pi_d(F)^(p^r) != 0 mod pfrak.             (11)
```

This is precisely the whole-layer half of THM-2022: first obtain a nonzero
characteristic-zero seed, then choose a good finite place, then amplify the
whole tied layer. It is important that the seed comes before reduction.
THM-1630 records that a naive positive-characteristic version of the
Duistermaat--van der Kallen seed theorem is false.

## 6. Sharp bad-prime boundary

The condition `p` not dividing `N` cannot be removed. If `N=p`, then in
`F_p[C_p]`

```text
h=u-1 != 0,            h^p=u^p-1=0.                      (12)
```

The group algebra is nonreduced, no division by `N` is available in (1), and
Frobenius need not preserve nonzero layers. In particular, the apex prime
`7` is a bad prime for the period-`14` LRC deck. Any period-14 Frobenius
argument must use another prime or retain the nonsemisimple 7-primary
filtration as extra data.

## 7. Exact LRC scope

For LRC, (5) proves preservation but not production or interpretation:

```text
nonzero base safe/avoidance layer
  + endpoint-labelled implication to a strict lonely time
  + good-prime Frobenius preservation
  => a usable seed-and-amplify step.                       (13)
```

Only the third line is supplied here. The distinction is visible in the
stored HYP-2979 audit. The AP row and the row `{1,...,11,13,84}` both have
the same nonzero danger-count energy

```text
E_14(N)=216,
```

but their primitive phase states differ, `((danger,boundary)=(0,2))` versus
`((1,2))`, and their weak-safe energies are `6` versus `0`. A nonzero danger
layer therefore does not imply a strict safe phase. Endpoint ownership,
open-versus-boundary status, or a proved dual inequality remains essential.

HYP-3036's positive primitive safe-residue count at some `q<=13` is already a
direct witness; it needs no amplification. The live use of this theorem is
instead to make any future familywise seed certificate stable under every
coprime dilation without losing collisions inside its exact-period shell.

## 8. Verification and assumption challenge

The stable-twist script checks primitive-support and Ramanujan kernels for
`2<=N<=24`, every tested coprime prime through `23`, five moment levels, and
a non-invariant negative control. The independent idempotent script checks all
`111` pairs `(N,d)` with `d|N`, `N<=30`, against
every prime `p<=47` coprime to `N`: `1479` idempotent/Frobenius cases and
`4437` each of projector commutation, nonzero preservation, and Ramanujan
energy dilation. It also checks five witnesses of (12).

Tournament vertices are proof carriers rather than runners. The observable
ranks retention of the strict LRC predicate, endpoint owners, a nonzero whole
layer, Frobenius stability, linearity, and compression. The resulting path is

```text
endpoint-labelled primitive packet
  > primitive safe-residue count
  > primitive-period idempotent
  > Ramanujan shell energy
  > raw Ramanujan trace
  > raw denominator blockedness.
```

The quotient preserves seed-and-amplify obligations but destroys actual arc
geometry. The challenged assumption is explicit: arithmetic stability of a
projector is not the same as an LRC witness theorem.
