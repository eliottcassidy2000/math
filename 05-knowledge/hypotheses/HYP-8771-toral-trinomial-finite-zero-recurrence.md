# HYP-8771 -- finite zero recurrence for primitive toral trinomials

**Status: OPEN for higher primitive charges; announced externally for the
symmetric Legendre case.** Owner: codex-2026-07-21.

For coprime positive `p0,q0`, put `r=p0+q0` and

```text
A_m(kappa)=sum_(0<=k<=m/r)
  m!/((q0*k)!(p0*k)!(m-r*k)!) * kappa^k.             (1)
```

## Conjecture

For every fixed `kappa in C`, the set

```text
{m>=1 : A_m(kappa)=0}                                (2)
```

is finite.

This is exactly the toral factor on THM-2021's proportional central-resonance
slice `h=kappa*b^r`. If (2) is finite, nullity forces `L(b^m)=0` eventually;
EMP applied to a power of `b` gives `b=0`, and then one charged endpoint dies.
Thus HYP-8771 proves NC2 on the full proportional central hypersurface for the
chosen primitive charge pair.

## What is already proved

1. Every `A_m` has constant term one. Hence the union of all exceptional
   `kappa` is countable, and THM-2021 proves NC2 off that set.
2. At `p0=q0=1`,

   ```text
   A_m(kappa)=(1-4kappa)^(m/2)
              P_m((1-4kappa)^(-1/2)),
   ```

   so every exceptional `kappa` is negative real. The February 2026 IAS
   announcement for Mangoubi--Kadets--Weller Weiser states precisely that a
   fixed nonzero point is a zero of only finitely many Legendre polynomials.
   Subject to that announced external theorem, HYP-8771 is settled for `(1,1)`.
3. Exact computation finds no common polynomial factor at all among 1,379 level
   pairs: `(1,1)` through `m=36`, `(1,2)` through 30, `(1,3)` through 26, and
   `(2,3)` through 22. Pairwise coprimality is stronger than finite recurrence,
   but this is evidence only.

## Why three-term descent is not the statement

For `(1,1)` the Legendre recurrence forbids **consecutive** common zeros. An
infinite set of isolated zero levels would still block the simple EMP deduction.
The needed predicate is tail finiteness, not adjacency. This is the same
controlled-forgetting warning as THM-2005: projecting the full zero profile to
nearest-neighbour data discards the theorem predicate.

For higher charges the return sequence has holonomic order `p0+q0`, consistent
with THM-1670. There is no reason to expect a three-term orthogonal descent.

## Two proposed proof routes

### 1. Finite-place separation

If `A_m(kappa)=0` even once, then `kappa` is algebraic. At a nonarchimedean place,
the channel valuations are

```text
v( m!/((q0k)!(p0k)!(m-rk)!) ) + k v(kappa).          (3)
```

Legendre's factorial formula turns (3) into a carry/digit-sum statistic. A unique
minimum prevents cancellation. THM-2020 develops exactly this finite-place
certificate for general Gaussian channels. The target is to prove that for a
fixed algebraic `kappa`, some place has a unique minimum in (3) for every
sufficiently large `m`, or after splitting `m` into finitely many residue classes.

### 2. Algebraic-function singularities

The ordinary generating function is the constant term

```text
sum_(m>=0) A_m(kappa)t^m
 = CT_u 1/(1-t(1+x u^p0+y u^(-q0))),
   x^q0 y^p0=kappa.                                  (4)
```

It is algebraic. A zero-recurrence theorem for coefficients of this one-parameter
family, classifying the root-of-unity collision cases of its dominant
singularities, would also prove (2). This is the higher-charge analogue of the
spectral-geometric argument announced for Legendre polynomials.

## Tournament Analysis and challenged vertex assumption

Candidate vertices considered: return channels, individual roots, moment levels,
and proof obligations. Moment levels are the faithful choice because (2) is a
property of the level set. The pairwise observable is `deg gcd(A_m,A_n)`; the
chronological tournament edge is flipped on a shared-root event, and chronological
order is the tie Hamiltonian path. The four exact scans have transitive score
histograms, no flips, no directed triangles, singleton SCCs, and one Hamiltonian
path. The quotient detects recurrence but forgets which root recurs and therefore
cannot by itself prove tail finiteness.

The challenged assumption is that the natural vertices are channels. They are
not: channel sets grow with `m` and a channel-wise tournament destroys the fixed
`kappa` recurrence relation. Moment levels are the correct carriers.

## Cross-links

THM-2021 (exact NC2 reduction and symmetric Legendre transform); concurrent
THM-2018 (independent symmetric factorization); THM-2020 (finite-place channel
separation); THM-1670 (higher recurrence order); HYP-8766 (finite resonance
band); HYP-8769 (Sheffer no-common-zero); THM-2005 (profile versus scalar
endpoint). Script/output:
`gmc2_proportional_legendre_finite_recurrence_thm2021.py/.out`.
