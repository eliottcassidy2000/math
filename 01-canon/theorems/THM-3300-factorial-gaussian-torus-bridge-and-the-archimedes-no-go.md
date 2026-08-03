---
id: THM-3300
title: "Factorial-Gaussian torus bridge and the Archimedes no-go"
status: >
  PROVED + VERIFIED-EXACT (bridge, closed form, no-go) + FINITE-EXACT
  (cyclic exclusion through degree three); INDEPENDENT IMMUTABLE AUDIT PENDING.
  The factorial functional IS a Gaussian moment functional:
  `L_fac(x^alpha) = alpha! = E[prod |z_i|^(2 alpha_i)]` for independent
  standard complex Gaussians, and the moment map pushes the uniform measure of
  `S^(2n-1)` to the uniform measure of `Delta_(n-1)`.  Hence the homogeneous
  Factorial-Conjecture class HFC(n) is exactly the `U(1)^n`-INVARIANT part of
  the Gaussian-moment problem in `2n` real variables.  THM-3290's Archimedes
  mechanism therefore CANNOT be transported: its annihilator is the projection
  onto torus weight zero, which is the identity on the factorial class, and its
  engine is a coefficient extraction that exists only because the test function
  has a pole in the angular variable.  What replaces it is a finite group, and
  that is exactly why the port fails: a U(1) character has infinite order, so a
  nonzero weight kills EVERY `m>=1`, while the simplex's full affine symmetry
  group is only `S_n`, whose characters have finite order `e`, so the cyclic
  route always leaves `m in e Z` alive.  For `n=3` the surviving obligation is
  `<g^(3k)>=0`, and the cyclic-eigenvector families of degree 1, 2 and 3 on
  `Delta_2` are excluded outright.
audit: >
  The exact companion verifies the bridge on 63 monomial rows for `n=1,2,3`
  with a hostile control pinning the covariance normalization, the moment-map
  identity on 43 rows, and the closed form
  `<(sum omega^j s_j)^m> = (n-1)! m!/(m+n-1)! [n|m]` for `m=1..12`.  It then
  excludes the degree-1 cyclic line explicitly, the degree-2 plane by an exact
  `Q(omega)` gcd of `<g^3>` and `<g^6>` plus its missing projective point, and
  the degree-3 space by a mod-p resultant certificate (`p=10^9+9`, `p=1 mod
  3`) with interpolated resultants of degree 18 and 27 whose gcd is constant.
  The omitted line `c0=0` is separately parameterized and has constant moment
  gcds, with both endpoints nonzero.  Positive controls confirm that all gcd
  routines detect a shared root when one exists.  Normal and `-O` replay are
  byte-identical.  Independent immutable audit is pending.
source: death-star-fc-archimedes-port-2026-08-03
depends_on: []
related:
  - THM-3290-archimedes-flatness-and-the-gmc3-gvc3-counterexample-family
  - THM-3018-factorial-conjecture-as-a-simplex-moment-problem
  - THM-2022-gmc2-frobenius-lowest-balanced-face
script: 04-computation/factorial_gaussian_torus_bridge_thm3300.py
output: 05-knowledge/results/factorial_gaussian_torus_bridge_thm3300.out
script_sha256: f0384ba08f857491d28fba27c048bb9bb8e7da65277014bbbd66abbd91c13f2c
output_sha256: 708e2493a7d4a941a12d8bcbba7c3e71490672f5f3cda179eee5b73e74e8d31d
hash_basis: LF-normalized bytes
---

# THM-3300 -- factorial-Gaussian torus bridge and the Archimedes no-go

**PROVED + VERIFIED-EXACT + FINITE-EXACT; INDEPENDENT IMMUTABLE AUDIT
PENDING.**

THM-3290 produced a Gaussian-moment counterexample by a spherical-average
collapse.  The natural question is whether the same device attacks the
Factorial Conjecture, whose functional THM-3018 reduces to a simplex moment
problem.  This file answers it: the two lanes are the *same* lane, and the
answer is a structural no, with the exact reason.

## 1. The bridge

Let `L_fac : C[x_1,...,x_n] -> C` be the factorial functional,
`L_fac(x^alpha)=alpha!`.  Let `z_1,...,z_n` be independent standard complex
Gaussians normalized by `E|z_i|^2=1`, i.e. `z_i=u_i+i v_i` with all `u_i,v_i`
independent real `N(0,1/2)`.

**Theorem 1.**  For every multi-index `alpha`,

```text
L_fac(x^alpha) = alpha! = E[ prod_i |z_i|^(2 alpha_i) ].            (1)
```

Equivalently, writing `q(z)=(|z_1|^2,...,|z_n|^2)`,

```text
L_fac(f) = L_Gauss^(2n)( f o q )   for every f in C[x_1,...,x_n],   (2)
```

where `L_Gauss^(2n)(g)=(exp(Delta/4)g)(0)` is the Gaussian moment functional on
`R^(2n)` with covariance `1/2` per real coordinate.

*Proof.*  `|z_i|^2` is `Exp(1)`-distributed, because the sum of squares of two
independent `N(0,1/2)` variables is exponential with mean 1; so
`E|z_i|^(2k)=k!` and independence gives `(1)`.  Both sides of `(2)` are linear,
so `(1)` suffices.  QED

The covariance is forced, not chosen: the companion checks that covariance `1`
gives `E|z|^2=2 != 1` and breaks `(1)`.

**Theorem 2 (moment map).**  For `f` homogeneous,

```text
< f o q >_{S^(2n-1)}  =  < f >_{Delta_(n-1)},                       (3)
```

both averages against the uniform probability measure.  Concretely both sides
equal `(n-1)! alpha! / (n-1+|alpha|)!` on `f=x^alpha`.

So `Delta_(n-1) = S^(2n-1)/U(1)^n` as measure spaces, and:

```text
HFC(n)  is exactly the U(1)^n-invariant part of the Gaussian
        moment problem in 2n real variables.                        (4)
```

That is the first payload.  The Factorial lane and the GMC/NC2 lane are not
neighbours; the first is a *subalgebra* of the second.

## 2. Why the Archimedes mechanism cannot be transported

THM-3290's engine has two parts, and the invariance kills both.

**(a) The annihilator degenerates to the identity.**  THM-3290 works on the
sphere by expanding a polynomial as a Laurent polynomial in an angular variable
and taking its constant term -- that is, projecting onto torus weight zero.  By
`(4)` the factorial test functions are *already* weight zero.  The projection is
the identity on them and annihilates nothing; by `(3)` the spherical average
simply is the simplex average, with no intermediate structure to exploit.

**(b) The coefficient extraction requires a pole.**  The step that did the work
in THM-3290 was `[w^(k-delta)]`, available only because the restricted
polynomial had the shape `A(1-t^2A^(2nu))^2 / z^2` -- a genuine negative power
of the angular variable.  A `U(1)^n`-invariant polynomial has no angular
variable at all, hence no pole, hence no coefficient functional that kills low
degrees.  Note also that THM-3290's ambient real dimension is `3`, which is
odd; it is not of the form `2n`, so it is not in the image of `(2)` even before
invariance is considered.

## 3. What replaces it, and the exact reason the replacement is weaker

The general principle behind both mechanisms is:

> if a compact group `G` preserves the measure and `g` is a `G`-eigenvector
> with character `chi`, then `<g^m> = chi^m <g^m>`, so `<g^m>=0` whenever
> `chi^m != 1`.

**Theorem 3 (torsion dichotomy).**  The set of exponents this argument fails to
kill is `{m >= 1 : chi^m = 1}`, the order subgroup of `chi`.

- For `G=U(1)` the character group is `Z`, which is **torsion free**: a nonzero
  weight has infinite order, so *every* `m>=1` is killed.  This is the Gaussian
  side, and it is why GMC admits counterexamples of this type.
- The affine automorphism group of a simplex is exactly the vertex permutation
  group `S_n`, which is **finite** of exponent `e`.  Every character has order
  dividing `e`, so the surviving set `e Z` is always nonempty.  No symmetry
  argument on `Delta_(n-1)` can kill all `m`.

For `n=3` and the cyclic character of order 3 this is completely explicit:

**Theorem 4.**  With `omega=e^(2 pi i/3)` and `g = s_1 + omega s_2 + omega^2
s_3` on `Delta_(n-1)` (the degree-one cyclic eigenvector),

```text
< g^m > = (n-1)! m! / (m+n-1)! * [ n divides m ].                   (5)
```

*Proof.*  Multinomially, `<g^m> = ((n-1)! m!/(m+n-1)!) * sum_{|a|=m}
omega^(sum_j j a_j)`, and the inner sum is
`[t^m] prod_(j=0)^(n-1) (1-omega^j t)^(-1) = [t^m] (1-t^n)^(-1)`.  QED

So the cyclic route reduces the obligation "`<g^m>=0` for all `m>=1`" to
"`<g^(3k)>=0` for all `k>=1`" -- a real reduction, but never a solution.

## 4. The surviving cyclic families are empty through degree three

**FINITE-EXACT.**  Let `g` range over the `omega`-eigenspace of homogeneous
forms of degree `d` in the barycentric coordinates of `Delta_2`.

```text
d=1: eigenspace is a line;   <g^3> = 1/10 != 0.                     (6)
d=2: eigenspace is a plane;  gcd_a( <g^3>, <g^6> ) is constant.     (7)
d=3: eigenspace has dim 3;   gcd_a( Res_b(<g^3>,<g^6>),
                                    Res_b(<g^3>,<g^9>) ) is constant. (8)
```

Hence no member of these families satisfies all of the displayed surviving
conditions: degrees one and two already fail within the first two, while the
degree-three certificate uses the first three.  `(7)` is computed in exact
`Q(omega)` on the chart `g=B0+a B1`; its missing point `[0:1]` has
`(<B1^3>,<B1^6>)=(1/11340,1/43783740)`, so the projective line is complete.
For `(8)`, the chart `g=B0+a B1+b B2` uses a mod-`p` certificate at
`p=10^9+9` (`p=1 mod 3`, so `omega` exists in `F_p`), with the two resultants
recovered by interpolation in degrees 18 and 27.  On its missing projective
line `c0=0`, write `g=B1+t B2`; the gcds of `M3` with `M6` and with `M9` are
both constant modulo `p`, and both projective endpoints are nonzero already
for `M3` and `M6`.  Thus both projective charts are empty.  All denominators
are prime to `p`, and the input/resultant degrees are retained modulo `p`.
After primitive denominator clearing, a characteristic-zero common factor
would therefore retain positive degree modulo `p`; the constant modular gcd
is a sound exclusion certificate.  Every gcd routine carries a positive
control confirming it reports a shared root when one exists.

## 5. Scope

Proved: Theorems 1--4 and identity `(5)`.  Finite-exact: `(6)`--`(8)`.

**Not proved.**  Nothing here proves or refutes `FC(n)`, `HFC(n)`, `SFC`, or
`JC`.  The no-go in section 2 is about *one* mechanism -- THM-3290's -- and does
not exclude other routes to an HFC counterexample; THM-3018's own Laplace
closure remains AUDIT-REQUIRED, so HFC is open on both sides.  The exclusions in
section 4 cover cyclic eigenvectors of degree at most three on `Delta_2` only,
not general polynomials and not `n>3`.  `(4)` is an identification of classes,
not an implication between the conjectures: GMC and FC have different
conclusions (Mathieu-subspace versus `f=0`), so neither is formally derived from
the other here.

**What it changes for the Factorial lane.**  Any HFC(n) counterexample is a
`U(1)^n`-invariant Gaussian-moment counterexample in `2n` real variables.  That
is a usable search criterion: it says to look for torus-invariant members inside
the known GMC counterexample families, and it explains in advance why the known
ones -- built from nonzero-weight vectors like `rho + z^2` -- cannot be among
them.

Run

```text
python 04-computation/factorial_gaussian_torus_bridge_thm3300.py
python -O 04-computation/factorial_gaussian_torus_bridge_thm3300.py
```

and compare LF-normalized bytes with the declared output.  Exact rational,
cyclotomic and modular arithmetic only; no floating point, random sampling,
imported executable, or assertion-sensitive test.

**QED** for sections 1--3; section 4 is finite-exact in its stated universe.
