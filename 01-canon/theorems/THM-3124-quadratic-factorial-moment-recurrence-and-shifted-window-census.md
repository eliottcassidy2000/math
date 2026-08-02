---
id: THM-3124
title: "Quadratic factorial-moment recurrence and shifted-window census"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For an exact three-slot polynomial q(t)=a+bt+ct^2, its factorial moments
  satisfy a division-free second-order recurrence.  Three consecutive zero
  moments beginning at r force the unique resonance b/a=-1/(r+2).  Hence a
  fixed exact quadratic has at most one candidate bad window, and none at all
  unless -a/b is an integer at least three.  At the resonant ratios, a
  two-prime exact polynomial-gcd census excludes every window start
  1<=r<=200.  The recurrence proof is symbolic; the bounded census is
  FINITE-EXACT.
audit: >
  An independent derivation checked every integration-by-parts boundary sign,
  the resonance and automatic third zero, and the primitive modular lifting
  argument.  A separate SymPy finite-field engine replayed all 200 starts at
  both primes with no degree loss or common factor.  Normal, optimized, and
  stored transcripts and declared hashes agree.  The O(N) recurrence cost is
  a direct consequence, not a fitted extrapolation.
source: root/frontier-synthesis-2026-08-02
depends_on: []
related:
  - THM-2173-sparse-projective-factorial-moment-floor
  - THM-2824-arbitrary-three-slot-factorial-moment-divisibility-and-atomic-orientation-boundary
  - THM-2836-sfc3-arbitrary-support-shifted-window-census
script: 04-computation/factorial_quadratic_recurrence_census_thm3124.py
output: 05-knowledge/results/factorial_quadratic_recurrence_census_thm3124.out
script_sha256: 575fd6a8ae2bbedb4a6991a66382d2772033ea666424a1b9b57244f01deb1b7c
output_sha256: 2383aa94f5ae03d5a4009fe596b3e4cc675b557020c77828e15ece45b885e05a
hash_basis: LF-normalized bytes
---

# THM-3124 -- quadratic factorial-moment recurrence and shifted-window census

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

This theorem concerns the original one-variable factorial functional

```text
L : C[t] -> C,                 L(t^j)=j!,                     (1)
```

and the exact three-slot support `{0,1,2}`.  In the notation of
Edo--van den Essen, this is a three-monomial slice of `SFC(1)`: the `1`
indexes the ambient variable dimension, while the window length `3` comes
from the number of nonzero monomials.  It is not the full ambient
three-variable conjecture `SFC(3)`.

The new mechanism is a boundary-forced recurrence.  It collapses a putative
three-moment window to one rational resonance in the coefficient ratio.  An
exact modular census then closes the first two hundred resonant windows.

## 1. Exact recurrence

Let

```text
q(t)=a+bt+ct^2,                  M_n=L(q^n),
Delta=b^2-4ac,                  M_0=1.                       (2)
```

Then for every integer `n>=1`,

```text
M_(n+1)
 =a^n(a+(n+1)b)
  +2(n+1)(2n+1)c M_n
  +n(n+1)Delta M_(n-1).                                      (3)
```

This identity is division-free and is valid for arbitrary complex
`a,b,c`, including a repeated root of `q`.

Together with `M_0=1` and `M_1=a+b+2c`, `(3)` computes the first `N`
factorial moments in `O(N)` ring operations and constant recurrence state.
This is the structural closed form for the sequence; it replaces repeated
multinomial expansion and is not inferred from a finite fit.

### Proof

The factorial functional is the exponential moment

```text
M_n=int_0^infinity q(t)^n e^(-t) dt.                          (4)
```

Put

```text
J_n=int_0^infinity q'(t)q(t)^n e^(-t) dt.                    (5)
```

Integrating the derivative of `q^(n+1)e^(-t)` and using exponential decay
at infinity gives

```text
-a^(n+1)=(n+1)J_n-M_(n+1),
M_(n+1)=a^(n+1)+(n+1)J_n.                                    (6)
```

The quadratic identities

```text
q''=2c,                    (q')^2=4cq+Delta                  (7)
```

make a second integration by parts close on two preceding moments:

```text
-b a^n
 =int_0^infinity (q' q^n e^(-t))' dt
 =2c M_n+n int_0^infinity (q')^2 q^(n-1)e^(-t)dt-J_n
 =2c(2n+1)M_n+n Delta M_(n-1)-J_n.                           (8)
```

Therefore

```text
J_n=b a^n+2c(2n+1)M_n+n Delta M_(n-1).                       (9)
```

Substitution of `(9)` into `(6)` is exactly `(3)`.  All boundary terms used
above are explicit: the first is `q(0)^(n+1)=a^(n+1)`, and the second is
`q'(0)q(0)^n=b a^n`.  This boundary forcing is why three consecutive zeros
cannot occur at a generic coefficient ratio.

## 2. Resonance reduction

Assume from now on that

```text
abc != 0,                                                     (10)
```

so `q` has exactly the three monomials on `{0,1,2}`.  Suppose a shifted
three-moment window beginning at `r>=1` is bad:

```text
M_r=M_(r+1)=M_(r+2)=0.                                       (11)
```

Apply `(3)` with `n=r+1`.  Every moment term on its right vanishes, leaving

```text
0=a^(r+1)(a+(r+2)b).                                         (12)
```

Since `a!=0`, a necessary condition is

```text
b/a=-1/(r+2),                    equivalently -a/b=r+2.       (13)
```

Consequently:

1. a fixed exact quadratic has at most one candidate bad window;
2. if `-a/b` is not an integer at least `3`, every three-moment shifted
   window contains a nonzero moment;
3. at the resonance `(13)`, the third equation in `(11)` is automatic once
   the first two hold.

Normalize by the nonzero scalar `a` and put `u=c/a`.  At the candidate start
`r`, write

```text
Q_(r,u)(t)=1-t/(r+2)+u t^2,
P_(n,r)(u)=L(Q_(r,u)^n).                                     (14)
```

Then `(11)` is equivalent to the common-root problem

```text
P_(r,r)(u)=P_(r+1,r)(u)=0,                                   (15)
```

because `(3)` supplies `P_(r+2,r)(u)=0`.  Direct multinomial expansion gives
the exact coefficient formula

```text
[u^j]P_(n,r)
 =binom(n,j) sum_(ell=0)^(n-j)
    binom(n-j,ell)(-1/(r+2))^ell (2j+ell)!,
                  0<=j<=n.                                  (16)
```

In particular `deg_u P_(n,r)=n`, with leading coefficient `(2n)!`.

## 3. Finite-exact closure through window start 200

For every `1<=r<=200`, the companion constructs the two polynomials in
`(15)` from `(16)` and runs the Euclidean algorithm over each of

```text
F_1000003,                     F_1000033.                     (17)
```

Both moduli are prime.  They exceed `2(r+1)`, do not divide `r+2`, and hence
preserve every denominator and both leading coefficients `(2r)!` and
`(2r+2)!`.  For every one of the `200` candidate starts, both reductions
have gcd `1`.

This is a rigorous one-way certificate over `C`, not merely a search.  If
the rational polynomials in `(15)` had a nonconstant common factor, clear
their powers of `r+2` and take primitive integral representatives.  Since
the leading coefficients remain nonzero modulo either prime, every factor
and cofactor retains its positive degree after reduction.  A common factor
would therefore remain common modulo that prime, contradicting the computed
gcd `1`.  One prime already proves emptiness; the second is an independent
prime replay of the same exact engine.

Thus:

```text
For every exact q=a+bt+ct^2 and every 1<=r<=200,
at least one of L(q^r), L(q^(r+1)), L(q^(r+2)) is nonzero.     (18)
```

As orthogonal controls, SymPy over `Q` returns gcd degree zero for
`1<=r<=7`; the symbolic expansion verifies `(3)` identically for
`1<=n<=6`; and the modular Euclidean routine returns the intended positive
control

```text
gcd(x^2+2x+1, x^2+3x+2)=x+1.                                (19)
```

## 4. Boundary and nonclaims

The analytic recurrence `(3)` and its resonance consequence `(13)` hold at
all orders.  The gcd closure is finite-exact only for `r<=200`.  At a
resonant start `r>=201`, the common-root problem `(15)` remains open here.

No claim is made for:

- supports other than the exact univariate support `{0,1,2}`;
- translated supports `{p,p+1,p+2}` with `p>0`, because multiplying by
  `t^p` tilts the factorial functional differently in different powers;
- missing-slot cases with one of `a,b,c` equal to zero;
- all of `SFC(1)`, ambient `SFC(3)`, or the full three-variable `FC(3)`;
- a closed form for the possible common roots at resonant starts beyond the
  certified range; or
- any transfer from the false three-dimensional Gaussian moment conjecture.

This result complements THM-2824: that theorem closes the first window for
every arbitrary three-slot support, while `(3)` controls all shifted windows
for one support up to a single coefficient resonance.  THM-2173 explains why
three equations, rather than two, are the sharp generic detection scale.

## 5. Reproduction

Run

```text
python3 04-computation/factorial_quadratic_recurrence_census_thm3124.py
python3 -O 04-computation/factorial_quadratic_recurrence_census_thm3124.py
```

Both executions must byte-match
`05-knowledge/results/factorial_quadratic_recurrence_census_thm3124.out`
after LF normalization.  The transcript records the symbolic recurrence
checks, rational gcd controls, both exact prime censuses, degree preservation,
and both positive controls.

**End of proof.**
