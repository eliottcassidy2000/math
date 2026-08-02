---
id: THM-3131
title: "Prime-resonance Newton-slope separation for quadratic factorial windows"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.  For an exact
  quadratic q(t)=a+bt+ct^2, three consecutive zero factorial moments beginning
  at r force the THM-3124 resonance -a/b=r+2.  If r+2 is prime, the two
  remaining resonant moment polynomials have disjoint p-adic Newton slopes and
  hence no common root.  Thus every possible bad exact-quadratic window has
  composite r+2.  The argument is uniform in the prime and is not a finite
  census.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
script: 04-computation/factorial_prime_resonance_newton_slopes_thm3131.py
output: 05-knowledge/results/factorial_prime_resonance_newton_slopes_thm3131.out
script_sha256: 793cbe9510f8cf20b12faf30f2e4bcba3b2a6c4522a178c2c94c0162b85a2c1c
output_sha256: 1ad97b44f8ca402022d3b885fcd1860766396f7075005d62fb3153d2925e4ffd
hash_basis: LF-normalized bytes
---

# THM-3131 -- prime-resonance Newton-slope separation for quadratic factorial windows

**PROVED CANDIDATE + VERIFIED-EXACT; INDEPENDENT AUDIT PENDING.**

Let

```text
L(t^k)=k!,                    q(t)=a+bt+ct^2,                 (1)
```

with `abc!=0`.  If, for some `r>=1`,

```text
L(q^r)=L(q^(r+1))=L(q^(r+2))=0,                              (2)
```

then `r+2` is composite.  Equivalently, every prime-resonant shifted
three-moment window contains a nonzero moment.

This gives an unbounded family beyond THM-3124's finite-exact range: whenever
`p=r+2` is prime, the start `r=p-2` is closed by one uniform argument.

## 1. Reduction to an integral polynomial pair

THM-3124 proves that `(2)` forces

```text
b/a=-1/(r+2).                                                  (3)
```

After dividing `q` by `a` and writing `u=c/a`, the only remaining question is
whether

```text
P_r(u)=L((1-t/(r+2)+u t^2)^r),
P_(r+1)(u)=L((1-t/(r+2)+u t^2)^(r+1))                         (4)
```

have a common complex root.

Assume `p=r+2` is prime and make the invertible change `v=pu`.  Put

```text
A_n(v)=p^n P_n(v/p)=L((p-t+v t^2)^n) in Z[v].                 (5)
```

Thus `(4)` has a common root if and only if `A_(p-2)` and `A_(p-1)` do.
Writing

```text
A_n(v)=sum_(j=0)^n a_(n,j)v^j,                               (6)
```

direct expansion gives

```text
a_(n,j)=binom(n,j) sum_(ell=0)^(n-j)
  binom(n-j,ell)p^(n-j-ell)(-1)^ell(2j+ell)!.                 (7)
```

Modulo `p`, only the terminal summand `ell=n-j` remains:

```text
a_(n,j) == binom(n,j)(-1)^(n-j)(n+j)!             (mod p).   (8)
```

## 2. The two exact Newton polygons

Use coefficient index on the horizontal axis and `p`-adic valuation on the
vertical axis.

For `n=p-2`, `(8)` says that `a_(p-2,0)` and `a_(p-2,1)` are units, while all
coefficients with `j>=2` are divisible by `p`.  The leading coefficient is

```text
a_(p-2,p-2)=(2p-4)!.                                         (9)
```

For `p>=5`, its valuation is exactly one.  Every intervening nonzero
coefficient has integral valuation at least one, hence lies strictly above
the segment from `(1,0)` to `(p-2,1)`.  Therefore

```text
NP_p(A_(p-2)) : (0,0)--(1,0)--(p-2,1),
slopes            0             1/(p-3).                     (10)
```

For the special prime `p=3`, `A_(p-2)=A_1` has both coefficients units, so
its sole slope is `0`.

For `n=p-1`, `(8)` makes the constant coefficient a unit and every
positive-degree coefficient divisible by `p`.  Its leading coefficient is

```text
a_(p-1,p-1)=(2p-2)!,                                         (11)
```

again of valuation exactly one.  All intermediate points lie strictly above
the endpoint segment, so for every odd prime `p`

```text
NP_p(A_(p-1)) : (0,0)--(p-1,1),
sole slope                 1/(p-1).                           (12)
```

The slope sets in `(10)` and `(12)` are disjoint.  This is also true at
`p=3`, where they are `{0}` and `{1/2}`.

## 3. Why slope separation excludes a common factor

For completeness, the Newton product lemma says that over a valued field the
multiset of lower Newton slopes of a product is the union of the slope
multisets of its factors.  One quick proof uses the weighted Gauss valuation

```text
w_s(sum_j c_j v^j)=min_j(v_p(c_j)+sj),                        (13)
```

which is additive under multiplication; its breakpoints record the lower
polygon.  In particular, after splitting over an algebraic closure of
`Q_p`, a linear factor `v-alpha` contributes slope `-v_p(alpha)`.

If `A_(p-2)` and `A_(p-1)` had a nonconstant common factor over `Q`, the same
factor would remain nonconstant over `Q_p` and would contribute at least one
Newton slope to both polynomials.  Equations `(10)`--`(12)` make this
impossible.  The pair is therefore coprime over `Q`, so it has no common
complex root.  This contradicts `(4)` and proves the claim.  QED.

## 4. Connection contract and hostile boundary

The source object is THM-3124's resonant common-root pair `(4)`.  The map is
the invertible coefficient rescaling `u -> v=pu` followed by the `p`-adic
valuation polygon.  It preserves the existence of a common algebraic factor;
it destroys residue and angular information about individual roots.  The
needed sidecars are precisely primality of `p`, the two unit coefficients,
and the exact leading valuations `(9)` and `(11)`.

Primality is not cosmetic.  At composite `d=r+2`, the same construction can
produce overlapping slopes:

```text
d=4, p=2:  slopes(A_2)={1},       slopes(A_3)={1};
d=6, p=3:  slopes(A_4)={0,1/3},   slopes(A_5)={1/3,1}.        (14)
```

These are hostile controls on the observer, not counterexamples to the
factorial statement: THM-3124 already excludes these small windows by exact
gcd.  For composite resonances the slope map has lost too much information,
so residual polynomials or a different place are still needed.

## 5. Scope

Combining this theorem with THM-3124, any still-open bad exact-quadratic
window must satisfy

```text
r>=201,             r+2 composite,             b/a=-1/(r+2). (15)
```

No claim is made here for translated three-slot supports, arbitrary
univariate polynomials, ambient `SFC(3)`, or full `FC(3)`.

The companion is an exact verifier, not the proof source.  It independently
constructs small bivariate powers, checks `(7)` and `(8)`, verifies the stated
polygons for every odd prime through `223` (including new starts above `200`),
and records both composite hostile controls.  Run

```text
python3 04-computation/factorial_prime_resonance_newton_slopes_thm3131.py
python3 -O 04-computation/factorial_prime_resonance_newton_slopes_thm3131.py
```

and compare byte-for-byte with the declared output.

**End of proof candidate.**
