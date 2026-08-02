---
id: THM-3142
title: "Prime-power-predecessor Newton irreducibility and composite-window census"
status: >
  PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  If the resonance predecessor d-1 is a prime power p^k, the degree-(d-1)
  resonant factorial polynomial has one primitive p-adic Newton slope and is
  irreducible over Q_p; it cannot share a factor with its degree-(d-2) mate.
  By contrast, for every even d and every n, the divisor-place 2-adic polygon
  of L((d-t+vt^2)^n) has sole slope one, proving this observer is uniformly
  blind.  Together with THM-3124 and THM-3131, an exact two-prime gcd census
  closes every exact-quadratic resonant window through d=1000 (r=998).  Any residual window
  has r>=999, d=r+2 composite, and d-1 not a prime power.  The wider place
  patterns seen through d<=220 are finite-exact observations only.
audit: >
  An independent hostile derivation rechecked the prime-power binomial
  valuation, odd-prime and binary digit-sum bounds, primitive slope
  denominator, Q_p factor-degree implication, uniform even-d slope-one
  theorem, recurrence/direct construction agreement, and degree-preserving
  modular certificate.  It replayed all 676 composite d=203..1000 at the two
  declared primes and an independent third prime, with planted common-factor,
  lost-leading-degree, and nonprimitive-slope controls.  Canonical normal,
  optimized, and stored transcripts match exactly.
source: root/frontier-synthesis-2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3131-prime-resonance-newton-slope-separation
related:
  - THM-3138-left-factorial-resonance-newton-slope-separation
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
script: 04-computation/factorial_prime_power_predecessor_newton_census_thm3142.py
output: 05-knowledge/results/factorial_prime_power_predecessor_newton_census_thm3142.out
script_sha256: 656b2a98dc275112845728cede18de44c4a506faa64196f14ce193a579ee57d3
output_sha256: 967217609727bbaae83272521b7c49f595cc096f64768c2a323608fb93b719cd
hash_basis: LF-normalized bytes
---

# THM-3142 -- prime-power-predecessor Newton irreducibility and composite-window census

**PROVED + FINITE-EXACT + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^m)=m!,
A_n^(d)(v)=L((d-t+vt^2)^n) in Z[v].                         (1)
```

THM-3124 shows that three zero moments of an exact quadratic beginning at
`r` force `d=r+2` and reduce to a common root of

```text
A_(d-2)^(d),                   A_(d-1)^(d).                  (2)
```

This theorem supplies a uniform adjacent-prime-power exit, proves a global
blindness boundary for the most tempting divisor place, and extends the
finite exact closure.

## 1. Prime-power-predecessor theorem

> If `N=d-1=p^k` for a prime `p` and `k>=1`, then
> `A_N^(N+1)` is irreducible over `Q_p`.  In particular,
> `A_(N-1)^(N+1)` and `A_N^(N+1)` are coprime over `Q`.

Write

```text
B(v)=A_N^(N+1)(v)=sum_(j=0)^N b_jv^j.                       (3)
```

Direct expansion gives

```text
b_j=binom(N,j) sum_(ell=0)^(N-j)
    binom(N-j,ell)(N+1)^(N-j-ell)(-1)^ell(2j+ell)!,          (4)
```

and hence

```text
b_j=binom(N,j)(2j)! Z_j,                       Z_j in Z.      (5)
```

For `0<j<N`, the standard prime-power binomial valuation is

```text
v_p(binomial(N,j))=k-v_p(j).                                (6)
```

Consequently

```text
v_p(b_j)>=k-v_p(j)+v_p((2j)!).                              (7)
```

The constant coefficient is a unit modulo `p`: in `(4)` at `j=0`, the
`ell=0` term is `(N+1)^N==1`; every `0<ell<N` term contains
`binomial(N,ell)`, and the `ell=N` term contains `N!`.  The leading
coefficient and its valuation are

```text
b_N=(2N)!,                       e=v_p((2N)!).               (8)
```

### Odd `p`

Put `h=v_p(j)` and `K=k-h`.  Since `j=p^h u` with `1<=u<p^K`,

```text
s_p(2j)=s_p(2u)<=(p-1)K.                                   (9)
```

Legendre's formula and `(7)` give

```text
v_p(b_j)>=K+(2j-s_p(2j))/(p-1)>=2j/(p-1).                  (10)
```

But

```text
e=2(N-1)/(p-1),
ej/N=2j/(p-1)-2j/(N(p-1)).                                 (11)
```

Thus every internal coefficient lies strictly above the endpoint segment,
by at least `2j/(N(p-1))`.

### `p=2`

Here `s_2(2j)=s_2(j)` and `s_2(j)<=K`.  Equations `(7)` and Legendre's
formula give

```text
v_2(b_j)>=K+2j-s_2(j)>=2j.                                 (12)
```

Since `e=2N-1`, the endpoint line has height

```text
ej/N=2j-j/N,                                                (13)
```

so every internal point is again strictly above it.

In both cases the lower Newton polygon is the single segment

```text
NP_p(B):(0,0)--(N,e).                                       (14)
```

Moreover `p` does not divide `e`: for odd `p`, `e==2 (mod p)`, while for
`p=2`, `e` is odd.  Hence `gcd(N,e)=1`, and the slope `e/N` has reduced
denominator `N`.

Every irreducible `Q_p` factor of `B` must therefore have degree divisible
by `N`: its roots all have valuation denominator `N`, and the ramification
degree divides the factor degree.  Since `deg(B)=N`, `B` is irreducible over
`Q_p`.  Its degree-`N-1` mate in `(2)` cannot share a nonconstant factor.
Base change then proves coprimality over `Q`.  QED.

This strictly extends the exclusion conclusion of THM-3138 from `k=1` and
odd `p` to every prime power, including powers of two.  THM-3138 retains the
finer left-factorial description of the first polynomial's Newton chamber.

## 2. A uniform blind divisor place

The adjacent place succeeds because the obvious divisor place can be
intrinsically blind.  Let `d` be even and put `e_d=v_2(d)>=1`.  Expand a
term contributing to the coefficient of `v^j` in `(1)` by choosing `j`
quadratic terms, `ell` linear terms, and `k` constant terms, where

```text
j+ell+k=n.                                                   (15)
```

Relative to `v_2(n!)`, its exact valuation excess is

```text
j+s_2(j)+s_2(ell)+s_2(k)-s_2(2j+ell)+k(e_d-1).              (16)
```

Binary digit subadditivity gives

```text
s_2(2j+ell)<=s_2(j)+s_2(ell),                               (17)
```

so `(16)` is at least `j`.  At `j=0`, the choice `k=0`, `ell=n` is the
unique term of minimal valuation; hence the constant coefficient has exact
valuation `v_2(n!)`.  At `j=n`, the coefficient is `(2n)!`, whose valuation
is exactly `v_2(n!)+n`.  Every intervening point lies on or above the joining
line.  Therefore, for every `n>=1`,

```text
NP_2(A_n^(d)) has sole slope 1.                              (18)
```

In particular, the `p=2` divisor-place polygons of the two polynomials in
`(2)` always overlap completely when `d` is even.  Any shared rescaling only
translates both polygons and cannot separate them.  This is a proved
observer no-go, not evidence for a common algebraic factor.

## 3. Finite-exact place scouts

The canonical companion makes the following bounded observations.  They
motivate the adjacent-place theorem but are not extrapolated beyond their
declared universes.

- Among all `172` composite `4<=d<=220`, all `365` prime divisor places
  `p|d` have at least one shared Newton slope.
- At the `303` predecessor places `p|d-1` in the same universe, separation
  occurs in exactly `58` cases, precisely when `d-1` is a power of that
  `p`.
- Testing every prime `p<=2d` for every composite `d<=160` gives `4,643`
  place tests and `47` separations, again exactly the prime-power
  predecessors in this finite universe.
- Exact rational gcds are one for
  `d=4,6,8,9,10,12,15,16,18,20`, despite divisor-place slope overlap.

The cases `d=4`, `d=6`, and `d=132` are retained as explicit shared-slope
hostiles.  They show why slope overlap is observer loss, not a common-factor
certificate.

## 4. Finite-exact closure through `d=1000`

THM-3124 already excludes every start `1<=r<=200`, equivalently every
resonance `3<=d<=202`.  For every composite

```text
203<=d<=1000,                                                (19)
```

the companion constructs the pair `(2)` by the division-free THM-3124
recurrence and computes its gcd over both

```text
F_1000003,                         F_1000033.                 (20)
```

There are exactly `676` such composite `d`, and every computed gcd has
degree zero at both primes.  Both moduli are prime, exceed `2(d-1)`, and do
not divide `d`; therefore the integral rescaling and both leading factorial
coefficients retain their degrees.  If the pair over `Q` had a nonconstant
common factor, primitive factorization and reduction would leave a
positive-degree common factor at either prime, contradicting the exact gcd.

The recurrence agrees with the direct integer formula for `4<=d<=15`.
The independent audit replays the entire census at both declared primes and
also at `1000037`; it detects a planted preserved-degree common factor and
demonstrates the failure of an otherwise identical certificate when the
leading degree is lost.

Combining THM-3131's prime exit, the prime-power-predecessor exit, and this
census, any
still-open bad exact-quadratic window must satisfy

```text
r>=999,
d=r+2 is composite,
d-1 is not a prime power,
b/a=-1/d.                                                    (21)
```

The finite statement is `r<=998`, not `r<=1000`.

THM-3143 subsequently excludes the additional lane `d-2` prime by resolving
its shared Frobenius face with one exact Euclidean cancellation.  It further
forces `d-2` composite but is not a dependency of the theorem proved here.

## 5. Reproduction and audit

Run

```text
python3 04-computation/factorial_prime_power_predecessor_newton_census_thm3142.py
python3 -O 04-computation/factorial_prime_power_predecessor_newton_census_thm3142.py
```

and compare byte-for-byte with the declared output.  The independent audit
does not import this companion.  Its script and byte-identical
normal/optimized output hashes are

```text
3aeb3cadcd0494ee00d3875e69fe9c9a9373e99355f66487fda509c4aadbd564
d7c898c078a9b25f01e3b205bedf4f8abc56e2bebd51108c123dc0a21a55e0c4
```

## 6. Connection contract and scope

The source is THM-3124's resonant polynomial pair.  The map is the invertible
integral rescaling followed by the Newton polygon at a prime dividing
`d-1`, rather than `d`.  It preserves common factors and forgets root phases
and residual polynomials.  The load-bearing sidecars are the prime-power
binomial valuation, digit-sum bounds, and primitive endpoint slope.  The
even-`d` theorem identifies an entire observer kernel and explains why the
place must sometimes change.

These are exact `{0,1,2}`-support results inside ambient `SFC(1)`.  They do
not settle arbitrary-support `SFC(1)`, ambient `SFC(3)`, or the full
three-variable Factorial Conjecture `FC(3)`.

**End of proof.**
