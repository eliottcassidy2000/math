---
id: THM-3153
title: "Four-step prime resonance second Euclidean-Newton separation"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  For every odd prime p, no exact quadratic has three zero factorial moments
  beginning at r=p+2.  Two full linear Euclidean cancellations separate all
  positive Newton slopes.  THM-3148's fixed offset-four resultant excludes
  the remaining unit roots, with exact lawful treatment of every arithmetic
  wall and the genuine p=4547 residual exception.
audit: >
  An independent hostile audit rederived the resonance normalization, both
  quotient signs, the offset-four residual resultant, every Lucas/factorial
  valuation tier, the midpoint transfer, both high identities, generic and
  raised-endpoint polygons, the lawful A/R rather than spurious A/S unit-face
  test, the p=4547 prime-window fallback, and the five-composite consequence.
  Fresh normal and optimized replays match the stored transcript and declared
  LF-normalized hashes exactly.
source: root/multiscale-newton-flag/2026-08-02
depends_on:
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
  - THM-3143-two-step-prime-resonance-euclidean-newton-separation
  - THM-3148-fixed-offset-frobenius-endpoint-resultant-classification
related:
  - THM-3131-prime-resonance-newton-slope-separation
  - THM-3138-left-factorial-resonance-newton-slope-separation
  - THM-3146-three-step-prime-resonance-full-euclidean-newton-separation
script: 04-computation/factorial_four_step_prime_resonance_second_euclidean_thm3153.py
output: 05-knowledge/results/factorial_four_step_prime_resonance_second_euclidean_thm3153.out
script_sha256: d98a87d6bf56b45646a0342e41e5fc5ecbd49e753f6a2efe3f9f3db2bbc342c6
output_sha256: 5d5d3ea66e84f3238988fbcfad917e9f8584d754f7737901767f3e1f3a3300af
hash_basis: LF-normalized bytes
---

# THM-3153 -- four-step prime resonance second Euclidean-Newton separation

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

Let

```text
L(t^j)=j!,                    q(t)=a+bt+ct^2,                (1)
```

with `abc!=0`.  For every odd prime `p`, put

```text
r=p+2,                        d=r+2=p+4.                     (2)
```

Then

```text
L(q^r), L(q^(r+1)), L(q^(r+2))                              (3)
```

cannot all vanish.

The new feature is that one full Euclidean quotient still leaves the same
positive Newton slope.  A second full quotient breaks that slope.  Its only
remaining overlap is at height zero, where THM-3148 supplies the exact fixed
resultant.

## 1. Resonant pair and its first polygon

By THM-3124, `(3)` forces `b/a=-1/d`.  Divide by `a`, put `v=d(c/a)`, and
define

```text
A_n(v)=L((d-t+v t^2)^n).                                    (4)
```

It is enough to prove coprimality of

```text
A=A_(p+2),                         B=A_(p+3).                (5)
```

The coefficient formula is

```text
[v^j]A_n
=binom(n,j) sum_(ell=0)^(n-j)
 binom(n-j,ell)d^(n-j-ell)(-1)^ell(2j+ell)!.                (6)
```

Write `p=2m+1`.  Lucas reduction and the factorial valuation in `(6)` give

```text
v_p(A_j)=0                           for j=0,1,2,
v_p(A_j)>=1                          for 3<=j<=m,
v_p(A_j)>=2                          for m+1<=j<=p+2.        (7)
```

The last coefficient is `(2p+4)!`, of valuation two.  Hence, for `p>=7`,

```text
NP_p(A):(0,0)--(2,0)--(p+2,2),       slopes 0, 2/p.         (8)
```

The small primes will be discharged in Section 6.

## 2. First full Euclidean quotient

The linear quotient of `B` by `A` is

```text
C v+q,
C=(2p+6)(2p+5),
q=-2(p+3)(p+4)/(2p+3).                                     (9)
```

Clear the `p`-adic-unit denominator and define

```text
R=(2p+3)B-
 [(2p+3)C v-2(p+3)(p+4)]A.                                (10)
```

Then `deg R<=p+1`.  Its leading coefficients are

```text
R_(p+1)=-2(p+1)(p+2)(p+3)(24p+61)(2p)!,                   (11)

R_p=2p(p+1)(p+2)(p+3)
    (24p^2-67p-337)(2p-2)!.                                (12)
```

Modulo `p`, THM-3148's offset-four profiles give

```text
A == 96v^2+16v+40,
B == 2880v^3-288v^2+96v+136,
R == 1368-2928v=24(57-122v).                       (mod p) (13)
```

Thus the first quotient lowers the height-zero degree from three to one,
but for generic `p` it still has polygon

```text
(0,0)--(1,0)--(p+1,2),                  slopes 0,2/p.       (14)
```

The positive slope `2/p` still coincides with `(8)`.  This is why stopping
after the first leading cancellation cannot prove the theorem.

## 3. Second full Euclidean quotient

The quotient of `A` by `R` is `q_1 v+q_0`, where

```text
q_1=-2(2p+1)(2p+3)/[(p+3)(24p+61)],

q_0=4(p+4)(2p+1)(28p+67)
    /[(p+3)(2p-1)(24p+61)^2].                              (15)
```

Put

```text
D =(p+3)(2p-1)(24p+61)^2,
N_1=-2(2p+1)(2p+3)(2p-1)(24p+61),
N_0=4(p+4)(2p+1)(28p+67),                                 (16)

S=D A-(N_1 v+N_0)R.                                       (17)
```

The two top coefficients vanish exactly, so `deg S<=p`.  Every common
factor of `A,B` divides `R` and `S`; clearing `D` loses no common factor.

The low residual is the fixed linear polynomial

```text
S ==168(14640v-11387)                              (mod p). (18)
```

Its constant and linear coordinates are both units except at

```text
p in {7,59,61,193}.                                        (19)
```

Set

```text
k=m+2=(p+3)/2.                                             (20)
```

The valuation tiers following from `(6),(10),(17)` are

```text
v_p(R_j)>=1 for j>=2,        v_p(R_j)>=2 for j>=m+2,
v_p(S_j)>=1 for j>=2,        v_p(S_j)>=2 for j>=m+3.        (21)
```

At the missing transition index, only the shifted `R_(m+1)` term survives
after division by `p`.  Exact simplification gives

```text
S_k/p ==
 65880(-1)^(m-2)4^(m+3)/[m(m-1)(m-2)]              (mod p). (22)
```

This is nonzero for every `p>=7` except `p=61`.

The two high coefficients are

```text
S_p=4p(p-1)(p+1)(p+2)(p+3)(2p+3)
    E_4(p)(2p-4)!,                                         (23)

E_4(p)=256p^4-29696p^3-193568p^2-406368p-276169,

S_(p-1)=-4p(p-1)(p-2)(p+1)(p+2)(p+3)(2p+3)
         E_5(p)(2p-6)!,                                    (24)

E_5(p)=256p^5-32000p^4+97504p^3+1530816p^2
       +3882095p+2891889.
```

Since

```text
276169=277*997,                                             (25)
```

the endpoint in `(23)` has valuation two except at `p=277,997`, where it
has valuation three.  At those two primes, `(24)` has valuation exactly two.
The factors `7,11,13,107` of `2891889` only raise a nonvertex in the generic
polygon.

## 4. Generic and raised-endpoint polygons

For every `p>=7` outside `(19)` and `{277,997}`, the low units `(18)`, the
midpoint unit `(22)`, the bounds `(21)`, and endpoint `(23)` force

```text
NP_p(S):(0,0)--(1,0)--(k,1)--(p,2),                        (26)

slopes       0,       2/(p+1), 2/(p-3).                    (27)
```

At `p=277` and `p=997`, `(24)` replaces the raised endpoint segment:

```text
NP_p(S):(0,0)--(1,0)--(k,1)--(p-1,2)--(p,3),               (28)

slopes       0,       2/(p+1), 2/(p-5), 1.                 (29)
```

Every positive slope in `(27)` and `(29)` is different from `2/p`, the
only positive slope of `A`.

The large `p=997` polygon is certified without giant integers: exact
arithmetic modulo `p^5` determines every valuation below five, while all
displayed vertices have height at most three.  This is a rigorous truncated
valuation computation, not floating point or a fitted pattern.

## 5. Height-zero roots and the fixed resultant

The only generic slope shared by `(8)` and `(26)` or `(28)` is zero.  A
common unit root of `A,B` would reduce to a common nonzero root of the first
two residuals in `(13)`.  Their exact resultant is

```text
Res_v(96v^2+16v+40,1368-2928v)
=586672128
=2^11*3^2*7*4547.                                          (30)
```

This is precisely `delta_4` from THM-3148.  Under `p>6`, the only residual
exceptions are therefore

```text
p=7,4547.                                                   (31)
```

Thus every prime governed by `(26)` or `(28)`, except possibly `4547`, has
neither a common positive-slope root nor a common unit root.  Hence `A` and
`B` are coprime there.

Notice that the extra factor `18089` which can appear in
`Res(A mod p,S mod p)` is irrelevant: the second quotient can introduce a
spurious common residual root through its linear multiplier.  A genuine
common factor must also divide `R`, so the lawful residual test is the
three-row condition represented by `(30)`.  This is the first failed
implication of a naive pairwise-second-remainder argument.

## 6. Arithmetic walls and completion

The primes

```text
p=3,5,7,59,61,193                                           (32)
```

all have `r=p+2<=200`, so THM-3124 closes their resonant windows exactly.
The remaining genuine low-face exception is

```text
p=4547,                      r=p+2=4549.                    (33)
```

The integer `4549` is prime.  Therefore THM-3143, applied with prime
`r=4549` and `d=r+2=4551`, closes exactly this same window.

All odd primes are now covered, proving `(3)`.

Combining this theorem with THM-3131, THM-3138, THM-3143, THM-3146, and the
finite range of THM-3124, any still-open bad exact-quadratic resonance must
satisfy

```text
r>=201,
d=r+2, d-1, d-2, d-3, and d-4 are all composite.            (34)
```

This is a five-composite residual, one step stronger than THM-3146.

## 7. Exact companion

Run

```text
python 04-computation/factorial_four_step_prime_resonance_second_euclidean_thm3153.py
python -O 04-computation/factorial_four_step_prime_resonance_second_euclidean_thm3153.py
```

and compare byte-for-byte with

```text
05-knowledge/results/factorial_four_step_prime_resonance_second_euclidean_thm3153.out.
```

The companion independently expands the original bivariate moments at
`p=7`; checks both exact quotients, `(22)--(24)`, and full polygons at

```text
p=7,11,13,59,61,107,193,277;                               (35)
```

checks `(28)` at `p=997` modulo `p^5`; and scans all `1,228` odd primes below
`10,000` for every closed exceptional constant.  The exact controls include
the three low-coordinate walls, both raised endpoints, all raised
penultimate nonvertices, the resultant exceptions, and the primality of
`4549`.

## 8. Scope

The source is THM-3124's resonant exact `{0,1,2}` quadratic pair.  The maps
are the two exact Euclidean rows `(10),(17)`.  They preserve common factors
but not a distinguished root, phase, or coefficient sign.  THM-3148 supplies
only the zero-slope sidecar; `(21)--(29)` are independently necessary to
dispose of positive slopes.

No arbitrary fixed offset, arbitrary support, arbitrary radial polynomial,
all of `NC2`, or the Gaussian Moment Conjecture in two dimensions is proved
here.  The pattern suggests an iterated Euclidean-Newton flag indexed by the
prime gap, but denominator and endpoint walls proliferate and no induction
is asserted.

**QED.**
