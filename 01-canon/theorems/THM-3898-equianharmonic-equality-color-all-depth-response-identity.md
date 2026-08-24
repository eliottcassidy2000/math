---
id: THM-3898
title: "Equianharmonic equality colors have an exact all-depth response identity"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
  independent hostile audit.  On every positive equal-y-degree seam of the
  THM-3881 residual, seven normalized macros give one exact all-depth color
  identity.  Their base arrival depths are
  0,n,n+2,2n+2,2n+4,3n+4,4n+4.  The canonical zero-lower-sidecar payment
  consequently requires the nonpolynomial response L^2/h at its first
  marked depth n; in particular its n=4 fourth response fails.  This closes
  only that canonical family.  Arbitrary sidecars, the full equality seam,
  a Keller atlas, and JC(2) remain OPEN.
source: >
  jc_zero_debt_lift / root-authorized takeover of the abandoned THM-3898
  reservation, post-THM-3905 response compression, 2026-08-23
audit: >
  SELF-AUDITED EXACT PROOF CANDIDATE.  The companion independently normalizes
  all seven residual macros, proves the universal identity with epsilon^n
  left symbolic, checks the affine arrival schedule, replays the canonical
  source at n=1,...,8, verifies the n=4 response and rational-x hostile, and
  retains the exact C3 color constants in 46 active gates.  Normal and
  optimized runs byte-match the frozen output.  Independent audit must
  recheck every normalization exponent, formal-square-root uniqueness at the
  marked depth, the gcd(h,L) obstruction, and the arbitrary-sidecar boundary.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
  - THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors
related:
  - THM-3902-nonzero-sidecar-equality-color-two-jet-response
  - THM-3905-nonzero-sidecar-equality-color-third-response
script: 04-computation/jc2_equianharmonic_equality_color_all_depth_response_thm3898.py
output: 05-knowledge/results/jc2_equianharmonic_equality_color_all_depth_response_thm3898.out
script_sha256: 1e7f8050e63b1c2116c0e6988d8f13b4c21e3245de8fcc640c961ff5db956e91
output_sha256: 0e25c7134fa4b8be82c68929bec17cfab9928e77f0b525d2b5ca9c62af0e1ab5
semantic_sha256: 0fcd414b8ab768168444169f1d318aab15889d5199178711ebfc3491ab83c122
hash_basis: raw LF bytes
---

# THM-3898 -- one identity contains every equality-color response

**RESERVED / PROVISIONAL PROOF CANDIDATE + VERIFIED-EXACT; awaiting
independent hostile audit.**  Work over an algebraically closed field `k` of
characteristic zero.  In `D=k[x,y]` retain

```text
a=x+1,                  L=9x+4,
kappa=-15x^2-15x-4,     K=y^2+kappa,
P=aL^2,                 r=aT+Kf,
A=KT+aPf,               B=Pf^2-T^2.                      (1)
```

The THM-3881 residual is

```text
S=L^4+6AL^2f+6PL^2f+2r^2L^2f+8AB+6PB+3r^2B.             (2)
```

Suppose

```text
f!=0,                   S(T,f)=G^2,
deg_y(T)=deg_y(f)=n>=1.                                  (3)
```

The identity below does not use the source address, so it applies in
particular to the addressed THM-3881 universe.  It is exact at every response
depth; it is not a truncated jet calculation.

## 1. Coordinates at y-infinity

Put `epsilon=y^(-1)` and write in `k(x)((epsilon))`

```text
T=y^nU,                 f=y^nV,
G=y^(2n+2)Gamma,        H=Gamma/V.                        (4)
```

Here `U,V in k[x][epsilon]`, their constant terms are the nonzero leading
coefficients `u,v`, and `V` is therefore a unit after passing to `k(x)`.
Choose `d in k` with `d^2=-3` and define the two colors

```text
C_-=H-dU,                C_+=H+dU.                        (5)
```

Introduce the three exact normalized sidecars

```text
E=PV^2-U^2,
Q=kappa+aU/V,
Rhat=1+Q epsilon^2,
Ahat=(1+kappa epsilon^2)U+aPV epsilon^2.                  (6)
```

They are not arbitrary notation: directly from `(1),(4)`,

```text
r=y^(n+2)V Rhat,
A=y^(n+2)Ahat,
B=y^(2n)E.                                                 (7)
```

## 2. The exact all-depth identity

The full color product is

```text
C_-C_+
=3U^2+3 Rhat^2 E
 +2L^2 epsilon^n V Rhat^2
 +8 epsilon^(n+2) Ahat E/V^2
 +6L^2 epsilon^(2n+2) Ahat/V
 +6P epsilon^(2n+4) E/V^2
 +6PL^2 epsilon^(3n+4)/V
 +L^4 epsilon^(4n+4)/V^2.                                (8)
```

This is an identity in `k(x)((epsilon))`.  Every denominator is an explicitly
displayed power of `V`; the identity records that denominator debt rather
than hiding it.  Passing to `k(x)` does not by itself prove that any quotient
is polynomial in `x`.

To prove `(8)`, divide `(2)` by `y^(4n+4)V^2`.  Since `S=G^2`, the result is
`H^2`; adding `3U^2` gives `C_-C_+`.  Equations `(7)` send the seven macros
of `(2)` respectively to

| macro | normalized contribution | base arrival |
|---|---|---:|
| `3r^2B` | `3 Rhat^2E` | `0` |
| `2r^2L^2f` | `2L^2 epsilon^n V Rhat^2` | `n` |
| `8AB` | `8epsilon^(n+2)Ahat E/V^2` | `n+2` |
| `6AL^2f` | `6L^2epsilon^(2n+2)Ahat/V` | `2n+2` |
| `6PB` | `6Pepsilon^(2n+4)E/V^2` | `2n+4` |
| `6PL^2f` | `6PL^2epsilon^(3n+4)/V` | `3n+4` |
| `L^4` | `L^4epsilon^(4n+4)/V^2` | `4n+4` |

This proves both the identity and the base arrival schedule

```text
0, n, n+2, 2n+2, 2n+4, 3n+4, 4n+4.                      (9)
```

The factors `Rhat` and `Ahat` contain their visible shifts by two and four;
the coefficients of `U,V` carry the ordinary lower-jet responses.  Thus
`(9)` records when a new residual macro first enters, not every coefficient
that the macro subsequently produces.

Expanding `(8)` modulo `epsilon^4` gives THM-3905's third-response law.  The
new compression is that no later response can receive an unlisted macro.

## 3. The persistent response terminates at depth four

The part of `(8)` arriving at depth zero is

```text
3U^2+3(1+Qepsilon^2)^2E.                                 (10)
```

It has responses only at depths zero, two, and four.  In particular its
depth-four coefficient is exactly

```text
3Q^2E.                                                     (11)
```

This explains the cancellation seen below: after depth four, every genuinely
new obligation comes from one of the six marked macros in `(9)`, rather than
from an indefinitely continuing quartic background.

## 4. Canonical colors fail at their first marked source

Use the canonical leading payment from THM-3899:

```text
v=1,
h=(a+3L^2)/2,            u=(3L^2-a)/(2d),
h-du=a,                  h+du=3L^2,
h^2=3(P-u^2).                                               (12)
```

Take the address-compatible zero-lower-sidecar family

```text
f=y^n,                    T=uy^n.                          (13)
```

Then `U=u`, `V=1`, `E=h^2/3`, and `Rhat=1+Qepsilon^2` with
`Q=kappa+au`.  Before the first marked source at depth `n`, equation `(8)` is
the exact persistent square

```text
C_-C_+=h^2 Rhat^2+3u^2       modulo epsilon^n.            (14)
```

Choose the square-root sign whose leading term is `h`.  Formal square-root
uniqueness gives

```text
H=h Rhat                  modulo epsilon^n.               (15)
```

At depth `n`, the unique new macro is `2r^2L^2f`; its coefficient is
`2L^2`.  If `J_n` denotes the response added to `(15)`, coefficient
comparison gives

```text
2hJ_n=2L^2,              J_n=L^2/h.                       (16)
```

But

```text
h=(x+1+3(9x+4)^2)/2
```

is a nonunit quadratic and `gcd(h,L)=1`.  Hence `L^2/h` is not in `k[x]`.
Consequently the **specific canonical family `(13)` is never a polynomial
square**, for any `n>=1`.  This does not constrain arbitrary lower sidecars.

The rational response `J_n=L^2/h` does exist over `k(x)`.  Thus the failure
mechanism is polynomial descent at a new denominator divisor, not formal
incompatibility of the equianharmonic colors.

## 5. The new fourth-response boundary

The cases `n=1,2,3` of `(16)` recover the canonical hostiles in THM-3902 and
THM-3905.  The first new row is `n=4`.  Through depth three,

```text
H=h+hQepsilon^2+O(epsilon^4).                             (17)
```

At depth four the left side contributes

```text
h^2Q^2+2hJ_4,                                              (18)
```

while `(8),(11)` contribute

```text
h^2Q^2+2L^2.                                               (19)
```

The persistent terms cancel exactly, leaving `J_4=L^2/h`, which is not
polynomial.  For `n>=5`, the identity proves only that the same canonical
family has the positive fourth response

```text
H=hRhat+O(epsilon^5);                                      (20)
```

it then fails later, at its own marked depth `n`, by `(16)`.  No claim is
made that an arbitrary `n>=5` sidecar survives even four responses, much
less that it extends to a square.

## 6. Reservation lineage and the exact C3 constant field

The original THM-3898 file was an **empty, unproved reservation** for an
equianharmonic cubic-order carrier.  Its owner exited before a theorem was
canonized, and root explicitly reassigned the namespace to the stronger
all-depth identity `(8)`.  No cubic-order or affine-atlas claim is activated
here.

One exact constant-field observation from that exploration survives as
provenance.  If

```text
lambda^2-lambda+1=0,             d=2lambda-1,              (21)
```

then

```text
d^2=-3,
(3+d)/2=1+lambda,         (3-d)/2=2-lambda,
(1+lambda)(2-lambda)=3.                                  (22)
```

Thus the former carrier's two constants were exactly the two norm colors
already used in `(5),(12)`.  The uncanonized carrier construction and its
infinity geometry remain outside the theorem.

## 7. Scope and reproduction

Equation `(8)` is exact and all-depth.  Its canonical corollary closes only
the zero-lower-sidecar family `(13)`.  Polynomial choices of lower
coefficients can change every response beginning at depth one; THM-3902 and
THM-3905 exhibit such repairs.  Positive equality, strict-degree lanes, a
Keller atlas, and `JC(2)` remain **OPEN**.

Run

```bash
python3 04-computation/jc2_equianharmonic_equality_color_all_depth_response_thm3898.py
python3 -O 04-computation/jc2_equianharmonic_equality_color_all_depth_response_thm3898.py
```

Both streams must byte-match
`05-knowledge/results/jc2_equianharmonic_equality_color_all_depth_response_thm3898.out`.
The companion has 46 active gates.  It treats `q=epsilon^n` as an independent
symbol for the identity itself, then samples `n=1,...,8` only as a hostile
replay of the all-degree source argument.  The samples are controls; the
affine schedule `(9)` is the proof for every `n>=1`.
