---
id: THM-3901
title: "Nonzero-sidecar osculating response fan"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  In the complementary
  THM-3899 region deg_y(T)>deg_y(f), the complete leading fan and its exact
  first osculating response are classified, including every equality wall,
  the constant-sidecar boundary, and the deg_y(T)=deg_y(f)+2 cancellation
  stratum.  The formerly automatic first-response cell has a second exact
  response and a new divisibility tariff.  This is a necessary passport,
  not an existence theorem or an emptiness theorem.
source: tournament-jc-breakthrough / post-THM-3899 response scout, 2026-08-23
audit: >
  INDEPENDENTLY HOSTILE-AUDITED on 2026-08-23 by two audits.  They rebuilt
  every affine degree fan, the UFD wall passport, both osculating defects,
  the sharpened second-response tie, and a separate 49-cell model of the
  deg_y(T)=deg_y(f)+2 cancellation wall.  All boundaries and payable controls
  passed.  One omitted wall was restored to the summary list; it was already
  present and correct in the detailed proof.  The focused companion expands
  both exact defects and checks
  192 generic leading cells, 192 first-response cells, 160 second-response
  cells, all twelve sampled constant-sidecar cells, the corrected
  deg_y(T)=deg_y(f)+2 wall simplification, the exact r=0 address branch, and
  sharp leading payments in 1152 active gates.  Normal and optimized runs
  byte-match the frozen output.
depends_on:
  - THM-3881-cusp-ideal-residual-transport-rank-two-matrix-factorization
related:
  - THM-3888-f-zero-equianharmonic-jacobian-and-two-section-integrality
  - THM-3895-f-zero-quartic-covariant-and-high-y-degree-emptiness
  - THM-3897-f-zero-residual-all-degree-global-emptiness
  - THM-3899-nonzero-sidecar-y-degree-tariff-and-equianharmonic-equality-colors
  - THM-3904-nonzero-sidecar-constant-y-seam-emptiness-and-primitive-equality-colors
script: 04-computation/jc2_nonzero_sidecar_osculating_response_fan_thm3901.py
output: 05-knowledge/results/jc2_nonzero_sidecar_osculating_response_fan_thm3901.out
script_sha256: 69b04b7686ff3be0051cf1f2c81672892db49b0027f735be49ecf4a2f8446876
output_sha256: 4366aef6140e13fa501e619260c6515a40fb7f4da09d5bbed80cfc48209095f6
semantic_sha256: ad2fd5d83e8b53b75b9721915701085a71bdc44d694ec516ec1a516a28799f7d
hash_basis: raw LF bytes
---

# THM-3901 -- the nonzero-sidecar osculating response fan

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  In `D=k[x,y]` put

```text
a=x+1,                       L=9x+4,
F=15x^2+15x+4,              K=y^2-F,
P=aL^2,                     r=aT+Kf,
A=KT+aPf,                   B=Pf^2-T^2.                  (1)
```

Retain the THM-3881 residual and module address

```text
S=L^4+2(3A+3P+r^2)L^2f+(8A+6P+3r^2)B=G^2,
T(0,0)=4f(0,0).                                             (2)
```

Suppose `f!=0` and

```text
m=deg_y(T)>n=deg_y(f),      delta=m-n,
u=lc_y(T),                  v=lc_y(f).                    (3)
```

Fix `s in k` with `s^2=-3`.  The assertions below are necessary conditions
for `(2)`.  They do not assert that any displayed leading payment extends to
a polynomial square.

## 1. The fifteen-term support and the complete positive-degree top fan

Apart from the base term `L^4`, the exact expansion is

```text
 6aL^4f             12a^2L^4f^2          8a^3L^4f^3
 6KL^2Tf            12aKL^2Tf^2          6a^2KL^2Tf^3
-6aL^2T^2           -6a^2L^2T^2f         3a^3L^2T^2f^2
-8KT^3              -6aKT^3f            -3a^2T^4
 2K^2L^2f^3          3aK^2L^2f^4         -3K^2T^2f^2.     (4)
```

For `n>=1`, their y-degrees in the same order are

```text
n, 2n, 3n,
m+n+2, m+2n+2, m+3n+2,
2m, 2m+n, 2m+2n,
3m+2, 3m+n+2, 4m,
3n+4, 4n+4, 2m+2n+4.                                    (5)
```

Writing `w=au+v`, the complete top fan is

| regime | `deg_y(S)` | `lc_y(S)` |
|---|---:|---|
| `delta=1` | `4n+6` | `-3u^2v^2` |
| `delta=2`, `w!=0` | `4n+8=4m` | `-3u^2w^2` |
| `delta>=3` | `4m` | `-3a^2u^4` |

The first row is the unique `K^2T^2f^2` term.  The middle row is the exact
sum of `T^4`, `KT^3f`, and `K^2T^2f^2`.  The last row is the unique `T^4`
term.  Thus, after replacing `G` by `-G` if necessary, its leading term
agrees respectively with

```text
Q0=sTr,                                                     (6)
```

whose leading coefficients are `suv`, `suw`, and `sau^2`.  The other
square-root sign merely exchanges the two arms `G-Q0` and `G+Q0`.

The sole cancellation of this top fan is

```text
delta=2,                     w=au+v=0.                    (7)
```

It is treated in Section 4 rather than silently divided away.

## 2. Exact first response and its divisibility passports

The osculating subtraction has the exact two-arm identity

```text
(G-Q0)(G+Q0)=R,                                             (8)

R=S+3T^2r^2
 =L^4+6(A+P)L^2f+(8A+6P)B+r^2L^2f(2+3af).                 (9)
```

Its full expansion is

```text
R=3aK^2L^2f^4+2K^2L^2f^3
 +6a^2KL^2Tf^3+12aKL^2Tf^2+6KL^2Tf-8KT^3
 +8a^3L^4f^3+12a^2L^4f^2+6aL^4f+L^4
 +3a^3L^2T^2f^2-6a^2L^2T^2f-6aL^2T^2.                  (10)
```

Only four possible top responses survive the affine comparisons:

```text
term             degree                    leading coefficient
KT^3             3m+2                     -8u^3
T^2f^2           2m+2n                    3a^3L^2u^2v^2
KTf^3            m+3n+2                   6a^2L^2uv^3
K^2f^4           4n+4                     3aL^2v^4.       (11)
```

Put `W=G-Q0`.  Since `lc_y(G+Q0)` is twice the chosen leading coefficient,
comparison of `(10)` with `G+Q0` gives the following complete first-response
table for `n>=1`.  In a tie row, the displayed degree and leading coefficient
apply when the numerator is nonzero.  If it vanishes, the degree drops, but
the stated divisibility passport still follows from the vanishing identity.

| regime | top of `R` | `deg_y(W)` | response coefficient / necessary passport |
|---|---:|---:|---|
| `delta=1,n=1` | `-8u^3+3aL^2v^4` | `3` | `(-8u^3+3aL^2v^4)/(2suv)`; hence `u|aL^2v^4`, `v|u^3` |
| `delta=1,n>=2` | `3aL^2v^4` | `2n+1` | `-saL^2v^3/(2u)`; hence `u|aL^2v^3` |
| `delta=2,w!=0,n<4` | `-8u^3` | `n+4` | `-4u^2/(sw)`; hence `w|u^2` |
| `delta=2,w!=0,n=4` | `-8u^3+3aL^2v^2w^2` | `8` | quotient by `2suw`; hence `w|u^3` and `u|aL^2v^2w^2` |
| `delta=2,w!=0,n>4` | `3aL^2v^2w^2` | `2n` | `3aL^2v^2w/(2su)`; hence `u|aL^2v^2w` |
| `delta>=3,n<delta+2` | `-8u^3` | `m+2` | `4su/(3a)`; hence `a|u` |
| `delta>=3,n=delta+2` | `u^2(-8u+3a^3L^2v^2)` | `m+2=2n` | `(-8u+3a^3L^2v^2)/(2sa)`; hence `a|u` |
| `delta>=3,n>delta+2` | `3a^3L^2u^2v^2` | `2n` | `-sa^2L^2v^2/2`; automatic |

All divisibilities are in the UFD `k[x]`.  For example, in the second row
`W in k[x,y]` forces its leading coefficient to lie in `k[x]`, so
`u|aL^2v^3`.  No parity conclusion is available from `W`: it is a response
polynomial, not a square.

Some useful primitive corollaries are immediate.  If `gcd(u,v)=1`, then
`gcd(u,w)=1`; consequently

```text
delta=1,n=1:        v is a unit and u|aL^2;
delta=1,n>=2:       u|aL^2;
delta=2,n<4:        w is a unit;
delta=2,n=4:        w is a unit and u|aL^2;
delta=2,n>4:        u|aL^2.                              (12)
```

No primitivity is assumed in the theorem itself.

On the first tie wall, cancellation has an exact UFD parametrization:

```text
8u^3=3aL^2v^4

iff  u=c a^3L^2q^4,          v=d a^2Lq^3,
     8c^3=3d^4,              q in k[x]\{0}.               (13)
```

This follows prime by prime: away from `aL`, the exponents are `(4e,3e)`;
at `a` they are `(4e+3,3e+2)`, and at `L` they are `(4e+2,3e+1)`.

## 3. Second response in the formerly automatic cell

Only the final row of the first-response table is advanced a second step.
There the leading term of `W` is exactly
`-sa^2L^2v^2y^(2n)/2`, so define

```text
Q1=s(Tr-a^2L^2f^2/2).                                    (14)
```

The second defect is exact:

```text
S-Q1^2=
 3aK^2L^2f^4+2K^2L^2f^3
 +3a^2KL^2Tf^3+12aKL^2Tf^2+6KL^2Tf-8KT^3
 +(3/4)a^4L^4f^4+8a^3L^4f^3+12a^2L^4f^2+6aL^4f+L^4
 -6a^2L^2T^2f-6aL^2T^2.                                (15)
```

Assume here, and only here, that

```text
delta>=3,                     n>delta+2.                  (16)
```

Put `H=G-Q1`.  The second fan and its new tariffs are

| regime inside `(16)` | top of `S-Q1^2` | `deg_y(H)` | response coefficient / passport |
|---|---:|---:|---|
| `delta+2<n<2delta` | `-8u^3` | `m+2` | `4su/(3a)`; hence `a|u` |
| `n=2delta` | `u(-8u^2+3a^2L^2v^3)` | `m+2` | `(-8u^2+3a^2L^2v^3)/(2sau)` |
| `n>2delta` | `3a^2L^2uv^3` | `2n-delta+2` | `3aL^2v^3/(2su)`; hence `u|aL^2v^3` |

At the equality row, polynomiality gives the sharper exact passport

```text
a|u,             and, writing u=au1,          u1|L^2v^3. (17)
```

Indeed the response coefficient becomes

```text
(-8u1^2+3L^2v^3)/(2su1).                                 (18)
```

The same conclusions hold when a displayed tie numerator vanishes.  This
second response is not asserted outside `(16)`.

## 4. The `delta=2`, `w=0` cancellation stratum

Now assume `n>=1`, `delta=2`, and `w=0`.  Thus

```text
v=-au,                         deg_y(r)<=n+1.              (19)
```

If `r!=0`, put `t=deg_y(r)` and `rho=lc_y(r)`.  From the structural form
of `(2)`, the only possible new top terms are `8AB` and `3r^2B`:

| r-order | `deg_y(S)` | `lc_y(S)` | necessary consequence |
|---|---:|---|---|
| `r=0` or `2t<n+4` | `3n+8` | `-8u^3` | `n` is even and `u` is a square up to a unit |
| `2t=n+4` | `3n+8` | `-u^2(8u+3rho^2)` | `n` is even; norm law below |
| `2t>n+4` | `2t+2n+4` | `-3u^2rho^2` | square-colored; use `Q0=sTr` |

In the tie row, if the displayed coefficient is nonzero, `u|lc_y(G)` and
there is an `h in k[x]` with

```text
h^2=-(8u+3rho^2),
(h-srho)(h+srho)=-8u.                                    (20)
```

If it vanishes, then `rho^2=-8u/3`, so `u` is already square-colored and a
lower response is required.  The parity restriction comes from the degree
of the original square `G^2`; it is not imposed on either response arm.

In the last, r-dominant row, `(8)` applies again.  The top of `R` is now a
two-term fan between `8AB` and `3ar^2L^2f^2`.  Since
`deg_y(G+Q0)=n+2+t`, the corrected wall response is

| nested r-order | top of `R` | `deg_y(W)` | response coefficient / passport |
|---|---:|---:|---|
| `n+4<2t<n+8` | `-8u^3` | `2n+6-t` | `-4u^2/(srho)`; hence `rho|u^2` |
| `2t=n+8` | `-8u^3+3aL^2rho^2v^2` | `2n+6-t` | `u(-8u+3a^3L^2rho^2)/(2srho)`; hence `rho|u^2` |
| `2t>n+8` | `3aL^2rho^2v^2` | `t+n-2` | `3a^3L^2rho u/(2s)`; automatic |

The last column has used `v=-au`.  Retaining the unsimplified quotient would
hide that the final cell is automatic and would state only a vacuous
divisibility.

The exact `r=0` branch is particularly cheap and respects the origin
address.  Since `gcd(a,K)=1`, it has

```text
f=aq,                         T=-Kq,
C=a^3L^2-K^2,
S=L^4(1+6a^2q)+12aL^2q^2C+8q^3C^2.                      (21)
```

At `(x,y)=(0,0)`, `K=-4`, so `T=4q=4f`; hence the address does not remove
this branch.  More generally the address makes `r(0,0)=0`, but a pointwise
zero is not a polynomial divisibility and is not used in the leading fan.

## 5. Constant-sidecar boundary `n=0`

Let `f=v in k[x]\{0}`.  The three distinct boundary cells are

| `m` | top of `S` | exact consequence |
|---:|---|---|
| `1` | degree `6`, coefficient `-3u^2v^2` | `Q0` matches; `deg_y(W)=2`, `lc_y(W)=-4u^2/(sv)`, hence `v|u^2` |
| `2` | degree `8`, coefficient `-u^2(8u+3(au+v)^2)` | if nonzero, `h^2=-(8u+3(au+v)^2)` and `(h-s(au+v))(h+s(au+v))=-8u`; if zero, lower response required |
| `m>=3` | degree `4m`, coefficient `-3a^2u^4` | `Q0` matches; `deg_y(W)=m+2`, `lc_y(W)=4su/(3a)`, hence `a|u` |

This boundary is not covered by the positive-`n` threshold table and is
therefore stated separately.

## 6. Proof candidate and audit boundary

Expansion gives `(4)` and `(5)`.  Every entry of each fan is the maximum of
the displayed affine degree functions.  Their pairwise differences have the
following only equality walls:

```text
top S:             delta=2;
first R, delta=1:  n=1;
first R, delta=2:  n=4;
first R, delta>=3: n=delta+2;
second defect:     n=2delta;
wall S:            2t=n+4;
wall R:            2t=n+8.                               (22)
```

This proves the degree partitions once the corresponding leading sums are
collected.  Identities `(8)`, `(9)`, and `(15)` then turn each nonzero defect
coefficient into the quotient of two leading coefficients.  Because both
arms lie in `k[x,y]`, every such quotient lies in `k[x]`; Euclid's lemma in
the UFD `k[x]` gives the listed divisibilities.  On a zero tie coefficient,
the same divisibilities follow directly from the tie equation.  Section 4
repeats the comparison after the sole top cancellation and Section 5 audits
the excluded constant boundary.

The tariffs are genuinely payable at this layer.  Representative controls
include

```text
delta=1,n=1 tie:   u=216a^3L^2q^4, v=72a^2Lq^3;
delta=1,n>=2:      u=aL^2, v=1;
delta=2,n<4:       u=1, v=1-a, so w=1;
delta=2,n=4:       the same unit-w choice pays the exact quotient;
delta=2,n>4:       u=aL^2, v=1;
delta>=3 low/tie:  u=a, v=1;
Q1 high cell:      u=aL^2, v=1.                          (23)
```

For the first control, `8*216^3=3*72^4`.  These are hostile controls against
promoting a necessary response tariff to emptiness.  They pay only the
displayed leading equation, not all coefficients of `(2)`.

The theorem deliberately stops before any third response, before solving
the `delta=2` zero-norm walls, and before using a Keller atlas.  It proves
neither an actual survivor nor global emptiness in the `f!=0` lane.  Its
relationship to THM-3888/3895 is mechanistic: when `f=0`, `Q0=saT^2` is the
same equianharmonic osculating arm used in the proved quartic descent.  Here
the sidecar moves that arm through the response fan, and the extra
divisibility symbols are the information needed to continue an x-integral
descent.  THM-3897 closes only the separate `f=0` lane.

## 7. Reproduction

```bash
python3 04-computation/jc2_nonzero_sidecar_osculating_response_fan_thm3901.py
python3 -O 04-computation/jc2_nonzero_sidecar_osculating_response_fan_thm3901.py
```

Both runs must byte-match
`05-knowledge/results/jc2_nonzero_sidecar_osculating_response_fan_thm3901.out`.
The companion has 1152 active gates.  It is an exact identity/support audit
and a broad finite hostile control, not a substitute for independent review
of the all-degree UFD proof candidate.
