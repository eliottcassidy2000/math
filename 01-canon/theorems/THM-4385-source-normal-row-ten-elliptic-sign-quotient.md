---
id: THM-4385
title: "Source-normal row-ten elliptic sign quotient"
status: >
  PROVED GEOMETRIC COROLLARY OF THM-4380 + VERIFIED-EXACT + INDEPENDENTLY
  AUDITED. The reduced row-ten source projection in the fixed THM-4308
  residual-weight-at-most-twelve chart is a smooth elliptic curve with its
  four two-torsion points removed. Simultaneous sign change is elliptic
  negation, and r=eta/Phi is its degree-two quotient; the branch polynomial
  A(r)B(r) is squarefree of degree four. The THM-4380 row-eleven septimic
  avoids this branch divisor, explaining its fourteen reduced lifts, while
  its row-twelve quadratic is disjoint. This is geometry internal to one
  finite source chart. Chart or seam entry, stability at higher residual
  weights, a Keller pair, JC(2), and DC(2) remain OPEN.
source: root / JC2 continuation session, 2026-09-03
audit: >
  PASS, 21 exact checks plus an independent proof and implementation audit.
  Exact rational
  Groebner and infinity-gcd certificates prove projective smoothness, with
  an independent good-reduction route modulo 7. The verifier checks the
  boundary factorization, sign action, quadratic function-field extension,
  squarefree branch quartic, and the K/J divisor separations. Normal,
  optimized, and hash-seeded runs byte-match the frozen LF transcript.
depends_on:
  - THM-4380-source-normal-weight-twelve-row-twelve-extinction
related:
  - THM-4377-reciprocal-source-mutation-and-boundary-jet-cokernel
  - THM-4378-bilateral-packet-quotient-reciprocal-eigenlattice
script: 04-computation/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.py
output: 05-knowledge/results/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.out
script_sha256: 780481b1610e281846450edf3b2cbc38a2750c0d7cd36a213b7e90908ed2f4df
output_sha256: 2e68ebaa0c15e44347078adc166a328eb96e64b1664dadf5980dc3df1595aba8
hash_basis: raw LF bytes
---

# THM-4385 -- the row-ten carrier is elliptic before taking signs

**PROVED GEOMETRIC COROLLARY OF THM-4380 + VERIFIED-EXACT + INDEPENDENTLY
AUDITED.** This theorem identifies the geometry hidden by THM-4380's
ratio substitution. It changes neither that theorem's finite universe nor
its extinction conclusion.

Work over an algebraically closed field `k` of characteristic zero. Retain
THM-4380's fixed source-normal residual-weight-at-most-twelve chart and write

```text
P=Phi,                  e=eta.
```

## 1. The smooth cubic

THM-4380's `P!=0` row-ten compatibility equation is

```text
D(P,e)=7231154026500P^3+50541940696500P^2e
      +6793915500000Pe^2+353642000625e^3
      -631918028977864704P-91584545734393856e.          (1)
```

Let `Ebar` be the projective plane cubic

```text
Dbar(P,e,T)=7231154026500P^3+50541940696500P^2e
           +6793915500000Pe^2+353642000625e^3
           -631918028977864704PT^2
           -91584545734393856eT^2=0.                    (2)
```

This cubic is smooth. The exact affine singular ideal is

```text
(D, partial_P D, partial_e D)=(1) in Q[P,e].             (3)
```

At infinity, `P` cannot vanish. In the coordinate `r=e/P`, the binary cubic
is

```text
A(r)=353642000625r^3+6793915500000r^2
    +50541940696500r+7231154026500,                       (4)
```

and `gcd(A,A')=1`. Thus the three geometric points at infinity are distinct
and nonsingular. The companion verifies `(3)--(4)` exactly over `Q`; its
separate unit-ideal and squarefree checks modulo 7 are a good-reduction
control, not a rational-point census.

The point

```text
O=[0:0:1]                                                   (5)
```

is rational and smooth. Hence `Ebar`, with origin `O`, is an elliptic curve.

## 2. The two THM-4380 strata form `Ebar` minus its two-torsion

On the boundary `P=0`, equation `(1)` factors as

```text
D(0,e)=19e F(e),
F(e)=18612736875e^2-4820239249178624.                    (6)
```

THM-4380 proves on reduced loci that its row-ten projected source locus is
`D=0` on `P!=0`, while on `P=0` it consists of the two reduced nonzero roots
of `F`. Equation `(6)` therefore supplies the missing gluing statement:
the complete reduced `(P,e)` projection is

```text
Ebar intersect {T!=0}, with O removed.                    (7)
```

Equivalently, it is `Ebar` minus `O` and the three points at infinity.

The involution

```text
sigma[P:e:T]=[-P:-e:T]                                    (8)
```

preserves `(2)`. In the affine chart its only fixed point is `O`; every point
at infinity is fixed projectively, and `(4)` says there are exactly three of
them. Since `sigma` is a nonidentity order-two automorphism fixing the chosen
origin, it is elliptic negation. Its four fixed points are exactly `Ebar[2]`.
Thus `(7)` becomes the intrinsic identity

```text
boxed: row-ten reduced source projection = Ebar minus Ebar[2]. (9)
```

In particular the two genuine `P=0` points are not a detached exceptional
stratum: they are an exchanged pair in the same smooth elliptic carrier.

## 3. The ratio coordinate is the sign quotient

On `P!=0`, put `r=e/P` and define

```text
B(r)=91584545734393856r+631918028977864704.              (10)
```

Direct substitution into `(1)` gives

```text
D(P,rP)=P(P^2A(r)-B(r)).                                  (11)
```

Moreover

```text
deg A=3,  deg B=1,  gcd(A,B)=1,
gcd(A B,(A B)')=1.                                       (12)
```

Consequently the elliptic function field is

```text
k(Ebar)=k(r)(P),                 P^2=B(r)/A(r),           (13)
```

and `sigma` fixes `r` and changes the sign of `P`. Thus

```text
q:Ebar -> P^1_r                                            (14)
```

is the degree-two quotient by `sigma`, with branch divisor the four simple
roots of `A B`. Equivalently, after `y=A(r)P`, its quartic model is

```text
y^2=A(r)B(r).                                             (15)
```

The root of `B` is the value of `r` at `O`, read from the tangent
`631918028977864704P+91584545734393856e=0`; the three roots of `A` are the
points at infinity. The two nonzero `P=0` boundary points of `(6)` lie over
`r=infinity` and are unramified. Therefore `(9)` is precisely the etale
double-cover locus of `(14)`. The ratio quotient loses only the sign sheet,
but that sheet is indispensable for recovering the source-point count.

## 4. Why seven becomes fourteen, and why row twelve is empty

Let `K(r)` be THM-4380's primitive row-eleven septimic and `J(r)` its
row-twelve quadratic. Exact Euclidean calculations give

```text
deg K=7,       gcd(K,K')=1,
gcd(K,A B)=1,  gcd(K,J)=1.                                (16)
```

The first two statements give seven distinct finite quotient points. The
third says all seven avoid the ramification divisor, so each has exactly two
distinct lifts under `(14)`. This recovers THM-4380's fourteen reduced
row-eleven source points as the etale pullback

```text
q^* V(K).                                                  (17)
```

The two points over `r=infinity` are the `P=0` pair already killed by
THM-4380's row-eleven residual. Finally `gcd(K,J)=1` says the row-twelve
divisor misses every point of `(17)`. This is the quotient-geometric content
of the literal mod-29 Bezout extinction in THM-4380.

## 5. Connection contract and stopping boundary

| field | contract |
|---|---|
| source | the reduced row-ten source projection in the fixed THM-4308 weight-at-most-twelve chart |
| target | the elliptic curve `Ebar` and quotient line `P^1_r` |
| map | projective completion followed by `q`, with `r=e/P` |
| preserved predicate | row-ten bracket and joint projected-depth compatibility; row-eleven and row-twelve quotient divisors |
| destroyed information | the sign sheet, terminal affine-eight fibre, chart ownership, higher residual weights, and global-map data |
| required sidecar | the squarefree branch divisor `A B` and the two-sheet lift `P^2=B/A` |
| cheapest next test | add the first omitted residual-weight-13 face and recompute whether its row-nine/ten equations remain tangent to `Ebar` |

The signed reciprocal quotient of THM-4378 is a useful warning, not a proved
identification: both quotients lose a two-sheet coordinate, but no map from
its packet two-glue to the present elliptic cover has been constructed.

The sharp next question is higher-weight stability. The first omitted face at
residual weight 13 can already enter rows seven and eight, so the elliptic
carrier may be a cap-specific object. Only if that hostile extension remains
tangent should one seek a global seam map. This theorem itself proves no chart
entry, no seam entry, no all-weight lift, no Keller pair, and no consequence
for `JC(2)` or `DC(2)`, which remain **OPEN**. **QED.**

## 6. Exact replay

The 21-check companion independently verifies the geometry extracted from
THM-4380:

```text
python3 -B 04-computation/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.py
python3 -O -B 04-computation/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.py
PYTHONHASHSEED=17 python3 -B 04-computation/jc2_source_normal_row10_elliptic_sign_quotient_thm4385.py
```

All three streams byte-match the stored LF output. Raw-byte hashes are

```text
script: 780481b1610e281846450edf3b2cbc38a2750c0d7cd36a213b7e90908ed2f4df
output: 2e68ebaa0c15e44347078adc166a328eb96e64b1664dadf5980dc3df1595aba8
```
