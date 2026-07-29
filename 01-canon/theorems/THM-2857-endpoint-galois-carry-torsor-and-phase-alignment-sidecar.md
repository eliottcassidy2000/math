---
id: THM-2857
title: "Endpoint Galois carry torsor and phase-alignment sidecar"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT.  The THM-2847 endpoint scalar has an
  exact free C13 Galois orbit over Q(zeta_91); its centered orbit is the
  faithful character 3, and its tenth powers remain point-separating with
  Fourier support F13 minus {7,10}.  External tensoring turns any first
  factorial exit m<13 into exactly m+1 carry channels.  The existing
  proved coefficient-support ancestry action is K0-linear and fixes the
  endpoint scalar, so a new semilinear clutch is necessary.  Final immutable
  audit is pending.
source: root/lrc-endpoint-galois-carry-torsor-2026-07-28
depends_on:
  - THM-2380-cross-word-charged-target-correlation-and-pair-twist-gate
  - THM-2611-principal-c13-bibundle-lift-torsor-and-holonomy-section-obstruction
  - THM-2835-q11-semantic-word-horn-and-bockstein-blind-support-no-go
  - THM-2847-q3-q11-transverse-endpoint-edge-and-e3-realization-horn
  - THM-2851-q3-q11-q7-bockstein-holonomy-and-realization-sidecar
related:
  - THM-2460-idempotent-semantic-word-copy-and-word-block-cosupport-boundary
  - THM-2788-physical-modular-odometer-versus-heisenberg-bockstein-extension
  - THM-2820-boolean-idempotent-rigidity-and-norm-top-cotangent-jet-no-go
  - THM-2855-shifted-positive-cone-transverse-asymptotic-family
  - THM-2858-six-ray-open-half-plane-moment-eight-gaussian-hostile
script: 04-computation/lrc14_endpoint_galois_carry_torsor_thm2857.py
output: 05-knowledge/results/lrc14_endpoint_galois_carry_torsor_thm2857.out
script_sha256: 0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c
output_sha256: ac1194c46db2cdf43c807ece781b63971c081cc5f9070964007fdecdc20f1583
hash_basis: LF-normalized bytes
---

# THM-2857 -- endpoint Galois carry torsor and phase-alignment sidecar

**PROVED CANDIDATE + VERIFIED-EXACT.**

THM-2847's apparently inconvenient coefficient field already contains a
sharp algebraic model of the missing first-carry fibre.  The endpoint scalar
has a free thirteen-element Galois orbit, and its centered orbit is one
faithful `C_13` character.  After the exact q3-to-q11 exponent alignment, the
thirteen tenth powers are still distinct.  Thus phase separation needs a
point-separating orbit, not all twelve nontrivial Fourier characters.

This is not an LRC(14) closure.  The theorem constructs a coefficient-field
torsor and proves that the current canonical carry action does not realize
it: that action is linear over the endpoint field and leaves `c` fixed.
It does not realize the tenth power by a lawful packet operation or supply
the endpoint-matched reference required for a physical pair twist.

## 1. The relative cyclotomic torsor

Let

```text
zeta=zeta_1183,          1183=7*13^2,
K=Q(zeta),               F=Q(zeta^13)=Q(zeta_91),
omega=zeta^91=zeta_13.
```

The relative Galois group itself is

```text
Gal(K/F)={sigma_r:r in F_13},
sigma_r(zeta)=zeta^(1+91r).                              (1)
```

Indeed, the units `1+91r` form a subgroup modulo `1183` because
`91^2=7*1183`, and they fix `zeta^13`.  Moreover

```text
[K:F]=phi(1183)/phi(91)=936/72=13,                       (2)
```

so these are all relative automorphisms.

THM-2847 gives the exact nonzero endpoint scalar

```text
c=zeta^624-zeta^510.
```

Put `A=zeta^624 in F` and `B=zeta^510`.  Since `624=13*48` and
`510=13*39+3`,

```text
c_r:=sigma_r(c)=A-B omega^(3r).                          (3)
```

The thirteen values are distinct, because multiplication by `3` permutes
`F_13` and `B` is nonzero.  Hence this orbit is a free `C_13`-set.  A
faithful transitive `C_13`-set has at least thirteen points, so its cardinality
exactly attains the first-carry invoice in THM-2851.  This is only an
algebraic realization of that invoice.

The mean, centered orbit, forward Galois difference, and norm are

```text
(1/13) sum_r c_r = A,
c_r-A              = -B omega^(3r),
c_(r+1)-c_r         = -B(omega^3-1)omega^(3r),
prod_r c_r          = A^13-B^13 !=0.                    (4)
```

Thus the orbit has normalized Fourier support `{0,3}`; both its centered
part and its forward difference have support exactly `{3}`.  Character `3`
is faithful because `13` is prime.

## 2. Exact q7-to-q11 coefficient alignment

Let `x` be the allocation-label generator of `Q[C_13]`.  The THM-2847 raw
endpoint edge is `c_r x^3`; the THM-2835/2851 semantic coefficients on the
q7/QAB and q11/QA legs are

```text
B_sem=449 x^7,             A_sem=449 x^11.
```

Since `3*10+7=11 mod13`,

```text
(c_r x^3)^10 B_sem = c_r^10 A_sem.                       (5)
```

This is an identity in the coefficient group ring.  It is not yet a
physical tenfold-current construction.

The values `c_r^10` remain pairwise distinct.  If `c_r^10=c_s^10`, then
`c_r/c_s` is a tenth root of unity in `K`.  The roots of unity of the
odd-conductor field `K` are `mu_2366`, and

```text
mu_2366 intersect mu_10={+1,-1}.
```

The plus sign gives `r=s`.  For the minus sign,

```text
2A=B(omega^(3r)+omega^(3s)).                              (6)
```

In the standard complex embedding, equality in the triangle inequality
forces `r=s` and `A/B=omega^(3r)`.  But
`A/B=zeta^114`, whose exponent is not divisible by `13`, whereas every
`omega^(3r)` has exponent divisible by `13`.  This is impossible.

There is a useful sharp spectral boundary.  The binomial expansion is

```text
c_r^10 =
  sum_(j=0)^10 (-1)^j binom(10,j) A^(10-j)B^j omega^(3jr). (7)
```

The eleven residues `3j mod13`, `0<=j<=10`, are distinct, and every
coefficient in `(7)` is nonzero.  Therefore

```text
supp_Fourier(r |-> c_r^10)
  ={0,1,2,3,4,5,6,8,9,11,12}
  =F_13\{7,10}.                                           (8)
```

The aligned scalar is injective on the carry torsor although two Fourier
channels vanish.  Consequently an argument demanding every carry character
would impose more than coefficient-level phase selection requires.

## 3. Factorial detection depth becomes carry bandwidth

There is a direct algebraic interface with the factorial/Gaussian-moment
frontier that does not require either related THM-2855 or THM-2858.  Let
`L(s^n)=n!`, and let `H` be any complex polynomial whose first nonzero
factorial moment has order `m<13`:

```text
L(H^j)=0  (1<=j<m),              L(H^m)!=0.              (9)
```

For `H_r=c_r H`, scalar homogeneity gives

```text
L(H_r^j)=c_r^j L(H^j).                                  (10)
```

The first exit remains `m`, and the binomial expansion of `c_r^m` gives
the exact carry Fourier support

```text
supp_Fourier(r |-> L(H_r^m))
  ={3j mod13:0<=j<=m}.                                  (11)
```

The residues in `(11)` are distinct, so the first exit has exactly `m+1`
carry channels.  For `m=4`,

```text
support={0,3,6,9,12},                                   (12)
```

and the thirteen values `c_r^4` are pairwise distinct by the same
`mu_2366 intersect mu_4={+1,-1}` argument used above.  Thus a fourth
factorial exit is already point-separating on the full carry torsor while
using only five Fourier channels.  For `m=5`, the support is
`{0,2,3,6,9,12}`.

This passes canonically to the standard complex-Gaussian construction.  If
`H=sh`, put

```text
P_r(Z,W)=W+Z c_r h(ZW),              W=conj(Z).
```

Charge balance gives

```text
E[P_r^(2j)]=binom(2j,j)c_r^j L(H^j),
E[P_r^(2j-1)]=0.                                         (13)
```

Hence a fourth factorial exit gives
`E[P_r^8]=70c_r^4L(H^4)`, and a fifth gives
`E[P_r^10]=252c_r^5L(H^5)`.  This preserves degree and all preceding
moment vanishings, but it is external complex scalar tensoring.  It need
not preserve a real positive-cone presentation and is not a physical LRC
packet product.  The result is an exact bridge between detection order and
carry bandwidth, not an identification of the two problems.

## 4. The current carry is scalar-trivial

THM-2851 proves on the comparison nerve

```text
L_9 L_8 = T L_4,                                        (14)
```

where `T` is one positive ancestry carry.  Suppose a future physical
construction provides a fully typed endpoint packet and an equivariance

```text
coefficient(T^r packet)=sigma_r(c),                      (15)
```

This equivariance is absent from the existing carrier.  THM-2851's proved
oriented object is the group-algebra mask

```text
h_L in Z[C_(13^5)],            partial_T h_L=(Z-1)h_L.
```

After extension to `K_0=Q(zeta_1183)`, `T` is multiplication by `Z` and is
`K_0`-linear.  Tensoring with the THM-2847 endpoint scalar therefore gives

```text
T(c h_L)=c Z h_L,              partial_T(c h_L)=c(Z-1)h_L. (16)
```

But

```text
sigma_1(c)=zeta^624-zeta^783 != zeta^624-zeta^510=c.      (17)
```

Thus the current canonical action is not `(15)`.  The first failed
implication is now exact: a free Galois orbit in the coefficient field does
not twist a coefficient-linear ancestry module.  In local-system language,
canon has the constant `K_0` local system on the ancestry circle, whereas
the endpoint torsor would require monodromy `sigma_1`.  The missing datum is
a semilinear clutch, not another scalar normalization.

Conditionally, if a future physical construction supplies that clutch while
preserving the target, semantic word, right-endpoint gauge, and phase
reference, then `(4)` detects the Bockstein loop by `omega^3`.  If it also
lawfully realizes the tenth-power coefficient in `(5)`, the distinct values
select a unique carry section.  Under the separate endpoint-matched
hypotheses of THM-2380, one fixed nonzero QA reference then recovers that
complex scalar.

None of those conditional hypotheses follows from the Galois calculation.
A physical carry acts on interval/root coordinates, whereas `sigma_1`
multiplies cyclotomic exponents.  Treating them as identical despite
`(16)--(17)` would manufacture the desired basepoint.

## 5. Exact verification

The companion pins the exact THM-2847 and THM-2851 scripts and transcripts.
It verifies:

1. all thirteen relative Galois automorphisms and the orbit formula `(3)`;
2. freeness, the symbolic supports in `(4)`, `(8)`, and `(11)`, and `(5)`;
3. the relative norm identity, Gaussian coefficients `70/252`, and the
   exact scalar-linear/Galois mismatch `(16)--(17)`;
4. an independent specialization in `F_4733` with an exact primitive
   `1183`rd root, separating all thirteen `c_r` and all thirteen `c_r^10`;
5. the fourth-, fifth-, and tenth-power supports and point-separation
   controls by exact finite-field transforms.

Normal, optimized, and stored-output replay are byte-identical after LF
normalization.  LF-normalized SHA-256:

```text
script  0bae59c9b1460f37e1879a81746154593cb0699ee13b3e5e800ba0af95ea5e4c
output  ac1194c46db2cdf43c807ece781b63971c081cc5f9070964007fdecdc20f1583
```

An independent hostile audit rederived the relative field, orbit,
tenth-power distinctness, exact support, alignment, and finite-field
certificate.  Final immutable-path replay remains pending before promotion.

## 6. Connection contract and stopping boundary

```text
source:
  THM-2847 c=zeta_1183^624-zeta_1183^510 and edge c*x^3;

target:
  THM-2851 first C13 ancestry-carry fibre and the
  THM-2835 q7/QAB -> q11/QA coefficient alignment;
  independently, any factorial/Gaussian first-exit witness of order m<13;

map:
  r |-> sigma_r(c)=A-B omega^(3r),
  (c_r*x^3)^10*(449*x^7)=c_r^10*(449*x^11),
  H |-> H_r=c_r H with first-exit spectrum {3j:0<=j<=m};

preserved algebraically:
  thirteen distinct section values, semantic coefficient 449,
  target-label alignment, faithful centered character, and relative norm;

not supplied physically:
  a semilinear T=sigma_1 clutch, a lawful tenfold endpoint product,
  common endpoint gauge,
  fixed source, positive/current realization, response basepoint,
  E3/complement truth, owner phase, or preservation of a positive cone;

cheapest decisive test:
  the current K0-linear action already fails;
  construct one genuinely semilinear carried endpoint coefficient,
  then test a lawful aligned power/reference pair by THM-2380.
```

No row is excluded and the LRC(14) ledger remains `165`.
