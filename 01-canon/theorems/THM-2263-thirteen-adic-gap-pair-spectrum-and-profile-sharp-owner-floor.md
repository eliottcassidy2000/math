---
id: THM-2263
title: "Thirteen-adic gap pair spectrum and profile-sharp owner floor"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENT DIRECT-INTERVAL REFEREE. At a
  prescribed positive 13-adic valuation gap d, the complete sharp
  danger-pair overlap interval depends only on the parity of d. Even gaps
  have upper defect 24 and lower defect -5/3; odd gaps interchange them.
  This explains THM-2255's exceptional ratio 1:169 and sharpens its
  exclusive-owner ledger profile by profile. Among the 150 strict
  first-depth-one rows, the unique worst pair ledger is (1,3,5), with
  total exclusive mass at least 15041431/197927730 and expiration image
  at least 15041431/70270200>1/7. Among the 15 repeated-first rows, the
  compatible pair maximum is uniquely at (1,1,5), with exclusive mass at
  least 5229541/197927730 and expiration image
  5229541/70270200>1/14. These are sharp for the stated overlap ledgers,
  but the missing post-expiration target inclusion remains; no profile is
  excluded and LRC(14) remains open.
source: codex-2026-07-25-gap-sensitive-owner-spectrum
depends_on:
  - THM-1166-seven-wall-fano-gcd-discrepancy
  - THM-2255-valuation-separated-pair-cap-and-exclusive-owner-mass
related:
  - THM-1191-four-comb-pair-floor-refuted-by-thirteen-adic-ladder
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2234-first-depth-one-private-two-owner-mass-and-two-step-expansion
  - THM-2246-depth-one-private-joint-two-step-fibre-cap
  - THM-2257-depth-three-common-core-169-image-sieve-exclusion
  - THM-2261-expiration-image-surjectivity-and-one-core-carrier-no-go
script:
  - 04-computation/lrc14_gap_sensitive_owner_spectrum_thm2263.py
  - 04-computation/lrc14_gap_sensitive_owner_spectrum_referee_thm2263.py
output:
  - 05-knowledge/results/lrc14_gap_sensitive_owner_spectrum_thm2263.out
  - 05-knowledge/results/lrc14_gap_sensitive_owner_spectrum_referee_thm2263.out
script_sha256:
  - 57c20994df57a5510868c113758ae42dc90beab3ef6d6061f28f67e537a6577e
  - 21c7668417dcdf3092dac1863f9a2674a8aaceb0b734fe1add7c40ab39b3fd80
output_sha256:
  - 3511eb4cd8ce5bd23037850152e3a4fe1dc9a6dda573a5b92b73d4d8dcb933e1
  - b013c51b834dd74afe10cb581fa7bd77f01c48b12086be64bf61cd0f4a76b7da
hash_basis: working-tree bytes (LF)
---

# THM-2263 -- the parity spectrum behind the 1:169 pair

Put

```text
D_s={x in R/Z:||sx||<1/14},
rho(s,t)=measure(D_s intersection D_t).               (1)
```

THM-2255 proves the global ramified cap

```text
nu_13(s)!=nu_13(t)  implies  rho(s,t)<=25/1183,       (2)
```

with equality at reduced ratio `1:169`. The present theorem resolves the
whole overlap interval at each prescribed valuation gap. The result is not
another finite-cell accident: it comes from a two-endpoint spectrum of the
folded defect.

## 1. Exact overlap interval at gap d

Assume, after possibly interchanging the two speeds,

```text
nu_13(t)-nu_13(s)=d>=1.                               (3)
```

Divide by `gcd(s,t)`. There are positive integers `a,k` such that

```text
s=a,                 t=13^d k,
gcd(a,k)=1,          13 does not divide ak.           (4)
```

Write `p=13^d`. Then the sharp bounds are:

```text
d even:
  1/49-5/(588p) <= rho(s,t) <= 1/49+6/(49p);

d odd:
  1/49-6/(49p) <= rho(s,t) <= 1/49+5/(588p).         (5)
```

The equality loci in reduced oriented coordinates `(a,k)` are:

```text
d even:
  upper equality  (a,k)=(1,1),
  lower equality  (a,k)=(1,12) or (12,1);

d odd:
  lower equality  (a,k)=(1,1),
  upper equality  (a,k)=(1,12) or (12,1).            (6)
```

Thus even gaps have upper reduced ratio `1:p`; odd upper equality consists
of the unordered ratio families

```text
{1,12p},                  {12,p}.                    (7)
```

In particular,

```text
d=1:
  1/91 <= rho <= 23/1092;

d=2:
  289/14196 <= rho <= 25/1183.                       (8)
```

Equation (8) both recovers and explains THM-2255's global maximum: the
largest positive ramified defect occurs at the smallest even gap.

### Proof of the gap law

THM-1166 gives

```text
rho(a,pk)
 =1/49+
  [F((a+pk) mod 14)-F((pk-a) mod 14)]/(196apk),

F(r)=r(14-r),                         0<=r<14.        (9)
```

The folded function is even modulo fourteen. Since

```text
13^d = 1 mod 14       if d is even,
13^d =-1 mod 14       if d is odd,                   (10)
```

define

```text
Delta(a,k)=F((a+k) mod 14)-F((k-a) mod 14).          (11)
```

Then the correction in (9) is respectively

```text
 Delta(a,k)/(196pak),          d even,
-Delta(a,k)/(196pak),          d odd.                (12)
```

The complete normalized folded spectrum is

```text
-5/3 <= Delta(a,k)/(ak) <= 24.                       (13)
```

For the upper endpoint, `|Delta|<=49`. If `ak>=3`, then

```text
Delta/(ak)<=49/3<24.
```

The only coprime cells with `ak<3` are `(1,1),(1,2),(2,1)`;
their quotients are `24,10,10`. Hence the upper endpoint is unique at
`(1,1)`.

For the lower endpoint, if `ak>=30`, then

```text
Delta/(ak)>=-49/30>-5/3.
```

The remaining finite universe

```text
ak<=29,  gcd(a,k)=1,  13 does not divide ak           (14)
```

has `75` ordered cells. Exact evaluation gives unique ordered minima

```text
(a,k)=(1,12),(12,1),            Delta/(ak)=-5/3.     (15)
```

This proves (13). Substitution into (12) proves (5)--(7).

## 2. The 150 strict profiles: a compatible sharp pair ledger

Retain THM-2255's positive residual and exclusive owner strata:

```text
A_0=C_H minus union_(i=1)^5 D_(q_i),

E_j=A_0 intersection D_(c_j)
     minus union_(k!=j)D_(c_k),

measure(A_0)>=delta_5:=961/6930.                     (16)
```

The owner-multiplicity identity gives

```text
sum_j measure(E_j)
 >=delta_5-sum_(i<j)rho(c_i,c_j).                   (17)
```

For a strict remaining profile write

```text
(lambda_1,lambda_2,lambda_3)=(1,b,c),
2<=b<c,                         5<=c<=19.            (18)
```

Define the positive gap correction

```text
g(d)=
  6/(49*13^d),                 d even,
  5/(588*13^d),                d odd.                (19)
```

Equations (5) and (17) give

```text
sum_j measure(E_j)
 >=delta_5-3/49
   -g(b-1)-g(c-b)-g(c-1).                            (20)
```

The exact maximum of the three correction terms over (18) is unique:

```text
(b-1,c-b)=(2,2),

max [g(b-1)+g(c-b)+g(c-1)]
 =2g(2)+g(4).                                       (21)
```

One can see (21) without scanning 150 rows. If both gaps are even, each is
at least two and their sum is at least four, giving the displayed bound
with equality only at `(2,2)`. If both are odd, their sum constraint forces
the two smallest possible gaps to be `1,3`; direct cross-multiplication
gives `g(1)+g(3)<2g(2)`. In the mixed case the sum is at least five and
either `g(2)+g(3)+g(5)` or `g(1)+g(4)+g(5)` bounds the correction; both are
strictly below `2g(2)+g(4)`.

Consequently the unique worst strict profile is

```text
(lambda_1,lambda_2,lambda_3)=(1,3,5),               (22)
```

and its pair-cap sum is

```text
2 rho_gap(2)+rho_gap(4)
 =12531/199927.                                      (23)
```

Thus every strict profile satisfies

```text
sum_j measure(E_j)
 >=15041431/197927730
  =88159/1171170+144/199927.                         (24)
```

Some labelled owner therefore has mass at least

```text
15041431/593783190.                                  (25)
```

The bound is sharp for the simultaneous pair ledger. At the worst profile,
all three pair equalities are compatible on the common-core ladder

```text
(c_1,c_2,c_3)=g(1,169,28561).                        (26)
```

This is a stronger boundary statement than three separately sharp pair
caps: the equality geometries actually coexist.

## 3. The 15 repeated-first profiles: compatibility matters

Now let

```text
(lambda_1,lambda_2,lambda_3)=(1,1,c),
5<=c<=19,                    d=c-1.                  (27)
```

THM-2255 bounded the same-valuation shallow pair by `1/14` and each cross
pair independently by `25/1183`. Those three independent maxima cannot
coexist. The compatibility correction is small but exact.

For distinct reduced speeds, the universal overlap maximum `1/14` occurs
only at ratio `1:2`. If the ratio is not `1:2`, then

```text
rho<=1/21.                                           (28)
```

Indeed (9) gives

```text
rho<=1/49+1/(4ab).
```

For `ab>=10` this is below `1/21`; the coprime bank `ab<=9` has second
maximum `1/21`, at ratios `1:3` and `2:3`.

Suppose the shallow ratio is `1:2`. After common scaling the three speeds
have the form

```text
(g,2g,pu),                 gcd(g,u)=1.               (29)
```

For even `d`, the sum of the two cross-pair corrections above `2/49` is at
most

```text
17/(98p),                                           (30)
```

with equality at normalized `(g,u)=(1,1)` or `(1,2)`. For odd `d`, the
corresponding maximum is

```text
8/(539p),                                           (31)
```

with equality at `(g,u)=(1,11)` or `(11,2)`.

Here is a short proof. Let

```text
q(a,k)=Delta(a,k)/(ak).
```

If `u` is odd, `196p` times the cross correction in (29) is

```text
q(g,u)+q(2g,u);
```

if `u` is even, it is

```text
q(g,u)+q(g,u/2).                                     (32)
```

For the positive endpoint, `|Delta|<=49` reduces the odd-`u` branch to
`gu<=2` and the even-`u` branch to `gu<=4`; exact evaluation gives maximum
`34` at `(1,1),(1,2)`. For the negative endpoint the same tail reduces to
`gu<=25` and `gu<=50`, respectively; exact evaluation gives minimum
`-32/11` at `(1,11),(11,2)`. Odd valuation gaps negate the folded
correction by (12), proving (30)--(31).

For every allowed `d=4,...,18`, the ratio-`1:2` cap exceeds the
non-`1:2` upper bound from (28), and the unique largest value occurs at
`d=4`. It is

```text
1/14+2/49+17/(98*13^4)
 =3206/28561.                                        (33)
```

The equality triples are the common dilates of

```text
(1,2,13^4),                    (1,2,2*13^4).         (34)
```

Equations (17) and (33) prove, uniformly over all repeated-first profiles,

```text
sum_j measure(E_j)
 >=5229541/197927730
  =14627/585585+577/399854.                          (35)
```

Some labelled owner has mass at least one third of (35).

## 4. Expiration images and the exact stopping boundary

THM-2255 proves that an exclusive stratum `E_j` attached to an owner of
valuation `lambda_j` satisfies

```text
measure(T^(lambda_j+1)(E_j))
 >=(169/20)measure(E_j),             T(x)=13x.       (36)
```

Combining (24)--(25) with (36), every strict profile has a labelled
expiration image of mass at least

```text
15041431/70270200
 =1/7+5002831/70270200.                              (37)
```

For the repeated-first profiles, (35)--(36) give

```text
5229541/70270200
 =1/14+210241/70270200
 <1/7.                                               (38)
```

Thus the strict branch remains above the one-comb threshold, and the
repeated branch now crosses the half-comb threshold. Neither numerical
crossing closes a row.

The first failed implication is unchanged from THM-2255. There is no proved
post-expiration inclusion putting the image in a named successor union.
The exact local image can be large while its support is redistributed
among unit-mask/carry sheets. THM-2261 records the corresponding
surjectivity/no-one-core stopping boundary.

What THM-2263 adds is the sharp arithmetic address of the obstruction:

```text
global ramified cap             -> even gap two;
strict worst compatible ladder  -> gaps (2,2,4);
repeated worst compatible triple-> (1,2,13^4);
missing coordinate              -> post-expiration sheet/carry target. (39)
```

No scalar profile is excluded here, and LRC(14) remains open.

## 5. Verification

Run

```bash
python3 04-computation/lrc14_gap_sensitive_owner_spectrum_thm2263.py
python3 -O 04-computation/lrc14_gap_sensitive_owner_spectrum_thm2263.py

python3 04-computation/lrc14_gap_sensitive_owner_spectrum_referee_thm2263.py
python3 -O 04-computation/lrc14_gap_sensitive_owner_spectrum_referee_thm2263.py
```

The primary companion checks the `75`-cell folded endpoint bank, every gap
used by the `165` profiles, the exact strict and repeated ledgers, all
fractions, the same-valuation second maximum, the cross-score tails, and
both compatible equality carriers.

The referee does not call the folded formula. It constructs the danger
combs as exact rational interval unions and intersects them directly. It
checks eight parity endpoints through gap four, the small same-valuation
bank, and both simultaneous equality triples. Normal and optimized
transcripts match their stored outputs byte for byte. QED.
