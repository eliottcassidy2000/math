---
id: THM-4101
title: "AP7 weight-seven gap-four second-moment absorption"
status: >
  PROVED + PROVED RELATIVE TO THM-2072/4092 + VERIFIED-EXACT +
  INDEPENDENTLY VERIFIED-EXACT. Four outliers of antipodal weight seven over
  the exact AP7 safe interval are absorbable when three are odd, one is even,
  an odd pair is {u,u+4}, and an explicit reciprocal tax is smaller than a
  forced pair-overlap credit. Uniformly, minimum outlier speed 264 suffices.
  This closes the associated eleven-core dyadic two-tail seams through common
  dilation. It is neither arbitrary weight-seven absorption nor LRC(14).
source: codex-arithmetic-boundary-breakthrough-20260825
depends_on:
  - THM-2072-fixed-owner-clock-bank-no-go-and-half-shift-certificate
  - THM-4092-parity-weighted-antipodal-k-comb-density-compiler
related:
  - THM-1155-threespeed-exhaustive-and-ceiling
  - THM-2061-lrc14-dyadic-two-tail-folded-seam
  - THM-4079-lrc14-antipodal-outlier-absorption-and-adaptive-clock
  - THM-4098-weight-seven-antipodal-scale-escape-and-missing-parity-rows
script: 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101.py
output: 05-knowledge/results/lrc_ap7_weight_seven_gap4_second_moment_thm4101.out
independent_audit_script: 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc_ap7_weight_seven_gap4_second_moment_thm4101_independent_audit.out
script_sha256: 9937f686321c7950fb6c9c907946506aaf37d7c9e7f118c5942074d7a92b74c0
output_sha256: ec138ce28df0e921d45c462ef19d81308e0725863e85479b0f340033303bc038
independent_audit_script_sha256: 34cd587c7bb1dad0a50d0d624f06ef21d64479952340c42f65dcb2060c332745
independent_audit_output_sha256: 06703acf7534a918813f7bcf863155442e03f8011e1ff3778d074621a171e743
hash_basis: raw LF bytes
---

# THM-4101 -- AP7 weight-seven gap-four absorption

Put `delta=1/14` and retain THM-4092's exact two-phase AP7 interval

```text
D={1,...,7},
J=[4/35,13/98],                  L=|J|=9/490,
J subset G_D^+-.                                          (1)
```

For a speed `v`, let

```text
U_v={theta:min(||v theta||,||v(theta+1/2)||)<delta},
omega(v)=1 if v is even, and 2 if v is odd.                (2)
```

THM-4092 proves

```text
mu(J intersect U_v)
 <= omega(v)L/7+(7-omega(v))/(49v).                       (3)
```

When four outliers have total weight seven, summing `(3)` leaves no positive
multiple of `L`. The missing coordinate is pair overlap. A gap of four is
resonant with the unique internal split `theta=1/8` of the AP7 interval and
forces enough overlap to restore a strict certificate.

## 1. A four-comb second-moment inequality

Let `V` contain four speeds and put, pointwise on `J`,

```text
F(theta)=sum_(v in V) 1_(U_v)(theta).                      (4)
```

For every integer `r` from zero through four,

```text
1_(r>0) <= r-(1/2) binom(r,2).                             (5)
```

The right sides are respectively

```text
0, 1, 3/2, 3/2, 1.                                       (6)
```

Integrating `(5)` at `r=F(theta)` gives

```text
mu(J intersect union_(v in V) U_v)
 <= S_1-(1/2)S_2,                                         (7)

S_1=sum_v mu(J intersect U_v),
S_2=sum_{v<w} mu(J intersect U_v intersect U_w).
```

This inequality is specific to at most four combs: at multiplicity four its
right side returns to one, while at multiplicity five it is zero and the
majorization already fails. That arity boundary is load-bearing.

## 2. The gap-four overlap packet

Let `u` be odd and put `n=2u`. On the two pieces

```text
J_-=[4/35,1/8],       |J_-|=3/280,
J_+=[1/8,13/98],      |J_+|=3/392,                        (8)
```

write

```text
epsilon(theta)=8theta-1.
```

Then

```text
epsilon(J_-)=[-3/35,0],
epsilon(J_+)=[0,3/49].                                    (9)
```

Because both `u` and `u+4` are odd, their two-phase danger tests in the
`n theta` coordinate are

```text
U_u:     ||n theta||<1/7,
U_(u+4): ||n theta+epsilon(theta)||<1/7.                  (10)
```

Consequently the following fixed arcs are contained in the pair overlap:

```text
n theta mod 1 in (-2/35,1/7)   on J_-,   duty rho_-=1/5,
n theta mod 1 in (-1/7,4/49)   on J_+,   duty rho_+=11/49. (11)
```

For example, on `J_-` the first arc lies in `(-1/7,1/7)`, and adding any
epsilon from `[-3/35,0]` still lies strictly between `-1/7` and `1/7`.
The `J_+` verification is identical using the second line of `(9)`.

We need the lower, rather than upper, interval-discrepancy direction.

> **Lemma 2.1 (periodic-arc lower bound).** If a period-`1/n` set occupies one
> arc of duty `rho` in each period, then on every interval `I` of length
> `ell`,
>
> ```text
> mu(I intersect A)>=rho ell-rho(1-rho)/n.                (12)
> ```

Write `ell=m/n+r`, `0<=r<1/n`. Complete periods contribute exactly their
duty. In the residual interval the least occupied length is

```text
max(0,r-(1-rho)/n).
```

Its deficit from `rho r` is maximized at `r=(1-rho)/n` and equals
`rho(1-rho)/n`, proving `(12)`.

Apply `(12)` separately to `(11)`, recalling `n=2u`. We obtain

```text
boxed:
mu(J intersect U_u intersect U_(u+4)) >= P-E/u,           (13)

P=(1/5)(3/280)+(11/49)(3/392)=927/240100,
E=((1/5)(4/5)+(11/49)(38/49))/2=10027/60025.
```

This is a literal containment proof, not an asymptotic decorrelation claim.

## 3. Exact absorption criterion

Let `V` consist of four distinct speeds, exactly one even and three odd, and
suppose

```text
{u,u+4} subset V,                  u odd.                 (14)
```

Thus

```text
sum_(v in V) omega(v)=1+2+2+2=7.                          (15)
```

Define the exact first-moment boundary tax

```text
T(V)=sum_(v in V) (7-omega(v))/(49v).                     (16)
```

Summing `(3)` and using `(15)` gives `S_1<=L+T(V)`. In `(7)`, retain only
the selected pair from `(13)`. Therefore

```text
mu(J intersect union_(v in V) U_v)
 <= L+T(V)-P/2+E/(2u).                                   (17)
```

We have proved the primary criterion:

> **Theorem 3.1 (gap-four absorption).** Under `(14)--(16)`, if
>
> ```text
> boxed: T(V)+10027/(120050u) < 927/480200,               (18)
> ```
>
> then
>
> ```text
> J minus union_(v in V) U_v is nonempty.                 (19)
> ```

Every survivor in `(19)` is safe at both antipodal phases for all eleven
speeds `D union V`.

## 4. The uniform height 264 corollary

Assume in addition

```text
min(V)>=264.                                               (20)
```

Let `e` be the even speed and `w` the third odd speed outside the selected
pair. Then `e>=264` and `u>=265`.

If `u=265`, distinctness forces `w>=267`; monotonicity reduces `(18)` to

```text
6/(49*264)+(5/49)(1/265+1/267+1/269)+E/(2*265) < P/2.
```

Its exact margin is

```text
489784241/100536614409000>0.                              (21)
```

If `u>=267`, the worst remaining placement is

```text
e=264,       w=265,       {u,u+4}={267,271}.
```

The exact margin is

```text
203219143/20256819706200>0.                               (22)
```

Thus `(18)` holds in every case:

```text
boxed:
min(V)>=264 and (14) imply G_({1,...,7} union V)^+- nonempty. (23)
```

The integer `264` is a uniform sufficient height for this mechanism, not a
claimed sharp global threshold. At the adjacent floor, the bank

```text
{263,264,265,267}
```

misses criterion `(18)` but has an actual survivor. It is a criterion near
miss, not a counterexample and not evidence that `264` is optimal.

## 5. Endpoint clock and seam consequence

The survivor set in `(19)` has positive measure by the strict inequality
`(17)`, hence contains a nondegenerate closed component. An endpoint is an
endpoint of `J` or of an outlier tooth, so it has an even presentation

```text
theta=r/N,             N<=max(98,14 max(V)).               (24)
```

For the boundary-positive bank

```text
V_0={264,265,267,269},                                    (25)
```

the independent wall audit gives

```text
survivor measure =128698793/23448773040,
theta=214/1869,        N=3738.                            (26)
```

The two phases `theta,theta+1/2` are both in the weak-safe set of the
eleven-speed core `C={1,...,7} union V`. THM-2072 therefore closes the strict
dyadic seam: for all distinct positive odd `x,y`,

```text
boxed: 2C union {x,y} is 1/14-lonely.                     (27)
```

As in THM-4098, the certificate transports through every common dilation
`d>=1`: at `theta/d`, an odd `d` preserves the half phase and an even `d`
collapses it to the safe zero phase. Hence

```text
boxed: 2dC union {x,y} is 1/14-lonely.                    (28)
```

## 6. Failure controls and exact audits

The THM-4092 bank

```text
{8,9,11,13}                                               (29)
```

contains the gap-four pair `{9,13}` but literally covers `J`. It fails
`(18)`, showing that the resonance alone is insufficient; the reciprocal
boundary tax is essential.

Primary reproduction:

```bash
python3 -B 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101.py
python3 -B -O 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101.py
```

The primary exact referee replays `997` odd gap-four pairs `9<=u<=2001`, the
periodic lower bound, literal containment, both monotonicity boundary cases,
the hostile `(29)`, two positive controls, and the height-263 criterion near
miss.

Independent reproduction:

```bash
python3 -B 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101_independent_audit.py
python3 -B -O 04-computation/lrc_ap7_weight_seven_gap4_second_moment_thm4101_independent_audit.py
```

The independent referee never merges comb intervals. It constructs all
literal original-theta walls and integrates predicates cell by cell, checking
`497` pairs, `94,715` containment points, the hostile walls, and the endpoint
in `(26)`. Both normal/optimized pairs byte-match their frozen outputs. All
arithmetic is rational, all gates remain active under optimization, and no
script uses floating point.

The connection ledger is

```text
source:       four literal antipodal danger combs on the AP7 interval
observable:   multiplicity F(theta), not merely union support
operation:    integrate a pointwise second-moment majorant
preserved:    exact union upper bound and strict wall conventions
destroyed:    all pair labels except the selected gap-four pair
sidecar:      epsilon=8theta-1 and the split theta=1/8
recovery:     two frozen phase arcs plus lower interval discrepancy
hostile:      {8,9,11,13}
target:       an antipodal safe pair and adaptive even owner clock.          (30)
```

This theorem handles only AP7 banks with the declared pair resonance. It does
not absorb every weight-seven bank, improve THM-4098's all-shape dilation
threshold, produce an arbitrary-core supplier, or prove LRC(14). **QED.**
