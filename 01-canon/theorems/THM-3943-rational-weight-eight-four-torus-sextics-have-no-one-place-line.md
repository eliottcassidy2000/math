---
id: THM-3943
title: "Rational simple weight-eight four-torus sextics have no one-place line"
status: >
  PROVED + CITED CLASSIFICATION + VERIFIED-EXACT + INDEPENDENTLY
  HOSTILE-AUDITED. Over C, Degtyarev's classification says that simple
  irreducible weight-eight sextics have seven singularity
  configurations and four torus structures. Exactly three configurations are
  rational: E6+A5+4A2, 2A5+4A2, and A5+6A2+A1. In Degtyarev's complete
  trigonal-section models, no projective line pulls back to a divisor
  supported at one normalization point. The rigid row is killed by a missing
  S5T coefficient in its binary-sextic line system; the two moving rows are
  killed by exact quadratic-normalization norm equations, including all
  finite and infinity seams. Thus multiple torus structures do not supply a
  one-place affine Cardano branch in any complex rational simple weight-eight
  family.
  The non-simple J2,3+3A2 and J2,0+4A2 families are not classified here.
source: jc_zero_debt_lift / Degtyarev--Oka weight-eight sextic reframe, 2026-08-24
audit: >
  The 35-gate assertion-free companion verifies
  Degtyarev's trigonal parametrization, the rationality ledger, all three
  double-cover factorizations, the rigid basepoint-free binary sextics and
  pure-power rank tests, both moving-family line norms, their coefficient
  contradictions, conic degenerations, and finite/infinity address packets.
  The cited configuration and four-torus classification is kept separate from
  the repo-derived exact obstruction. An independent hostile audit
  reconstructed the line norms, checked multiplicities both off and on the
  conic branch locus, audited both points over infinity and the gamma=0
  seams, verified every excluded endpoint against the primary section atlas,
  and reproduced the normal/-O/frozen output and all hashes.
related:
  - THM-3879-rational-torus-sextic-c3-packet-one-place-tradeoff
  - THM-3925-fivefold-conic-contact-torus-sextic-one-place-fold
  - THM-3928-split-affine-conic-one-place-fold-degree-barrier
  - THM-3932-infinity-component-linear-conic-torus-sextic-fold-classification
  - THM-3945-nonsimple-weight-eight-j-sextics-have-no-one-place-line
script: 04-computation/jc2_rational_weight8_four_torus_sextic_line_thm3943.py
output: 05-knowledge/results/jc2_rational_weight8_four_torus_sextic_line_thm3943.out
script_sha256: b7ffa43f5f0ef30f2e4af5280e1a892a1164b7b54ecb3d676fd0c66bfca772d5
output_sha256: 70047c45202139904d65b95f516165502a7181a82260c89748a82729dd6fd673
semantic_sha256: 6df55931aa9ef0f9dde3eb4d6aaf92028f618a8acf0014f10f835b81f63837fd
hash_basis: raw LF bytes
---

# THM-3943 -- four torus structures still cannot buy one affine place

**PROVED + CITED CLASSIFICATION + VERIFIED-EXACT + INDEPENDENTLY
HOSTILE-AUDITED.** The exhaustive classification statement is over
`k=C`, the field of Degtyarev's cited complex-plane classification.  The
displayed parametrizations, factorizations, and line-support contradictions
below are algebraic and remain valid after base change to any algebraically
closed characteristic-zero field, but no such base-change argument for the
exhaustiveness of the classified families is asserted here.

The imported input is [Degtyarev, *Irreducible plane sextics with large
fundamental groups*, §§1.1 and 3.2--3.6](https://arxiv.org/abs/0712.2290).
It classifies the simple irreducible weight-eight sextics into the seven
configurations

```text
E6+A5+4A2,  E6+6A2,  2A5+4A2,
A5+6A2+A1,  A5+6A2,  8A2+A1,  8A2,                      (1)
```

and constructs their four torus structures from a single four-cuspidal
trigonal quotient. The classification and the number four are **CITED**, not
reproved by the companion.

A plane sextic has arithmetic genus ten. Using

```text
delta(E6)=delta(A5)=3,       delta(A2)=delta(A1)=1,        (2)
```

the rational configurations in `(1)` are exactly

```text
E6+A5+4A2,        2A5+4A2,        A5+6A2+A1.             (3)
```

The theorem proves that for every complex curve in each family `(3)`, and
every projective line `ell`, the pullback to the normalization has at least two
support points. Equivalently, no affine chart cut out by a line makes the
normalization an affine line. Since a line section has degree six, the exact
test is

```text
nu^*ell has one support  iff  its binary sextic is lambda*L^6. (4)
```

This is a curve-and-line obstruction. It does not depend on which of the four
torus structures is selected.

## 1. The universal four-cuspidal quotient

Use affine coordinates `(x,v)` on the Hirzebruch quotient. Degtyarev's
four-cuspidal trigonal curve is

```text
f(x,v)=4v^3-(24x^3+3)v+8x^6+20x^3-1=0                  (5)
```

with normalization

```text
x=3u/(u^3+2),
v=-(u^6-20u^3-8)/(2(u^3+2)^2).                          (6)
```

For a quadratic section `v=s(x)`, the corresponding plane sextic is

```text
B_s: f(x,y^2+s(x))=0.                                    (7)
```

Thus its normalization is the quadratic cover of the `u`-line

```text
y^2=v(u)-s(x(u)).                                         (8)
```

The special sections giving `(3)` make `(8)` rational. More importantly,
their factorizations retain enough information to audit every line.

## 2. The rigid E6+A5+4A2 row

For the ordered tangent-cusp pair used in Degtyarev's model,

```text
s(x)=2x^2-1/2.                                            (9)
```

Substitution into `(8)` gives

```text
v(u)-s(x(u))=6(u-1)^2(2u+1)/(u^3+2)^2.                  (10)
```

Put `2u+1=w^2`. After a harmless nonzero scaling of coordinates, the
projective normalization is

```text
[W:T] |-> [X:Y:Z],
X=12(W^2-T^2)T^4,
Y=4 sqrt(6) W(W^2-3T^2)T^3,
Z=W^6-3W^4T^2+3W^2T^4+15T^6.                            (11)
```

The three binary sextics in `(11)` have no common zero. Every line pullback

```text
alpha X+beta Y+gamma Z                                  (12)
```

has zero coefficient at `W^5T`. If `(12)=lambda(pW+qT)^6`, then

```text
6 lambda p^5q=0.                                         (13)
```

The line is nonzero and the map is nondegenerate, so `lambda!=0`. Hence
`p=0` or `q=0`.

If `p=0`, `(12)` would be a pure `T^6`. Its odd coefficients first force
`beta=0`; its `W^6` coefficient then forces `gamma=0`; and its `W^2T^4`
coefficient forces `alpha=0`. If `q=0`, the same odd coefficients force
`beta=0`, the `W^4T^2` coefficient forces `gamma=0`, and `W^2T^4` forces
`alpha=0`. Both endpoints give the zero line. This proves the rigid no-go.

## 3. The 2A5+4A2 family

The complete section line through the two selected quotient cusps is

```text
s_a(x)=a x^2+(2-a)x-1/2.                                 (14)
```

Equation `(8)` factors as

```text
v(u)-s_a(x(u))
 =3(u-1)^2 Q_a(u)/(u^3+2)^2,
Q_a=(a-2)u^2+2au+2.                                      (15)
```

Thus, writing `z^2=3Q_a`, the normalization map is

```text
[X:Y:Z]=[3u:(u-1)z:u^3+2].                               (16)
```

The parameters realizing the advertised irreducible configuration form the
open subset specified by Degtyarev. The visible degeneration ledger is

```text
disc(Q_a)=4(a^2-2a+4);
Q_a(1)=3a;                   lc(Q_a)=a-2.                 (17)
```

The square-discriminant seam splits the quadratic pullback.  The endpoint
`a=2` is literally the rigid row `(9)`; the uniform norm proof below also
includes `a=0`, so no appeal to its cusp-swapped rigid interpretation is
needed.

Let `ell=alpha X+beta Y+gamma Z`. First suppose `gamma!=0` and normalize it
to one. The norm down to `k(u)` of its numerator is

```text
N_a=(u^3+2+3alpha u)^2-B(u-1)^2Q_a,      B=3beta^2.      (18)
```

If all six zeros lie at a finite normalization point over `u=r`, then

```text
N_a=(u-r)^6.                                              (19)
```

This remains correct when the point is fixed by the conic involution: local
order doubles upstairs and `u-r` has order two. Comparing coefficients in
`(18)--(19)` gives

```text
[u^5]: r=0,
[u^3]: 4(1-B)=0,
[1]:   4-2B=0.                                           (20)
```

The middle equation gives `B=1`, while the last then reads `2=0`.

If `gamma=0` and `a!=2`, the smooth conic closure has two points at infinity.
Homogenizing `(16)` to degree three gives

```text
[3UW^2:(U-W)VW:U^3+2W^3].                                (21)
```

Both points `W=0` map to `[0:0:1]`, and every `gamma=0` line contains both
addresses. The remaining value `a=2` is the rigid row. This closes the whole
`2A5+4A2` family, including its line-at-infinity seam.

## 4. The A5+6A2+A1 family

Let `tau` be the normalization parameter of the selected smooth tangency
point. The section tangent there and passing through the selected cusp is

```text
s_tau(x)=a_tau x^2+b_tau x+c_tau,
a_tau=2(tau^3-3tau-1)/[(tau-1)(tau+2)^2],
b_tau=6tau(tau+1)/[(tau-1)(tau+2)^2],
c_tau=-(tau^3+3tau^2+8)/[2(tau-1)(tau+2)^2].              (22)
```

The advertised smooth-tangency row excludes the cusp collisions and the
singular denominators. The exact radicand is

```text
v(u)-s_tau(x(u))
=6(tau-u)^2(u-1)^2 Q_tau(u)
 /[(tau-1)(tau+2)^2(u^3+2)^2],
Q_tau=u^2+2(tau+1)u+tau+3,
disc(Q_tau)=4(tau-1)(tau+2).                              (23)
```

After absorbing the nonzero scalar in `(23)` into the second plane
coordinate and writing `z^2=Q_tau`, the normalization map is

```text
[X:Y:Z]=[3u:(u-tau)(u-1)z:u^3+2].                        (24)
```

### 4.1 Finite support

For a line with `gamma!=0`, normalize `gamma=1`. Its norm is

```text
N_tau=(u^3+2+3alpha u)^2
      -B(u-tau)^2(u-1)^2Q_tau.                            (25)
```

If its sole support is finite, write

```text
N_tau=lambda(u-r)^6,             lambda!=0.               (26)
```

The `u^5` coefficient gives `r=0`. The difference between half the `u^3`
equation and the constant equation is `2(B-1)=0`; hence `B=1`. But the
leading coefficient in `(25)` is `lambda=1-B=0`, contradicting `(26)`.

If `gamma=0` and `beta!=0`, the norm is

```text
9alpha^2u^2-B(u-tau)^2(u-1)^2Q_tau.                      (27)
```

Its leading coefficient is nonzero and its quintic coefficient is zero, so
a finite pure sixth power again has `r=0`. The `u^4`, `u^3`, and constant
conditions, divided by the nonzero `B`, are

```text
3tau(tau+1)=0,
-tau^3-3tau^2+2=0,
-tau^2(tau+3)=0.                                         (28)
```

The first leaves `tau=0,-1`; the second removes zero, and the third removes
`-1`. If also `beta=0`, the line is `X=0`, which contains both infinity
addresses described next.

### 4.2 Infinity support

Homogenizing `(24)` gives

```text
[3UW^2:(U-tau W)(U-W)V:U^3+2W^3],
V^2=U^2+2(tau+1)UW+(tau+3)W^2.                           (29)
```

At `W=0` there are two points `V=+U,-U`, mapping to
`[0:+1:1]` and `[0:-1:1]` after scaling. A `gamma=0,beta!=0` line misses both;
`beta=0` gives `X=0` and contains both. Thus only a nonzero-`gamma` line could
have sole support at exactly one infinity point.

For such a line, absence of finite zeros makes `(25)` a nonzero constant.
The degree-six coefficient gives `B=1`; the degree-four coefficient gives
`2alpha+tau^2+tau=0`; and the degree-three coefficient becomes

```text
-2(tau-1)(tau+2)^2=0.                                    (30)
```

These are precisely the excluded parameters where `(22)--(23)` degenerates.
Hence the infinity seam is empty on the advertised family.

## 5. Consequence and exact boundary

Sections 2--4 exhaust the rational **simple** configurations in the cited
weight-eight list. Therefore no rational simple weight-eight sextic admits a
projective line whose pullback is a sixth power. In particular, having four
torus structures does not improve THM-3879's two-place barrier to one place;
the extra Cardano presentations do not alter the line linear system.

The two non-simple weight-eight configurations

```text
J2,3+3A2,             J2,0+4A2                           (31)
```

are outside the cited simple list and are not dependencies of this theorem.
They are closed separately by
[THM-3945](THM-3945-nonsimple-weight-eight-j-sextics-have-no-one-place-line.md).
No claim is made about non-line boundaries, birational target changes, a
Keller atlas, or JC(2).

Run

```bash
python3 04-computation/jc2_rational_weight8_four_torus_sextic_line_thm3943.py
python3 -O 04-computation/jc2_rational_weight8_four_torus_sextic_line_thm3943.py
```

The companion independently replays all displayed identities and hostile
endpoint controls; it does not use Python assertions.
