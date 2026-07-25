---
id: THM-2232
title: "Same-core signed-eigen Markov dual exclusion"
status: >
  PROVED + VERIFIED-EXACT + HOSTILE-AUDITED. In either high-first scalar
  residue (4,6,8) or (5,7,9), the three normalized blocker coefficients
  cannot share one thirteen-unit core. More generally, an exact
  signed-transfer dual gives a quantitative event-level and atom-level
  stability obstruction to approaching any common-core model. The unsigned
  common-core three-checkpoint mass is 916159/4826809, above the residual
  floor, but the signed capacities are respectively about 0.045006 and
  0.044771, both far below 961/6930. This is a conditional equality-carrier
  exclusion and stability theorem; it closes no unrestricted valuation
  profile and does not prove LRC(14).
source: audit-bellman-frontier-2026-07-25-same-core-signed-dual
depends_on:
  - THM-2198-scalar-five-plus-three-image-pump-and-first-depth-exclusion
  - THM-2222-scalar-transfer-parity-tower-and-four-checkpoint-survivor-reduction
related:
  - THM-2218-labelled-guard-hole-fourier-and-signed-lift-energy
  - THM-2227-sharp-parity-three-checkpoint-bellman-profile-exclusion
  - THM-2229-unit-time-positive-set-bellman-profile-exclusion
script: 04-computation/lrc14_same_core_signed_eigen_dual_thm2232.py
output: 05-knowledge/results/lrc14_same_core_signed_eigen_dual_thm2232.out
script_sha256: acd1199d9f860762cdd7b81b2ce317ebca1208787bd70c451edfd20fadd55e55
output_sha256: e607373141e7a4744cc943b5556375404154a50728838df3876e70cd2f2436a4
hash_basis: working-tree bytes (LF)
---

# THM-2232 -- the same-core signed-eigen obstruction

THM-2229 leaves two high-first profiles:

```text
(4,6,8), (5,7,9).                                   (1)
```

Their common geometric-chain model passes the unsigned three-checkpoint
test. It nevertheless cannot coexist with the original signed five-unit
residual. This theorem isolates that incompatibility and records the precise
inverse/stability statement still needed for unrestricted cores.

## 1. Setup and scope

On `T=R/Z`, put

```text
D_b={x:||bx||<1/14},       C_H={x:||Hx||>1/7}.       (2)
```

Work in THM-2198's scalar five-unit/three-blocker branch. Thus

```text
R=1_(C_H)-sum_(i=1)^5 1_(D_(q_i)),                  (3)
```

where `H,q_1,...,q_5` are thirteen-units, and

```text
p:=integral R_+ >= delta_5:=961/6930.               (4)
```

Since `R` takes values in `{-5,...,1}`, its positive part is the indicator
of `A_+={R>0}` and

```text
0<=R_+<=1,        0<=R_-<=5.                        (5)
```

Also `measure(C_H)=5/7` and `measure(D_(q_i))=1/7`, so

```text
integral R=0,       integral R_+=integral R_-=p.     (6)
```

Fix `a in {4,5}` and write a general blocker triple with the relevant
profile as

```text
c_j=13^(a+2(j-1))u_j,       13 does not divide u_j,
j=1,2,3.                                             (7)
```

The same-core case means

```text
u_1=u_2=u_3=u.                                      (8)
```

The theorem proves that (8) is impossible. It does **not** infer (8) from
the valuation profile.

## 2. Positive and negative parity carriers

For a thirteen-unit `v`, define the stationary danger process

```text
X_t^(v)(x)=1_(D_v)(13^t x).                          (9)
```

For a triple `bold(u)=(u_1,u_2,u_3)`, let

```text
P_a(bold(u))
 = intersection_(0<=k<=a, k even)
     {OR_(j=1)^3 X_(a+2(j-1)-k)^(u_j)=1},

N_a(bold(u))
 = intersection_(1<=k<=a, k odd)
     {OR_(j=1)^3 X_(a+2(j-1)-k)^(u_j)=1}.           (10)
```

THM-2222's transfer-parity inclusions, applied at every available divided
blocker checkpoint, say

```text
support(R_+) subset P_a(bold(u)),
support(R_-) subset N_a(bold(u))                    (11)
```

almost everywhere. Thus the positive and negative parts of the *same*
residual must occupy the even and odd carriers respectively. Keeping both
signs is the new information absent from an unsigned checkpoint bound.

## 3. A general signed-transfer dual lemma

Let `S(x)=13x mod 1` and let `L` be its unnormalized transfer operator.
THM-2222 proves

```text
L^t R=(-1)^t R.                                     (12)
```

Transfer duality applied to (9) therefore gives

```text
integral R X_t^(v)
 =(-1/13)^t integral R X_0^(v).                     (13)
```

Put `rho=-1/13` and center every time coordinate along this eigenline:

```text
Y_t^(v)=X_t^(v)-rho^t X_0^(v).                      (14)
```

Equations (6) and (13) imply

```text
integral R q=0                                     (15)
```

for every finite affine combination

```text
q=c+sum_t alpha_t Y_t^(v).                          (16)
```

This produces the following reusable inequality. If `P,N` support `R_+`
and `R_-` as in (11), then

```text
p
 =integral [R_+(1-q)+R_-q]
 <=integral [
      1_P (1-q)_+ + 5 1_N q_+
   ].                                               (17)
```

The pointwise step is hostile-safe even when `q` lies outside `[0,1]`:
for `0<=A<=1_P` and `0<=B<=5 1_N`,

```text
A(1-q)<=1_P(1-q)_+,
Bq<=5 1_Nq_+.                                      (18)
```

No independence between `R`, the blocker carriers, or time coordinates is
used in (17).

## 4. Exact one-core Markov law

For one unit core, root counting gives

```text
P(X_t=1 | X_(t+1)=b)=(2-b)/13,       b in {0,1}.    (19)
```

The stationary atom law is `(6/7,1/7)`. Detailed balance turns (19) into
the same forward and backward transition matrix

```text
       next 0   next 1
0       11/13     2/13
1       12/13     1/13.                             (20)
```

The conditional in (19) depends on the terminal point only through its
danger bit. Backward root induction therefore factors every finite word,
so (20) is the exact Haar law of `(X_t^(v))`; it is not a Markov
approximation.

## 5. Two explicit dual certificates

In the common-core model abbreviate `X_t=X_t^(u)`,
`P_a=P_a(u,u,u)`, and `N_a=N_a(u,u,u)`. Use

```text
q_4=2-Y_1-2Y_3-2Y_5-Y_7

   =2-X_1-2X_3-2X_5-X_7
      -(4884270/62748517)X_0,                       (21)

q_5=2-Y_2-2Y_4-2Y_6-Y_8

   =2-X_2-2X_4-2X_6-X_8
      +(4884270/815730721)X_0.                      (22)
```

The displayed `X_0` corrections are load-bearing. Without them, the
coefficient of `integral R X_0` in `integral Rq` would be respectively

```text
4884270/62748517,       -4884270/815730721,          (23)
```

not zero.

Enumerating the exact Markov words through times eight and nine gives

```text
C_4
 :=E[1_(P_4)(1-q_4)_+ + 5 1_(N_4)(q_4)_+]
  =2303649491556761/51185893014090757
  =0.04500555438...,

C_5
 :=E[1_(P_5)(1-q_5)_+ + 5 1_(N_5)(q_5)_+]
  =2711005672365842843/60552911435669365531
  =0.04477085590....                                  (24)
```

Both inequalities are strict:

```text
delta_5-C_4
 =33225352210052863747/354718238587648946010 >0,

delta_5-C_5
 =5629154082883281339043/59947382321312671875690 >0. (25)
```

Equations (4), (17), and (24) contradict each other. Therefore neither
profile in (1) can occur with a common normalized blocker core. QED.

## 6. Why the unsigned geometric chain survives

For both values of `a`, stationarity and (20) give

```text
measure(P_a(u,u,u))
 =916159/4826809
 =delta_5+1710418421/33449786370
 =0.1898063503... > delta_5.                        (26)
```

This is exactly the three-window, no-`000` geometric-chain carrier already
seen in THM-2222. Thus the common core is not excluded by positive-set
occupancy. It is excluded because the original cover simultaneously
requires:

```text
positive residual mass on P_a,
negative residual mass on N_a,
and the signed correlations (13).                  (27)
```

The surviving unsigned object is therefore a false equality model for the
full signed cover, not evidence that the signed branch is sharp there.

## 7. Quantitative stability dichotomy

The dual also states exactly what an unrestricted survivor must destroy.
Fix any thirteen-unit reference core `v` and abbreviate

```text
P=P_a(bold(u)),       N=N_a(bold(u)),
P^v=P_a(v,v,v),       N^v=N_a(v,v,v).               (28)
```

Let

```text
M_a^+=||(1-q_a^(v))_+||_infinity,
M_a^-=||(q_a^(v))_+||_infinity.                    (29)
```

Comparing (17) with the common-core expectation yields the proved necessary
condition

```text
delta_5-C_a
 <=M_a^+ measure(P triangle P^v)
   +5M_a^- measure(N triangle N^v).                 (30)
```

The exact ranges of the duals are

```text
q_4 in [-255878338/62748517, 2],
q_5 in [-4, 1636345712/815730721],                  (31)
```

so

```text
(M_4^+,M_4^-)=(318626855/62748517, 2),
(M_5^+,M_5^-)=(5, 1636345712/815730721).            (32)
```

There is also an atom-level corollary. Put

```text
d(s,v)=measure(D_s triangle D_v).                   (33)
```

Boolean union/intersection and invariance of Haar measure under `S^t`
give

```text
measure(P triangle P^v)
 <=3 sum_(j=1)^3 d(u_j,v),

measure(N triangle N^v)
 <=o_a sum_(j=1)^3 d(u_j,v),                        (34)

o_4=2,       o_5=3.
```

Consequently every unrestricted survivor must obey, for every reference
unit core `v`,

```text
sum_j d(u_j,v)
 >=33225352210052863747/12497970889120926859650
   =0.0026584597...                    if a=4,

sum_j d(u_j,v)
 >=5629154082883281339043/2703016319465029128110550
   =0.0020825454...                    if a=5.       (35)
```

Equation (35) is the equality/stability dichotomy: exact or sufficiently
close common-core geometry is impossible. What remains open is an inverse
step showing that a general high-mass three-core Bellman survivor must
approach such a common-core model in the event metric of (30), or a
different signed certificate that works without alignment.

## 8. Connection and loss ledger

```text
source:
  the original signed five-unit residual and both parity carriers;

map:
  project one reference danger process onto the transfer eigenline,
  then use its centered time coordinates as an affine dual;

preserved:
  the sign of the residual, every even and odd checkpoint support,
  the exact -1/13 correlation phase, and the complete one-core word law;

destroyed:
  cross-core phase, carry and owner labels, relations among distinct
  normalized cores, and any inverse theorem from Bellman mass to alignment;

needed sidecar:
  an event-level alignment/stability theorem, a three-core Fourier/carry
  certificate, or a direct optimization over several centered eigenlines;

cheapest next hostile probes:
  triples with two equal cores, small multiplicative relations among the
  cores, and profiles saturating THM-2229 while violating (30) maximally.
                                                               (36)
```

## 9. Exact audit and scope

Run

```bash
python3 04-computation/lrc14_same_core_signed_eigen_dual_thm2232.py
python3 -O 04-computation/lrc14_same_core_signed_eigen_dual_thm2232.py
```

The companion uses exact rational arithmetic and checks:

```text
the stationary law, detailed balance, and reverse root law;
direct Fraction and independent integer-weight Markov evaluations;
all 4,536 and 8,512 endpoints of enlarged pointwise residual boxes;
the centered-correlation cancellation and hostile uncentered defects;
the exact capacities, margins, dual ranges, and stability constants;
the unsigned common-core mass above the target;
explicit witnesses showing that omitting positive parts in (17) fails;
ordinary/optimized/stored transcript identity.      (37)
```

This theorem excludes two **same-core subfamilies**, but zero unrestricted
valuation profiles. The current scalar ledger remains `240`, and LRC(14)
remains open.
