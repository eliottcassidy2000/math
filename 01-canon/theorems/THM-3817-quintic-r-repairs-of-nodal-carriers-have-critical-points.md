---
id: THM-3817
title: "Quintic R-repairs of nodal carriers have critical points"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the c=1
  cubic pseudo-plane, every canonical
  nodal carrier A=e^2-z/3+r sum_(i=0)^5 b_i e^i has a critical point.  In the
  genuine quintic case, a finite weight cover turns the top logarithmic-
  remainder equations into a triangular one-parameter resonance.  One later
  row splits into C_4=0 or a cubic S(C_4^2)=0; the first branch dies directly
  and the second dies by an exact T-resultant coprimality certificate.  This
  closes pure-r profiles through degree five, but not degree at least six or
  mixed corrections.
source: jc_sparse_direct_search / quintic logarithmic-divisibility lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (root, 2026-08-23).  The audit independently
  checked the multiplicity-safe boundary implication, the b_5-localized
  quotient, the surjective square-root cover, every Laurent weight and
  primitive back-substitution, the invariant universe and sole T!=0
  saturation, the four triangular pivots, the retained C_4=0 branch, and the
  T-resultant/coprimality contradiction on S(C_4^2)=0.  It also rechecked the
  finite-root reconstruction and all dangerous u-values.  The deterministic
  companion derives the
  universal Hamiltonian reduction and resultant, degree and leading
  coefficient, exact quotient denominator, surjective b_5=s^2 weight cover,
  twelve invariant remainder rows, a 149-element exact Groebner basis, four
  triangular normal forms, the row-sixteen branch factor, the C_4=0 exit,
  the T-resultant factorization, three exact coprimality certificates, hostile
  sparse/repeated-root controls, and finite source reconstruction.  Normal and
  optimized runs match the frozen transcript exactly after LF normalization,
  all raw-LF artifact hashes match, and the 1,448 active gates replay without
  assertions disabled.  A second hostile pass manually reconstructed the
  Sylvester division, weight rows, and post-basis exits; its stronger alternate
  F5B basis computation was boundedly stopped inside reduction and is recorded
  as INCOMPLETE, not used as proof evidence.
depends_on:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3813-quartic-r-repairs-of-nodal-carriers-have-critical-points
related:
  - THM-3807-trinomial-and-full-cubic-r-repairs-have-critical-points
  - THM-3800-sharp-torus-escaping-nodal-carrier-has-fourteen-critical-points
script: 04-computation/jc2_cubic_pseudoplane_quintic_r_repair_thm3817.py
output: 05-knowledge/results/jc2_cubic_pseudoplane_quintic_r_repair_thm3817.out
script_sha256: 9df9627bfd2ae8be0ffdad2d8ff8f917fb2f35fda95fa4dbbfe5fb313f647110
output_sha256: fd988b916a79e682a3c31d2c007d5821bb742fc30320019a1b14382504300c80
semantic_sha256: 2fa1efbadbfeff8ff24257948412e1db1f461851c15434949cae63fbd5c76f69
hash_basis: raw LF bytes
---

# THM-3817 -- every quintic r-repair remains critical

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Work over an
algebraically closed field `k` of characteristic zero.  On the `c=1` cubic
pseudo-plane put

```text
Y=Spec k[r,z,e]/(r^2e-z^3+r),
{r,z}=3r^2,       {r,e}=9z^2,       {z,e}=3+6re.       (1)
```

For arbitrary `b_0,...,b_5 in k`, let

```text
g(e)=sum_(i=0)^5 b_i e^i,
A=e^2-z/3+r g(e).                                      (2)
```

Then `A` has a critical point on `Y`, and therefore no regular Darboux mate.
The case `b_5=0` is THM-3813.  It remains to prove

```text
b_5 != 0.                                               (3)
```

## 1. Universal residual and multiplicity-safe boundary test

Put

```text
u=re,                         K=1+2u,
P=g u^2-K(2e^3+u e g'),
Q=e^2K^3-729g^3u^2(1+u)^2.                            (4)
```

As in THM-3813, the Hamiltonian components are

```text
{A,r}=r^2-9z^2(2e+rg'),
{A,z}=3gr^2-3(1+2re)(2e+rg'),
{A,e}=9gz^2-(1+2re),                                  (5)
```

and exact elimination, with independent symbols `G,D` for `g,g'`, gives

```text
Res_u(P,Q)=G^3e^4 H_univ(e,G,D).                       (6)
```

After the quintic substitution, write the residual as `H(e)`.  Then

```text
deg H=27,                       LC(H)=34012224 b_5^5.  (7)
```

Suppose every root of `H` lay on `V(e g)`.  Factoring with arbitrary
multiplicities over `k` gives exactly as in THM-3813

```text
H divides e g H'.                                     (8)
```

This uses no squarefreeness of `g`.  Divide `e g H'` by `H` over the
parameter fraction field.  Its quotient has degree five, leading coefficient
`27b_5`, and exact global denominator

```text
512b_5^4.                                              (9)
```

Let

```text
R=512b_5^4(e g H'-qH)=sum_(j=0)^26 r_j e^j,            (10)
```

where `q` is that quotient.  Under `(3)`, divisibility `(8)` is equivalent to
the vanishing of every coefficient of `R`.

## 2. The finite weight cover and top ideal

The quintic leading coefficient has mod-seven grading weight `-2`.  Since
`k` is algebraically closed, choose `s!=0` with

```text
b_5=s^2.                                               (11)
```

Set

```text
C_i=b_i s^(3-i)  (0<=i<=4),              T=s^7.       (12)
```

Equivalently,

```text
b_0=C_0/s^3, b_1=C_1/s^2, b_2=C_2/s,
b_3=C_3, b_4=C_4s, b_5=s^2.                           (13)
```

This is a surjective finite cover of the entire locus `b_5!=0`, not a generic
normalization; in particular `T!=0`.  For `26>=j>=15`, the pullback of `r_j`
is a nonzero rational scalar times

```text
s^(j-12) F_j(C_0,C_1,C_2,C_3,C_4,T),                  (14)
```

where `F_j` is the signed integer primitive part.  Thus the powers are
`14,13,...,3`, and no coefficient equation is lost.  If `(8)` held, then

```text
F_26=F_25=...=F_15=0.                                 (15)
```

Work in

```text
Q[C_0,C_1,C_2,C_3,C_4,T]                              (16)
```

with graded reverse lexicographic order on the displayed variables, and put

```text
I=(F_26,F_25,F_24,F_23,F_22).                         (17)
```

Exact Buchberger reduction gives a 149-element Groebner basis for `I`.  The
companion reconstructs every generator from `(6),(9)--(14)`, checks primitive
back-substitution coefficient by coefficient, and freezes the full invariant
packet and basis by hashes.

## 3. Four rows force a one-parameter resonance

Modulo `I`, the next four rows are

```text
NF(F_21)=-512T L_21,
L_21=210C_3-83C_4^2+678,

NF(F_20)=-128T L_20,
L_20=6060C_2+5322C_3C_4-2567C_4^3+25238C_4,

NF(F_19)=-128T L_19,
L_19=9920C_1+16317C_2C_4+4548C_3^2+436C_3C_4^2
     +16256C_3-2201C_4^4+30018C_4^2-9680,

NF(F_18)=-(256/153)T L_18,
L_18=4088250C_0+1340840C_1C_4+693900C_2C_3
     +1400082C_2C_4^2+666198C_2+732067C_3^2C_4
     -642798C_3C_4^3+2909169C_3C_4
     +1283409C_4^3-1975219C_4.                       (18)
```

Because `T!=0`, these are triangular with nonzero constant pivots.  They
force uniquely

```text
C_3=(83C_4^2-678)/210,

C_2=26C_4(156C_4^2-2711)/53025,

C_1=(13008795C_4^4-322587262C_4^2+2738675802)
    /1841028000,

C_0=C_4(148537264577C_4^4-4086004213596C_4^2
        -24178143226221)/282246852037500.             (19)
```

No profile coefficient was divided out in `(18)--(19)`.

## 4. The last two branches are impossible

Substitute `(19)` into the normal forms of the remaining rows and take their
integer primitive numerators in `Q[C_4,T]`; call them `P_j` for row `j`.
Put `x=C_4^2`.  The row-sixteen numerator factors exactly as

```text
P_16=C_4 T S(x),                                      (20)

S(x)=133257064516197275229295x^3
    -3208348664856734112057306x^2
    +349275415597792003907999514x
    -100800589473950881879245888.                     (21)
```

There are two branches.

* If `C_4=0`, row seventeen specializes to

  ```text
  P_17=49208702469888916532372952734335038991872656250000 T,
                                                               (22)
  ```

  which is nonzero because `T!=0`.

* Suppose `C_4!=0`.  Then `(20)` forces `S(x)=0`.  Both `P_17` and `P_15`
  are affine-linear in `T`.  Their exact resultant, eliminating **only** `T`,
  is

  ```text
  Res_T(P_17,P_15)
   =-3363550544979296875 C_4 Q_4(x)Q_8(x)Q_10(x),      (23)
  ```

  where

  ```text
  Q_4 =148537264577x^2-4086004213596x-24178143226221,

  Q_8 =22172090448333650297319377595200x^4
       -782476464930749584609729836216645x^3
       +43455557119423145069424296532772242x^2
       -196699263766307709411772113438437964x
       +637608548591074605888432844835874192,

  Q_10=709557651508549958433446867x^5
       -38269012170154817355067751208x^4
       +413899168413250905311964034530x^3
       +3370561945779184898686058483676x^2
       -51355055679882251411436522510600x
       +80149644817765957262829141790560.              (24)
  ```

  Exact Euclidean gcds in `Q[x]` give

  ```text
  gcd(S,Q_4)=gcd(S,Q_8)=gcd(S,Q_10)=1.                 (25)
  ```

  As independently cheap nonvanishing certificates, the corresponding
  integer resultants reduce modulo `103` to

  ```text
  Res_x(S,Q_4)=20,  Res_x(S,Q_8)=37,
  Res_x(S,Q_10)=37                 (mod 103).           (26)
  ```

  Thus none of `Q_4,Q_8,Q_10` can vanish at a root of `S` in characteristic
  zero.  Since `C_4!=0`, `(23)` is nonzero, contradicting the common equations
  `P_17=P_15=0`.

Both branches contradict `(15)`.  Hence `(8)` is impossible, and `H` has a
root `eta` with

```text
eta g(eta) != 0.                                      (27)
```

The only saturation in the argument is `T=s^7!=0`, inherited from `b_5!=0`.
The `C_4=0` seam is retained and killed explicitly.  Repeated roots of `g`,
sparse quintics, and both choices of `s` are all included.

## 5. The residual root is an actual critical point

At `e=eta`, `(6)--(7)` make the `u`-resultant vanish.  The leading coefficient
of the quartic `Q` is

```text
LC_u(Q)=-729g(eta)^3 !=0,                             (28)
```

so the common projective root is finite even if `P` drops degree.  Moreover,

```text
Q(eta,0)=eta^2,
Q(eta,-1)=-eta^2,
Q(eta,-1/2)=-729g(eta)^3/16,                          (29)
```

excluding `u=0,-1,-1/2`.  For a common root `u_0`, put

```text
r_0=u_0/eta,
z_0=9g(eta)u_0(1+u_0)/(eta(1+2u_0)).                 (30)
```

The exact identities

```text
z_0^2-K_0/(9g(eta)) = -Q/(9g(eta)eta^2K_0^2),
r_0^2eta-z_0^3+r_0 = u_0(1+u_0)Q/(eta^3K_0^3),
{A,e}|_0 = -Q/(eta^2K_0^2),
{A,z}|_0 = 3P/eta^2                                  (31)
```

put the point on `Y` and kill the last two Hamiltonian components.  The
surface Casimir identity

```text
(1+2re){A,r}-3z^2{A,z}+r^2{A,e}=0                    (32)
```

then kills `{A,r}` because `1+2u_0!=0`.  This is a genuine critical point.

Together with THM-3813, this theorem closes every pure `r g(e)` carrier
with `deg g<=5`.  It does **not** address degree at least six, mixed
`z^2h(e)+r g(e)` corrections, another pseudo-plane arm profile, or arbitrary
planar Keller maps.  The exact companion has 1448 active gates; normal and
optimized runs match the frozen transcript exactly after LF normalization.
**QED.**
