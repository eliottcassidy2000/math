---
id: THM-3861
title: "Cubic normal-strip Keller pairs are automorphisms"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over every
  characteristic-zero field, a polynomial Keller pair in k[s,z] whose two
  transverse z-degrees are at most three is a polynomial automorphism.  After
  a constant target change, every genuinely cubic survivor has one z-linear
  coordinate C=b+beta*z and the other is
  rho*C^3+d*C^2+e*C-(lambda/beta)*s+a0.  The apparent (3,2) Kummer branch
  factors through a rational profile X=beta/h; the constant Jacobian forces X
  to have a finite pole when h is nonconstant, and its cubic principal part
  then makes the arm coefficient nonpolynomial.  The h-constant edge is also
  contradictory.  Cubic polynomial strip depth therefore cannot realize the
  self-identifying Russell arm.  Rational and infinite formal strip expansions
  and arbitrary global Darboux pairs remain open.
source: root / higher_normal_lift / first live cubic normal-depth lane, 2026-08-23
audit: >
  INDEPENDENT HOSTILE PROOF/REPLAY AUDIT PASSED.  The independent companion
  reconstructs the six buckets by coefficient convolution, both SL2 target
  charts and the quadratic handoff, the beta=0 contradiction, all (3,1)
  quotient integrations and the two-sided inverse, and the (3,2) transformed
  equations and factored E0.  Its 622 optimization-safe gates include 361
  valuation profiles, all seven Kummer lattice points, a nonclosed-Q
  scalar-root-free control, local DVR and constant-h edges, degenerate (3,1)
  controls, and the sharp rational bracket -3/16 hostile.  Normal and
  optimized executions of both companions byte-match their frozen LF outputs.
  A second independent 17,524-gate audit rebuilds the proof from the six
  buckets without using the first audit.  Its determinant-one target shear
  gives a different local-DVR exclusion of every nonconstant `(3,2)` factor;
  the constant edge reduces to `J=(KB+d)B'`.  It checks 14,064 valuation
  profiles, both target charts, every zero edge, the explicit inverse and the
  sharp rational hostile.  Normal and optimized replay again byte-match.
depends_on:
  - THM-3856-quadratic-normal-strip-keller-pairs-are-automorphisms
related:
  - THM-3785-linear-higher-pole-russell-pseudoplane-maximal-observable
  - THM-3843-russell-arm-birational-immersion-and-forced-self-identification
  - THM-3846-formal-arm-darboux-lift-and-algebraization-gate
  - THM-3860-russell-higher-normal-rational-lifts-and-vertical-pole-barrier
script: 04-computation/jc2_cubic_normal_strip_keller_thm3861.py
output: 05-knowledge/results/jc2_cubic_normal_strip_keller_thm3861.out
script_sha256: 056f177d0586362fb0a1fa3daa3e57e77be0af6ac105a58fe856cb1a3ba4182e
output_sha256: 51b5c1e05e67c6043228ae3656486afb6e5a7be417ab78ec4a4f0b9269fb16a0
semantic_sha256: 372e956d295f82c3e21540e368d5b9ad932fdb3e954a36a97f0c86e346de4b7f
independent_script: 04-computation/jc2_cubic_normal_strip_keller_independent_audit_thm3861.py
independent_output: 05-knowledge/results/jc2_cubic_normal_strip_keller_independent_audit_thm3861.out
independent_script_sha256: 3c658af035d67b6fa7d230a3985147e5ee29e39476c3cb94b845c191fe1cca44
independent_output_sha256: 858e6f8aca4a38a24724a29aa0e9df30161109337a11c6edf2792757027443d9
independent_semantic_sha256: 7862efeee48d97b47b7c6c55fab9a3ad83d2d137815b8b71a3f2eab10e5724b5
second_independent_script: 04-computation/jc2_cubic_normal_strip_keller_second_hostile_audit_thm3861.py
second_independent_output: 05-knowledge/results/jc2_cubic_normal_strip_keller_second_hostile_audit_thm3861.out
second_independent_script_sha256: d5adf63924d63f326c2e25767e01f70bb732b977783375f391d41812fd4e75d6
second_independent_output_sha256: e9f4acf4f82118729ba2ebdc4067cbc199d732b7d8c38153706a08cbb56f1e82
second_independent_semantic_sha256: e286e5276ec64f1e06a4f8aa187edd2cbef8fdfefb14acf47da602d95d0eb44c
hash_basis: raw LF bytes
---

# THM-3861 -- cubic transverse depth is still triangular

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  Let `k` be a
field of characteristic zero.  Suppose

```text
A,C in k[s,z],                  deg_z A,deg_z C <=3,             (1)
J_(z,s)(A,C)=A_zC_s-A_sC_z=lambda in k*.                        (2)
```

Then `(A,C)` is a polynomial automorphism.  More precisely, after a
constant target `SL_2(k)` change, either the pair has transverse degree at
most two and THM-3856 applies, or there are

```text
rho,beta in k*,                 d,e,a_0 in k,       b in k[s]  (3)
```

such that

```text
C=b(s)+beta z,
A=rho C^3+d C^2+e C-(lambda/beta)s+a_0.                       (4)
```

The inverse in the genuinely cubic case is explicit:

```text
s=(beta/lambda)[rho C^3+d C^2+e C+a_0-A],
z=[C-b(s)]/beta.                                               (5)
```

The only new top shape not visibly covered by `(4)` has normalized
transverse degrees `(3,2)`.  The proof below shows that this shape is empty
in `k[s,z]`.  It also records a rational-coefficient pair in `k(s)[z]`
that satisfies every bucket and fails at exactly one cubic principal part;
thus polynomiality in the arm coordinate is load-bearing.

## 1. The six exact coefficient buckets

Write

```text
A=a+alpha z+u z^2+p z^3,
C=b+beta z+v z^2+q z^3,                                      (6)
```

where all eight coefficients lie in `k[s]`.  Expanding `(2)` gives exactly

```text
z^0: alpha b'-a'beta=lambda,                                  (7)
z^1: alpha beta'-alpha'beta+2(ub'-a'v)=0,                     (8)
z^2: alpha v'+2u beta'+3pb'-2alpha'v-u'beta-3a'q=0,           (9)
z^3: alpha q'+2uv'+3p beta'-3alpha'q-2u'v-p'beta=0,          (10)
z^4: 2u q'+3p v'-3u'q-2p'v=0,                               (11)
z^5: 3(pq'-p'q)=0.                                           (12)
```

These are identities in the ordinary polynomial ring; there is no
completion, reduction, or saturation.

If `p=q=0`, THM-3856 proves the result, including every lower-degree edge.
Suppose `(p,q)!=(0,0)`.  Equation `(12)` says that the nonzero rational
ratio of `p,q` has derivative zero.  The constant field of `k(s)` is `k`,
so

```text
(p,q)=h(s)(P,Q),                    (P,Q) in k^2\{0}.          (13)
```

A constant target `SL_2(k)` matrix sends `(P,Q)` to a nonzero multiple of
`(1,0)` and preserves `(2)`.  Absorb that multiple into `p`.  Hence every
genuinely cubic branch may be normalized as

```text
q=0,                              p!=0.                         (14)
```

The remaining top equation is

```text
3pv'-2p'v=0.                                               (15)
```

It splits the proof into the `(3,1)` branch `v=0` and the apparent `(3,2)`
branch `v!=0`.

## 2. Complete `(3,1)` classification and inverse

Let `v=0`.  If `beta=0`, equation `(9)` becomes `3pb'=0`, so `b'=0`;
then `(7)` contradicts `lambda!=0`.  Thus `beta!=0`.  Equation `(10)` is

```text
3p beta'-p'beta=0,                                           (16)
```

or `(p/beta^3)'=0` in `k(s)`.  Hence

```text
p=rho beta^3,                         rho in k*.                (17)
```

Successively, equations `(9)` and `(8)` become

```text
(u/beta^2)'=3rho b',
(alpha/beta)'=2(3rho b+d)b',                                (18)
```

and therefore

```text
u=beta^2(3rho b+d),
alpha=beta(3rho b^2+2db+e),                 d,e in k.          (19)
```

Finally `(7)` factors as

```text
beta[rho b^3+d b^2+e b-a]'=lambda.                          (20)
```

A product of polynomials is a nonzero scalar only if both factors are
units.  Thus `beta in k*`, and integration gives

```text
a=rho b^3+d b^2+e b-(lambda/beta)s+a_0.                     (21)
```

Expanding `(4)` in powers of `z` reproduces `(17),(19),(21)`, while `(5)`
recovers first `s` and then `z` polynomially.  This proves both the
classification and automorphism claim in the entire `(3,1)` branch,
including the cases in which some of `u,alpha,b'` vanish.

## 3. Kummer normalization of the apparent `(3,2)` branch

Now suppose `v!=0`.  At an irreducible prime `pi`, write the local orders of
`p,v` as `A_pi,B_pi`.  Characteristic zero makes `pi'` a unit in the DVR
`k[s]_(pi)`, so the lowest-order coefficient of `(15)` is proportional to
`3B_pi-2A_pi`.  Hence

```text
2 ord_pi(p)=3 ord_pi(v)                                     (22)
```

for every irreducible `pi in k[s]`.  Since the only units of `k[s]` are
`k*`, unique factorization gives

```text
p=P h^3,                         v=V h^2,
P,V in k*,                       0!=h in k[s].                (23)
```

This includes constant `h`; no algebraic closure or extraction of scalar
roots is used.  In `k(s)` put

```text
X=beta/h,              K=3P/(2V^2),             M=3P/(2V).   (24)
```

Equation `(10)` can be written

```text
-2v^2(u/v)'+3P h^4 X'=0,                                   (25)
```

and hence

```text
u/v=KX+d,                              d in k.                 (26)
```

Write `alpha=hY`.  Substitution into `(9)` cancels every `h'` term and
gives

```text
2Y'=(KX+2d)X'+(3P/V)b'.                                    (27)
```

Therefore

```text
alpha/h=(K/4)X^2+dX+Mb+e,                 e in k.             (28)
```

Substituting `(26),(28)` into `(8)` and integrating once gives the exact
arm coefficient

```text
a= -P X^3/(16V^3)+3P bX/(4V^2)+eX/(2V)+db+a_0,
a_0 in k.                                                     (29)
```

No polynomiality has yet been used in `(24)--(29)`; these identities hold
in the rational function field and are exhaustive because each integration
uses the characteristic-zero constant field `k`.

Define

```text
T=b-X^2/(4V).                                                (30)
```

A final substitution into the constant bucket `(7)` factors it completely:

```text
h(MT+e)T'=lambda.                                           (31)
```

Equations `(23),(26),(28)--(31)` are the complete normalized `(3,2)`
packet.  The obstruction is now local and all-degree.

## 4. Polynomiality makes the `(3,2)` packet empty

First suppose `h` is nonconstant, and choose any irreducible factor `pi` of
`h`.  Work in the DVR

```text
R=k[s]_(pi).                                                  (32)
```

The ordinary derivative preserves `R`: if `f/g in R` with `g` a `pi`-unit,
then

```text
(f/g)'=(f'g-fg')/g^2 in R.                                  (33)
```

If `X` belonged to `R`, then `b,T,MT+e,T'` would all belong to `R`.
Because `h` is in the maximal ideal, the left side of `(31)` would also be
in the maximal ideal (or would vanish), contradicting the unit `lambda`.
Consequently

```text
ord_pi(X)=-m<0.                                              (34)
```

But `(29)` then has a unique most singular summand:

```text
ord_pi[-P X^3/(16V^3)]=-3m,
ord_pi(bX),ord_pi(eX)>=-m,
ord_pi(db),ord_pi(a_0)>=0.                                  (35)
```

The coefficient of `X^3` is nonzero, so cancellation at order `-3m` is
impossible.  Thus `a` has a pole at `pi`, contrary to `a in k[s]`.

It remains to treat the exceptional case `h in k*`.  Then `X,T,MT+e` are
polynomials, and `(31)` is a product of polynomials equal to a nonzero
scalar.  Both `MT+e` and `T'` must therefore be nonzero constants.  Since
`M!=0`, the first condition makes `T` constant, while the second says it is
nonconstant.  This contradiction closes the constant-`h` edge.

Hence the normalized `(3,2)` branch is empty.  Together with Section 2 and
THM-3856 for `(p,q)=(0,0)`, this exhausts every degree and zero-component
case in `(1)`.

## 5. A sharp rational hostile

The polynomiality step cannot be omitted.  In `k(s)[z]`, set

```text
C=s^4z+s^10z^2,
A=-1/(16s^3)+(3/8)s^3z+(3/2)s^9z^2+s^15z^3.                 (36)
```

Direct differentiation gives

```text
J_(z,s)(A,C)=-3/16.                                         (37)
```

This is `(23)--(31)` with

```text
h=s^5,        P=V=1,        X=1/s,        b=d=e=0,
T=-1/(4s^2).                                                   (38)
```

It satisfies all six buckets, including

```text
h(MT)T'=-3/16.                                               (39)
```

The only failed requirement is exactly the pole
`-X^3/16=-1/(16s^3)` in `(29)`.  Thus `(36)` is a positive control for the
Kummer integrations and a hostile against treating rational arm
coefficients as polynomial ones.  It is not a polynomial Keller map and
has no Jacobian-conjecture consequence.

## 6. Russell boundary

THM-3846 identifies the completed Russell arm as

```text
Bhat=k[s][[z]],                         {z,s}=1.                (40)
```

The restriction to `z=0` of a polynomial automorphism of `A2_(s,z)` is
injective.  THM-3843 forces the arm restriction of every hypothetical
global Russell Darboux pair to be noninjective.  Therefore its completed
coordinates cannot both lie in `k[s,z]` with transverse degree at most
three.

This conclusion is one depth stronger than THM-3856, but its boundary is
the same: a global element of `B` can have an infinite canonical `z`-series,
and THM-3860 exhibits rational formal flexibility outside the polynomial
strip.  No arbitrary higher-normal obstruction, global Darboux pair, or
planar Jacobian-conjecture conclusion follows.  **QED.**

## 7. Reproduction

```bash
python3 -B 04-computation/jc2_cubic_normal_strip_keller_thm3861.py
python3 -B -O 04-computation/jc2_cubic_normal_strip_keller_thm3861.py
python3 -B 04-computation/jc2_cubic_normal_strip_keller_independent_audit_thm3861.py
python3 -B -O 04-computation/jc2_cubic_normal_strip_keller_independent_audit_thm3861.py
python3 -B 04-computation/jc2_cubic_normal_strip_keller_second_hostile_audit_thm3861.py
python3 -B -O 04-computation/jc2_cubic_normal_strip_keller_second_hostile_audit_thm3861.py
```

Each primary execution must byte-match
`05-knowledge/results/jc2_cubic_normal_strip_keller_thm3861.out`; the two
independent pairs must byte-match their corresponding independent-audit
outputs.
