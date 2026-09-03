---
id: THM-4397
title: "Rank-two Poisson counterexample symplectic gauge equivalence"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Christopher D. Long's
  arXiv:2608.23777v1 rank-two Poisson map and THM-2044 are polynomially
  symplectically right-left equivalent over Q. The source equivalence is an
  explicit Hamiltonian momentum translation and the target equivalence is a
  linear symplectic quarter-turn. Both lift to exact A_2 Weyl automorphisms,
  so finite A_2 quantizability and polynomial termination are gauge-invariant.
  This gives no A_2 lift, DC(2), planar JC(2), or de-stabilization result.
source: root + independent referee / JC2 and August-preprint continuation session, 2026-09-03
depends_on:
  - THM-2044-explicit-rank-two-poisson-counterexample-by-symplectic-suspension
related:
  - THM-2045-the-smooth-factorized-R-family-has-no-planar-jacobian-mate
  - THM-2049-the-DC2-Ore-boundary-correction-complex-is-acyclic
primary_source: https://arxiv.org/abs/2608.23777v1
primary_script: 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397.py
primary_output: 05-knowledge/results/poisson_rank2_long_thm2044_gauge_equivalence_thm4397.out
primary_script_sha256: d72d527b24a2fb6b24408eadf20bb7e9e7e8810e12040f396f9300934a03e8ba
primary_output_sha256: ef2a2e8dfd1b790ee690037e21c678dc06a9176d87563efd360d2258c8033f9e
independent_audit_script: 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397_independent_audit.py
independent_audit_output: 05-knowledge/results/poisson_rank2_long_thm2044_gauge_equivalence_thm4397_independent_audit.out
independent_audit_script_sha256: 01d94738c60d6cdfa203cb9631da69b329a7929c5398eaaa3056e4c0b2b46979
independent_audit_output_sha256: 3d0b82065be6bb1b10346a4a5b4e5810033ae27629dabd98549ccf4635bf3a2a
hash_basis: raw LF bytes
audit: >
  PASS. The primary expands both four-variable maps and verifies their exact
  conjugacy, source and target symplectic identities, inverse, complexity
  profiles, Weyl lift, and three-point fibre transport. The independent audit
  imports no primary or THM-2044 code; it reconstructs both maps, uses the two
  adapted-coordinate increment identities, directly verifies all six paper
  brackets, recomputes the determinant-one core and reduced Groebner fibre,
  and checks the opposite source orientation. Normal and optimized runs
  byte-match both frozen LF outputs. No finite-field inference is used.
---

# THM-4397 -- rank-two Poisson counterexample symplectic gauge equivalence

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. THIS IDENTIFIES TWO
FOUR-VARIABLE POISSON COUNTEREXAMPLE PRESENTATIONS. IT DOES NOT PRODUCE AN
`A_2` ENDOMORPHISM, PROVE `DC(2)`, OR TOUCH PLANAR `JC(2)`.**

## 1. Statement

Use the convention

```text
{p,x}={z,q}=1
```

on `Q[x,q,p,z]`.  Let

```text
Phi_L=(R_L,T_L,D_L,S_L)
```

be the point map in Christopher D. Long,
[*An Explicit Counterexample to the Rank-Two Poisson Conjecture*](https://arxiv.org/abs/2608.23777v1),
and let

```text
Phi_R=(R_R,T_R,D_R,S_R)
```

be the map in THM-2044.  Put `s=xq` and define

```text
F=q^3(63s^2+318s+601)/840.                              (1)
```

The source map

```text
U_F:(x,q,p,z) |-> (x,q,p+F_x,z+F_q)                    (2)
```

and target map

```text
tau:(r,t,d,w) |-> (r,w,d,-t)                           (3)
```

are polynomial symplectomorphisms over `Q`, and

```text
boxed: Phi_R = tau o Phi_L o U_F.                       (4)
```

Thus the two maps lie in one polynomial symplectic right-left equivalence
class.  In the reverse orientation, put `A=U_F^(-1)=U_(-F)` and

```text
kappa:(r,t,d,w) |-> (r,-w,d,t).
```

Then equivalently

```text
Phi_L = kappa o Phi_R o A.                              (5)
```

## 2. The adapted-coordinate calculation

Both constructions begin with

```text
L =3x^2p+2(1-3xq)z,
D0=((1+3xq)/2)p-3q^2z,
R =x(2-3xq).                                            (6)
```

Long uses

```text
beta=L-9q^2,
y_L=q-x beta/3,                                         (7)
```

whereas THM-2044 uses

```text
G(s)=252s^3+1008s^2+1379s+659,
g=-q^2G(s)/140,
ell=L+g,
y_R=q-x ell/3.                                          (8)
```

Differentiating `(1)` gives

```text
F_x=q^4(21xq+53)/140,
F_q=q^2(105x^2q^2+424xq+601)/280.                       (9)
```

The first exact increment identity is

```text
3x^2 F_x+2(1-3xq)F_q=g+9q^2.                           (10)
```

Therefore

```text
beta o U_F=ell,        y_L o U_F=y_R.                  (11)
```

Regard the third adapted coordinate `e` as independent.  If `H_L(x,y,e)` is
Long's correction and `H_R(x,q,e)` is THM-2044's correction, direct expansion
gives the second identity

```text
H_R-H_L=-q^4(18x^2q^2+78xq+125)/20.                   (12)
```

The same polynomial is the `D0` increment under `(2)`:

```text
((1+3xq)/2)F_x-3q^2F_q
 =-q^4(18x^2q^2+78xq+125)/20.                          (13)
```

Equations `(11)--(13)` prove

```text
D_L o U_F=D_R.                                          (14)
```

The remaining three coordinates are the same determinant-one core in
`(x,y,e)`:

```text
Pi=(1+xy)^3e+y^2(1+xy)(4+3xy),
Sigma=y+3x(1+xy)^2e+3xy^2(4+3xy),
Rcore=2x-3x^2y-x^3e.                                   (15)
```

Long orders it as `(Rcore,-Pi/2,Sigma)` and THM-2044 as
`(Rcore,Sigma,Pi/2)`.  The quarter-turn `(3)` and `(11),(14)` now give every
coordinate of `(4)`.

The generating coefficient can also be recovered rather than guessed.  If
`F=q^3 k(s)`, equation `(10)` is

```text
s(2-3s)k'(s)+6(1-3s)k(s)
 =(601-1379s-1008s^2-252s^3)/140,                      (16)
```

whose polynomial solution is `k=(63s^2+318s+601)/840`.

## 3. Symplectic and Weyl gauges

Because `(2)` fixes `x,q` and translates the two momenta by a gradient,

```text
d x wedge d(p+F_x)+d q wedge d(z+F_q)
 =d x wedge dp+d q wedge dz.                           (17)
```

The mixed term vanishes by `F_xq=F_qx`; `U_(-F)` is the polynomial inverse.
The exact matrix identity `J_tau^T Omega J_tau=Omega` proves the same for
`(3)`.

The identical calculation works in the second Weyl algebra with
`[p,x]=[z,q]=1`.  Multiplication by `F_x,F_q` gives

```text
[p+F_x,z+F_q]=partial_x(F_q)-partial_q(F_x)=0,          (18)
```

and all other generator relations are immediate.  The linear target
quarter-turn is also a Weyl automorphism.  Hence any finite polynomial `A_2`
quantization of one classical map transports to one of the other, and any
proof of nontermination invariant under Weyl automorphisms transports back.
This equivalence does not assert that either quantization exists.

## 4. Complexity and fibre transport

Direct expansion gives

```text
                  nonzero terms          total degrees
Phi_L             (2,47,139,22)          (3,15,23,11)
Phi_R             (2,35,246,78)          (3,21,48,30).  (19)
```

Thus Long's gauge is substantially smaller in the high-degree coordinates,
even though it is not a new counterexample class.  It is the better concrete
representative for future coupled-`D` Weyl computations.

Long's three source points are

```text
(0,0,1/24,-1/8),
(1,2/3,247/96,-89/64),
(-1,-2/3,247/96,-89/64).                               (20)
```

Applying `U_F^(-1)` sends them exactly to THM-2044's points

```text
(0,0,1/24,-1/8),
(1,2/3,224839/90720,-173417/60480),
(-1,-2/3,224839/90720,-173417/60480).                  (21)
```

The target quarter-turn sends `(0,1/8,0,0)` to `(0,0,0,-1/8)`.  Independently,
the reduced lexicographic core fibre has basis

```text
e-(27x^2-1)/4,
y+3x/2,
x(x-1)(x+1),                                           (22)
```

so the fibre is reduced and contains exactly three points.

## 5. Consequences and stopping boundary

1. The new preprint and THM-2044 use the same THM-1300 three-dimensional core
   and the same noninvertibility mechanism.  Complexity is gauge-dependent.
2. Long's smaller representative can reduce symbolic support in the open
   coupled-`D` quantization computation.  It cannot change the existence or
   finite-termination answer because the gauge itself is an exact Weyl
   automorphism.
3. THM-2049's abstract `(Sigma,Pi/2)` Ore boundary pair becomes Long's
   `(S,-T)` after the target rotation, so its graded acyclicity and open
   termination gate are unchanged.
4. The first coordinate remains `R=x(2-3xq)`.  THM-2045 proves that it has no
   polynomial planar Jacobian mate.  Stabilization is therefore essential for
   this route, and planar `JC(2)` remains open.
5. Long's Appendix-B construction is in `A_4`, not `A_2`; it neither proves
   `DC(2)` nor supplies the missing rank-preserving quantization.

## 6. Exact audit and reproduction

The primary certificate expands both full four-variable maps and proves
`(4)` coefficientwise.  It checks the source and target symplectic identities,
the polynomial inverse, the Weyl mixed relation, profiles `(19)`, and every
point in `(20)--(21)`.

The independent audit imports neither the primary nor THM-2044 code.  It uses
the reverse orientation `(5)`, separately reconstructs `(10)--(14)`, directly
checks Long's six Poisson identities, recomputes the determinant of `(15)` and
the Groebner basis `(22)`, and transports source and target fibres.

Replay from the repository root:

```text
python3 -B 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397.py
python3 -B -O 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397.py
python3 -B 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397_independent_audit.py
python3 -B -O 04-computation/poisson_rank2_long_thm2044_gauge_equivalence_thm4397_independent_audit.py
```

Normal and optimized executions byte-match the frozen LF outputs.  No
finite-field specialization, numerical tolerance, or unproved inverse is
used.
