---
id: THM-3215
title: "Arbitrary-degree root-jet Hamiltonian, affine-dihedral holonomy, and p-fold carry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Every polynomial has a universal trace-zero two-jet Hamiltonian whose square
  is scalar.  At a selected common simple root these Hamiltonians normalize to
  affine reflections, arbitrary words have one explicit alternating shear,
  and genuine pair transitions form a flat dlog cocycle.  The pair shear is
  the next Pluecker jet; in residue characteristic p a primitive shear has
  exact projective order p and returns as a nonzero first p-adic carry.
source: root/multiscale-newton-flag/low-child-flag-extension/2026-08-02
audit: >
  The assertion-independent exact companion pins THM-3206; proves the
  universal jet square and factorial conjugacy symbolically; checks 48
  arbitrary-degree simple-root polynomials, reflection words through length
  eight, the strict transition cocycle, 2,612 coordinate changes, 2,450
  factorial common-root pairs, 2,184 nontrivial exact-order-p pairs, and
  24,660 longer words.  It verifies the heterogeneous scalar four-word,
  five exact p-adic carries, and sharp higher-jet, multiple-root, and infinity
  boundaries.  An independent immutable audit rederived the conjugacy,
  word order/sign, coordinate law, flat cocycle, Pluecker orientation,
  transition-germ and original-frame signs, exact projective order and carry,
  and every stated scope boundary.  Fresh normal and `-O` replay byte-match
  the stored output and both LF-normalized hashes.
depends_on:
  - THM-3206-heterogeneous-factorial-exterior-reflection-groupoid-and-fixed-plane-holonomy
related:
  - THM-3175-first-witt-hensel-wedge-and-infinitesimal-pluecker-covariance
  - THM-3178-squarefree-resultant-tangent-cone-and-first-witt-norm
  - THM-3192-reciprocal-coefficient-jet-transfer-and-z-adic-pluecker-return
  - THM-3204-parabolic-continuant-single-gate-and-jacobi-smith-obstruction
  - THM-3214-two-jet-pseudo-division-locality-and-catalan-sharpness
script: 04-computation/arbitrary_degree_root_jet_affine_dihedral_carry_thm3215.py
output: 05-knowledge/results/arbitrary_degree_root_jet_affine_dihedral_carry_thm3215.out
script_sha256: 1799ad3b3c296bdb2fc41b98abb54eed4ae40937a02ca7519c5b2c3eee2fb4ff
output_sha256: 6c430edb048e5bd859c712068db07d9cf18867244e165795cac57ff3968c35e7
hash_basis: LF-normalized bytes
---

# THM-3215 -- arbitrary-degree root-jet Hamiltonian, affine-dihedral holonomy, and p-fold carry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-3206 identifies every quadratic factorial exterior block with a
projective reflection, but its common-root wall was described only through a
pair resultant.  The local mechanism is neither quadratic nor factorial.  It
is the two-jet representation of an arbitrary function at a selected root.

This theorem separates three levels which must not be blended:

```text
arbitrary reflection word:       one affine-dihedral alternating shear;
composable root-chart transport: a flat exact dlog cocycle;
global polynomial problem:       still needs a root branch and higher jets. (1)
```

## 1. The universal two-jet Hamiltonian

Let `A` be a commutative ring, let `f` be a polynomial or formal function in
one variable, and fix a point `a`.  Derivatives below are evaluated at `a`.
Define

```text
J_a(f)=[ f'    f'' ]
       [-2f   -f' ].                                      (2)
```

Direct multiplication gives the universal identity

```text
J_a(f)^2=(f'^2-2ff'')I_2.                                (3)
```

No degree bound and no division is used.  If `f(a)=0` and `f'(a)` is a unit,
put

```text
lambda_f=f'(a),              kappa_f=f''(a)/f'(a),
R(kappa)=[1  kappa]
         [0   -1  ].                                      (4)
```

Then

```text
J_a(f)=lambda_f R(kappa_f),       R(kappa_f)^2=I_2.       (5)
```

Thus every arbitrary-degree polynomial with a selected simple root supplies
a projective reflection.  The degree-independent statement is local: only
the value, first derivative, and second derivative at that root survive.

## 2. Exact identification with the factorial reflection

For THM-3206's reduced factorial quadratic

```text
q(t)=vt^2-t+s
```

and its internal reflection

```text
F(s,v)=[-(1+2v)              2v]
       [-2(s+v+1)          2v+1],                         (6)
```

put

```text
P_x=[1  0]
    [1+x 1].                                               (7)
```

For every `x`, whether or not it is a root,

```text
P_x^(-1)F(s,v)P_x
 =[ q'(x)     q''(x)]
  [-2q(x)    -q'(x)]
 =J_x(q).                                                  (8)
```

Moreover

```text
q'(x)^2-2q(x)q''(x)=1-4sv=Delta.                          (9)
```

Equation `(8)` is the missing literal map between THM-3206's exterior
reflection and a coefficient jet.  At a root `x`, `Delta=q'(x)^2`; hence the
off-discriminant hypothesis is exactly simplicity of the selected root.
The split maps `iota,pi` of THM-3206 then lift every identity below to its
fixed exterior carrier without changing the internal `2x2` calculation.

## 3. Complete common-root word law

Let `K` be a field of characteristic different from two.  Suppose
`f_1,...,f_N` have the same selected simple root `a`.  Write

```text
lambda_i=f_i'(a),             kappa_i=f_i''(a)/f_i'(a),
Lambda_N=product_(i=1)^N lambda_i,
Theta_N=sum_(i=1)^N (-1)^(i+1) kappa_i.                  (10)
```

Then the ordered word is

```text
J_a(f_N)...J_a(f_1)
 =Lambda_N [1  Theta_N]
           [0  (-1)^N].                                  (11)
```

Consequently:

* for odd `N`, the normalized word is again the reflection `R(Theta_N)` and
  is never scalar;
* for even `N`, the normalized word is the translation

  ```text
  T(Theta_N)=[1 Theta_N;0 1],                             (12)
  ```

  and it is scalar exactly when `Theta_N=0`;
* translations add and reflections reverse them:

  ```text
  T(u)T(v)=T(u+v),       R(kappa)T(u)R(kappa)=T(-u).       (13)
  ```

Thus the common-simple-root stratum of the reflection word lies in the
affine-dihedral group `G_a semidirect C_2`.

### Proof

Equation `(5)` reduces the assertion to the matrices `R(kappa_i)`.  The two
identities

```text
R(v)R(u)=T(u-v),            R(w)T(u)=R(u+w)               (14)
```

give `(11)` by induction.  Equations `(12),(13)` and the scalar
classification follow immediately.  Characteristic two is excluded because
the diagonal characters `1,-1` then collide. QED.

## 4. The actual root-chart cocycle is flat

An arbitrary reflection word and a composable chart transition are not the
same object.  Define the directed transition from the `i` chart to the `j`
chart by

```text
T_(j<-i)=R(kappa_j)R(kappa_i)=T(kappa_i-kappa_j).          (15)
```

Then every composable triangle telescopes strictly:

```text
T_(k<-j)T_(j<-i)=T_(k<-i).                                (16)
```

There is a coordinate-free reason.  In a local coordinate `x`, put

```text
theta_ij=(kappa_i-kappa_j)dx
        =d log(f_i'/f_j').                                (17)
```

If `y=phi(x)`, with `c=phi'(a)` a unit and `e=phi''(a)`, the chain rule gives

```text
kappa_i^(y)=kappa_i^(x)/c-e/c^2.                          (18)
```

The inhomogeneous term is common to every `i`, so

```text
(kappa_i^(y)-kappa_j^(y))dy
 =(kappa_i^(x)-kappa_j^(x))dx.                            (19)
```

Therefore `theta_ij` is an exact additive one-cocycle on a fixed common-root
branch, and its local triangle holonomy is zero.  The same difference is
unchanged if every `f_i` is multiplied by one common local unit: each
`kappa_i` acquires the same `2h'/h` connection term.

This is the holotopy distinction in `(1)`.  A nonzero even word in `(11)` is
a legitimate ordered block transport, as in THM-3206, but it is not by itself
a nontrivial Cech loop of the pair transitions `(15)`.  Genuine global
root-chart holonomy must enter through root permutation, branch monodromy,
a multiple-root seam, an infinity chart, or another label which prevents the
telescoping in `(16)`.

## 5. The shear is the next Pluecker jet

For two functions `f,g` with a common simple root, define

```text
Omega_2(f,g)=f''g'-g''f'.                                 (20)
```

In their common root frame,

```text
J_a(g)J_a(f)=mu I_2+Omega_2(f,g)E_12,
mu=f'g',                    E_12^2=0.                     (21)
```

Hence

```text
kappa_f-kappa_g=Omega_2(f,g)/(f'g').                      (22)
```

This is literally the next oriented Pluecker coordinate after the common
value row has vanished:

```text
Omega_2(f,g)=det [ f''  g'' ].                            (20a)
                     [ f'   g'  ]
```

The displayed orientation and the use of raw derivatives are intentional;
they avoid a sign or factorial ambiguity with coefficient-normalized
`P_1` conventions.  This coordinate is complementary to THM-3175/3178's
first-Witt wedge.  Their wedge is normal to the simple resultant divisor and
detects whether a residual common root lifts; `Omega_2` is tangent to the
exact common-root divisor and detects whether the two simple-root germs have
different second jets.

There is also an intrinsic transition-germ interpretation.  Since `f'(a)` is
a unit, let

```text
h=g composed with f^(-1)
```

be the formal coordinate change near zero.  Then

```text
h'(0)=g'/f',
h''(0)=(g''f'-g'f'')/f'^3,
kappa_f-kappa_g=-f' h''(0)/h'(0).                         (23)
```

Thus `(21)` is precisely the quadratic coefficient of the formal transition
germ.  Higher coefficients of that germ are the missing sidecar when
`Omega_2=0`.

### Quadratic specialization

Suppose `q_i(t)=v_i t^2-t+s_i` share the affine root `x`, so

```text
s_i=x-v_i x^2,                lambda_i=2v_i x-1.          (24)
```

Then

```text
Omega_2(q_1,q_2)=2(v_2-v_1).                              (25)
```

Returning from the root frame `(7)` to THM-3206's `F` frame gives

```text
F_2F_1=lambda_1lambda_2 I_2+N,

N=2(v_1-v_2)[ x+1       -1   ]
              [(x+1)^2  -(x+1)],          N^2=0.          (26)
```

The matrix in brackets has rank one over every field.  Hence distinct points
of the normalized quadratic common-root pencil give a nontrivial parabolic
pair, while equal points give a scalar pair, recovering and refining
THM-3206's resultant-zero classification.

## 6. Exact p-fold return and first carry

Let `O` be a DVR of mixed characteristic `(0,p)`, with `p` odd.  Suppose
`f,g in O[t]` have an exact common root `a in O`, both derivatives are units,
and `Omega_2(f,g)` is a unit.  Put `Q=J_a(g)J_a(f)` and `mu=f'g'`.  From
`(21)` and `E_12^2=0`, for every `r>=1`,

```text
Q^r=mu^r I_2+r mu^(r-1)Omega_2 E_12.                     (27)
```

Therefore the reduction of `[Q]` in `PGL_2(k)` has exact order `p`, and

```text
Q^p ==mu^p I_2                                      (mod p),

(Q^p-mu^p I_2)/p
 ==mu^(p-1)Omega_2 E_12                            (mod p). (28)
```

The second line is nonzero and rank one.  The mod-`p` scalar reset therefore
does not disappear: it returns one layer higher as a primitive first carry.
This is the pair-holonomy analogue of THM-3204's parabolic period-`p`
signature, without asserting that two individual factorial continuant steps
are simultaneously parabolic.

For the exact integer quadratics

```text
q_1=-t+1,                     q_2=2t^2-t-1,               (29)
```

the common root is `1`, `(mu,Omega_2)=(-3,4)`, and `(28)` is primitive for
every prime `p>=5`.  The companion verifies the full integer powers at
`p=5,7,11,13,17`, rather than merely reducing a formal binomial identity.

## 7. Sharp boundaries and cancellation hostiles

### 7.1 Longer heterogeneous words can return scalarly

THM-3206's two-reflection scalar criterion does not extend to a longer word.
Over `F_7`, take common root `x=1` and the four rows

```text
(s,v,lambda,kappa)
 =(1,0,6,0), (6,2,3,6), (5,3,5,4), (3,5,2,5).             (30)
```

All four parameter pairs are distinct and off the discriminant, but

```text
Theta_4=0-6+4-5=0 mod 7,
F_4F_3F_2F_1=5I_2.                                       (31)
```

Thus the alternating shear is the complete invariant on one common-root
word, and it can cancel globally.

### 7.2 Two jets do not identify an arbitrary polynomial

At `a=0`,

```text
f(t)=t,                         g(t)=t+t^3                (32)
```

are distinct and have the same value, first derivative, and second
derivative.  Hence `J_0(f)=J_0(g)` and their pair return is scalar.  This is
the sharp arbitrary-degree hostile: `(2)--(23)` remove the degree restriction
from the local root mechanism, but higher jets are genuinely necessary for
arbitrary coefficient recovery.

### 7.3 Multiple roots are not reflections

For `f(t)=t^2` at zero,

```text
J_0(f)=[0 2;0 0],               rank=1,       J_0(f)^2=0. (33)
```

The normalization `(4)` is undefined.  Multiple roots require a nilpotent
or higher-osculating chart; they are not a limiting simple-root reflection
inside the stated group law.

### 7.4 Infinity requires a homogeneous trivialization

The arbitrary-degree theorem is stated at a selected finite affine root.
At infinity, polynomials of different degrees are sections of different line
bundles, and `f''/f'` has no canonical meaning until a reciprocal coordinate
and a common homogeneous trivialization are chosen.  THM-3206's binary
quadratics do have such a fixed chart.  When `v_1=v_2=0`, their common root is
infinity and

```text
F(s_2,0)F(s_1,0)
 =I_2+[0 0;2(s_2-s_1) 0].                                (34)
```

Equation `(34)` is a valid separate boundary chart, not permission to apply
the affine arbitrary-degree formula at infinity without new data.

## 8. Consequence and scope

The proved connection contract is

```text
source:       arbitrary polynomial germs at one selected common simple root;
map:          f |-> J_a(f), then pairwise Pluecker shear Omega_2;
target:       flat affine root transitions and a primitive p-fold carry;
preserved:    simple-root branch, two jets, coordinate-invariant dlog class;
destroyed:    all third and higher jets, root choice, global branch labels;
needed next:  a canonical root selector plus the next formal-transition jet. (35)
```

This is the first mechanism in the factorial reflection lane whose local
statement is genuinely independent of polynomial degree.  It can therefore
be applied to arbitrary radial coefficient polynomials **after** a common
simple root has been selected and typed.  It does not prove that such a root
exists, survives the moment sum, or is compatible across faces.  No
arbitrary-support `NC(2)`, Gaussian Moment conjecture, Euclidean-depth
staircase, or `LRC(14)` conclusion follows.

## 9. Exact evidence

Run

```text
python 04-computation/arbitrary_degree_root_jet_affine_dihedral_carry_thm3215.py
python -O 04-computation/arbitrary_degree_root_jet_affine_dihedral_carry_thm3215.py
```

and compare LF-normalized bytes with the declared output.  The companion is
assertion-independent and uses exact symbolic or integer arithmetic only.
Besides the symbolic identities, it exhausts `2,450` admissible common-root
quadratic pairs and `24,660` longer words over the stated finite fields,
checks exact projective order on all `2,184` nontrivial pairs, and carries the
four sharp boundaries above.  There is no floating point or randomness.

**QED.**
