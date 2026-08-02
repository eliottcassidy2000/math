---
id: THM-3172
title: "Shear-invariant differential-owner filtration and transverse recurrence"
status: >
  PROVED + VERIFIED-EXACT.  For a marked finite separable inverse pair over
  k(u,v), the filtration obtained by adjoining all iterated target
  derivatives is unchanged by every polynomial target shear fixing u.  If
  the pair is owned by one polynomial A2 Keller map, every filtration stage
  is exactly k[x,y].  Consequently a shear representative v-H(u)=Phi(z)
  with deg(Phi)>1 forces a nonconstant unit at the first derivative stage and
  cannot have such an owner.  This rejects the punctured cubic component of
  the THM-3167 hostile at one derivative, while rejecting its diagonal and
  nondiagonal members equally and therefore not deciding the inverse-
  different line gate.  Separately, an explicit polynomial recurrence
  computes every transverse root derivative.  Under THM-2241's complex
  monic polynomial-Keller hypotheses it converts the sharp response-tail
  criterion into one inverse-minimal-polynomial divisibility test.  None of
  this proves JC(2), excludes a genuine A4/S4 branch, supplies abstract Faber
  chart entry, or raises the THM-3151 degree floor.
source: root/frontier-synthesis/jc-owner-beach/2026-08-02
audit: >
  The exact companion checks the target-shear derivation change, seven
  nonlinear polynomial automorphisms and 49 implicit derivatives, a
  nonlinear sheared automorphism, the THM-3167 cubic one-jet unit, six
  independently differentiated cubic levels, 42 power-branch closed-form
  values, the polynomial Gate-II blind spot, the Laurent finite-tail near
  miss, and the compositional derivative-unit identity.  Ordinary and
  optimized executions are byte-identical to the stored output.  The proofs
  are algebraic and do not depend on the computation.
depends_on:
  - THM-2230-planar-jacobian-response-fiber-and-exact-target-shear-quotient
  - THM-2241-monic-transverse-response-depth-and-resultant-nonproper-quotient
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-3167-inverse-different-three-gate-target-shear-descent-and-full-marked-jet-no-go
related:
  - THM-2102-power-free-weight-face-and-first-defect-descent
  - THM-2240-dc2-grade-response-gauge-is-not-a-continuation-state
  - THM-3068-c3-escape-inverse-pole-ledger-and-reciprocal-cofactor-shift
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3166-falling-factorial-order-join-path-colour-transform
script: 04-computation/jc_differential_owner_filtration_thm3172.py
output: 05-knowledge/results/jc_differential_owner_filtration_thm3172.out
script_sha256: d3a15e0f435a25302234017a248341ab266ddea4c30ca5ff31fa47cc0eee4d62
output_sha256: 4265da18737d11d1537ec2f160a82f45cfc8b1fc235316951eb3d5b5c1a9f0c9
hash_basis: LF-normalized bytes
---

# THM-3172 -- shear-invariant differential owners and transverse recurrence

**PROVED + VERIFIED-EXACT.**

THM-3167 separates three predicates: diagonal inverse-Jacobian response over
the target field, membership of the diagonal scalar in the geometric constant
field, and ownership by one connected polynomial `A2` source.  It also proves
that bounded marked local jets cannot decide the first predicate.  The present
theorem addresses the third predicate instead.  Its carrier is not another
branch scalar: it is the marked order together with the target derivations
that a polynomial owner must preserve.

The closest positive mechanism is THM-2241's finite Hamiltonian response tail.
The corrected near miss is MISTAKE-301: the trace--Liouville class is exact on
every polynomial source, and THM-3068 shows that exact trace data can persist
on a punctured Laurent hostile.  The missing sidecar is the global regular
owner ring, including its connectedness and unit group.

## 1. Differential-owner filtration

Let `k` be a characteristic-zero field and put

```text
R=k[u,v],                         K=k(u,v).              (1)
```

Let `A` be a finite separable `K`-algebra and let `x,y in A` be a marked
pair.  The base derivations extend uniquely to `A`; their extensions commute.
Define

```text
B_0=k[u,v,x,y] subset A,
B_(n+1)=B_n[partial_u(a),partial_v(a):a in B_n],
B_infinity=union_(n>=0) B_n.                            (2)
```

If the marked pair is presented in THM-2621 form

```text
A=K[T]/(f),                  x=T mod f,
y=b(x),                                                    (3)
```

with `f` monic separable, the first stage is explicit.  All coefficient
derivatives below are taken at fixed `T`:

```text
x_i=-f_i/f_T,
y_i=b_i+b_T x_i=(f_T b_i-b_T f_i)/f_T,
i in {u,v}.                                               (4)
```

Thus `B_1` is computable from the complete marked inverse pair without
choosing roots in a splitting field.

## 2. Exact target-shear covariance

Let `H in k[U]` and perform the target shear

```text
u'=u,                            v'=v+H(u).              (5)
```

At fixed new target coordinates,

```text
partial_(v')=partial_v,
partial_(u')=partial_u-H'(u)partial_v.                   (6)
```

Conversely,

```text
partial_u=partial_(u')+H'(u)partial_(v').                (7)
```

Also `k[u',v']=k[u,v]`.  Hence the modules generated by the two target
derivations over the base polynomial ring are identical.  Starting with the
same `B_0`, induction in `(2)` gives

```text
B_n'=B_n for every n,                  B_infinity'=B_infinity. (8)
```

This is equality of rings, not merely equality of valuations or associated
graded pieces.

## 3. Polynomial ownership freezes the filtration

Suppose now that the marked pair comes from one polynomial map

```text
u=P(x,y),                       v=Q(x,y),
P,Q in k[x,y],                  Jac(P,Q)=kappa in k*.    (9)
```

Inverting its Jacobian matrix gives polynomial derivations

```text
partial_u=(Q_y/kappa)partial_x-(Q_x/kappa)partial_y,
partial_v=-(P_y/kappa)partial_x+(P_x/kappa)partial_y.   (10)
```

Because `B_0` contains `x,y` and `u,v` are polynomials in them,

```text
B_0=k[x,y].                                             (11)
```

Equations `(10)` preserve this ring, so induction gives

```text
B_n=B_infinity=k[x,y]                    for every n.   (12)
```

In particular every polynomial Keller owner filtration is a domain, has no
idempotents except `0,1`, and has unit group `k*`.

## 4. One-jet compositional gate

Suppose a candidate marked owner satisfies

```text
v-H(u)=Phi(z),                                           (13)
```

where `z in B_0` is nonconstant and `Phi in k[T]`.  Differentiating at fixed
`u` gives

```text
Phi'(z) partial_v(z)=1.                                 (14)
```

Thus `Phi'(z)` is a unit of `B_1`.  If the candidate has a polynomial Keller
owner, `(12)` forces this unit to lie in `k*`.  A nonconstant element of
`k[x,y]` is transcendental over `k`, so in characteristic zero this implies

```text
deg(Phi)=1.                                             (15)
```

Equivalently, for every polynomial Keller pair and every `H`, the target-shear
representative `S=Q-H(P)` is not a nontrivial polynomial composite.  Directly,

```text
Jac(P,S)=kappa=Phi'(z)Jac(P,z),                          (16)
```

and both factors on the right must be units of `k[x,y]`.

This is an exact full-polynomial owner gate.  It is not a claim that a
power-free top face proves invertibility; THM-2102's first-defect quotient and
future layers remain necessary in that filtration problem.

## 5. The THM-3167 hostile fails the owner gate at one derivative

On THM-3167's cubic component the forward map is

```text
u=-rY/(3X^2),                       t=X^3,
Jac_(X,Y)(u,t)=r.                                       (17)
```

The split inverse polynomial is

```text
f(T)=(T-u)(T^3-t).                                      (18)
```

Using `(4)` on the cubic factor gives

```text
partial_t(X)=1/(3X^2),
X^(-1)=3X partial_t(X) in B_1.                          (19)
```

The first differential owner ring therefore contains the nonconstant unit
`X`.  It cannot equal a polynomial `A2` owner.  Under `t'=t+H(u)`, equation
`(6)` gives `partial_(t')=partial_t`, so the obstruction is shear invariant.

The full fixed-plus-cubic marked algebra fails even earlier.  With
`x=(u,X)` in the product algebra,

```text
x-u=(0,X-u),                       x^3-t=(u^3-t,0),
(x-u)(x^3-t)=0,                                        (20)
```

and both displayed factors are nonzero generically.  Its marked order is not
a domain and hence cannot be a connected `A2` owner.

These obstructions hold for every `r`, including THM-3167's diagonal control
`r=1`.  They do not decide its Gate I and do not weaken the all-fixed-depth
marked-jet no-go: both members of that hostile pair fail Gate III equally.

## 6. Transverse inverse-spectral recurrence

Let `f(T;u,v)` be monic and separable and put `xi=T mod f`.  Define

```text
Lambda_f=f_T partial_v-f_v partial_T,
A_1=-f_v,
A_(n+1)=f_T Lambda_f(A_n)
          -(2n-1)A_n Lambda_f(f_T).                    (21)
```

Then for every `n>=1`,

```text
partial_v^n(xi)=A_n(xi)/f_T(xi)^(2n-1).                (22)
```

Indeed, implicit differentiation proves `(22)` for `n=1`.  For any
coefficient polynomial `h`, total differentiation along `f=0` is

```text
partial_v(h(xi))=Lambda_f(h)(xi)/f_T(xi).              (23)
```

If `(22)` holds at `n`, applying `(23)` to
`A_n f_T^(-(2n-1))` gives denominator `f_T^(2n+1)` and exactly the numerator
in `(21)`.  This proves the recurrence by induction.

For the shear `(5)`, `partial_(v')=partial_v`; applying the same substitution
to `f` and every `A_n` proves exact shear covariance of `(21)--(22)`.

## 7. Precisely scoped THM-2241 consequence

Work now over `C`.  Let

```text
P,Q in C[x,y],                   Jac(P,Q)=kappa in C*,
P monic in y of degree d>0,
u=P,                             v=Q,                    (24)
```

and let `f in C(u,v)[T]` be the inverse minimal polynomial of the chosen
source coordinate `x`.  Put `D_P=Jac(P,-)`.  On the source function field,

```text
D_P=kappa partial_v.                                      (25)
```

Combining `(22)` with THM-2241 gives

```text
(P,Q) is a polynomial automorphism
  iff f divides A_(d+1) in C(u,v)[T].                    (26)
```

On the automorphism branch THM-2241 also gives `D_P^d(x)!=0`, hence

```text
f does not divide A_d.                                  (27)
```

The complex field, polynomial Keller pair, monicity, positive `y`-degree,
and inverse-minimal-polynomial hypotheses in `(24)--(26)` are load bearing.
The congruence is a decision reformulation, not a proof that it always holds.

## 8. Falling-factorial controls and failure boundaries

For the nonlinear automorphisms

```text
P=x+(y+x^2)^m,                     Q=y+x^2,
x=u-v^m,                                                (28)
```

the recurrence gives

```text
partial_v^n(x)=-(m)_(falling n)v^(m-n),       1<=n<=m,
partial_v^(m+1)(x)=0.                                  (29)
```

For the punctured power branch `t=X^m`, `m>1`, it gives

```text
partial_t^n(X)=(1/m)_(falling n)X^(1-mn),              (30)
```

which never vanishes in characteristic zero.  At the `X=0` valuation the
orders form the arithmetic ray `1-mn`; its first negative point is the unit
obstruction `(19)`.  The shared falling-factorial syntax with THM-3166 is a
sequence lens only, not a tournament-to-Jacobian object map.  Tropicalizing
`(21)` can also lose cancellation at ties, so any Newton implementation must
retain leading residues as in THM-3161's Newton--Euclidean ledgers.

Two controls delimit `(26)`:

1. the polynomial map `u=X^4,v=XY` has nonconstant Jacobian `4u` and
   `partial_v^n X=0`, while `partial_uX=1/(4X^3)` exposes `X^(-1)`; the
   transverse tail alone does not replace the constant-field gate;
2. the Laurent pair `h=xy^2,s=-1/y` has `Jac(h,s)=1` and
   `D_h^3(x)=0`, but it lacks a polynomial `A2` mate and the monic properness
   package required by THM-2241.

The viable next test is therefore to build `B_1` from `(4)` for an actual
quartic or heptic source-lift candidate, retain every conductor/Jelonek
denominator, and print a certified idempotent, nonconstant unit, or surviving
owner order.  Before an affine source lift, the abstract Faber passport does
not define `f,b`, so this theorem alone yields no new passport exclusion.

## Exact reproduction

Run

```text
python3 04-computation/jc_differential_owner_filtration_thm3172.py
python3 -O 04-computation/jc_differential_owner_filtration_thm3172.py
```

Both modes must reproduce the stored output byte for byte.  The companion
uses explicit `require` checks, so optimized mode executes every gate.  QED.
