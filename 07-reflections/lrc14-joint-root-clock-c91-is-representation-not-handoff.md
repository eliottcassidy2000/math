---
source: codex-2026-07-27-joint-root-clock-c91-audit
status: noncanonical reflection / proved algebra plus open handoff test
tags: [lrc14, c91, tournament, cayley, sawtooth, clock, ancestry]
---

# The joint root--clock C91 object is a representation, not yet a handoff

This reflection does not promote a theorem or alter the LRC(14) proof graph.
It records an exact representation-level synthesis of
[THM-2508](../01-canon/theorems/THM-2508-affine-cut-bundle-covariance-and-carry-permutation.md),
[THM-2532](../01-canon/theorems/THM-2532-cyclic-tournament-cayley-algebra-and-chi7-even-quotient.md),
and
[THM-2535](../01-canon/theorems/THM-2535-boundary-tooth-clock-intertwiner-and-neutral-collapse.md),
then isolates why it does not yet compose with
[THM-2349](../01-canon/theorems/THM-2349-first-depth-one-delayed-shallow-restart.md).

Scope is explicit:

```text
PROVED here:       finite-group algebra, eta-cut identity, carry hostile,
                   conductor and equivariant-map no-go.
CITED canonical:  cut covariance/nonvanishing, odd-cycle Cayley inverse,
                   scheduler cells, and the marked 91-unit Fourier edge.
OPEN:              same-trajectory clock transfer, inherited source/deep
                   ancestry, and semantic arrival at the empty endpoint.
```

## 1. The four conductors explain the existing counts

Put

```text
G=F_13 x F_7 congruent Z/91Z.
```

Over `Q`, its regular representation splits exactly as

```text
Q[G]
 = Q
   direct-sum Q(zeta_7)
   direct-sum Q(zeta_13)
   direct-sum Q(zeta_91),                                  (1)

dimensions: 1, 6, 12, 72.
```

For a table `d(h,r)`, the THM-2508 row-zero law

```text
sum_r d(h,r)=0 for every h                                  (2)
```

kills the trivial and conductor-13 summands.  The common kernel of the full
cut bank is the six-dimensional pure-clock module `Q(zeta_7)`.  Therefore the
`72` mixed modes retained by the cut characters are not merely a matching
count: they are exactly the primitive conductor-91 module
`Q(zeta_91)`, whose dimension is `phi(91)=72`.

For fixed THM-2535 cut parameters `(sigma,a)`, the boundary-relative table
`Q_(sigma,a)(kappa,t)` has Fourier law

```text
Qhat(alpha,beta)=M_(sigma,a)(alpha,beta)dtilde(alpha,beta), (3)
```

where the cited geometric multiplier is nonzero whenever
`alpha,beta!=0`.  Thus `Q_(sigma,a)` is an isomorphism on the primitive
conductor-91 summand.

## 2. The exact 91-cycle Cayley and sawtooth calculus

For `tau in F_13^*` and `eta in F_7^*`, translation by

```text
g=(tau,eta)
```

is one cycle of length `91`.  With `P=P_g`, define

```text
C_91=sum_(q=1)^90 (-1)^(q+1)P^q,
B_91=1/91 sum_(q=1)^90 (2q-91)P^q.                           (4)
```

The universal odd-cycle calculation in THM-2532 gives

```text
(I+P)C_91=P-I,
C_91 B_91=B_91 C_91=I-J/91.                                 (5)
```

Because (3) is Fourier diagonal, on `Q(zeta_91)` one has the exact
intertwiners

```text
Q_(sigma,a) C_91=C_91 Q_(sigma,a),
Q_(sigma,a) B_91=B_91 Q_(sigma,a).                           (6)
```

This is the valid joint root--clock tournament statement: the full cut table
supports a lossless `91`-carousel derivative and sawtooth inverse on every
primitive mixed mode.

The two old tournaments are exact odd-fibre quotients.  If `pi_13` sums the
clock coordinate and `pi_7` sums the root coordinate, grouping the
coefficients in (4) by residue proves

```text
pi_13 C_91=C_13 pi_13,       pi_13 B_91=B_13 pi_13,
pi_7  C_91=C_7  pi_7,        pi_7  B_91=B_7  pi_7.            (7)
```

For a nonzero quotient residue, an odd number of alternating signs sums to
the corresponding sign; over the zero residue, the remaining even number
sums to zero.  The affine-ramp coefficients in `B_91` group to the smaller
affine ramp by the same direct sum.  Thus `C_13` and `C_7` are literally the
two majority quotients of one root--clock carousel, not merely analogies.

## 3. General clock step: the local output is exact

Let `s->t` be a root chord and put `tau=t-s`.  For any `eta!=0`, define

```text
D_kappa(h,r)
 =1_(h=s)[x_kappa 1_(r=kappa)
          -x_(kappa+eta)1_(r=kappa+eta)].                    (8)
```

Choose cut scale `eta^(-1)` and cut intercept `-eta^(-1)kappa`.  The two
clock labels become teeth zero and one, so the cut output is exactly

```text
R_kappa=x_kappa delta_s-x_(kappa+eta)delta_t,

sum_kappa R_kappa
 =M(delta_s-delta_t),       M=sum_kappa x_kappa.              (9)
```

Under `r->Br+C`, the data transform as

```text
eta->B eta,
kappa->B kappa+C,
eta^(-1)->eta^(-1)B^(-1),                                  (10)
```

which is precisely the THM-2508/2535 affine cut law.  The construction is
not special to the step `eta=1`.

## 4. Unequal endpoint mass is incidence, not an edge current

Equation (8) is not the boundary of one weighted directed edge unless

```text
x_kappa=x_(kappa+eta).                                      (11)
```

Its augmentation is `x_kappa-x_(kappa+eta)`.  The whole labelled family
cancels before charting only because every vertex mass occurs once with each
sign:

```text
sum_kappa D_kappa=0.                                        (12)
```

Calling each component an edge current would therefore overstate the result.
It is a two-end incidence component with unequal endpoint weights.  The
algebraic replacement

```text
x_kappa(delta_kappa-delta_(kappa+eta))                       (13)
```

is a genuine edge, but it falsely places the Boolean event `X_kappa` in the
successor scheduler cell.  Mixing gives asymptotic equality of the seven
masses, not the exact same-event transport needed by (13).

For a genuine deterministic `eta`-cycle using all the mass, conservation at
every clock vertex forces (11) for every `kappa`, hence one common weight.
If only a positive subcurrent is needed, a common
`0<w<=min_kappa x_kappa` is measure-theoretically feasible, but a physical
proof must still exhibit the same trajectory at both endpoints.

## 5. The seven moving charts are not one 91-cycle

Fix one chart origin `kappa_0` and put

```text
j(r)=rep(eta^(-1)(r-kappa_0)) in {0,...,6}.                  (14)
```

For seven same-root clock pairs, the fixed chart produces six root steps
`tau`; at the wrap it produces `-6tau`, not `tau`.  The correct `91`-cycle
requires thirteen root layers.  For `q=7ell+j`, use source cells

```text
u_(7ell+j)=(s+7ell tau, kappa_0+j eta),
0<=ell<=12, 0<=j<=6.                                        (15)
```

The fixed cut sends them to

```text
z_q=(s+q tau,kappa_0+q eta),                                (16)
```

which is the genuine `(tau,eta)` cycle.  In source coordinates its carry is

```text
(h,r_j)->(h,r_(j+1))          for j<6,
(h,r_6)->(h+7tau,r_0)         at the wrap.                   (17)
```

THM-2535 has seven scheduler intersections at one selected root, not the
thirteen root layers in (15).  Moreover, conjugating by this cut shear makes
(17) fail to preserve the original per-root row-zero subspace.  Thus there
are two honest but different structures:

```text
ordinary joint translation:  preserves row-zero, but is not the local
                             same-root clock incidence;

cut-conjugate carry:          realizes the local cut geometry, but is not
                             an internal symmetry of the lawful defect.   (18)
```

## 6. The conductor no-go and the THM-2349 type mismatch

The THM-2535 output `delta_s-delta_t` lies in the conductor-13 root module.
But

```text
Hom_G(Q(zeta_91),Q(zeta_13))=0.                              (19)
```

Hence no translation-equivariant linear combination of `C_91`, `B_91`, or
the primitive cut intertwiner can manufacture the root dipole from the
lawful primitive defect.  The dipole is possible precisely because the
transported charts are summed after, rather than before, evaluation.  It is
chart holonomy, not a hidden consequence of the `91`-cycle spectral algebra.

[THM-2349](../01-canon/theorems/THM-2349-first-depth-one-delayed-shallow-restart.md)
provides a genuine ordinary-frequency edge

```text
Y-X=m c_3,                  gcd(m,91)=1.                      (20)
```

Arithmetic compatibility is exact:

```text
(Z/91Z)^* congruent F_13^* x F_7^*,                          (21)
```

so `m` can reparameterize a retained joint cycle by `P_g->P_g^m`.  But `m`
is an edge multiplier/automorphism, not yet a physical root--clock step.
The endpoints `X,Y` have the same prescribed root character because
`c_3/13=0 mod 13`.  THM-2349 does not prescribe `m` and proves neither

```text
m mod 13=t-s,                  m mod 7=eta                    (22)
```

nor any same-event map from its marked Fourier edge to a scheduler handoff.
Affine scaling cannot repair a failed equality in (22): both candidate
steps scale together.  This agrees with
[THM-2327](../01-canon/theorems/THM-2327-two-colour-marked-unit-c3-triangle.md),
whose two-colour triangle is deliberately undirected and whose word mark is
lost in a bare tournament shadow.

## 7. Minimal repair and the next decisive test

The smallest missing object is not another Fourier mode.  It is one positive
Boolean orbit segment that visits the seven scheduler cells in cyclic order.
For a positive selector/scheduler event

```text
X_kappa=E_s intersection P_kappa^L,
```

test separated delays `N_0<...<N_6` and

```text
Y_(kappa_0,eta)
 = intersection_(j=0)^6 T^(-N_j)X_(kappa_0+j eta).           (23)
```

The precise next test is:

1. On one explicit full-support row from
   [THM-2517](../01-canon/theorems/THM-2517-cubic-anchored-spectrum-flt3-and-three-time-boolean-lift.md),
   compute or prove `measure(Y_(kappa_0,eta))>0` for a lacunary delay bank.
2. Use the one common mass `y=measure(Y)` on all seven consecutive temporal
   edges and verify that the generalized cuts give
   `7y(delta_s-delta_t)` without unequal endpoint weights.
3. Check whether the two endpoints retain the inherited
   [THM-2449](../01-canon/theorems/THM-2449-coprime-owner-anova-and-delta-replica-boundary.md)
   source/deep coordinate and whether the negative endpoint is an actual
   later Boolean arrival.

A positive first two steps would close the clock same-ancestry seam.  Failure
of the third would isolate arrival semantics as the remaining obstruction.
If even (23) cannot be made positive with the required labels, the scheduler
clock and selector ancestry are genuinely incompatible, and the `91`-cycle
route should be stopped rather than enlarged.

LRC(14) remains **OPEN**.
