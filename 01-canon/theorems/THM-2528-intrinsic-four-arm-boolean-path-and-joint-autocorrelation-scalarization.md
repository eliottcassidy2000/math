---
id: THM-2528
title: "Intrinsic four-arm Boolean path and joint autocorrelation scalarization"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.  Before the live
  A_tau H_tau odd bank is collapsed to coefficients of sizes 1,3,5,7, its
  two convolution stages have coefficients only +/-1.  Retaining those
  stages as two genuine intermediate predecessor addresses realizes every
  odd coordinate as the difference of two disjoint Boolean unions on the
  ordered four-arm depth-one fibre product, with no copy tag.  At the fixed
  THM-2527 coordinate -4 tau the positive union is strictly larger whenever
  the late owner sees a nonconstant root fibre; one matched pair of lawful
  collision intersections already has a strict imbalance.  The shallow,
  deep, terminal-word, and late-owner sidecars remain on the same ancestry
  object.  Independently, two-dimensional autocorrelation of the complete
  THM-2521 pulled-back K_14 potential scalarizes all 72 mixed Hilbert modes:
  its normalized Fourier coefficients are ||d_tilde||^2/91>0.  Applying the
  guard Hilbert operator in the root displacement makes a joint-converse-odd
  scalar table without losing a mode.  Doing the same only after retaining
  the full THM-2508 affine cut bank gives 5,184 indexed scalar norm
  certificates.  Forgetting the intermediate addresses merges the Boolean
  path components, while the THM-2527 mask layers are converse-invariant;
  neither construction identifies source time with arrival time or emits a
  signed scalar-cover current.  No row exclusion or LRC(14) proof follows.
source: codex-2026-07-27-intrinsic-boolean-path-joint-autocorrelation
depends_on:
  - THM-2508-affine-cut-bundle-covariance-and-carry-permutation
  - THM-2518-perron-inverse-branch-owner-word-cospan-recovery
  - THM-2521-k13-drift-k14-potential-module-bridge
  - THM-2522-intrinsic-collision-depth-toothpick-descent-and-late-owner-decoupling
  - THM-2526-affine-skew-orientation-gauge-boundary
  - THM-2527-owner-weighted-all-mode-odd-bank-and-boolean-cut-coordinate
related:
  - THM-2471-owner-first-collision-weighted-root-service-and-temporal-atom-boundary
  - THM-2478-delayed-owner-handoff-graft-and-deep-sheet-rebase-boundary
script: 04-computation/lrc14_intrinsic_path_joint_autocorrelation_thm2528.py
output: 05-knowledge/results/lrc14_intrinsic_path_joint_autocorrelation_thm2528.out
script_sha256: 335dadbc6f2837de8e0e8705ec4ab7fd3c651cd77a7d3102dc89ddf9e24884ac
output_sha256: d296c095e94c86107f5f723bd940517ade505f483d151bf6ff6e0140aaa753ec
hash_basis: working-tree bytes (LF)
---

# THM-2528 -- keep the operator path, then scalarize the complete table

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.**

THM-2527 closes one Booleanization seam: on the live Boolean sheet the
fixed coordinate `-4 tau` of the owner-weighted odd bank is a positive sum
of full-root-mask events.  Two sharper structural questions remain.

```text
Can the original signed collision summands themselves be separated
without adding a formal multiplicity tag?

Can THM-2521's 72/5,184 Hilbert-valued mixed modes be made into
nonzero scalar modes without integrating a phase that may cancel?          (1)
```

Both answers are yes, with different carriers.  The first answer is the
ordered four-arm inverse-branch path already allowed by THM-2518.  The
second is the two-dimensional Wiener--Khinchin autocorrelation of the
whole potential or cut table.  Neither answer is a same-time semantic
owner current; their exact quotient losses are recorded below.

## 1. Boolean depth-one collision notation

Let `Z` be the future base with its probability measure, let

```text
e_r(z) in {0,1},                         r in F_13,             (2)
```

be the thirteen predecessor occupancies, and let `g(z) in {0,1}` be a
common future owner or owner--word factor.  In the live application,

```text
e_r(z)=F((z+r)/13),
g(z)=G(13^R z).                                                (3)
```

Put

```text
c_u(z)=sum_r e_r(z)e_(r+u)(z),
B_u=1/13 integral_Z g(z)c_u(z)dz,
Bbar=1/13 sum_u B_u,
b_u=B_u-Bbar.                                                 (4)
```

Thus `B` is exactly the THM-2527 owner-weighted depth-one collision
profile.  It is even, and positive weighted replica defect means

```text
b_0>0
 iff g gives positive measure to nonconstant root fibres.                    (5)
```

Fix the guard-transported slope `tau`.  Write the two circulant kernels as

```text
(A_tau x)_t=sum_(v!=0)a_tau(v)x_(t+v),
(H_tau x)_t=sum_(h!=0)eta_tau(h)x_(t+h),                       (6)
```

where, for `s=1,...,6`,

```text
a_tau(+/-2tau s)=-chi_7(s),
eta_tau(-2tau s)=+1,
eta_tau( 2tau s)=-1.                                         (7)
```

Every nonzero residue occurs exactly once in each kernel.  In particular,

```text
a_tau(-v)=a_tau(v),               eta_tau(-h)=-eta_tau(h),
a_tau(v),eta_tau(h) in {+1,-1}.                                (8)
```

THM-2527's owner-weighted odd bank is

```text
O_tau=13A_tau H_tau b=13A_tau H_tau B.                         (9)
```

The second equality uses that both factors kill constants.

## 2. The uncollapsed operator is an intrinsic four-arm Boolean path

For fixed output displacement `t`, form the ordered root path

```text
r_0=r,
r_1=r+t,
r_2=r+t+v,
r_3=r+t+v+h,                 v,h in F_13^*.                  (10)
```

Its four physical predecessors are

```text
x_i=(z+r_i)/13,                         i=0,1,2,3.             (11)
```

All four have the same future `13x_i=z`.  The tuple in (10) recovers
`r,v,h` uniquely: this is a genuine component of the ordered four-fold
fibre product, not a collision term plus a newly invented copy label.

Give

```text
Z x F_13 x F_13^* x F_13^*                                  (12)
```

the product of base measure and uniform counting measure.  Define two
Boolean subsets

```text
P_t={g(z)=1, e_(r_0)(z)=e_(r_3)(z)=1,
     a_tau(v)eta_tau(h)=+1},

N_t={g(z)=1, e_(r_0)(z)=e_(r_3)(z)=1,
     a_tau(v)eta_tau(h)=-1}.                                  (13)
```

They are disjoint unions of literal Boolean collision intersections on
distinct ancestry components.  Expanding (9), then using (4), gives

```text
O_tau(t)
 =sum_(r,v,h)a_tau(v)eta_tau(h)
    integral_Z g e_r e_(r+t+v+h)

 =13*12^2 [measure(P_t)-measure(N_t)].                         (14)
```

The factor `13*12^2=1872` is only the normalization of the three discrete
address coordinates.  Equation (14) is the desired Boolean split before
the convolution paths with the same endpoint are combined into the
weights `1,3,5,7` of THM-2526 equation (44).

THM-2527 proves, whenever (5) holds,

```text
O_tau(-4tau)>0.                                               (15)
```

Consequently the fixed, predeclared path coordinate satisfies

```text
measure(P_(-4tau))>measure(N_(-4tau)).                         (16)
```

This is a strict inequality between two disjoint Boolean unions on an
already-lawful ancestry space.  No positivity was inferred from one
arbitrarily chosen signed term, and no extra Bernoulli or colour tag was
inserted.

## 3. One matched Boolean wedge already has the strict sign

Let

```text
H_tau^+={h in F_13^*:eta_tau(h)=+1}.                          (17)
```

Pair the paths with increments `h` and `-h`.  Equation (9) becomes

```text
O_tau(t)
 =13 sum_(v!=0) sum_(h in H_tau^+)
    a_tau(v)[B_(t+v+h)-B_(t+v-h)].                            (18)
```

Therefore (15) forces some `(v,h)` with

```text
a_tau(v)[B_(-4tau+v+h)-B_(-4tau+v-h)]>0.                     (19)
```

This can be refined to one absolute predecessor sheet.  Put

```text
E_(r,u)={z:g(z)=1,
             e_r(z)=e_(r+u)(z)=1}.                            (20)
```

Since

```text
B_u=1/13 sum_r measure(E_(r,u)),                              (21)
```

equation (19) gives one `r` for which the same strict inequality holds
between the two corresponding Boolean component measures.  Their endpoint
addresses differ by `2h!=0`.  Thus the global path imbalance is not hiding
behind a cancellation in every fixed root sheet.

One of the two gaps in (19) can be zero.  This is sharp: a centred-delta
root mask has all six nonzero antipodal correlations equal, so every
strict comparison uses the identity component.  No theorem may uniformly
upgrade (19) to two distinct first-collision edges without an additional
live-mask hypothesis.

## 4. Relation to the positive THM-2527 mask layers

For a Boolean mask `e`, THM-2527 defines

```text
Psi_tau(e)=(A_tau H_tau c(e))_(-4tau)                         (22)
```

and proves

```text
Psi_tau(e) in {0,...,98},
Psi_tau(e)=0 iff e is empty or full,
42Psi_tau(e)>=n(13-n),              n=sum_r e_r.              (23)
```

The relative floor has a transparent discrete mechanism.  On a mixed
mask the integer gap gives `Psi>=1`, while

```text
n(13-n)<=42.                                                  (24)
```

Thus (23)'s `13/42` drift floor is exactly the integer gap plus the largest
possible thirteen-point cut.  The threshold decomposition

```text
Psi_tau=sum_(j=1)^98 1_(Psi_tau>=j)                           (25)
```

is a Boolean function on the full thirteen-arm fibre.

That positive factorization and the path factorization (13) remember
different information.  If `e^c=1-e`, then

```text
c_u(e^c)=13-2n+c_u(e).                                       (26)
```

Since `A_tau H_tau` kills constants, `Psi_tau(e^c)=Psi_tau(e)`.
Root translation and reflection also leave every `c_u` unchanged, and
simultaneously replacing `tau` by `-tau` transports (22) to the same
positive score.  Hence every threshold layer in (25) is converse-blind.

The four-arm split instead keeps the two signed operator increments.  It
therefore retains the odd path ledger in (16), at the cost of living on the
natural-extension address space.  The two constructions should not be
identified:

```text
positive full-mask layer:      Boolean downstairs, converse-blind;
unequal path unions P_t,N_t:   Boolean upstairs, address-dependent.          (27)
```

## 5. Shallow, deep, word, and late-owner sidecars survive the lift

For every live THM-2349 row, THM-2522 supplies

```text
F=1_(E_j)(x)1_Q(13^k x),             k>=2,
c_j=13a,                    13 does not divide a,
nu_13(c_3)>1.                                                  (28)
```

Every root displacement `w/13` obeys

```text
c_j(x+w/13)-c_jx=aw in Z,
c_3(x+w/13)-c_3x in Z.                                       (29)
```

Thus the old shallow septimal and deep thirteenth phase labels are
unchanged on all four addresses in (10), not merely recovered after a
sum.  The terminal word is also common because

```text
13^k(x+w/13)=13^kx                         mod 1.              (30)
```

Finally `g(z)=G(13^Rz)` is literally the same late-owner factor on every
arm.  Equations (13), (19), and (20) can therefore be refined by all these
retained labels without altering their inequalities.

The intermediate addresses `r_1,r_2` are operator-path witnesses.  No
packet is asserted to occur there.  In particular, their presence does not
turn the first arm into a source-time atom and the fourth into an
arrival-time atom; THM-2471's temporal adjoint boundary remains in force.

## 6. Two-dimensional autocorrelation scalarizes all 72 potential modes

There is an orthogonal scalarization of the complete THM-2521 potential.
Let

```text
mathcal H=L^2(g(z)dz),
d(h,r)=d_p(h,r) in mathcal H,          (h,r) in F_13 x F_7,   (31)
```

be the pulled-back `K_14` potential in one retained affine chart.  Use the
unnormalized table transform

```text
d_tilde(alpha,beta)
 =sum_(h,r)d(h,r)zeta_13^(-alpha h)zeta_7^(-beta r).          (32)
```

Define its scalar two-dimensional autocorrelation by

```text
C(t,s)=sum_(h,r)<d(h,r),d(h+t,r+s)>_mathcal H.                (33)
```

Normalize only the Fourier transform:

```text
C_hat(alpha,beta)
 =1/91 sum_(t,s)C(t,s)
        zeta_13^(-alpha t)zeta_7^(-beta s).                   (34)
```

Finite Wiener--Khinchin gives the exact identity

```text
C_hat(alpha,beta)=1/91 ||d_tilde(alpha,beta)||_mathcal H^2.   (35)
```

The factor is `1/91`, not `1/91^2`, because the autocorrelation in (33)
is an unnormalized sum while the transform in (34) is normalized.

THM-2521 equation (36), in chart `(tau,a,c)`, says

```text
d_tilde(alpha,beta)
 =zeta_7^(beta a^(-1)c)
  L(alpha tau,beta a^(-1))p_hat(alpha).                       (36)
```

Moreover the normalized collision coefficient satisfies

```text
B_hat(alpha)=1/169 ||p_hat(alpha)||_mathcal H^2.              (37)
```

Combining (35)--(37) gives the scalar formula

```text
C_hat(alpha,beta)
 =13/7 |L(alpha tau,beta a^(-1))|^2 B_hat(alpha)>0
                         for alpha,beta!=0.                   (38)
```

The last inequality uses THM-2527's owner-weighted all-colour conclusion.
Thus all `12*6=72` formerly Hilbert-valued mixed modes have nonzero scalar
norm certificates on one table.

## 7. The guard Hilbert transform makes the scalar table joint-odd

Apply the guard-selected cyclic Hilbert operator only in the root
displacement:

```text
Q(t,s)=(H_(tau_H)^(t) C)(t,s).                               (39)
```

Its nontrivial root multiplier is

```text
m_(tau_H)(alpha)
 =(zeta_13^(alpha tau_H)-1)
   /(zeta_13^(alpha tau_H)+1)!=0.                             (40)
```

Therefore

```text
Q_hat(alpha,beta)
 =m_(tau_H)(alpha)C_hat(alpha,beta)!=0
                                      for alpha,beta!=0.      (41)
```

The symmetry type is exact.  Autocorrelation gives

```text
C(-t,-s)=C(t,s),                                             (42)
```

and the odd root kernel then gives

```text
Q(-t,-s)=-Q(t,s).                                            (43)
```

This is **joint-converse oddness**.  In general neither
`C(-t,s)=C(t,s)` nor `Q(-t,s)=-Q(t,s)` holds.  Applying `H` has not
manufactured a pure root orientation after forgetting the septimal
displacement; it has changed the joint-even autocorrelation into a
joint-odd scalar bank while preserving every mixed mode.

## 8. The full affine cut bank scalarizes all 5,184 indexed modes

The same construction reaches the larger THM-2508 bank only after its cut
coordinate is retained.  For each

```text
sigma in F_13^*,                  a in F_7^*,                 (44)
```

let

```text
R_(sigma,a,c)^d(v),               (v,c) in F_13 x F_7,       (45)
```

be the full affine-cut table of `d`.  Define

```text
C_(sigma,a)(t,s)
 =sum_(v,c)<R_(sigma,a,c)^d(v),
             R_(sigma,a,c+s)^d(v+t)>_mathcal H.              (46)
```

If `Psi_(sigma,a)(alpha,beta)` is THM-2508's transform of (45), the same
calculation gives

```text
C_hat_(sigma,a)(alpha,beta)
 =1/91 ||Psi_(sigma,a)(alpha,beta)||^2.                       (47)
```

THM-2508's factorization is

```text
Psi_(sigma,a)(alpha,beta)
 =K(alpha sigma,beta)d_tilde(alpha,-beta a),                  (48)
```

and both factors are nonzero on primitive mixed modes.  Hence (47) is
strictly positive for all

```text
(sigma,a,alpha,beta)
 in F_13^* x F_7^* x F_13^* x F_7^*,                         (49)
```

exactly

```text
12*6*12*6=5,184                                              (50)
```

indexed scalar norm certificates.  Applying `H_(tau_H)` in `t` again
makes each chart table joint-odd and loses none of these modes.

The order of operations is load-bearing.  Summing away the affine cut
translation takes `beta=0`, recreating THM-2508's pure-root kernel.  The
`5,184` statement is a family over the retained `(sigma,a)` charts, not a
claim that one cut-free scalar has `5,184` independent coordinates.

## 9. Boolean provenance and the exact remaining quotient

The potential scalarization has honest Boolean-before-scalar provenance.
On a live root fibre,

```text
13p_i=sum_(j in F_13)(e_i-e_j).                              (51)
```

Every entry of `d_p`, every affine-cut entry, and hence every product in
(33) or (46) expands after clearing powers of thirteen into a signed
integer sum of lawful two-arm Boolean collision intersections.  This proves
that the scalar norm certificates have not introduced arbitrary real
weights.

It does **not** give them the stronger disjoint-union interpretation of
(13).  Potential centering, `K_14` incidence, and the affine-cut sum can
send several operator routes to the same physical endpoint pair.  Keeping
all those incidence routes would give a larger Boolean path space, but
forgetting them recreates signed multiplicity.  The theorem makes the
no-external-tag claim only for the explicit two-stage `A_tau H_tau` path,
where the two intermediate physical addresses recover the route uniquely.

The three available scalar/Boolean objects now have a sharp ledger.

1. THM-2527's positive mask layers are Boolean on the future base but are
   invariant under converse.
2. The unequal sets `P_t,N_t` are Boolean and chiral on the retained
   four-arm ancestry path, but their projections overlap after the path is
   forgotten.
3. The joint autocorrelations are genuine nonzero scalar tables with all
   `72/5,184` modes, but they are quadratic norm shadows and retain only
   joint-converse oddness after `H` is applied.

None identifies `P_omega(x)` with `P_omega(T^Kx)`, turns an intermediate
address into a semantic source or arrival event, or couples the norm
certificate to THM-2365's signed scalar-cover current.  No live row is
excluded and LRC(14) remains open.

## 10. Exact dependency-free referee

Run

```bash
python3 04-computation/lrc14_intrinsic_path_joint_autocorrelation_thm2528.py
python3 -O 04-computation/lrc14_intrinsic_path_joint_autocorrelation_thm2528.py
```

Both runs reproduce

```text
05-knowledge/results/lrc14_intrinsic_path_joint_autocorrelation_thm2528.out
```

byte-for-byte.  The exact companion:

- checks (14) at every coordinate of all `8,192` Boolean root masks;
- checks the fixed path imbalance, all `8,190` matched-pair extractions,
  `1,872` distinct four-arm components, and every live sidecar shift;
- independently reproduces THM-2527's integer range, sharp floor, and
  `52` equality masks;
- checks the two-dimensional Wiener--Khinchin identity, normalized
  `1/91` factor, joint-even/joint-odd laws, and all `72` potential modes;
  and
- checks, on an independent centred-delta control, all `5,184` affine-cut
  coefficients, norm coefficients, and odd multipliers in the exact
  order-`91` subgroup of `F_547^*`.

Nonvanishing modulo `547` is an exact certificate: `2^6 mod 547` has order
`91`, so a zero cyclotomic integer would map to zero under this evaluation.
The universal Wiener--Khinchin and path identities are the finite reindexing
proofs above; the modular calculation is an independent arithmetic control.
An independent line audit rederived the `1,872` path normalization and
fixed-root wedge extraction, checked every retained sidecar and quotient
claim, verified the `1/169`, `13/7`, and `1/91` Fourier normalizations, and
confirmed that the scalar bank is joint-odd rather than root-only odd.  It
also reproduced normal, optimized, and stored output byte-for-byte and
checked both recorded hashes.
**QED.**
