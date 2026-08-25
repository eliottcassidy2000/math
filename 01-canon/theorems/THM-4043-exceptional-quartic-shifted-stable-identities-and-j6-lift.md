---
id: THM-4043
title: "Exceptional-quartic shifted stable identities and actual J6 lift"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Multiplication
  of an arbitrary target two-form by the stable target coordinate w shifts
  the universal THM-3683 order-six identity down the coefficient ladder.
  On the complete zero-fourth parabola this gives, in particular,
  Lambda(J_3)+R_4(J_1)=0 and
  Lambda(J_5)+R_4(J_3)+R_2(J_1)=0.  The THM-3688 actual lift therefore
  continues stagewise through J_5; the unshifted exceptional-root identity
  then continues it through J_6.  This holds over the common irreducible
  quartic field and hence at all four embeddings.  No coherent all-order
  series, degree control, algebraization, global pair, Keller map, or JC(2)
  conclusion is asserted.
source: jc2-double-zero-rebuild-20260824 / stable-identity shift, 2026-08-24
audit: >
  PASS -- two independent audits checked the pullback coefficient shift,
  target-ring typing, finite truncation, negative-index convention, and the
  F_6/F_7 stage indices.  An independent exact QQ(r) reconstruction of the
  complete 18-by-168 retained odd window recovered rank 15, nullity three,
  a one-dimensional J_5 projection, and the shifted R_2/R_4/Lambda
  representative column by column; the off-parabola Q_6 hostile has rank 16
  and no J_5-active relation.  Normal
  and optimized production executions byte-match the frozen transcript.
depends_on:
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3688-russell-cylinder-exceptional-quartic-actual-j1-j2-lift
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
related:
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-4039-exceptional-quartic-j3-lift-and-adjacent-gate-rigidity
  - THM-4046-exceptional-quartic-j7-lift-and-j8-obstruction
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_shifted_stable_j6_lift_thm4043.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_shifted_stable_j6_lift_thm4043.out
script_sha256: ef76a4aefe213c63ff4ce40d97fb57f6a9cf1b6ea90a7a4032b57c6e9c462de3
output_sha256: 7566341e1eff4664ddd19ae2332455b10375d7ae730d19da530097f2a50055a9
hash_basis: raw LF bytes
---

# THM-4043 -- the exceptional actual lift reaches `J_6`

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  The exact
`J_3` cancellation of THM-4039 is part of a shorter structural mechanism:
the already-proved universal order-six identity can be shifted by powers of
the target coordinate `w`.  The same mechanism clears `J_5`, and the
unshifted identity then clears `J_6` at the exceptional quartic roots.

All rings below have characteristic zero.

## 1. The universal identity and its stable shifts

On THM-3677's complete zero-fourth parabola, let `Phi_r` be the quadratic
Russell-cylinder fold.  Its stable target coordinate is literally the source
parameter:

```text
Phi_r^*(w)=t.                                             (1)
```

For an arbitrary target two-form `Omega`, write

```text
Phi_r^*(Omega)=sum_(n>=0) t^n J_n(x) dx^dt.              (2)
```

THM-3683 proves the universal identity

```text
0=Lambda(J_6)+R_4(J_4)+R_2(J_2)+R_0(J_0),               (3)
```

where `Lambda,R_4,R_2,R_0` are the explicit retained Taylor functionals in
that theorem.  Their denominators are fixed integers, so `(3)` specializes
at every value of the parabola parameter, including all four exceptional
quartic roots.

For `0<=k<=6`, multiplication by `w^k` is an operation inside the same
target two-form module.  Equation `(1)` gives

```text
Phi_r^*(w^k Omega)
  =sum_(n>=k) t^n J_(n-k)(x) dx^dt.                     (4)
```

Put `J_m=0` for `m<0`.  Applying `(3)` to `w^k Omega` proves the entire
shifted ladder

```text
0=Lambda(J_(6-k))+R_4(J_(4-k))
                    +R_2(J_(2-k))+R_0(J_(-k))          (5)
```

for every arbitrary target two-form.  No closedness or decomposability
hypothesis is used.  The odd cases verified independently by the companion
are

```text
k=1:  0=Lambda(J_5)+R_4(J_3)+R_2(J_1),                 (6)
k=3:  0=Lambda(J_3)+R_4(J_1),                          (7)
k=5:  0=Lambda(J_1).                                   (8)
```

The production companion independently rebuilds the complete odd retained
window

```text
J_1 through x-degree 2:       9 rows,
J_3 through x-degree 1:       6 rows,
J_5 values:                   3 rows,                  (9)
```

against all

```text
3 binom(8,3)=168                                     (10)
```

target two-form monomials of total target degree at most five.  Exact
elimination over `Q(r)` gives rank `15` and nullity `3`.  The projection of
the left cokernel to the `J_5` block is one-dimensional.  The shifted
representative `(6)`, normalized to `Lambda` on `J_5`, has `R_4` on `J_3`
and `R_2` on `J_1`.  Every one of the 168 column residuals vanishes
identically as a polynomial in `r`.  Thus the computation verifies `(6)` at
every specialization, rather than only generically.  The off-parabola `Q_6`
window has rank `16` and no `J_5`-active relation, a hostile showing that the
parabola hypothesis is load-bearing.

## 2. The `J_3` and `J_4` stages are structural

THM-3688 gives actual target coefficients through `F_3,G_3` with

```text
J_0=1,                         J_1=J_2=0.               (11)
```

Take the finite actual target pair with provisional `F_4=G_4=0` and apply
`(7)` to its two-form `dF^# wedge dG^#`.  Equation `(11)` gives

```text
Lambda(J_3)=0.                                          (12)
```

This supplies a structural explanation of THM-4039's much larger retained-
jet and full-polynomial exact cancellation.  By THM-3737,
`L_0(S^2)=ker Lambda`, so `(12)` supplies actual `F_4,G_4` with `J_3=0`.

Now set `F_5=G_5=0` provisionally.  The `k=2` case of `(5)` is

```text
0=Lambda(J_4)+R_4(J_2)+R_2(J_0).                       (13)
```

The two value coefficients of `R_2` sum to zero and every other term of
`R_2` is a positive derivative, so `R_2(1)=0`.  Equations `(11)` and `(13)`
therefore give `Lambda(J_4)=0`.  THM-3737 supplies actual `F_5,G_5` with
`J_4=0`.  Thus the shift ladder independently reaches `J_4`; it does not
replace THM-4039's frozen-certificate audit or its solution-rigidity results.

## 3. The actual lift through `J_5`

Start with the stagewise actual lift through `J_4` supplied by Section 2.
Set the next provisional coefficients `F_6=G_6=0` and apply `(6)` to the
resulting finite actual target pair.  Its earlier identities give

```text
J_1=J_3=0,
```

so

```text
Lambda(J_5)=0.                                          (14)
```

The new order-six target coefficients enter the fifth source equation as

```text
J_5=D_5+6L_0(F_6,G_6).                                 (15)
```

Equations `(14)`, `(15)`, and THM-3737 therefore supply actual restrictions,
and hence actual target representatives, `F_6,G_6` for which

```text
J_5=0.                                                  (16)
```

This proof does not select or bound `F_6,G_6`; it needs only the exact,
cutoff-free image theorem.

## 4. The exceptional root clears `J_6`

Let `alpha` be the class of the debt-quartic parameter in the common field
`K`.  By definition `F(alpha)=0 in K`, and therefore

```text
R_0(1)=-256 F(alpha)/1594323=0.                         (17)
```

This field identity then holds after each of the four complex embeddings.

After choosing `(16)`, set `F_7=G_7=0` provisionally and apply the unshifted
identity `(3)` to that finite actual pair.  Since

```text
J_0=1,                         J_2=J_4=0,               (18)
```

equations `(3)`, `(17)`, and `(18)` give

```text
Lambda(J_6)=0.                                          (19)
```

The new order-seven target coefficients enter as

```text
J_6=D_6+7L_0(F_7,G_7).                                 (20)
```

Another application of THM-3737 supplies actual `F_7,G_7` that kill the
complete debt.  Combining the stages proves

```text
boxed: J_0=1,              J_1=J_2=J_3=J_4=J_5=J_6=0. (21)
```

All operations were performed over the irreducible quartic field, so `(21)`
holds uniformly after all four complex embeddings.

## 5. Finite typing and strict boundary

For each application, use only the finite polynomial in `w` containing the
already chosen target coefficients and set the next coefficient to zero.
Coefficient `J_n` depends on finitely many such terms, so higher powers of
`w` cannot affect the asserted equation.  Multiplication by `w^k` in `(4)`
is therefore an exact target-ring operation, not division by `w`, a formal
Laurent step, or an untyped source substitution.

The shift argument only moves the fixed order-six identity downward.  It
does not produce an identity at `J_7` or any higher unshifted order.  This
theorem proves neither

- a compatible all-order sequence or degree bounds for the stagewise lifts;
- convergence or algebraization of a formal sequence;
- a positive global pair on the Russell cylinder;
- a polynomial Keller map or its noninjectivity; nor
- a counterexample to `JC(2)`, which remains open.

Downstream THM-4046 proves a new order-eight identity, reaches `J_7`, and
shows that no actual order-nine choice can clear `J_8`; it does not
retroactively add an order-eight identity to this theorem.  **QED.**
