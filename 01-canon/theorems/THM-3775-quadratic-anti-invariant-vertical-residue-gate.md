---
id: THM-3775
title: "Quadratic anti-invariant vertical residue gate"
status: >
  PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.  For a smooth
  planar polynomial with complete rational Hamiltonian torsor and a quadratic
  generic-fibre model, the canonical response primitive is anti-invariant
  under the sheet involution.  At every reduced, component-exhaustive
  specialization, a nonzero anti-invariant leading pole coefficient cannot
  be canceled by a target shear, because target shears occupy only the
  diagonal trivial component representation.  Globally, if every pole fibre
  has such a model, a polynomial mate exists exactly when the canonical
  anti-invariant primitive is already polynomial.  The result recovers the
  vertical obstructions in THM-3758, THM-3765, and THM-3772, but makes no
  unproved higher-cyclic claim and has no canonical force before audit.
source: root + jc_quartic_c3_construct / 2026-08-23
audit: >
  SELF-AUDITED PROOF CANDIDATE.  The exact companion verifies the universal
  anti-commutation signs in a generic quadratic algebra, the diagonal/sign
  component decomposition, split and nonsplit controls, the exhaustive
  special-fibre factorizations and signs for THM-3765/3772, and THM-3758's
  separate two-component sign vector.  Normal and optimized runs byte-match
  the frozen transcript; independent hostile audit remains open.
depends_on:
  - THM-3770-vertical-principal-part-equalizer-and-log-canonical-dressing-gate
related:
  - THM-3758-quadratic-radial-carrier-rational-exact-split-fibre-nonentry
  - THM-3765-normalized-three-consecutive-charge-radial-keller-classification
  - THM-3772-variable-flank-three-charge-rational-exact-classification
script: 04-computation/jc2_quadratic_anti_invariant_vertical_residue_gate_thm3775.py
output: 05-knowledge/results/jc2_quadratic_anti_invariant_vertical_residue_gate_thm3775.out
script_sha256: 5082248e31d66f2c331951347b6cec90edd15988ac8f84c754ce77e2c46ec2ce
output_sha256: 3c2088e7b3147059afe99456a179fb0f1980a25d540a570d8668178434ac9975
semantic_sha256: ddc682a54eff9c8e47357b07ee76f34910538db51d2de96a50fbbdcdddbfcb4a
hash_basis: raw LF bytes
---

# THM-3775 -- target shears see the trivial sheet representation only

**PROVISIONAL PROOF CANDIDATE UNDER AUDIT + VERIFIED-EXACT.**  The opposite
principal coefficients in the recent rational-exact near misses are not
coincidental signs.  They are the nontrivial character of a quadratic sheet
involution.  A target shear is constant on every component over one target
value, so it lives in the diagonal trivial line and cannot touch that
character.

The theorem has two layers.  The first is the universal component quotient
behind THM-3770.  The second gives exact hypotheses under which a generic
quadratic involution specializes to that quotient without losing a boundary
component.

Work over an algebraically closed field `k` of characteristic zero.  Put

```text
A=k[x,y],                    K=k(x,y),
J(F,G)=F_xG_y-F_yG_x.                                    (1)
```

## 1. The principal coefficient is a component representation

Let `Q in A` have no critical point, let `lambda in k`, and write

```text
q=Q-lambda=prod_(i in I) f_i.                           (2)
```

The factors are distinct because a repeated component would make `dQ`
vanish there.  Let `D_i` be `f_i=0`.  If a rational function `P` has pole
order at most `r` on every `D_i`, define its leading component coefficient

```text
a_i=bar(q^rP) in k(D_i),                                (3)
```

with `a_i=0` when the pole order is smaller.  The coefficient packet is

```text
a(P;lambda,r)=(a_i)_i in C_lambda:=prod_i k(D_i).       (4)
```

A target shear with leading Laurent term `h q^(-r)`, `h in k`, adds the
single diagonal vector

```text
h 1=(h,...,h) in k 1 subset C_lambda.                  (5)
```

Therefore the order-`r` poles can all be raised if and only if `(4)` is one
common scalar.  Equivalently, its class

```text
[a] in C_lambda/(k 1)                                  (6)
```

must vanish.  This is the one-stage form of THM-3770's descending DVR
equalizer, and it is proved directly by multiplying `P-hq^(-r)` by `q^r`
and reducing on every `D_i`.

When all `a_i` are scalars, `(4)` lies in the permutation module `k^I`.
Target shears occupy only its all-ones trivial line; every nonzero projection
to the quotient is an obstruction.  In particular, on two components the
packet

```text
(a,-a),                       a!=0,                    (7)
```

is the sign representation and cannot be equalized.  No geometric
involution is needed for this first layer: the packet itself is the retained
sidecar.

## 2. Exact quadratic specialization hypotheses

Now let `Q in A` be smooth and suppose

```text
J(P_0,Q)=c in k*,              ker_K J(-,Q)=k(Q).       (8)
```

Assume its generic function field has a quadratic presentation

```text
K=k(Q)(z,Y),                  Y^2=Delta(z,Q),           (9)
J(Q,z)=M(z,Q)Y,              M in k(Q)(z)*.            (10)
```

Here `[K:k(Q)(z)]=2`; equivalently `Delta` is nonsquare in `k(Q)(z)`.

Let `sigma` fix `k(Q)(z)` and send `Y` to `-Y`.  At a target value `lambda`
call `(9)` **reduced and component-exhaustive** when all of the following
hold.

1. `Delta` has a `q=Q-lambda` specialization
   `Delta_lambda in k(z)*`; in particular it is not zero and the specialized
   quadratic algebra is reduced.
2. The functions `z,Y` are regular at the generic point of every component
   of `Q=lambda`.
3. Specialization induces an isomorphism of total quotient algebras

   ```text
   prod_i k(D_i)
     ~= E_lambda:=k(z)[Y]/(Y^2-Delta_lambda).           (11)
   ```

Condition `(11)` is the load-bearing no-leak hypothesis.  If
`Delta_lambda` is nonsquare, `E_lambda` is one quadratic field.  If
`Delta_lambda=s(z)^2`, it is `k(z) x k(z)` with sign branches
`Y=+s,-s`.  It excludes a collapsed affine component, an omitted branch,
an extra source component, and the nonreduced specialization
`Delta_lambda=0`.

For the global statement below, require reduced component-exhaustiveness at
every target value over which the canonical primitive defined in Section 3
has a pole.  No condition is imposed at irrelevant fibres.

## 3. The canonical primitive is anti-invariant

Let

```text
D=J(-,Q).                                               (12)
```

Equations `(9),(10)` give

```text
D(z)=-MY,                    D(Y)=-(M/2)Delta_z.        (13)
```

Thus

```text
sigma D=-D sigma.                                      (14)
```

Since `D(P_0)=c`, equation `(14)` gives `D(sigma P_0)=-c`.  Hence

```text
P_+=(P_0+sigma P_0)/2 in ker D=k(Q),
P_-=(P_0-sigma P_0)/2,       D(P_-)=c.                 (15)
```

Subtracting the target function `P_+` leaves the unique anti-invariant
representative of the complete rational response torsor.  Every
anti-invariant element of `(9)` has the form

```text
P_-=v(z,Q)Y,                    v in k(Q)(z).           (16)
```

In particular, every rational mate has the same `P_-`; target shears alter
only its invariant part.

Suppose `(16)` has a vertical pole over `lambda`, of order `r>=1`.  In the
`q`-adic field `k(z)((q))`, write

```text
v=q^(-r)v_(-r)(z)+higher terms,      v_(-r)!=0.         (17)
```

Reduced component-exhaustiveness makes `Y` a unit at the vertical generic
points, so the leading packet `(4)` identifies through `(11)` with

```text
v_(-r)(z)Y in E_lambda.                                (18)
```

It is nonzero and anti-invariant.  It cannot be a diagonal scalar.  There
are two equivalent exact proofs.

- If `Delta_lambda` is nonsquare, `(18)` lies in one quadratic field.  Were
  it a scalar `a in k`, applying `sigma` would give both `a` and `-a`, hence
  `a=0`, contradicting `(18)!=0`.
- If `Delta_lambda=s^2`, its two component values are

  ```text
  (v_(-r)s,-v_(-r)s).                                  (19)
  ```

  Equality to one scalar again forces both entries to vanish, contrary to
  `v_(-r)s!=0`.

By Section 1, no target shear can remove this leading pole.  A pole of
`P_-` on a divisor where `Q` is nonconstant is horizontal and cannot be
removed by a target function either.  Since `A` is normal and is the
intersection of its height-one DVRs, the global conclusion is

```text
some rational mate lies in A
iff the canonical anti-invariant mate P_- already lies in A.             (20)
```

The reverse implication is immediate because `D(P_-)=c`.  This proves the
quadratic anti-invariant vertical residue gate.  **QED conditional on audit
promotion.**

## 4. Three exact corollaries and one nonexample boundary

### 4.1 THM-3765

On the nonlinear rational-exact boundary of THM-3765,

```text
Q=h+pT+X(1+gT/2)^2,               g,p!=0,
Y=X-T[p+(g^2/4)XT],                P_-=cY/[g(Q-lambda_0)],
lambda_0=h-2p/g.                                      (21)
```

The exceptional fibre is exactly

```text
Q-lambda_0=(1+gT/2)[2p/g+X(1+gT/2)].                  (22)
```

Both components are rational sign branches of the quadratic model, and
`Y` takes the values `+2p/g,-2p/g`.  The leading packet is the nonzero sign
vector

```text
(+2cp/g^2,-2cp/g^2).                                  (23)
```

Thus `(20)` recovers its polynomial no-go without repeating the degree
descent.

### 4.2 THM-3772

The two mixed variable-flank orientations of THM-3772 merely rescale the
two factors in `(22)`.  Their quadratic model is component-exhaustive and
has the same `Y` values and the same packet `(23)`, with component labels
exchanged in the dual orientation.  This recovers the entire nonlinear mixed
wall.  THM-3772's pure-flank wall is caught by Section 1's unequal packet,
not by the quadratic specialization, because its generic cover has split.

### 4.3 THM-3758

For THM-3758,

```text
Q=XA(z)+X^2B(z),              Y=A+2XB,
A=a_0+a_1z,
B=gamma+beta a_0z+(beta a_1/2)z^2.                   (24)
```

Its affine source does **not** satisfy component-exhaustiveness: on the
generic quadratic model the `Y=+A(z)` branch has `X=0`, but the affine source
relation `z=XT` collapses that branch to `z=0`.  This is exactly why the
hypothesis in `(11)` cannot be omitted.  Nevertheless Section 1 applies
directly.  The primitive

```text
P_0=-c/[2Q] [z+Y/(a_1+2beta Q)]                       (25)
```

has scalar leading coefficients

```text
X=0:             -c a_0/(2a_1),
A+XB=0:          +c a_0/(2a_1).                       (26)
```

They form the nonzero sign packet `(7)`, recovering THM-3758's vertical
no-go while also exhibiting the precise degeneration that blocks a naive
generic-involution proof.

The same warning applies beyond degree two.  Section 1 always says that
target shears occupy the diagonal line, but it does **not** prove that a
cubic, cyclic, resolvent, or quartic primitive lies in a nontrivial
isotypic summand.  Such a claim requires its own trace/character calculation,
extension of the action to every pole valuation, and a component-exhaustive
specialization theorem.  No higher-cover conclusion is asserted here.

## 5. Exact controls and next use

Reproduce with

```bash
python3 -B 04-computation/jc2_quadratic_anti_invariant_vertical_residue_gate_thm3775.py
python3 -B -O 04-computation/jc2_quadratic_anti_invariant_vertical_residue_gate_thm3775.py
```

The assertion-free companion checks `(13),(14)` for generic coefficient
profiles, eight component-permutation quotients, eight split/nonsplit
hostiles, the exact exceptional factorizations and signs in THM-3765 and
both THM-3772 orientations, and THM-3758's separate collapsed-branch sign
packet.  These controls protect signs and boundary typing; Sections 1--4 are
the all-degree proof.

The reusable search rule is now exact: after integrating the generic time
form, compute the leading pole packet in the component quotient `(6)` before
attempting any larger coefficient ansatz.  In quadratic models, normalize
to `(16)` and ask first whether the sheet involution survives at the pole
fibre.  In the quartic `C3` lane, the analogous open task is not to assume a
sign: it is to calculate the trace-zero/V4-resolvent isotypic packet and prove
that the affine source retains it.
