---
id: THM-4046
title: "Exceptional-quartic actual J7 lift, sharp J8 obstruction, and critical-displacement closure"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  Over the
  exceptional irreducible quartic field, the complete retained order-eight
  arbitrary-two-form matrix has rank 40 in 45 retained coordinates.  The
  retained relation space (the dual cokernel) has dimension five; its
  J_8-inactive subspace has dimension four, so the quotient is
  one-dimensional.  Its normalized representative extends the THM-3683
  identity and has a nonzero constant response kappa.  The promoted
  THM-4043 actual stagewise lift therefore continues through J_7 on the
  quadratic exceptional fold,
  but no choice of actual F_9,G_9 extends it through J_8.  Character dilation
  handles monomial H, and target-compatible formal monomialization reduces
  general H to that case; together they show that, for each exceptional
  embedding and every 0!=H in t^2 C[t], no arbitrary target two-form has
  nonzero constant pullback in this Russell family.  No coherent infinite
  series, degree control, other fold or compiler, global Keller pair, or
  JC(2) conclusion is asserted.
source: jc2-double-zero-rebuild-20260824 / retained order-eight closure, 2026-08-24
audit: >
  PASS -- an independent SymPy algebraic-field/local-Taylor reconstruction,
  importing no production code, recovered the complete 45-by-495 universe,
  rank 40, nullity five, one-dimensional J_8 projection, 35-entry active
  relation, Lambda block, all four positive value-block cancellations, the
  same kappa coordinates and nonzero norm, and the identical canonical row
  hash.  It also embedded the full four-dimensional order-six cokernel.
  At the named exact and modular hostile controls, moving off the exceptional
  quartic or off the zero-fourth parabola removes the J_8-active class.
  Normal and optimized executions
  byte-match the frozen transcripts, and both exact implementations contain
  no Python Assert node.  A separate audit checked the w-shift, sign, stage indices,
  formal-conjugacy typing, and strict scope.
depends_on:
  - THM-3673-russell-cylinder-monomial-ramification-debt-dilation
  - THM-3675-russell-cylinder-critical-fold-formal-conjugacy-closure
  - THM-3683-russell-cylinder-sixth-debt-quartic-on-the-zero-fourth-parabola
  - THM-3737-russell-cylinder-exceptional-quartic-jacobian-image-hyperplane
  - THM-4043-exceptional-quartic-shifted-stable-identities-and-j6-lift
related:
  - THM-3677-russell-cylinder-degree-eight-fourth-debt-parabola
  - THM-3688-russell-cylinder-exceptional-quartic-actual-j1-j2-lift
  - THM-4039-exceptional-quartic-j3-lift-and-adjacent-gate-rigidity
  - THM-4054-exceptional-affine-simple-zero-retained-packet
script: 04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046.py
output: 05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046.out
script_sha256: c5032d3c066a31207a4f09583be012c2437e6650e23c2ccae6db542592dd060c
output_sha256: ae333c2baf1ab222d6b7f40fb400bac5159fbda6c9aca8a0b9a0d391ff10039b
hash_basis: raw LF bytes
---

# THM-4046 -- the exceptional lift reaches `J_7` and stops at `J_8`

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  A complete
order-eight retained calculation supplies the first unavoidable scalar debt
after the exceptional quartic's long sequence of cancellations.  The same
identity also closes every nonzero critical vertical displacement over the
four exceptional folds against nonzero constant pullback.

All rings in the lift statement have characteristic zero.

## 1. The complete order-eight retained universe

Let

```text
K=Q[alpha]/(72783360 alpha^4-77822208 alpha^3-28419741 alpha^2
                                      +7849770 alpha-1276420),       (1)
```

Using THM-3683's regular local target coordinates

```text
y=c/3,                       z=e+3,                       w
```

at the common retained target point, let `Phi_alpha` be the corresponding
quadratic exceptional fold on the complete zero-fourth parabola, with stable
target coordinate pulling back as `w=t`.  For arbitrary coefficient germs
`A,B,C`, write the arbitrary target two-form

```text
Omega=A dy^dz+B dy^dw+C dz^dw,
Phi_alpha^*(Omega)=sum_(n>=0)t^n J_n(x) dx^dt.             (2)
```

No closedness or decomposability condition is imposed on `Omega`.

Retain the Taylor coefficients at `x=-1,0,1` in the window

```text
J_0 through x-order 4:       15 coordinates,
J_2 through x-order 3:       12 coordinates,
J_4 through x-order 2:        9 coordinates,
J_6 through x-order 1:        6 coordinates,
J_8 values:                   3 coordinates.              (3)
```

The three wedge slots in `(2)` and all target monomials

```text
y^a z^b w^c,                    a+b+c<=8,                 (4)
```

give exactly

```text
3 binom(11,3)=495                                             (5)
```

coefficient-jet generators.  Coordinates are expanded through total order
nine before differentiation and only then truncated at density order eight,
so the top face is present.  Regularity and vanishing at the retained point
make all higher coefficient monomials invisible in this window.  Thus
`(3)--(5)` are a complete 45-by-495 retained matrix, not a sample or ansatz.

Exact elimination over `K` gives

```text
rank=40,            retained-relation dimension=5.         (6)
```

The lower 42-coordinate matrix has rank `38`, so its retained-relation space
has dimension four.  Extending those relations by zero on the `J_8` block
embeds them in the five-dimensional full retained-relation space;
restriction of the latter to the three `J_8` coordinates has rank one.  In
the canonical RREF gauge whose free coordinates are the `x=1` value rows of
`J_0,J_2,J_4,J_6,J_8`, set the first four free coordinates to zero and
normalize the last block to

```text
Lambda(P)=5P(-1)/18-P(0)+13P(1)/18.                       (7)
```

The resulting representative has 35 nonzero entries and annihilates all
495 columns exactly.  Its canonical serialization has SHA-256

```text
4efb72eafbb35333908ecdadf9ecbb1a861b14bb3615442343516c65b0c9f59b. (8)
```

Its `J_6,J_4,J_2` blocks are literally the `R_4,R_2,R_0` blocks of
THM-3683, reduced in `K`.  Calling the new `J_0` block `S_0`, the identity is

```text
0=Lambda(J_8)+R_4(J_6)+R_2(J_4)+R_0(J_2)+S_0(J_0).       (9)
```

The 12 nonzero entries of `S_0` in the 15-coordinate `J_0` block and all
predecessor coordinates are printed in the frozen transcript and pinned by
`(8)`.

## 2. The constant response is nonzero at all four embeddings

Define

```text
kappa:=S_0(1)
 =-5183766767360/3^19
   +(931699873280/3^16) alpha
   -(460399208960/3^14) alpha^2
   -(183606968320/3^13) alpha^3.                         (10)
```

Multiplication by `kappa` on the power basis of `K` has determinant

```text
N_(K/Q)(kappa)
 =8836893677589250971329876479284734525440
  /721180557365797959546519 !=0.                         (11)
```

The companion also verifies an exact inverse.  Hence `kappa` is nonzero in
`K` and remains nonzero under every one of its four complex embeddings.
The positive stable blocks have zero value sums, while their derivative
terms kill constants, so

```text
R_0(1)=R_2(1)=R_4(1)=Lambda(1)=0.                        (12)
```

## 3. The quadratic exceptional lift reaches `J_7`

THM-4043 supplies actual stagewise target coefficients over `K` with

```text
J_0=1,                    J_1=J_2=J_3=J_4=J_5=J_6=0.    (13)
```

Set the next coefficients `F_8=G_8=0` provisionally and denote the resulting
seventh coefficient by `D_7`.  Let `Omega=dF^# wedge dG^#` be the two-form of
this finite provisional actual pair.  Apply `(9)` to `w Omega`.  Since
`Phi_alpha^*(w)=t`, its coefficient sequence is shifted by one, and negative
indices vanish.  Therefore

```text
0=Lambda(D_7)+R_4(J_5)+R_2(J_3)+R_0(J_1).               (14)
```

Equation `(13)` gives `Lambda(D_7)=0`.  THM-3737 proves

```text
L_0(S^2)=ker Lambda,                                    (15)
```

where `S` is the actual restriction algebra.  The new order-eight target
coefficients enter the seventh equation as

```text
J_7=D_7+8L_0(F_8,G_8).                                  (16)
```

Equations `(14)--(16)` therefore supply actual restrictions, and hence
actual target representatives, `F_8,G_8` with `J_7=0`.  Combining with
`(13)` proves an actual stagewise lift

```text
boxed: J_0=1,                 J_1=J_2=...=J_7=0.         (17)
```

This construction is over the common field `K`, so `(17)` holds uniformly
after all four embeddings.

## 4. The `J_8` obstruction is sharp

After choosing `(16)`, set `F_9=G_9=0` provisionally and call the resulting
eighth coefficient `D_8`.  Apply the unshifted identity `(9)`.  Its positive
even inputs vanish by `(13)`, while `J_0=1`; consequently

```text
Lambda(D_8)=-kappa !=0.                                 (18)
```

The sign in `(18)` is forced by moving `S_0(1)=kappa` across `(9)`.  Any
order-nine correction has the form

```text
J_8=D_8+9L_0(F_9,G_9).                                  (19)
```

But `(15)` gives `Lambda(L_0(F_9,G_9))=0` for every actual choice.  Thus
`(18)--(19)` prove that no actual `F_9,G_9` can extend the quadratic
exceptional lift through `J_8`.  This obstruction is independent of the
choices made within earlier solution fibres: it uses only the identities in
`(13)` and the nonzero field element `(10)`.

## 5. Every nonzero critical displacement is closed in this family

This section is a separate obstruction statement, not an all-`H` lift.
First take

```text
H(t)=eta t^k,                    eta!=0, k>=2.            (20)
```

The character-decimation map of THM-3673 transports `(9)` to

```text
0=Lambda(J_(4k))+eta R_4(J_(3k))+eta^2 R_2(J_(2k))
                   +eta^3 R_0(J_k)+eta^4 S_0(J_0).      (21)
```

If the pullback density of an arbitrary target two-form were a nonzero
constant `c`, every positive coefficient in `(21)` would vanish and the
remaining equation would be

```text
0=eta^4 c kappa,                                         (22)
```

contradicting `(10)--(11)`.

Now let `0!=H in t^2 C[t]` be arbitrary, put `k=ord_0(H)`, and let `eta` be
its leading coefficient.  THM-3675's target-compatible formal
monomialization starts with

```text
phi(t)=t(H(t)/(eta t^k))^(1/k),             W=phi(w),
```

and an inverse coordinate `t=psi(u)` with

```text
psi(0)=0,              psi'(0)=1,
H(psi(u))=eta u^k.                                      (23)
```

Along the source `W=u`, and transporting a hypothetical target two-form
through this completed-target change gives another arbitrary formal target
two-form.  Identity `(9)` depends only on its finite order-eight jet, so it
applies to the transported form.

A constant density `c dx^dt` becomes

```text
c psi'(u) dx^du.                                        (24)
```

Every coefficient of `(24)` is constant in `x`.  Derivatives therefore
vanish, and `(12)` kills every positive stable value block.  The constant
term of `psi'` is one, so `(21)` again reduces exactly to `(22)`.  Hence,
for each of the four exceptional embeddings and each

```text
0!=H in t^2 C[t],                                       (25)
```

no arbitrary target two-form has nonzero constant pullback on this Russell
compiler.  A fortiori, no actual target pair has nonzero constant source
Jacobian on any fold in `(25)`.

## 6. Exact and hostile verification

The production companion performs exact power-basis elimination over `K`.
It verifies a nonzero rank-40 minor, every one of the 495 relation residuals,
the four-dimensional `J_8`-inactive relation subspace, the literal
predecessor blocks, `(10)--(12)`, and an inverse for `kappa`.  Over the split
prime `137`,
the four exceptional roots are

```text
44,82,92,134,
```

all four matrices have rank `40`, and the corresponding constant responses
are respectively

```text
105,71,89,8.                                             (26)
```

The exact hostile controls are

```text
zero-fourth parabola, r=0 off the quartic:
  rank 41, nullity 4, lower rank 38, J_8 projection 0;

Q_6 off the parabola:
  rank 42, nullity 3, lower rank 39, J_8 projection 0.   (27)
```

At the exceptional parameter with the parabola coefficient shifted by one,
a separate mod-137 hostile also has rank `42`, lower rank `39`, and no
`J_8` projection.  The independently constructed SymPy companion uses local
sparse Taylor arithmetic and differentiation before truncation; it imports
no production code, explicitly reconstructs and embeds the full order-six
cokernel, and reproduces the canonical hash `(8)` and norm `(11)`.  Its
frozen files are

```text
04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046_independent_audit.py
  sha256 e59e60fc6f431d61fc79321592a425ffee7c1a719f31f37cf682dbeb9d799a09
05-knowledge/results/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046_independent_audit.out
  sha256 aa2b4f1a182036b8934a3f5436dce8ab7f249bccb7b4bdb95acb9bbb0d550ef6
```

Normal and optimized executions of both exact companions are byte-identical
to their stored outputs.  Neither script contains a Python `Assert` node.

Reproduce with

```bash
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046.py
python3 -B 04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046_independent_audit.py
python3 -O -B 04-computation/jc2_russell_cylinder_exceptional_quartic_order8_closure_thm4046_independent_audit.py
```

## 7. Strict boundary

The lift statement proves only finite, stagewise actual coefficients through
`J_7` for the quadratic exceptional fold `H=t^2`; it proves their extension
through `J_8` impossible.  The closure statement proves only the displayed
nonzero-constant-pullback obstruction for the four conjugate exceptional
collision polynomials obtained from the embeddings of `Q_alpha` and
`0!=H in t^2 C[t]`.  It does not cover

- `H=0` or the affine boundary `H'(0)!=0`;
- another collision polynomial, source chart, or compiler;
- an arbitrary planar polynomial map;
- a coherent infinite lift, convergence, algebraization, or degree bounds;
- global polynomial pairs or Keller maps outside this restricted Russell
  family; or
- a counterexample to `JC(2)`, which remains **OPEN**.

This closes the exceptional critical-displacement lane inside the stated
Russell compiler and no larger universe.  **QED.**
