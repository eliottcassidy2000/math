---
id: THM-2615
title: "Physical-diagonal toric kernel and bare-endpoint dipole Radon invoice"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  Restriction from the independent endpoint torus T^E to the physical scalar
  curve D_w(x)=(w_e x) has group-algebra kernel generated exactly by the
  relation binomials chi_a-1, a in Lambda(w).  THM-2309's two explicit target
  relations show that this kernel contains every class of the lawful
  F13^2 target quotient.  Strictly positive product-torus functions can
  therefore agree on the entire physical curve, every inverse root, and
  every temporal itinerary while having different nonzero target-residue
  support.  The smallest natural repair retains independent present and
  conjugated bare-endpoint dipole translations in a 13 by 13 square M(r,s).
  Its lawful diagonal
  H(t)=M(t,t) has Fourier transform equal to the exact Radon line sum
  Hhat(q)=sum_(lambda-nu=q) Mhat(lambda,nu).  Two explicit nonnegative
  permutation squares have identical row/column marginals and base axes but
  diagonals with respectively zero and all twelve primitive colours, proving
  that one-sided data cannot infer this target.  The theorem types the missing
  carrier and decisive test for the present/bare-endpoint pair.  It does not
  identify the bare endpoint with a chronological future paired shift,
  construct a positive chronological lift of the square, prove a nonzero
  diagonal current, exclude a scalar row, or prove LRC(14).
source: root-long-frontiers-2026-07-28-diagonal-kernel
depends_on:
  - THM-2309-owner-aligned-pivot-packets-and-visible-height-separation
  - THM-2334-relation-residue-current-and-character-twist-pushforward
  - THM-2350-owner-pivot-dual-dipole-normal-form
related:
  - THM-2502-endpoint-boolean-newton-carry-tournament-and-dipole-boundary
  - THM-2605-inverse-root-dipole-connection-and-mixed-square-invoice
  - THM-2610-chronological-paired-slice-marked-triangle-graft-and-action-axis-boundary
  - THM-2613-canonical-root-diagonal-opposite-shift-section
script: 04-computation/lrc14_physical_diagonal_toric_kernel_thm2615.py
output: 05-knowledge/results/lrc14_physical_diagonal_toric_kernel_thm2615.out
script_sha256: ab63609e100895bcfc32de9f9d07fbda2183a01ca618c765a49bfae49e21a9d5
output_sha256: d1da2843febc8f51a059a27c5c14d1d95ee69be9a826d22e11c4cf4da7dfe893
hash_basis: LF-normalized bytes
---

# THM-2615 -- the old target residue lives off the physical diagonal

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2610 retains an old marked coefficient beside every future paired-shift
colour, and THM-2613 canonically identifies the physical future root with its
local paired-shift state.  The remaining temptation is to identify that future
shift with the conjugated old bare endpoint and thereby read the THM-2334
relation residue from this increasingly complete physical data.

That is impossible in principle.  The relation lattice is exactly the kernel
of restriction from the independent coordinate torus to the one-dimensional
physical orbit.  The kernel meets every target-residue class.  Recovering the
old target therefore requires an off-diagonal lift with two independent
coordinate actions, not another refinement of the physical root itinerary.

For the actual THM-2334 present/bare-endpoint pair, the exact finite consumer
is simple: form a `13 by 13` independent response square and sum its Fourier
transform along the lines `lambda-nu=q`.  This theorem proves both the no-go
and that lawful Radon invoice.  It does **not** prove that a chronological
future paired shift is the square's second axis.

## 1. The exact physical-diagonal kernel

Let `E` be a finite coordinate set, let

```text
w=(w_e)_(e in E) in Z^E,

Lambda(w)={a in Z^E:a.w=0},

g=gcd(|w_e|:e in E).                                     (1)
```

On the product torus `X=T^E`, write

```text
chi_u(Z)=exp(2 pi i u.Z),                  u in Z^E,       (2)
```

and define the physical scalar curve

```text
D_w:T->X,                  D_w(x)=(w_e x)_(e in E).       (3)
```

Restriction of Laurent characters is the group-algebra map

```text
R_w:C[Z^E]->C[gZ],

R_w(chi_u)=t^(u.w).                                      (4)
```

Its kernel is exactly

```text
ker R_w=<chi_a-1:a in Lambda(w)>.                         (5)
```

Indeed, the homomorphism `Z^E->Z`, `u->u.w`, has kernel `Lambda(w)` and
image `gZ`.  Hence

```text
Z^E/Lambda(w) isomorphic to gZ.                           (6)
```

Passing to complex group algebras, the quotient by the ideal in (5) is
`C[gZ]`, and the induced map is the identity on that quotient.  This proves
(5), including the reverse inclusion.

Pointwise, every `a in Lambda(w)` satisfies

```text
chi_a(D_w(x))=exp(2 pi i x(a.w))=1.                       (7)
```

Replacing `x` by `(z+h)/13`, by `13^N x`, or by any point in a digit
itinerary never leaves `D_w(T)`.  Thus every physical inverse root and every
positive temporal refinement still kills the whole ideal (5).

## 2. Every target-residue class occurs in the kernel

Use the live THM-2309 scalar types.  Let `u_0` be the omitted unit and let
`a,b` be the two target blockers.  The exact relations

```text
t_a=w_a e_(u_0)-w_(u_0)e_a,

t_b=w_b e_(u_0)-w_(u_0)e_b                              (8)
```

belong to `Lambda(w)`.  In THM-2350's canonical target quotient basis,

```text
pi(r)=(r_a-r_(k_a), r_b-r_(k_b)) mod 13.                 (9)
```

Therefore

```text
pi(t_a)=(-w_(u_0),0),             pi(t_b)=(0,-w_(u_0)).  (10)
```

Since `w_(u_0)` is a thirteen-unit, (10) proves

```text
pi(Lambda(w))=F_13^2.                                    (11)
```

The loss survives positivity.  Choose `a in Lambda(w)` with `pi(a)!=0` and
any rational `epsilon>0`.  On `X`, put

```text
F_0(Z)=1,

F_1(Z)=1+epsilon(2-chi_a(Z)-chi_(-a)(Z))
      =1+4epsilon sin^2(pi a.Z).                          (12)
```

Both functions are real and strictly positive, with `F_1>=1`.  Equation (7)
gives

```text
F_0(D_w(x))=F_1(D_w(x))=1                 for every x.    (13)
```

Yet their Laurent supports project respectively to

```text
{0},                         {0,pi(a),-pi(a)}.             (14)
```

Thus even strictly positive lifts with identical complete physical data can
carry different nonzero target residues.  This is a same-data hostile proving
non-identifiability, not merely the absence of a proposed formula.  It does
not assert that either hostile is a canonical THM-2334 current; it proves that
no reconstruction using only physical restriction can determine such a
current's residue decomposition.

## 3. The independent dipole square

Fix one lawful target dipole

```text
eta=e_a-e_(k_a),                 delta(u)=u_a-u_(k_a),

zeta=exp(2 pi i/13).                                      (15)
```

Let an absolutely convergent joint lifted endpoint series on `X x X` be

```text
mathcal F(Z,Z')=sum_(u,v)c(u,v)chi_u(Z)chi_v(Z').         (16)
```

Base-point phases may be absorbed into `c(u,v)`.  Define the independent
present/bare-endpoint response

```text
M(r,s)=mathcal F(r eta/13,-s eta/13)

      =sum_(u,v)c(u,v)zeta^(r delta(u)-s delta(v)),
                                      r,s in F_13.        (17)
```

The minus sign types the conjugated bare/right endpoint in the THM-2334
left-minus-right current.  It is not a declaration that `s` is the
chronological future local action.  Use normalized transform

```text
Mhat(lambda,nu)
 =1/169 sum_(r,s)M(r,s)zeta^(-r lambda+s nu).             (18)
```

Character orthogonality gives

```text
Mhat(lambda,nu)
 =sum_(delta(u)=lambda, delta(v)=nu)c(u,v).                (19)
```

Thus the two axes retain the left and right dipole residues separately.

## 4. The lawful diagonal is an exact Radon line sum

The lawful THM-2334 left-minus-right target action uses the same translation
at the present endpoint and its conjugated bare/right endpoint.  Its response
is the diagonal

```text
H(t)=M(t,t).                                               (20)
```

Taking its normalized one-dimensional transform and using (18)--(19) gives

```text
Hhat(q)
 =1/13 sum_t H(t)zeta^(-tq)

 =sum_(lambda-nu=q) Mhat(lambda,nu).                       (21)
```

This is the exact finite Radon transform of the mixed square.  If the deepest
axis contributes charge `m`, replace `H(t)` by `zeta^(mt)H(t)`.  Then

```text
Hhat_new(q)=Hhat_old(q-m)
           =sum_(lambda-nu=q-m)Mhat(lambda,nu),           (22)
```

or equivalently the selected charge equation is `lambda-nu+m=q`.

Finite Parseval supplies the cheapest decisive test once a lawful square has
been populated:

```text
sum_(q!=0)|Hhat(q)|^2
 =1/13 sum_t |H(t)-mean(H)|^2.                            (23)
```

Consequently the thirteen diagonal evaluations decide nonzero target residue
exactly.  Separate nonzero spectra on the two axes do not imply that their
line sums are nonzero.

## 5. Marginals and both base axes still do not determine the diagonal

The loss is sharp even for nonnegative `0/1` squares.  On `F_13^2`, define

```text
P(r,s)=1_(s=r+1),

Q(r,s)=1_(s=sigma(r)),

sigma=(0 1 12),                   sigma fixes 2,...,11.   (24)
```

Both are permutation matrices, so all row and column marginals equal one.
They also agree on the complete present base row `r=0` and complete bare/right
base column `s=0`:

```text
P(0,s)=Q(0,s)=1_(s=1),

P(r,0)=Q(r,0)=1_(r=12).                                  (25)
```

Their lawful diagonals are nevertheless

```text
P(t,t)=0,

Q(t,t)=1_(t in {2,...,11}).                              (26)
```

The latter is a nonempty proper consecutive block, so each of its twelve
primitive Fourier coefficients is nonzero.  Thus row marginals, column
marginals, a complete present-only profile, and a complete bare-endpoint-only
profile can all agree while the target residue changes from zero to full
primitive spectrum.

The `13 by 13` square is therefore the smallest **natural product carrier**
which retains two independent `F_13` translation coordinates and on which the
diagonal (20) is defined.  This is not an information-theoretic lower bound
against arbitrary compressed encodings.

## 6. Exact consequence and frontier

The theorem separates three objects which had repeatedly shared the label
`F_13`:

```text
physical root/local shift:        supplied canonically by THM-2613;

old relation-coordinate action:  killed by physical restriction (5);

left-minus-right target residue:  the diagonal Radon sum (21).    (27)
```

THM-2610 proves coexistence of an old marked coefficient and all future local
shift colours, but it neither populates `M(r,s)` under the independent present
and bare-endpoint actions nor identifies its future action with the bare/right
axis.  THM-2613 fixes the physical root-to-future-shift section, but
(12)--(14) prove that physical restriction cannot recover the old relation
action.

The next constructive target is therefore exact: retain the independently
translated present and conjugated bare endpoints of the THM-2334 current on a
positive chronological carrier, keep their full square before diagonal
restriction, and prove that (21) is nonzero.  If a chronological future paired
shift is to replace the bare/right axis, a separate lawful intertwiner must be
constructed and audited.  A root table, one-axis row, or separate pair of
primitive spectra cannot substitute for either construction.

No physical mixed square or nonzero Radon sum is proved here.  No scalar row
is excluded; the ledger remains `165` and LRC(14) remains open.

## 7. Exact companion

Run

```text
python3 04-computation/lrc14_physical_diagonal_toric_kernel_thm2615.py
python3 -O 04-computation/lrc14_physical_diagonal_toric_kernel_thm2615.py
```

Both modes byte-match the stored transcript.  Over `F_53`, with primitive
thirteenth root `16`, the dependency-free referee checks:

- all `2,197` delta-basis instances of the Radon identity (21);
- all `28,561` deepest-charge dephasings;
- every row/column/base-axis equality in (25) and all twelve primitive
  diagonal colours in (26);
- all `169` target residues generated by the explicit relations (8); and
- collapse of the physical exponents of the positive hostile.

The exact group-algebra kernel (5) and real positivity in (12) are analytic
proofs above, not extrapolations from the finite companion.  An independent
hostile audit rederived the kernel signs and quotient surjectivity, checked the
DFT normalization, Radon sign, Parseval identity, deep-charge shift, and `P/Q`
hostile, and caught the original draft's illicit naming of the second axis as
chronological future.  The repaired theorem types it as the conjugated old
bare endpoint throughout.  QED.
