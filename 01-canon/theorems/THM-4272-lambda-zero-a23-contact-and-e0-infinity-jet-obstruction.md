---
id: THM-4272
title: "Lambda-zero A23 contact, E0 infinity jets, and the fat-incidence wall collar"
status: >
  PROVED RELATIVE TO THM-4230/4241/4259/4260/4268 + VERIFIED-EXACT +
  INDEPENDENTLY AUDITED. Over B*=Spec C[U,W,Z,(UZD)^-1], without inverting
  Lambda=U+W+Z, the genus-seven component stays smooth and its contact with
  the rational component is finite flat of rank twelve. On W=0, Lambda=0,
  the contact is the single fat point 12Q, the two-component union has an
  A23 singularity, and every nonconstant map to E0 has ramification index
  in {1,2,4,7}, so it cannot be constant on the fat contact. Together with
  THM-4260/4268, the honest degree-34/42
  contact-collapse incidence has a non-effective Zariski collar across the
  entire UZ!=0, W=0 divisor. Raw Keller-map descent at the A23 contact,
  wall-candidate exclusion, exact-M=12 entry, JC(2), and DC(2) remain open.
source: codex-longer-frontiers-20260827
depends_on:
  - THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze
  - THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing
  - THM-4259-w0-explicit-hlambda-normalization-and-glue-dictionary
  - THM-4260-w0-canonical-node-reciprocal-denominator-attachment-exclusion
  - THM-4268-relative-abelian-map-incidence-properness-and-w0-collar
related:
  - THM-4265-w0-reduced-wall-factor-and-transverse-jacobian-reduction
external: >
  The relative-Hom graph construction and local finite-presentation input
  are Stacks Project Tags 0D1B/0D1C. The Neron extension step uses
  Bosch--Lutkebohmert--Raynaud, Neron Models (1990), Propositions 1.2.8 and
  1.4.4, as already audited in THM-4268.
script: 04-computation/jc23_lambda_zero_a23_fat_incidence_collar_thm4272.py
output: 05-knowledge/results/jc23_lambda_zero_a23_fat_incidence_collar_thm4272.out
script_sha256: a584753945f66aa70c5ff1672c209a6a6408b5f5459d286d92e491dc5e131b14
output_sha256: bf9e216dac041cc51af2dee53e79dbb257d68123308d86e8b4a442f595e5cf54
independent_script: 04-computation/jc23_lambda_zero_a23_fat_incidence_collar_independent_audit_thm4272.py
independent_output: 05-knowledge/results/jc23_lambda_zero_a23_fat_incidence_collar_independent_audit_thm4272.out
independent_script_sha256: 93b675e645c84c7fc88dc814811a0259ca79340a7fe012d74d4aa382ee5ba07f
independent_output_sha256: ad65372698afbfc5a0db9eced050d92138b50cde93365d5f33fa51b91fc86794
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The SymPy path checks the torus syzygy, C-only edge gates,
  finite-flat contact algebra, A23 local branch, hidden-differential
  resultants, and E0 vanishing orders. A standard-library path independently
  reconstructs the sparse syzygy, binomial discriminants, quotient-ring
  norms, local exponent data, and hostile controls. Ordinary, optimized, and
  fixed-hash-seed replays are byte-identical. The relative-Hom/Neron/Rosati
  argument is the load-bearing geometric proof of the incidence collar.
---

# THM-4272 -- `A_23` contact, `E_0` infinity jets, and the fat-incidence wall collar

**PROVED RELATIVE TO THM-4230/4241/4259/4260/4268 + VERIFIED-EXACT +
INDEPENDENTLY AUDITED. `JC(2)` REMAINS OPEN.**

This theorem has two deliberately separated outputs:

- **(A) PROVED:** the honest scheme-theoretic fixed-degree incidence which
  collapses the full fat contact has a non-effective Zariski collar across
  the entire `U*Z!=0`, `W=0` divisor, including `U+Z=0`;
- **(B) OPEN:** no theorem yet makes the Keller-induced map descend from its
  resolved rational chain to a morphism on this raw model. Thus (A) is not a
  Keller wall-candidate, exact-`M=12`, `JC(2)`, or `DC(2)` exclusion.

## 1. Statement

Work over `C`. Retain THM-4230's exact-weight-twelve main face

```text
R: P=S^2,
C: 1-U P^6-W S^2P^5-Z S^4P^4=0.                    (1)
```

For Parts 1--3 assume

```text
W=0,                    U*Z!=0,                    U+Z=0. (2)
```

Choose `alpha^6=U`, `beta^4=Z`, put

```text
epsilon=beta^2/alpha^3,              epsilon^2=-1,    (3)
```

and identify the normalization of `C` with

```text
C_0:x^6+y^4=1,
x=alpha P,                         y=beta SP.          (4)
```

Let `Q_epsilon` be the point at infinity of `C_0` at which
`y^2/x^3=epsilon`.

> **Theorem.**
>
> 1. In the raw toric face, `R` and `C` meet scheme-theoretically in
>
>    ```text
>    12 Q_epsilon.                                      (5)
>    ```
>
>    Their union has a two-branch `A_23` singularity there, with
>    intersection multiplicity and delta invariant both equal to `12`.
>
> 2. The order-twelve source action fixes `Q_epsilon` and acts on a local
>    parameter `b` by `b -> zeta_12^(-1)b`. Consequently
>
>    ```text
>    O_(12Q_epsilon)=C[b]/(b^12)                       (6)
>    ```
>
>    contains every `C_12` character exactly once.
>
> 3. For every nonconstant morphism `m:C_0->E_0`, where
>
>    ```text
>    E_0:Y^2=X^3+1,
>    ```
>
>    one has
>
>    ```text
>    ord_(Q_epsilon)(m^*(dX/Y)) in {0,1,3,6},
>    e_(Q_epsilon)(m)             in {1,2,4,7}.         (7)
>    ```
>
>    In particular, `m|_(12Q_epsilon)` is not constant. More strongly, no
>    nonconstant `C_0->E_0` map is constant on `8Q_epsilon`.
>
> 4. Put
>
>    ```text
>    D=W^2-4UZ,                 Lambda=U+W+Z,
>    B*=Spec C[U,W,Z,(UZD)^(-1)].                    (30)
>    ```
>
>    Do **not** invert `Lambda`. The positive-genus component has a smooth
>    projective genus-seven toric closure `Cbar->B*`, and its intersection
>    with `R` is a finite flat rank-twelve divisor
>
>    ```text
>    A=Spec_(B*) O_(B*)[b]/(b^12-Lambda).             (31)
>    ```
>
>    For `d in {34,42}`, let `I_d^0` parametrize degree-`d` maps
>    `Cbar_b->E_0` whose restriction to `A_b` is the zero map, and let `I_d`
>    parametrize such maps together with an arbitrary common value on `A_b`.
>    Then `I_d^0->B*` is finite,
>
>    ```text
>    I_d ~= (E_0 x B*) x_(B*) I_d^0                  (32)
>    ```
>
>    is proper, and both incidences are empty over the complete divisor
>    `B*_0=V(W)`. Consequently the complement of their finite
>    scheme-theoretic images is a Zariski-open neighbourhood of `B*_0`.
>    This is the non-effective collar in (A).

The choice of roots in `(3)` chooses one of the two infinity points. For a
fixed choice the contact is one length-twelve branch, not
`6Q_+ + 6Q_-`. Changing the radical normalization can exchange the two
descriptions, but does not create a second physical contact.

## 2. Raw contact and analytic type

Near `S=infinity`, use the pre-division toric coordinates

```text
b=1/S,                         r=P/S^2.                (8)
```

Multiplication of the second equation in `(1)` by `b^12` gives the regular
closure

```text
Ctilde=b^12-U r^6-W r^5-Z r^4,       R:r=1.           (9)
```

Thus on the complete coefficient base the contact algebra is

```text
C[b]/(b^12-(U+W+Z)).                                  (10)
```

It is finite free of rank twelve. Under `(2)`, equation `(9)` becomes

```text
Ctilde=b^12-U r^4(r^2-1).                             (11)
```

At `Q=(b,r-1)=(0,0)`,

```text
partial_r Ctilde=-2U!=0.                              (12)
```

Hence `C` is smooth at `Q`; its branch has

```text
r-1=(1/(2U))b^12+O(b^24).                             (13)
```

The two smooth branches `R` and `C` therefore have intersection order
twelve. In completed local coordinates their union is

```text
q(q-phi(b))=0,                   ord_b(phi)=12.        (14)
```

After translating `q` by `phi/2` and absorbing a unit in `b`, this is
`y^2-b^24=0`, the `A_23` singularity. Its delta invariant is `12`.

The one-branch label is also intrinsic. Equations `(3)--(4)` give

```text
y^2/x^3=epsilon*S^2/P=epsilon/r.                      (15)
```

Thus `R:r=1` selects `Q_epsilon`. The other boundary point has `r=-1` and
ratio `-epsilon`, so it is not on `R`. Under the source action
`S->zeta_12 S`, one has `b->zeta_12^(-1)b`, proving `(6)`.

## 3. The complete `E_0` pullback-differential subspace

Let

```text
M=Hom(J(C_0),E_0),                  O=Z[omega].        (16)
```

THM-4241 proves that `[u,f,g,h]` is a full `O`-basis and that

```text
2h=v+omega^2 f+g.                                      (17)
```

Therefore the complex span of all pullback differentials is generated by
`u,v,f,g`. This is the span through the chosen CM embedding of `O`; it must
not be misstated as `M tensor_Z C`.

For the visible maps

```text
u=(-x^2,y^2),                    v=(-x^-2,epsilon_0 y^2/x^3),
```

direct differentiation gives, up to nonzero constants,

```text
u^*(dX/Y)=x y dx/y^3,            v^*(dX/Y)=y dx/y^3.  (18)
```

For THM-4259's hidden basis, equations (16)--(17) there give

```text
eta_f=c x^-5(A_0+B_0y^2)dy,
eta_g=c' x^-5(A_0-B_0y^2)dy.                          (19)
```

Using `6x^5dx+4y^3dy=0`, these span

```text
dx/y^3,                         y^2 dx/y^3.            (20)
```

Both coefficients are nonzero. For

```text
q(a)=a^4-2a^3-2a+1,
A_0=-a^3+2a^2+3,                B_0=a^3-2a^2-1,
```

the exact resultants are

```text
Res(q,A_0)=6,                    Res(q,B_0)=-2.        (21)
```

Consequently the **complete** `E_0` pullback-differential subspace is

```text
span_C {dx/y^3, y dx/y^3, x y dx/y^3, y^2 dx/y^3}.    (22)
```

This uses THM-4241's integral exhaustion. A rational isotypic decomposition
alone would not prove completeness; MISTAKE-521 is controlling.

## 4. Infinity orders and the jet obstruction

At either infinity point of `C_0`, `b` is a uniformizer and

```text
ord_b(x)=-2,                ord_b(y)=-3,
ord_b(dx)=-3.                                         (23)
```

For

```text
eta_(a,c)=x^a y^c dx/y^3,
```

one has

```text
ord_b eta_(a,c)=6-2a-3c.                              (24)
```

The four basis vectors in `(22)` therefore have distinct orders

```text
6,3,1,0.                                               (25)
```

Distinct leading orders cannot cancel. Every nonconstant curve morphism to
`E_0`, after a target translation, induces a nonzero element of `M`, so its
pullback of the invariant differential has one of the four orders in `(25)`.
In characteristic zero the local ramification index is one more than that
order, proving `(7)`.

If `m|_(12Q_epsilon)` were constant and `z` were a target uniformizer at
`m(Q_epsilon)`, then

```text
m^*z in (b^12),                   ord_b(m^*(dX/Y))>=11, (26)
```

contradicting `(25)`. The same argument already contradicts constancy on
`8Q_epsilon`.

## 5. The `B*` family and the honest fat-incidence collar

The strengthening in Part 4 does not require `Lambda` to be a unit. On the
torus, write

```text
sbr=WP+2ZS^2,
pbr=6UP^2+5WS^2P+4ZS^4.                              (33)
```

The two face derivatives and their exact syzygy are

```text
F_S=-2SP^4 sbr,                 F_P=-P^3 pbr,
2Z pbr-(4ZS^2+3WP)sbr=-3DP^2.                        (34)
```

Thus `D!=0` rules out a torus critical point. For `C` alone, the three
nontrivial boundary edge polynomials have discriminants

```text
disc(1-ZX^4)=-256Z^3,
disc(UX^2+WX+Z)=D,
disc(U-X^6)=46656U^5.                                (35)
```

Hence their gates are exactly `Z,D,U`. The factor `X-1` belongs to the
rational component `R`, and

```text
Res_X(X-1,UX^2+WX+Z)=U+W+Z=Lambda.                   (36)
```

Consequently `Lambda` measures transversality of `R` with `C` at that
boundary edge; it is not a smoothness gate for `C` itself. Closing `C` in a
fixed smooth toric subdivision of its Newton fan, equations `(34)--(35)` and
the relative Jacobian criterion give a smooth projective family
`Cbar->B*`. Pick's formula gives genus seven on every fibre.

The toric chart `(8)--(10)` is valid without dividing by `Lambda`. It gives

```text
A=Cbar intersect R
 =Spec_(B*) O_(B*)[b]/(b^12-Lambda).                 (37)
```

Because the equation is monic in `b`, `A->B*` is finite faithfully flat of
rank twelve. It is étale where `Lambda!=0`; at `Lambda=0` it is the honest
nonreduced length-twelve contact. Reducedness and twelve separate sections
are not used below.

We now rerun, over `B*`, the relative-Hom argument audited in THM-4268.
The fixed-Hilbert-polynomial graph space makes the degree-`d` map functor
separated and of finite type. Requiring

```text
phi|A=0,       or more generally       phi|A=a        (38)
```

is a closed equalizer condition because `A` is finite flat and the abelian
target is separated. Define the unbased incidence using pairs `(phi,a)`;
faithful flatness of `A` makes `a` unique, and translation by `-a` proves
`(32)`.

For properness, let a DVR `R_DVR` map to `B*`. The source
`Cbar_(R_DVR)` is smooth and projective, and `E_0 x Spec(R_DVR)` is an
abelian scheme, hence the Néron model of its generic fibre. The Néron mapping
property extends a generic map uniquely over `R_DVR`. The pullback of a
fixed ample target bundle has locally constant fibre degree, so degree `d`
is preserved. Finally, the generic and zero restrictions on
`A_(R_DVR)` agree everywhere: the generic fibre is schematically dense by
flatness, and the target is separated. This proves the valuative criterion
for `I_d^0`, hence for `I_d` by `(32)`.

Every geometric fibre of `I_d^0` is finite. Choose a support point of `A_b`
as Abel--Jacobi basepoint. A based map determines an element of

```text
Hom(J(Cbar_b),E_0),                                   (39)
```

a finitely generated free lattice, and its degree is a positive-definite
Rosati quadratic form. Only finitely many lattice vectors have fixed value
`d`. Thus `I_d^0` is quasi-finite; proper plus quasi-finite makes it finite.

It remains to empty the complete `W=0` fibre. There `D=-4UZ` is already a
unit and `Lambda=U+Z`.

- If `Lambda!=0`, THM-4260 excludes every degree-`34` and degree-`42` map
  which collapses the twelve reduced contacts.
- If `Lambda=0`, Parts 1--3 identify `A_b` with one `12Q_epsilon`. A
  degree-`34` or degree-`42` map is nonconstant, while Part 3 says every
  nonconstant `C_0->E_0` map is already nonconstant on `8Q_epsilon`.
  It therefore cannot be constant on `A_b`.

The finite-type fibre has no geometric points in either case and is thus
empty as a scheme. If `Z_d` is the finite scheme-theoretic image of
`I_d^0->B*`, then

```text
Omega*=B* \ (Z_34 union Z_42)                         (40)
```

is a Zariski-open neighbourhood of the **entire** divisor `V(W)<=B*`,
crossing `Lambda=0`. No effective equation or analytic radius for this
collar is asserted.

## 6. Sharp hostile and application firewall

The `E_0`-isotypic restriction is load-bearing. Put

```text
s=y^2/x^3.
```

At `Q_epsilon`,

```text
s^2+1=x^-6,                     ord_b(s-epsilon)=12.  (27)
```

Therefore the canonical differential

```text
x^3 dx/y^3-epsilon^-1 y^2 dx/y^3                     (28)
```

has exact order twelve at `Q_epsilon`. It mixes a differential outside the
`E_0` pullback subspace with one inside it. Thus an argument using only genus,
the full canonical system, or the Riemann--Hurwitz total would be false.

Parts 1--4 prove the following unconditional statement about the honest raw
incidence:

```text
no nonconstant C_0->E_0 morphism can agree scheme-theoretically
with a constant R-map on the raw contact 12Q_epsilon.                 (29)
```

Section 5 promotes this to a proper/finite collar for the honest
scheme-theoretic contact-collapse incidence. It does **not** prove that the
Keller-induced rational face map is an object of that incidence at
`Lambda=0`. The present planar-Jacobian reduction extends that map only after
resolving indeterminacy. Resolution of the `A_23` contact produces rational
chain data and guarantees pointwise agreement at the resolved nodes; it need
not descend to equality on the raw nonreduced intersection. Cleared
equations or a resolved graph are not automatically the saturated raw graph
(MISTAKE-455).

A future wall exclusion needs one of:

1. regularity of the Keller face map on the raw toric model;
2. saturation showing that its resolved graph descends through `(14)`; or
3. a contact-order preservation theorem retaining the length-twelve jet.

No such result is asserted here. The exact status split is:

```text
PROVED (A): the honest degree-34/42 fat-contact incidence has a Zariski
            collar across the full UZ!=0, W=0 divisor;
OPEN   (B): the resolved Keller response descends into that raw incidence.
```

In particular this theorem proves no:

- isolated `U+Z=0` Keller wall-candidate exclusion;
- statement on `U=0` or `Z=0`, where the `C_0` scaling degenerates;
- off-`W=0` result or transverse denominator lift;
- exact-`M=12` entry statement; or
- `JC(2)` or `DC(2)` conclusion.

## 7. Dependencies, corrections, and reproduction

Load-bearing dependencies:

- `THM-4230-exact-weight-twelve-cyclotomic-prym-planar-jacobian-squeeze`:
  face equations, `C_0` normalization, infinity geometry;
- `THM-4241-w0-hidden-cm-saturation-and-visible-hidden-index-four-gluing`:
  complete integral `Hom(J(C_0),E_0)` basis;
- `THM-4259-w0-explicit-hlambda-normalization-and-glue-dictionary`:
  exact hidden pullback differentials;
- `THM-4260-w0-canonical-node-reciprocal-denominator-attachment-exclusion`:
  complete degree-`34/42` contact-collapse emptiness on the
  `W=0`, `Lambda!=0` gate.
- `THM-4268-relative-abelian-map-incidence-properness-and-w0-collar`:
  the graph-Hilbert/Néron/Rosati properness method. Its published base
  inverted `Lambda`; equations `(33)--(37)` are the new audit which permits
  the same argument over `B*`.

The external geometric inputs are the relative-Hom graph construction
(Stacks Project Tags `0D1B/0D1C`) and the Néron mapping property for abelian
schemes over DVRs (Bosch--Lütkebohmert--Raynaud, *Néron Models*, Props.
`1.2.8` and `1.4.4`), exactly as recorded and audited in THM-4268.

Controlling corrections/firewalls: MISTAKE-521, MISTAKE-455, MISTAKE-527,
and MISTAKE-509.

Reproduction from the repository root:

```bash
python3 -B \
  04-computation/jc23_lambda_zero_a23_fat_incidence_collar_thm4272.py
python3 -B -O \
  04-computation/jc23_lambda_zero_a23_fat_incidence_collar_thm4272.py
python3 -B \
  04-computation/jc23_lambda_zero_a23_fat_incidence_collar_independent_audit_thm4272.py
python3 -B -O \
  04-computation/jc23_lambda_zero_a23_fat_incidence_collar_independent_audit_thm4272.py
```

The matching outputs and hashes are recorded in the frontmatter. The theorem
proof itself has no known internal gap relative to the listed dependencies
and standard geometric inputs. The missing raw-model descent is an explicit
**application gap**, not a suppressed step in the honest-incidence theorem.
