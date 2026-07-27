---
id: THM-2557
title: "Doubly centered seven-step transfer, root-moment index, and the Hall-cone boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED for the abstract
  lattice, moment, cone, and Hall controls; CONDITIONAL APPLICATION inheriting
  THM-2550(B)'s current PROVED-CANDIDATE exact-input status.  On the integral doubly centered 7-by-13
  interaction lattice, every nonzero-slope seven-step transfer is injective
  with exact image the kernel of a six-dimensional mod-13 first-root-moment
  map; its Smith form is 1^66 direct-sum 13^6.  Consequently the nonzero
  THM-2550(B) typed-row interactions survive every such horizontal transfer,
  but remain signed zero-margin currents and cannot themselves be
  nonnegative THM-2545 ancestry couplings.  No actual word-resolved
  head/later-root graph is constructed, so no Hall deficit, scalar-row
  exclusion, or LRC(14) contradiction is claimed.  The typed row is not a
  scalar cover.  THM-2550(A)'s different k=2 drift packet is not composed
  with its large-clock Part (B) table here.
source: codex-2026-07-27-typed-row-transport
depends_on:
  - THM-2512-lawful-interaction-cut-bundle-transplant-and-replica-dichotomy
  - THM-2545-word-stratified-hall-arrival-criterion-and-owner-word-transportation-hostile
  - THM-2548-seven-step-c91-transfer-and-full-norm-separation
  - THM-2549-future-pullback-target-neutrality-and-cemetery-hall-boundary
  - THM-2550-canonical-typed-row-double-nondegeneracy
  - THM-2555-natural-extension-sheet-charge-and-future-digit-boundary
related:
  - THM-2538-anchored-transverse-gain-and-common-ancestry-arrival-boundary
  - THM-2551-horizontal-transfer-transverse-projector-bicomplex-boundary
  - THM-2552-flat-common-base-gain-versus-unjoinable-wall-spectrum
  - THM-2554-translation-quotient-root-displacement-and-endpoint-swap-parity
script: 04-computation/lrc14_doubly_centered_transfer_hall_boundary_thm2557.py
output: 05-knowledge/results/lrc14_doubly_centered_transfer_hall_boundary_thm2557.out
script_sha256: 776c7af64d09d697bc02a38fe270c7f68691f8eca4788eef51e7e1f5f0d8117f
output_sha256: 177b1a82d5db160a89ce6057b7186a38312969e71e5d3f2bac9766e8d63c5d68
hash_basis: working-tree bytes (LF)
---

# THM-2557 -- the transferred interaction has an exact moment invoice

THM-2512 turns a lawful `7 x 13` response table into a doubly centered
signed interaction.  THM-2550(B)'s candidate exact computation reports that
both interactions extracted from the canonical typed row are nonzero.
THM-2548 proves that the
seven-step `C91` transfer preserves every root-charged rational mode.

There is a sharper integral statement on exactly the interaction lattice.
The transfer is injective there, but its image is not primitive.  Six
first-root moments modulo thirteen are the complete obstruction.  This
simultaneously explains why the signal survives and why that survival does
not enter THM-2545: doubly centered interactions live in a signed linear
space whose nonnegative cone is only the origin.

The resulting typed verdict is exact:

```text
nonzero lawful ANOVA interaction
  -- seven-step horizontal transfer --> nonzero signed current
  -- unsigned atom pushforward ------> zero or still signed, never a
                                       nonzero Hall coupling;

owner clock ell and response co-shift s
  -/-> selected-head root h and genuinely later active root b.             (1)
```

The crossed arrow in (1) is absent data, not a vanishing theorem about an
already constructed compatibility graph.

## 1. The doubly centered lattice

Let

```text
G=F_7 x F_13
```

with clock coordinate `k` and root coordinate `r`.  Define

```text
V_00={f in Z^G:
      sum_r f(k,r)=0 for every k,
      sum_k f(k,r)=0 for every r}.                            (2)
```

It is a primitive lattice of rank

```text
(7-1)(13-1)=72.                                              (3)
```

An integral basis is given, for `0<=k<6` and `0<=r<12`, by the rectangles

```text
b_(k,r)=e_(k,r)-e_(k,12)-e_(6,r)+e_(6,12).                  (4)
```

For `a in F_13^*`, use the THM-2548 pullback convention

```text
(D_a f)(k,r)=sum_(j=0)^6 f(k-j,r-aj).                        (5)
```

Reversing the convention permutes the roots and clocks and changes none of
the image, index, or cone statements below.

For `g in V_00`, define its clockwise first root moments

```text
mu_k(g)=sum_(r=0)^12 r g(k,r)                    mod 13,     (6)
```

and put

```text
W={m in F_13^7:sum_k m_k=0}.                                (7)
```

The definition in (6) is independent of translating the chosen root
representatives: such a translation changes `mu_k` by a multiple of the
zero row margin in (2).

## 2. Exact transfer sequence and Smith form

For every `a!=0`, there is an exact sequence of abelian groups

```text
0 -> V_00 --D_a--> V_00 --mu--> W -> 0.                     (8)
```

Equivalently,

```text
D_a(V_00)=ker(mu),
V_00/D_a(V_00) is isomorphic to (Z/13 Z)^6,                  (9)
```

and the Smith normal form of the `72 x 72` restriction is

```text
1^66 direct-sum 13^6.                                      (10)
```

Thus THM-2548's primitive rank-85 image on the full `91`-lattice acquires a
real `13`-adic invoice only after intersecting both zero-margin conditions.
There is no contradiction: intersection with a primitive sublattice need
not preserve primitivity.

### Preservation and injectivity

Both margin families in (2) are preserved by (5).  THM-2548 identifies the
full integral kernel of `D_a` as

```text
f(k,r)=c_k,                  sum_k c_k=0.                    (11)
```

If (11) also lies in `V_00`, its row sum is `13c_k=0`, hence every `c_k=0`.
Therefore `D_a` is injective on `V_00`.

### The moment annihilator

If `g=D_a f`, then modulo thirteen, after substituting `s=r-aj`,

```text
mu_k(g)
 =sum_(j=0)^6 [mu_(k-j)(f)+aj sum_s f(k-j,s)]
 =sum_(ell=0)^6 mu_ell(f)
 =0.                                                        (12)
```

The middle equality uses the row margins and the last uses the column
margins.  Thus `D_a(V_00)` is contained in `ker(mu)`.

The moment map is onto `W`.  Indeed, the explicit basis element

```text
b_(k,0),                     0<=k<6,                         (13)
```

maps to `e_k-e_6` in `F_13^7`.  These six vectors form a basis of `W`.
Consequently

```text
[V_00:ker(mu)]=|W|=13^6.                                   (14)
```

### The determinant has exactly the same index

Let `xi` and `zeta` be primitive seventh and thirteenth roots.  The complex
span of `V_00` consists exactly of characters with both frequencies
nontrivial.  Up to reversing both frequencies, the eigenvalue of (5) on
`(beta,alpha)` is

```text
d_a(alpha,beta)=sum_(j=0)^6 (xi^beta zeta^(alpha a))^j,
                  alpha!=0, beta!=0.                        (15)
```

For fixed `alpha`, write `t=zeta^(alpha a)`.  Then

```text
product_(beta=1)^6 d_a(alpha,beta)
 =(1-t^7)^6 / product_(beta=1)^6(1-t xi^beta)
 =(1-t^7)^5(1-t).                                          (16)
```

As `alpha` runs through `F_13^*`, so do `alpha a` and `7 alpha a`.  Since

```text
product_(alpha=1)^12 (1-zeta^alpha)=Phi_13(1)=13,           (17)
```

(16) gives

```text
|det(D_a restricted to V_00)|=13^5*13=13^6.                (18)
```

The index in (18) equals (14), so the inclusion after (12) is equality.
This proves (8)--(10).  Notice that the six factors in (10) are `13`, not
higher powers: the quotient has already been identified with the elementary
abelian group `W`.

### Coprime rectangular form

The mechanism is not accidental to `7` and `13`.  For an integer `q>=2` and
a prime `p` coprime to `q`, on `C_q x F_p` with nonzero skew `a`, the
`q`-step transfer restricted to the doubly centered lattice has

```text
0 -> V_00 --D_a--> V_00 --mu--> {m in F_p^q:sum m=0} -> 0,
Smith(D_a|V_00)=1^((q-1)(p-2)) direct-sum (q-1 copies of p). (19)
```

The same proof applies: the full kernel is root-uniform clock augmentation,
the rectangle basis makes `mu` onto, and the Fourier product is `p^(q-1)`.
Coprimality is essential because it makes the `q`-step root translation a
generator of `F_p`.

## 3. Application to the canonical typed row

THM-2512 defines, for a rational response table `A_(ell,s)`,

```text
I_(ell,s)=A_(ell,s)-R_ell-C_s+bar A.                        (20)
```

After transposition if desired, (20) lies in `V_00 tensor Q`.  THM-2550(B)'s
candidate exact computation reports two such tables, `d_M` and `d_C`, on the
canonical typed row, both nonzero with all `91` entries nonzero.  Conditional
on that exact input, clearing denominators and applying (8) yields, for every
`a in F_13^*`,

```text
D_a d_M !=0,                       D_a d_C !=0.              (21)
```

Conditional on THM-2550(B)'s candidate input, this is an exact
no-cancellation composition with THM-2548.  It is a statement about the
labelled horizontal interaction table only.
The exact witnesses at `(ell,s)=(0,0)` are byte-anchored by the companion:

```text
d_M(0,0)=27141080175744172866190363
           /51745643666528683956098414400 !=0,

d_C(0,0)=665502894989984158227797
           /10720466392295341281600 !=0.                    (22)
```

No part of (21) uses THM-2550(A).  Its positive drift is computed on the
`k=2` owner packet, whereas (20)--(22) come from the sufficiently large
lawful clock class `k=6 mod 16380`.  They are not one packet or one clock.
In particular, this theorem does not infer THM-2466's supplied common
oriented root/service base or THM-2471's collision-root identification from
the drift computation.

Finally, this is the canonical **typed non-cover** row of THM-2309.  Nothing
in (21) promotes it to a scalar cover or removes any of the `165` surviving
rows.

## 4. Why the transferred signal cannot be a Hall table

Every nonzero element of `V_00 tensor Q` has both positive and negative
entries: a nonzero row has sum zero.  Consequently

```text
(V_00 tensor Q) intersect Q_(>=0)^G={0}.                    (23)
```

The currents in (21) are therefore necessarily signed.  THM-2545's joint
table, by contrast, is

```text
C^sigma_(h,b)=nu({omega:sigma(omega)=sigma,
                        head(omega)=h,
                        later-root(omega)=b}) >=0.           (24)
```

It is a measure pushforward on one common ancestry base.  ANOVA entries are
not atom masses, and changing signs or taking absolute values does not
preserve (20), (8), or a physical ancestry interpretation.

This remains true after any honest unsigned partition pushforward `P` that
preserves augmentation.  Since every `g in V_00` has total mass zero,

```text
Pg>=0  implies  sum Pg=sum g=0  implies Pg=0.               (25)
```

Thus no nonzero Hall coupling can be obtained merely by grouping the signed
cells of (21).  A nonlinear positive lift would be new mathematics and would
still have to construct the two maps in (24).

## 5. There is no typed-row Hall graph yet

The live coordinates must not be renamed cosmetically:

| object | first coordinate | second coordinate | sign/type |
|---|---|---|---|
| THM-2550(B) interaction | owner/source phase `ell in F_7` | intervention co-shift `s in F_13` | signed ANOVA current |
| THM-2548 transfer | clock chart `k in F_7` | root deck `r in F_13` | signed horizontal coefficient |
| THM-2545 Hall table | selected-head root `h in F_13` | genuinely later active root `b in F_13 union {partial}` | nonnegative common-base mass |

The typed row does retain the word `sigma={a}`.  It does **not** currently
provide atomwise maps from `(ell,s)` to `(h,b)`, nor a common rational-BV
refinement on which both `h` and `b` are defined.  Even the vertex sets have
different typed shapes: `F_7 x F_13` is not the ordered-root square
`F_13 x F_13`.

Therefore the question

```text
"does the actual word-resolved typed-row compatibility graph have a
 Hall deficit?"
```

is presently **undefined**, not answered negatively.  No such graph has yet
been constructed.  Fourier support, the moment quotient (8), and all `5,184`
nonzero cut coefficients describe linear observables; none supplies the
edge relation of THM-2545 equation (7).

## 6. The exact zero-arrival boundary

There are two useful controls, with different epistemic status.

### The currently proved-field completion: cemetery only

THM-2549 proves that every positive later field currently available on the
same base is old-target-neutral.  Its only equivariant partial selector is
the cemetery label.  For any source masses `p_h`, the induced table is

```text
C_(h,partial)=p_h,                    C_(h,b)=0 for b in F_13. (26)
```

After deleting the semantic diagonal, (26) is unchanged.  Every nonempty
source set has cemetery neighbourhood of mass `M>=p(S)`, so all Hall
inequalities hold and the arrival mass is zero.  This is the only lawful
completion furnished by the **currently proved fields**: it records absence
of a genuinely later root rather than inventing one.

### The smallest all-root symmetric cemetery-free hostile

If one asks only for an abstract cemetery-free table compatible with
simultaneous root translation, THM-2554's one-offset orbit is the smallest
nonempty equivariant support:

```text
C(h,b)=1_(b=h+1),                         h,b in F_13.       (27)
```

It has thirteen unit atoms, uniform margins, and zero diagonal.  On its
matching support, every Hall inequality is an equality because
`N(S)=S+1`.  Retaining all seven clock labels gives the analogous `91`-atom
control.  THM-2555 realizes (27) as the positive symbolic cylinder
`d_(L+1)=d_1+1` with every intervening sheet digit free.  This makes it an
exact natural-extension hostile, but still does **not** place the unique
THM-2461 target-active role on that cylinder or realize it on the typed row.
Without the translation requirement, the two-atom aligned/swap hostile of
THM-2545 is smaller still.

A diagonal-only support provides the positive control: after deleting the
diagonal, a positive singleton source has empty neighbourhood and hence an
immediate Hall deficit.

## 7. The precise missing coordinate

Write one natural-extension point schematically as

```text
y=0.d_1 d_2 ... d_L e_1 e_2 ...  in base 13.                (28)
```

THM-2555 proves that, for a unit role `y={wx}` with `x=(u+h)/13`,
the top old ancestry-sheet digit, corrected by its base carry, recovers `h`
exactly.  That is an **old-action ancestry root**.  The immediate future
digit

```text
e_1=floor(13 T^L y)                                         (29)
```

is a distinct coordinate.  With its later base carry it recovers the future
root `h_L=floor(13T^Lx)`, not `h`.  Current canon has no theorem deciding
whether the unique genuinely later target-active role on the selected packet
is semantically typed by the old-action ancestry chart or this future-action
chart, nor one fixing the invariant gauge between the relevant root torsors.

The missing datum can therefore be stated in either of two equivalent typed
forms:

1. an atomwise common-base map `b:Omega->F_13` for the genuinely later
   target-active role, together with the theorem selecting its old- or
   future-action chart, covariance, carry, and invariant gauge; or
2. after identifying the old head chart and the later chart as one root
   torsor, the displacement

```text
delta(omega)=b(omega)-h(omega) in F_13.                     (30)
```

Only after (30) is defined does THM-2554 reduce the Hall problem to the
zero-displacement orbit.  The first-moment quotient `mu` in (8) is not this
displacement: `mu` is a linear invariant of a signed `7 x 13` current,
whereas (30) is an atomwise relation between two semantic root maps.

## 8. Exact referee and stopping boundary

Run

```bash
python3 04-computation/lrc14_doubly_centered_transfer_hall_boundary_thm2557.py
python3 -O 04-computation/lrc14_doubly_centered_transfer_hall_boundary_thm2557.py
```

Both executions must byte-match

```text
05-knowledge/results/lrc14_doubly_centered_transfer_hall_boundary_thm2557.out.
```

The dependency-free referee constructs the rectangle basis and all twelve
`72 x 72` transfer matrices; verifies both margins, vanishing image moments,
rank `66` modulo thirteen and rank `72` modulo two for every nonzero slope;
computes all twelve exact determinants `13^6`; checks the six explicit
moment lifts and hence the Smith form; exhibits signed transferred controls;
byte-anchors, but does not independently recompute, the two THM-2550(B)
witnesses; exhausts all `8,192` root subsets
for both the cemetery and offset zero-arrival Hall controls; and detects the
diagonal singleton deficit.  It deliberately performs no large typed-row
recomputation and no Lean build.

An independent hostile audit rederived the pullback convention, moment
annihilator, six explicit lifts, determinant product (16), index comparison,
Smith form, and coprime rectangular generalization, and reran both normal and
optimized companions against the stored hashes.  It also enforced the
conditional status of the THM-2550(B) application: the present companion
byte-anchors those reported witnesses but does not recompute its large table.

The highest-leverage next experiment is consequently not another
nonvanishing census.  It is to atomize the unique target-active unit failure
together with the selected head on a common natural-extension sheet, decide
whether its semantic root is THM-2555's old carry root or future carry root,
and tabulate the actual word-resolved displacement support (30).  Support
confined to zero would force Hall arrival; an actual permitted nonzero
matching would exhibit the precise zero-arrival escape.  Until that semantic
object exists, LRC(14) remains open. **QED.**
