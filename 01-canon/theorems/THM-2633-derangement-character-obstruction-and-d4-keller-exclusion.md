---
id: THM-2633
title: "Affine inertia normal generation, point stabilizers, and D4 Keller exclusion"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  Let a complex
  polynomial Keller map of affine n-space have finite generic degree and
  geometric monodromy G on its generic
  sheets.  Every irreducible target divisor meets the open image in a dense
  open set, so every Jelonek component has at least one finite affine inverse
  branch and its local inertia fixes a sheet.  Purity and the trivial etale
  fundamental group of affine space then force the actual affine Jelonek
  inertia subgroups to normally generate G.  Hence every point stabilizer H
  normally generates G and surjects onto G^ab; restriction
  Hom(G,C_l)->Hom(H,C_l) is injective for every prime l.  Also every Jelonek
  component has 1<=k_D<=d-1 finite branches.  Consequently G cannot admit any
  nontrivial quotient killed by all sheet-fixing elements; in particular no
  nontrivial regular transitive group can occur, the source function-field
  extension contains no nontrivial Galois subextension of the target field,
  and no nonzero character can be supported on derangements.  In D4 the
  normal closure of the sheet-fixing elements is the proper source-deck
  kernel J, whose nonzero quotient coset consists of the two edge reflections
  and two four-cycles.  Hence D4 is impossible as the
  geometric monodromy of a polynomial Keller map in every dimension.  Among
  transitive quartic groups the same gate excludes C4 and V4, while it does
  not exclude A4 or S4.  Thus the existing planar and G1 degree-four
  frontiers reduce from D4/A4/S4 to A4/S4.  On character-support components,
  A4 is forced to have three-cycle inertia and exactly one finite branch,
  while S4 is forced to have transposition inertia and one or two finite
  branches.  No A4 or S4 exclusion, JC(2), general JC, or DC(2) follows.
source: inverse-spectral-residues-2026-07-27
related:
  - THM-2465-g1-exclusion-package-for-degree-four-twojet-keller
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2612-d4-deck-pole-tax-and-depressed-resolvent-gcd-gate
  - THM-2621-planar-degree-four-inverse-spectral-keller-congruence-and-sheet-defect-pole-ledger
  - THM-2627-d4-jelonek-quadratic-character-rank-and-component-gate
  - THM-2628-d4-opposite-pair-escape-and-deck-pole-census
script: 04-computation/jacobian_derangement_character_d4_exclusion_thm2633.py
output: 05-knowledge/results/jacobian_derangement_character_d4_exclusion_thm2633.out
script_sha256: 70535c0c85d949257b6a1af2a29358e83f603955fe6a5c75c95b6763971f2a83
output_sha256: e949156e96895a1c6ceb5ecd8a8ea50d1207070ae319c7a3ff47b4bcc7fb7dfb
hash_basis: working-tree bytes (LF)
---

# THM-2633 -- affine boundary inertia normally generates Keller monodromy

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**

The missing coordinate in the `D_4` boundary picture is not another quartic
coefficient or residue.  It is the image of the original affine morphism.
A Keller map is etale and therefore open, while affine space has no
nonconstant units.  Together these facts force a finite affine branch over
the generic point of every target divisor.  This turns the deck-character
support separation from a two-component normal form into a contradiction.

The mechanism is group-theoretic and applies in every dimension.  The
character obstruction is its abelian shadow.

## 1. General fixed-point normal-generation theorem

Let

```text
F:A^n_C -> A^n_C,                 n>=1,                  (1)
```

be a polynomial map with

```text
det JF in C*,

d=[C(x_1,...,x_n):C(F_1,...,F_n)]<infinity.              (2)
```

Let `Omega` be its `d` geometric generic inverse sheets and let

```text
G <= Sym(Omega)                                             (3)
```

be the transitive geometric monodromy group.  Let `A_F` be the affine Jelonek
set, let `I_D` be geometric inertia at each irreducible component `D`, and
define

```text
N_Jel=normal closure in G of {I_D:D component of A_F},     (4)

N_fix=normal closure in G of {g:Fix_Omega(g)!=empty}.     (4a)
```

Then the strongest form of the theorem is

```text
N_Jel=G,                    and therefore N_fix=G.         (5)
```

If `H` is the stabilizer of one sheet, an element fixes some sheet exactly
when it lies in a conjugate of `H`.  Hence `N_fix` is the normal closure of
`H`, and the point-stabilizer form of the primary theorem is

```text
normal closure_G(H)=G.                                  (5a)
```

There is an equivalent permutation-theoretic formulation which exposes the
lost quotient.  Identify `Omega=G/H` and put `N=normal closure_G(H)`.  The
`N`-orbits are indexed by `G/N`, and the induced action of `G/N` on those
orbits is regular.  Moreover every quotient of `G` which kills `H` factors
through `G/N`.  Thus `G/N` is the maximal regular block quotient of the
sheet action, and (5a) says precisely

```text
a Keller monodromy action has no nontrivial regular block quotient.       (5aa)
```

In particular, if the generic source extension is already Galois, then its
sheet action is regular and (5aa) forces degree one.  This consequence uses
the full normal-generation theorem; it need not be visible in abelianization.

Equivalently, write the target, source, and Galois-closure function fields as

```text
K=C(F_1,...,F_n) subset L=C(x_1,...,x_n) subset M,

G=Gal(M/K),                    H=Gal(M/L).                (5ab)
```

There is no intermediate field `K subsetneq E subseteq L` for
which `E/K` is Galois.  Indeed such an `E` would correspond to a proper
normal subgroup `Gal(M/E)` containing `H`, hence containing its normal
closure `G`, a contradiction.  Thus a Keller source extension has no
nontrivial Galois subextension.  The familiar exclusion of a nontrivial
Galois source extension is only the endpoint `E=L` of this stronger field
statement.

This gives the derangement-character obstruction as an immediate abelian
corollary.  Suppose that for some prime `ell` there is a nonzero character

```text
chi:G -> C_ell                                           (5b)
```

with

```text
chi(g)!=0  implies  Fix_Omega(g)=empty.                  (5c)
```

Then `chi` vanishes on every generator of `N_fix`, so `N_fix` lies in the
proper subgroup `ker(chi)`, contrary to (5).  Hence no such Keller map
exists.

In character language, this gives the following necessary condition for a
transitive group to occur as Keller monodromy:

> every nonzero prime-cyclic character is nonzero on at least one element
> that fixes a generic sheet.

Since a character is conjugacy invariant, (5c) is equivalent to
`chi|_H=0`.  The abelian shadow can therefore be stated in either of the
following equivalent forms:

```text
res_H^G:Hom(G,C_ell) -> Hom(H,C_ell) is injective
for every prime ell;

H -> G^ab is surjective.                                 (5d)
```

For the second equivalence, if the image of `H` in the finite abelian group
`G^ab` were proper, the nonzero quotient would have a `C_ell` quotient for
some prime `ell`; conversely such a character witnesses failure of
surjectivity.  Thus the derangement-character theorem and point-stabilizer
abelianization gate are the same abelian corollary.  The normal-generation
law (5a) is strictly stronger: it also forbids nonabelian quotients killed by
all point stabilizers.

The proof has two independent halves: an open-image branch lemma and a
purity/normal-generation lemma.  Kummer theory then localizes the abelian
corollary to a specific Jelonek component.

## 2. Open image forces a branch over every divisor

The Keller condition makes `F` etale.  In particular it is dominant,
quasi-finite, and open.  Let

```text
D=V(f) subset A^n_C                                      (6)
```

be any irreducible target divisor.  It cannot be wholly omitted by `F`.
Indeed, if `D` were disjoint from `F(A^n)`, then

```text
f o F in C[x_1,...,x_n]
```

would be nowhere zero.  By the Nullstellensatz it would be a unit, hence a
constant `c in C*`.  But dominance makes

```text
F^*:C[y_1,...,y_n] -> C[x_1,...,x_n]
```

injective, so `F^*(f-c)=0` would force `f=c`, a contradiction.

Therefore

```text
D intersection F(A^n) != empty.                          (7)
```

Because `F(A^n)` is open, the intersection in (7) is a nonempty relative
open subset of the irreducible divisor `D`, hence dense.  More explicitly,
the base change `F^(-1)(D)->D` is etale and quasi-finite.  It has finitely many
irreducible components, and the union of their images contains the dense open
set in (7), so at least one component dominates `D`.  At its generic
valuation, and then after passing to a geometric generic point and a
transverse strict-henselian DVR, an inverse sheet has a finite centre in the
affine source.  No global section over all of `D` is being asserted.

There is an independent algebraic proof of the same generic-hit statement.
Dominance makes `f o F` nonconstant.  Every irreducible component `E` of its
zero hypersurface has dimension `n-1` and maps into `D`.  Since an etale map
is quasi-finite,

```text
dim closure(F(E))=dim(E)=n-1.
```

The irreducibility of `D` forces `closure(F(E))=D`.  At the generic point of
`E`, etaleness again supplies the inertia-fixed inverse branch.  This
dimension proof cross-checks the open-image proof without relying on a
set-theoretic interpretation of the generic point.

Write

```text
S_D={generic inverse sheets with a finite affine centre over D},

k_D=|S_D|.                                                (8)
```

Then for every irreducible target divisor, and in particular for every
Jelonek component,

```text
k_D>=1.                                                   (9)
```

There is also a strict upper bound on a Jelonek component.  Factor `F` by
Zariski main as the open immersion `A^n -> Xbar` into the finite
normalization of the target in the source function field.  If all `d`
geometric points over the generic point of `D` belonged to the affine open,
the finite image of the closed complement would miss that generic point.
After shrinking around it, `F` would equal the finite normalization and
would be proper, contradicting that `D` is a nonproperness component.  Hence

```text
1<=k_D<=d-1             for every Jelonek component D.   (9a)
```

Let `I_D<=G` be geometric local inertia.  A finite branch passes through a
point where `F` is etale, so the local inverse extends uniquely across the
transverse disk.  Hence inertia fixes that branch pointwise:

```text
S_D subseteq Fix_Omega(I_D),

Fix_Omega(I_D) != empty.                                 (10)
```

This also proves the useful global side statement

```text
codim(A^n minus F(A^n))>=2.                              (11)
```

For `n=2`, a non-surjective Keller image can therefore miss only finitely
many target points, never a curve.

## 3. Boundary inertia normally generates

Let `A_F` be the nonproper-value set.  By the standard purity theorem for a
generically finite polynomial map between affine spaces, it is empty or a
pure hypersurface.  Write its reduced equation as

```text
A_F=V(f_1...f_c),

U=A^n minus A_F.                                         (12)
```

Over `U`, the restriction of `F` is finite etale.  Its connected Galois
closure is therefore a finite-etale `G`-cover.

Let `N_Jel` be the normal subgroup generated by the geometric inertia groups
`I_D` over all irreducible components of `A_F`.  The connected quotient cover
with group `G/N_Jel` has trivial inertia at every affine codimension-one
boundary valuation.  Normalize `A^n` in its function field.  The resulting
finite map is unramified in codimension one, so Zariski--Nagata purity makes
it finite etale everywhere over `A^n`.  But

```text
pi_1^et(A^n_C)=1.
```

Therefore the connected quotient is trivial and

```text
N_Jel=G.                                                 (12a)
```

Every `I_D` fixes the nonempty set `S_D` pointwise by (9)--(10).  Hence every
element used to generate `N_Jel` is sheet-fixing, so

```text
G=N_Jel subseteq N_fix subseteq G.
```

This proves the primary assertion (5).  If `A_F` is empty, the same argument
starts with `N_Jel=1` and says directly that `G=1`.

There are two scope guardrails.  First, the quotient cover is connected
because it is the quotient of the connected Galois closure by the normal
subgroup `N_Jel`; triviality is not being inferred from a disconnected finite
etale algebra.  Second, projective infinity is irrelevant to the purity step:
the normalization is being extended over the affine target `A^n`, and every
missing finite divisor is already a component of `A_F`.

### 3.1 Kummer localization of the character corollary

The stronger normal-generation proof does not need a character.  When a
prime-cyclic character is present, Kummer theory says exactly where its
abelian obstruction would have to occur.  Since

```text
O(U)=C[y_1,...,y_n,1/(f_1...f_c)]
```

is a UFD, `Pic(U)=0`, and its units modulo `ell`-th powers are freely
generated by the classes of the `f_i`.  The Kummer sequence gives

```text
H^1_et(U,mu_ell)
 =O(U)^*/O(U)^{*ell}
 =(Z/ell)^c.                                             (13)
```

Surjective monodromy `pi_1(U)->G` makes pullback injective on characters:

```text
Hom(G,C_ell) -> H^1_et(U,mu_ell).                        (14)
```

Thus the nonzero `chi` in (5b) has a nonzero valuation coordinate at some
component `D=V(f_i)`.  Equivalently, the restriction of `chi` to tame local
inertia at `D` is nonzero.  For a generator `g_D` of the cyclic inertia image,

```text
chi(g_D)!=0.                                             (15)
```

Projective infinity cannot absorb this coordinate.  Equivalently, if `chi`
were trivial on every finite codimension-one inertia group, purity would
extend its cyclic quotient etale over all of `A^n`; affine space has no
nontrivial connected finite-etale cover.  This is the same content as the
affine Kummer calculation (13).

Equation (10) says that `g_D` fixes a sheet, whereas (5c) says that it is a
derangement.  This recovers the character corollary componentwise.  If `c=0`,
(13)--(14) already give the contradiction before any component is chosen.

The roles of the hypotheses are now explicit:

- etaleness supplies openness and extension of every finite local inverse;
- the trivial unit group of affine space forbids an omitted divisor;
- purity forces the affine boundary inertia groups to normally generate all
  monodromy;
- Kummer theory locates every nonzero abelian character on an affine Jelonek
  component; and
- the derangement-support condition converts that parity into loss of all
  finite branches.

## 4. The `D_4` fixed-point closure is the deck kernel

Use the square conventions

```text
G=D4=<r,s | r^4=s^2=1, srs=r^(-1)>,

Omega={0,1,2,3},       r=(0 1 2 3),       s=(1 3),

z=r^2=(0 2)(1 3).                                      (16)
```

The point stabilizer is `H=<s>`, and

```text
J=N_G(H)=<s,z>.                                          (17)
```

The only nonidentity elements of `D_4` fixing a sheet are the two diagonal
reflections `s` and `zs`.  Their normal closure is

```text
N_fix=<s,zs>=J,                 |J|=4<8=|D4|.           (17a)
```

This already contradicts the primary normal-generation law (5).  The
two `J`-orbits are the two opposite vertex pairs, and `G/J=C2` swaps them
regularly.  Thus the obstruction is exactly the nontrivial regular block
quotient in (5aa).  The source-deck character is its abelian witness.

The quadratic source-deck intermediate has character

```text
chi_deck:G -> C2,             ker(chi_deck)=J.           (18)
```

Its nonzero set is

```text
G minus J={r,r^3,rs,r^3s}.                               (19)
```

The first two elements are four-cycles.  The last two are the edge
reflections `(0 1)(2 3)` and `(0 3)(1 2)`.  All four are fixed-point-free in
the root action.  Thus `chi_deck` satisfies (5c), and either form of the
general theorem gives

```text
no polynomial Keller map A^n_C -> A^n_C
has geometric monodromy D4, for any n.                  (20)
```

This proof uses the quadratic field `L^J/K`, not the different
matching/discriminant quadratic `L^(G intersect A4)/K`.  It does not require a
polynomial source quotient, a second deck transformation, a residue, or a
choice of primitive quartic coordinate.

## 5. Exact quartic group boundary

There are five transitive subgroups of `S_4` up to conjugacy.  Applying the
normal-closure theorem and its character shadow gives

| group | `|N_fix|` | character used | nonzero support | result |
|---|---:|---|---|---|
| `C4` | `1` | quotient by `<r^2>` | two four-cycles | excluded |
| `V4` | `1` | any nonzero `C2` character | two double transpositions | excluded |
| `D4` | `4` | source-deck character (18) | edge reflections and four-cycles | excluded |
| `A4` | `12` | quotient by the normal `V4` | three-cycles, each fixing one sheet | not excluded |
| `S4` | `24` | sign | transpositions as well as four-cycles | not excluded |

Thus the regular `C4,V4` lanes and the dihedral `D4` lane are removed by this
one gate.  (The solvable group `A4` is not removed.)  In the existing
degree-four reductions of THM-2598 and THM-2465, the live list

```text
D4, A4, S4
```

sharpens to

```text
A4, S4.                                                  (21)
```

The table is a boundary of the method, not evidence against `A_4` or `S_4`.
For `A_4`, every element detected by its `C_3` quotient is a three-cycle and
has a fixed sheet.  For `S_4`, the sign character detects transpositions,
which have two fixed sheets.  Their nonzero characters can therefore pay
Kummer parity on a component that still lies densely in the affine image.

Combining Kummer nontriviality with (9)--(10) gives the next exact typed
targets.  In an `A_4` lane, some character-support component must have
three-cycle inertia.  Its fixed set is a singleton, so

```text
A4 support component:       three-cycle inertia, k_D=1. (21a)
```

In an `S_4` lane, some sign-support component has odd inertia.  A four-cycle
is forbidden by (10), leaving a transposition; its fixed set has two sheets,
so

```text
S4 support component:       transposition inertia, k_D=1 or 2. (21b)
```

These are necessary component types, not exclusions.  They are the precise
residual boundary objects on which the degree-four search should now focus.

More generally, **every** nontrivial regular transitive group is excluded,
without any abelian-quotient hypothesis.  In a regular action the point
stabilizer is trivial, so its normal closure is trivial rather than `G`.
Equivalently, every nonidentity element is a derangement and `N_fix=1`.
The character obstruction alone sees only regular groups with a nonzero
prime-cyclic quotient; the primary normal-generation theorem also sees
perfect regular groups.

## 6. Why the earlier two-component atlas was a near miss

THM-2627 forces the `D_4` deck character to be nonzero on some Jelonek
component `B`.  The finite `D_4` census behind THM-2628 and the original
THM-2633 draft then gives

```text
chi_deck(I_B)=1  implies  k_B=0.                         (22)
```

That implication is correct.  What was missing was (9): an etale polynomial
map from affine space has `k_D>=1` on **every** target divisor.  Hence (22) is
already the contradiction; no second component, deck-pole ownership count,
resultant coefficient, or residue is needed.

The former two-component character matrices

```text
[0 1]          [0 1]
[1 0]    or    [1 1]
```

remain valid finite group atlases, but their deck-odd column necessarily has
`k=0`.  They destroy exactly the global open-image/trivial-unit coordinate
captured by Section 2 and therefore cannot be the boundary atlas of a Keller
map.  This identifies the first failed implication rather than merely
discarding the models.

The new proof also explains why the exact trace and opposite-pair residue
laws of THM-2621 did not close the lane: logarithmic residues are downstream
of a contradiction already visible in incidence between the affine image and
each boundary divisor.

## 7. Hostile controls and necessity of every input

### 7.1 Dropping etaleness destroys generic branch survival

The dominant polynomial map

```text
F_0(x,y)=(u,v)=(x,xy)                                   (23)
```

has nonproper divisor `D={u=0}` and

```text
F_0(A^2)={u!=0} union {(0,0)}.
```

Thus its image meets `D` but not densely; at the generic point of `D` one has
`k_D=0`.  There is no contradiction because

```text
Jac(F_0)=x,
```

so the image is not open and the finite inverse does not extend etale across
the special point.  This is the hostile to replacing “open image” by mere
dominance or by set-theoretic non-omission.

### 7.2 Dropping the affine-space unit condition permits an omitted divisor

The open immersion

```text
G_m x A^(n-1) -> A^n
```

is etale and open but omits `{y_1=0}`.  Its pullback `y_1` is a nonconstant
unit.  Hence the no-omitted-divisor lemma is a statement about the source
`A^n` (more generally, a source with only constant units), not about arbitrary
etale morphisms.

### 7.3 A character without derangement support is compatible with a branch

The strict-henselian cubic control

```text
R_fix=tT^3-T=T(tT^2-1)                                  (24)
```

has one finite branch `T=0` and two escaping branches `T=+-t^(-1/2)`.
Inertia is a transposition: its `S_3` sign is nonzero, but it fixes the finite
branch.  Thus Kummer parity plus `k_D>0` is not by itself contradictory.  The
support condition (5c) is load-bearing and fails exactly for `A_4` and `S_4`
in the quartic table.

### 7.4 The old `D_4` local rows are sharp but nonglobal

The all-escaping controls

```text
R_edge=(tT^2-1)(tT^2-2),

R_four=tT^4-1                                             (25)
```

have edge-reflection and four-cycle inertia after the opposite-pair labelling
of THM-2628.  Both specialize to degree zero, exactly as (22) predicts.  They
prove that the group/inertia census itself is sharp.  They are not global
Keller maps because a Jelonek component with `k=0` would violate Section 2.

Finally, `z -> z^2` shows why ramified finite maps are irrelevant controls:
its image is open and surjective, but it is not Keller and its character
ramification lies on a finite critical-value divisor rather than a Jelonek
component.

### 7.5 The affine target's trivial etale fundamental group is load-bearing

After restricting the same map to

```text
G_m -> G_m,                    z |-> z^2,
```

it becomes a connected finite-etale `C_2`-cover with regular monodromy.  Its
nonidentity element is a derangement and there is no finite boundary inertia
inside `G_m`.  The quotient survives because `pi_1^et(G_m)` is nontrivial;
the missing valuations `0` and infinity lie outside that target.  This is the
sharp hostile to replacing the affine target `A^n` by an arbitrary smooth
affine variety in the purity/normal-generation argument.

## 8. Exact verification and scope

Run

```text
python3 04-computation/jacobian_derangement_character_d4_exclusion_thm2633.py
python3 -O 04-computation/jacobian_derangement_character_d4_exclusion_thm2633.py
```

Both modes byte-match the stored transcript.  The companion constructs all
five transitive quartic groups, computes the fixed-point normal-closure orders
`1,1,4,12,24`, checks the displayed character kernels and fixed-point counts,
proves that the `C4`, `V4`, and `D4` character supports consist entirely of
derangements, and verifies that the `D4` normal closure is exactly `J`.  It
also exhibits the surviving fixed-point types in `A4/S4`, exhausts the `D4`
deck-odd survivor rows, and checks the exact local and non-etale hostiles
(23)--(25).

The theorem excludes one full geometric-monodromy lane in every dimension.
The point-stabilizer formulation (5a) and strict branch bound (9a) are formal
consequences of the same proof; the exact companion checks their quartic
character incarnation.
It does not exclude `A_4` or `S_4`, prove that a Keller map has degree one, or
settle `JC(2)`, the general Jacobian conjecture, or `DC(2)`.  The highest-
leverage remaining degree-four question is now to find an invariant of the
`A_4/S_4` inertia classes that have fixed sheets, rather than to continue
refining the impossible `D_4` deck-pole atlas.

**QED.**
