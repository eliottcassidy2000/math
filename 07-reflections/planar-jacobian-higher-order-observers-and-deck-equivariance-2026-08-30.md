# Planar Jacobian: higher-order observers, blowdown forms, and deck equivariance

**Session:** root/planar-jc-higher-order-20260830, 2026-08-30  
**Current truth:** [THM-4289](../01-canon/theorems/THM-4289-a23-blowdown-observer-kahler-dualizing-quotient.md)
is **PROVED FORMAL-LOCAL RELATIVE TO THM-4272/4279/4280/4284/4288 +
FINITE-EXACT**. [THM-4290](../01-canon/theorems/THM-4290-exact-weight-twelve-deck-equivariant-visible-quotient-exclusion.md)
is **PROVED RELATIVE TO THM-4012/4230 + FINITE-EXACT**. The complete
`U*Z*D*Lambda!=0` exact-`M=12` gate is closed. The coefficient walls, seam
entry, `JC(2)`, and `DC(2)` remain **OPEN**.

## Session verdict

The Ye--Xu higher-order-truth paper supplied a useful but deliberately narrow
research move: before trusting a compressed observer, test whether its fibres
are congruence classes for the next operation and pure for the target
predicate. On the `A_23` blowdown chain that principle becomes an exact linear
criterion,

```text
depth predicate P_j factors through q
iff ker(q) is contained in contact depth j.
```

This did not activate the paper's infinite-antichain/topos obstruction. The
fixed contact has only twelve order states. It instead led to a complete local
calculus distinguishing ambient Kahler forms from dualizing residues.

The wildcard then overtook the anchor. The same bookkeeping exposed a
coordinate that had been treated as a nuisance of stable reduction: the
12-fold deck action itself. Keeping that action, rather than quotienting it
away, forces every interior Keller specialization into a visible
degree-multiple-of-four quotient and contradicts the inherited degrees
`34/42`. This closes the complete exact-`M=12` gate, while preserving a sharp
firewall at `Lambda=0`.

## Inheritance pass

- **Closest proved mechanism:** THM-4288's terminal triple circuit and
  THM-4230's saturated genus-two quotient lattice.
- **Canonical hostile:** the actual degree-`34/42` branch map has formal-log
  order one and survives the resolved nodal tree, while failing the final raw
  blowdown.
- **Corrected near miss:** a saturated graph image is not a Cartesian etale
  pullback, and a regular dualizing/Poincare residue is not an ambient Kahler
  form.
- **Least-used sidecar:** the deck character of the scaled target differential.
  It records exactly which Hom eigenspace a genuine Keller specialization may
  occupy.

## Anchor / Niche / Wildcard portfolio

| lane | object | target | outcome |
|---|---|---|---|
| Anchor | exact-`M=12` Keller specialization | eliminate degree `34/42` | interior gate **CLOSED** by THM-4290 |
| Niche | `A_23` resolution triples and branch forms | identify the raw descent obstruction | exact `conductor/J_f` quotient **PROVED** by THM-4289 |
| Wildcard | `mu_12` deck action on scaled source and target | retain a symmetry discarded by coarse specialization | overtook anchor; forces visible quotient |

The local and cyclic outcomes are complementary. THM-4289 says exactly why
the invariant differential cannot be promoted from dualizing to ambient
Kahler by formal rhetoric. THM-4290 avoids that promotion entirely on the
interior gate by using equivariance before taking the lossy quotient.

## Live concept board

| concept | source / representation | invariant or predicate | operation | lost datum / restoring sidecar |
|---|---|---|---|---|
| blowdown depth | `D_s:bz(z-b^s)=0` and three branch logs | descent to the plane triple | contract next exceptional curve | coefficient of `b^s`; ambient circuit |
| sparse Hasse observer | channels `c_1,c_2,c_4,c_7` | `L in b^j` | pass to observer fibres | first invisible coefficient; kernel depth |
| differential carrier | `h*eta` in the dualizing module | ambient Kahler membership | normalize / take residue | class `[h] in conductor/J_f` |
| deck state | scaled `(sigma,S,P,X,Y)` | equivariance of the Keller map | `sigma->xi*sigma` | character if one retains only special equations |
| quotient degree | `C->B=C/<tau^6>->E_0` | degree modulo four | sixth iterate / quotient | saturated integral Hom lattice |
| wall specialization | raw `A_23` plus resolved Keller model | allocation of generic degree | stable reduction / blowup | multiplicities, extra components, response invoice |

Every pull changed another lane. The dualizing hostile made it unproductive to
seek a generic residue-extension slogan; this promoted the deck character.
The deck quotient closed the interior gate, which in turn changed the role of
the `A_23` work: it is now a wall diagnostic and a test for the missing
degree-allocation theorem, not the sole route through the interior.

## Result 1: the complete local blowdown calculus

For

```text
D_s=Spec k[[b,z]]/(b*z*(z-b^s)),
```

equal-value branch functions `(e,r,c)` descend exactly when

```text
c-r-e(b^s)+a in b^(s+1)k[[b]].
```

Regular normalized branch one-forms extend from the ambient plane exactly
when

```text
c-r-s*b^(s-1)e(b^s) in b^s k[[b]].
```

Differentiation identifies the two obstruction spaces in characteristic zero.
In dualizing coordinates,

```text
nu_*Omega_normalized=conductor*eta,
im Omega_D^tf=J_f*eta,
conductor/J_f ~= k[[b]]/(b^s).
```

The minimal hostile `bz*eta=(0,0,db)` is regular on all normalized branches
and has a regular numerator, but is not ambient Kahler. At the two-branch
`A_m` endpoint the quotient has length `m-1`; for `A_12` it recovers exactly
the eleven post-value characters.

This answers *why* the prior route stopped. Keller residue identities naturally
land on the dualizing side of the strict inclusion. The next local test is the
explicit ideal membership `h in (f_b,f_z)`, not another regularity adjective.

## Result 2: exact observer congruence, and where Ye--Xu stops

For `F^jV=V intersect b^jR/b^12`, a linear observer `q` preserves the descent
predicate `P_j` precisely when `ker(q) subset F^jV`. Here `P_(s+1)` is
descent to curve `D_s`. Applied to the two minimal integral channel triples
after complexification:

- `q_124` has kernel order seven, so it first fails at predicate depth eight,
  equivalently on curve `D_7`;
- `q_247` has kernel order one, so it first fails at predicate depth two,
  equivalently on curve `D_1`.

On the actual Eisenstein Hom lattice both are injective and preserve every
depth predicate. The exact order-state sets are twelve for arbitrary jets,
five for the complexified Hom span, four for the actual Hom lattice, and one
for the degree-`34/42` shells.

Thus there is no infinite orbit union outside a finitary certificate algebra
here. The source paper's higher-order theorem does not transfer. Its useful
contribution is the fibre-purity/congruence test, already isolated in the
[paper audit](../05-knowledge/reference/CORE-PAPERS-HEYTING-HIGHER-ORDER-TRUTH-2026-08-29.md).

## Result 3: deck equivariance closes the interior exact-M12 gate

With `Q=sigma^12`, the exact scalings are

```text
s=sigma^-1 S,        p=sigma^-2 P,
A=sigma^-4 X,        C_target=sigma^-6 Y.
```

The deck generator acts by

```text
(S,P)->(xi S,xi^2 P),
(X,Y)->(xi^4 X,xi^6 Y)=[-omega](X,Y).
```

Because the original Keller map is defined before base change, its scaled map
is exactly equivariant. Invariant graph closure and properness specialize this
identity to the unique positive-genus component:

```text
m o tau=[-omega] o m.
```

The sixth power makes `m` invariant under `tau^6:(S,P)->(-S,P)`, hence it
factors through the degree-two genus-two quotient `C->B`. THM-4230's
saturated integral lattice gives

```text
deg_C(m)=4(N(alpha)+N(beta)).
```

The sharper generator relation kills the first visible eigendirection by the
unit `omega^2+omega=-1`, so in fact `deg_C(m)=4N(beta)`. Either form contradicts
the only possible response degrees `34` and `42`.

This closes all hidden-Hom loci on `U*Z*D*Lambda!=0`. It does not assert that
hidden factors are absent; it proves a genuine Keller specialization cannot
use them because it must preserve the deck character.

## The differential-character reconciliation

On the good special target,

```text
[-omega]^*(dX/(2Y))=-omega*dX/(2Y).
```

The original invariant form is

```text
dA/(2C_target)=sigma^2*dX/(2Y).
```

The character-two factor makes the form deck-invariant, but it vanishes on
the special fibre. Therefore the useful nonzero limiting differential is only
a semi-invariant eigensection on the cover. This explains simultaneously:

1. why the raw ambient-Kahler descent sought in THM-4288/4289 is not automatic;
2. why taking the character seriously gives a global exclusion anyway; and
3. why field trace of the nontrivial character is expected to vanish rather
   than recover the missing raw form.

## Exact wall boundary

The intrinsic quotient statement survives on `U*Z*D!=0`, including
`Lambda=0`: every actual deck-equivariant morphism `C->E_0` has degree
divisible by four. The Keller exclusion does not yet survive there.

THM-4272 proves a smooth raw genus-seven component and an `A_23` contact, but
does not prove either of the following:

- that the resolved Keller restriction to `C` has degree `34` or `42`; or
- that `C` is the sole degree-carrying component after the required regular
  specialization.

The genus ledger is a useful hostile: the raw union has
`p_a=7+0+12-1=18`, so resolving the `A_23` contact cannot discard eleven genus
units without recording them in cycles, multiplicities, or additional stable
components. Until that inventory and response invoice are proved, extending
THM-4290 across the wall would be a type error.

A concurrent `THM-4291` reservation makes this distinction especially sharp.
Its messages propose an additional genus-five tail that may admit an abstract
degree-42 equivariant map, followed by a correction proposing that the Keller
differential vanishes on that tail. The theorem file is still **RESERVED /
UNPROVED EMPTY STUB**, so neither assertion is used here. The genuine research
signal is the missing coordinate: wall work must track the differential on
each component as well as its deck character and abstract Hom degree.

## Prioritized continuation

1. **Lambda wall degree-and-differential inventory:** construct an equivariant
   regular graph model and compute which components carry the generic degree
   and a nonzero Keller differential. If `C` carries `34/42`, THM-4290's
   intrinsic wall lemma closes the wall without raw fat-contact descent.
2. **Single Kähler class:** on the terminal chart, express the Keller residue
   as `h*eta` and test `[h] in conductor/J_f`. This remains the cheapest local
   alternative and gives a decisive positive or hostile output.
3. **Remaining coefficient walls:** analyze `U=0`, `Z=0`, and `D=0` with the
   deck action retained before normalization. The likely runner vertices are
   components, characters, multiplicities, and degree obligations—not contact
   points alone.
4. **Seam entry:** only after the wall atlas is complete should the reduction
   claim be promoted toward `JC(2)`.

## Reproduction

```bash
python3 -B 04-computation/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.py
python3 -B -O 04-computation/jc23_a23_blowdown_kahler_dualizing_observer_thm4289.py
python3 -B 04-computation/jc23_exact_weight12_deck_equivariance_thm4290.py
python3 -B -O 04-computation/jc23_exact_weight12_deck_equivariance_thm4290.py
```

Normal, optimized, and fixed-hash-seed outputs match for both audits. The
finite scripts check ranks and characters; the formal local-algebra,
equivariant graph specialization, quotient factorization, and saturated Hom
arguments live in the theorem files.
