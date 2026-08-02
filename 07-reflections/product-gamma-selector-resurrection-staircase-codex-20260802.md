# Product-Gamma selector resurrection staircase

**Status:** research reflection.  The depth-four/degree-12 cell is a proved
candidate in THM-3155.  The depth-five figures are finite-exact incoming work
and are not dependencies until separately scripted and audited.

## The first closed two-axis cell

For the support `(1,3)`, bank-`I2` selector problem, write

```text
C_D^(<=d)={laws on physical prefixes of depth <=d
           that are Hasse-positive in every degree 5,...,D}.
```

The exact proved-candidate pattern is

```text
C_11^(<=3)=empty,
C_11^(<=4) nonempty,
C_12^(<=4)=empty.
```

The positive corner has an explicit four-state law.  The two empty corners
have exact positive Farkas certificates.  This is stronger than a list of LP
outcomes: it is a closed birth/death cell in a monotone two-parameter family.

The monotonicities point in opposite directions:

- adjoining prefix states can only enlarge `C_D^(<=d)`;
- adjoining degree constraints can only shrink it.

Depth four crosses the first wall; degree 12 crosses the second.

## Why the cell is structurally small

The depth-three obstruction uses four upset rows.  The depth-four death adds
only the degree-12 singleton-complement row.  Every row factors through the
THM-3147 profile by partition length and singleton count.  Hence the active
dual carrier grows as

```text
four compressed rows  -->  the same four plus one rank tail.
```

The huge integer coefficients are coordinates of a ray in this small carrier,
not evidence that the obstruction intrinsically has huge dimension.  The
correct next algebraic object is the cone-valued pole-subtraction recurrence
on the endpoint jet.

## The next staircase step

A fresh exact depth-five law on seven states has denominator `10^9` and
numerators

```text
(1)             97185024
(1,1,1)         53307063
(1,1,1,1)      143473173
(1,1,2,2,2)    193810773
(1,2,2,3,3)     10572006
(2,2,2,3,3)    118925070
(5,5,6,7,8)    382726891.
```

They sum to `10^9` and have gcd one.  Exact enumeration reports strict
positivity on all `403,539` nontrivial upsets through degree 12.  Thus the next
positive corner is

```text
C_12^(<=5) nonempty.                                       (provisional)
```

That specific law fails in degree 13, but this does not yet prove
`C_13^(<=5)` empty.  The exact next dual search must include all physical
depth-five states, not only the seven support states of the positive law.

The depth-four five-row separator is crossed by exactly two depth-five states,
`(1,2,2,2,3)` and `(5,5,6,7,8)`.  As at depth three, the old obstruction dies
on a tiny new-state locus.  This is the finite-convex analogue of a bistellar
flip: the carrier remains almost valid, two cells change side, and a new
section appears.

## Possible global theorem and its obstruction

The physical pole multiset has only 16 letters.  Its complete submultiset
bank is finite:

```text
(4+1)(3+1)(2+1)^3(1+1)^3-1=4319 states.
```

Therefore this is not an infinite-dimensional escape problem.  If every
finite degree horizon were feasible on the full 4,319-state simplex, compactness
would produce one law satisfying every degree constraint.  Conversely, any
global failure has a finite Farkas certificate.

This makes two ambitious but decidable programmes visible:

1. continue the exact staircase until it stabilizes or dies before depth 16;
2. search directly for a structured law on the full state bank, using the
   finite staircase laws as moment constraints or warm starts.

Neither programme yet connects the averaged virtual-prefix current back to the
original product-Gamma response.  That realization map remains load-bearing.

## Product-law and stopping-law probes

The sparse laws found so far are arbitrary points of a simplex.  Two more
structured families should be tested before treating the coefficients as
accidental:

- **independent deletion:** each pole label is retained/deleted by a Bernoulli
  parameter, possibly depending only on pole value;
- **sequential hazard:** remove one remaining labelled pole at each step with
  a value-dependent transition weight and stop with a depth-dependent hazard.

For either family, the source is a small parameter space, the target is the
state-law simplex, and the preserved predicate is exact multiplicity validity.
The destroyed information is arbitrary correlation among deleted poles.  The
cheapest test is whether the depth-four and depth-five rational laws lie in the
closure of either family, or whether their low-order inclusion moments violate
the model identities.

If a structured family survives, it supplies the missing sequential sidecar.
If both fail, their first polynomial relation gives a useful no-go theorem.

## Updated concept board

1. **Two-axis selector barcode:** exact birth/death cells in `(D,d)`.
2. **Compressed dual jet:** every current active row uses only `(ell,m_1)`.
3. **Tiny crossing loci:** one depth-three crosser, two depth-four crossers,
   two depth-five crossers for the preceding dual.
4. **Finite full bank:** only 4,319 nonempty physical states exist.
5. **Structured stopping laws:** independent-deletion and hazard models are
   the cheapest route back toward a genuine process.
6. **Original-response map:** still absent; no selector theorem alone proves
   NC2 or GMC.

The creative lesson is to study the *motion of the dual face* and the *small
crossing locus* together.  Each alone looks unstable.  Their paired evolution
is rigid enough to suggest a recurrence, and finite enough to settle exactly.
