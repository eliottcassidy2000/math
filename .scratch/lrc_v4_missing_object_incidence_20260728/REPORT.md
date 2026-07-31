# Rank-587 missing-object incidence audit

**Status:** `FINITE-EXACT` scratch result.  This note makes no tracked canon
change and proves no LRC(14) row exclusion.  Its scope is the complete
THM-2825 cofiber-rooted path bank used by THM-2859.

## Verdict

The fourth `V_4` corner is not a previously overlooked collar vertex and is
not an integral signed combination of the 587 whole rooted paths.

There is a useful but lossy virtual identity

```text
X = R + M1 - M2.
```

It is the unique identity after retaining only dimension and the separate
semantic and carrier sign marginals.  It fails the mixed
semantic-by-carrier character by `4` at every root, hence by `2348` on the
rank-587 bank.  It is therefore not an object, a positive packet, or a typed
`V_4` globalization.

The smallest faithful description is

```text
X = Pi(R)
  = (R carrier and factor data)
      x_(root index)
    (M1 semantic data).
```

Abstractly this is one new joint parity tag.  It is not one uniform physical
mask: it refines into four factor-hole/semantic types, 42
physical-interval/factor-hole/semantic types, and 587 labelled atoms.
THM-2859's lower bound `dim X >= rank S = 587` is sharp.

## 1. Exact universe and inheritance

The audit reconstructs every nonempty THM-2825 cell and all 587
cofiber-rooted half-step paths:

```text
rooted paths                         587
common vertices on rooted paths  54,754
labelled vertices                 55,341
physical interval/weight rows      1,192
```

Every path has at least two common vertices after its right root.  In a
labelled cell, the forest paths are disjoint.

The right-root types, in `(factor holes, semantic live)` notation, are

```text
(E3,       live) 319
(E3+c2,    live)  37
(c2,       live) 217
(q1,       dead)  14.
```

The required parity-shifted copy has the same carrier and factor data with
the semantic bit reversed:

```text
(E3,       dead) 319
(E3+c2,    dead)  37
(c2,       dead) 217
(q1,       live)  14.
```

The two sets of decorated root types are disjoint.  Erasing labels leaves
37 physical root intervals and 42
`(physical interval, factor holes, semantic)` classes; no physical root
has conflicting semantics and none has the desired decorated type.

Nor can a later common vertex be `X`.  By construction

```text
common = source intersection pulled(target),
right  = pulled(target) minus common.
```

Thus every common vertex has source carrier present at twist zero, whereas
`X`, like `R`, needs the source-empty carrier.  Roots have the right carrier
but the wrong semantic bit; common vertices can have the opposite semantic
on alternating levels but have the wrong carrier.  This rules out the
entire existing rooted path, not only `M1` and `M2`.

## 2. Labelled whole-path incidence no-go

Index the rooted paths by `i=1,...,587`.  Let `A` be their labelled vertex
incidence matrix: column `i` is the characteristic vector of the whole path
from `R_i` through its common tail.  Let

```text
r = sum_i e_(R_i)
```

be the support selector of the required cofiber copy.  In

```text
A z = r,                                                   (1)
```

the unique row `R_i` forces `z_i=1`, while the unique first-common row
`M1_i` forces `z_i=0`.  This gives 587 independent local contradictions.
Equivalently,

```text
rank A = 587,       rank [A|r] = 588.                      (2)
```

The first exact witness is

```text
path       0
cell       (clock,sigma,target) = (1,0,3)
root index 0
R          (142004190428100,142004216872980,27581135604)
M1         (142004591508780,142004617953660,27581135604).
```

The same proof works for a prescribed signed/alternating path vector so
long as both its root and first-common coefficients are nonzero.  Permitting
a zero common tail is already permitting the missing root-supported
object.

## 3. Label-forgotten integral certificate

One might hope that physical interval collisions allow signed cancellation
after cell labels are forgotten.  They do not.

Let `A_phys` have the same 587 columns but rows indexed only by physical
interval/weight triples.  The following two rows have exactly the same seven
nonzero columns:

```text
u=(145418991337620,145419017782500,27580222516)
v=(145419392418300,145419418863180,27580222516)=u+h.
```

Their common column set is

```text
{265,269,273,277,281,285,289},
```

the root-index-zero paths in cells

```text
(clock,sigma,target)=(2,8,t),  t=3,4,5,6,7,8,9.
```

In every one of those paths, `u` is position zero and `v` is position one.
Consequently the two left-hand sides of the physical system are identical.
The labelled-multiplicity root target asks them to equal `7` and `0`:

```text
sum_(i in I) z_i = 7,
sum_(i in I) z_i = 0.                                    (3)
```

This is a direct integral contradiction.  If the physical target is
Boolean support rather than labelled multiplicity, the same rows demand
`1` and `0`.  Reducing the multiplicity system modulo two gives

```text
rank A_phys = 59,       rank [A_phys|r_phys] = 60,
```

with the two rows above as a complete inconsistency certificate.  Its
ordered-row digest is

```text
bcf9e54eb147b9cf33be8011ee110f2039ef3cdc035a387f6a1b5c17a5f34b79.
```

This witness also identifies information destroyed by the physical
quotient: `u` and `v` have identical path incidence but different
cofiber depth.  Any physical descent needs a root/common-depth sidecar; the
bare interval row cannot remember which member is the cofiber object.

## 4. Why the virtual signed formula almost works

Use relative bits `(semantic,carrier)`:

```text
R  =00,     M1=11,     M2=01,     X=10.
```

Forget the joint interaction and retain only

```text
q(s,c)=(1,(-1)^s,(-1)^c).
```

The three old columns have determinant `4`, so the coefficients expressing
`q(X)` are unique over characteristic zero:

```text
q(X)=q(R)+q(M1)-q(M2).                                   (4)
```

The omitted mixed character is

```text
chi_sc(s,c)=(-1)^(s+c).
```

On the right of `(4)` it has value `3`, while `chi_sc(X)=-1`.  The defect
is `4` per root.  Restoring this fourth coordinate makes the four-corner
Hadamard determinant `-16`: the old three corners have joint rank three and
adjoining `X` raises it to four.

Thus `(4)` is a genuine Grothendieck-style shadow in a quotient that forgets
semantic/carrier correlation.  Its negative coefficient and mixed-character
defect explain exactly why it cannot be promoted to a positive or faithfully
typed collar object.

## 5. Cheapest faithful enlargement

There are two equivalent invoices.

1. Adjoin the parity-tagged cofiber copy `X=Pi(R)`, retaining every root's
   carrier and factor annotations and importing the `M1` semantic bit over
   the same root index.
2. At bare incidence level, adjoin one common-tail object per root; then
   `whole path - common tail = root`.  These 587 independent tails merely
   repackage the same rank-587 filler and still need the parity tag for the
   typed problem.

No enlargement of smaller dimension can factor the rank-587 map `S`.
Scalar cocycles or signs only change coefficients on existing support; the
missing coordinate is the joint semantic/carrier bit itself.

## 6. Reproduction

```bash
python3 .scratch/lrc_v4_missing_object_incidence_20260728/missing_object_incidence.py
python3 -O .scratch/lrc_v4_missing_object_incidence_20260728/missing_object_incidence.py
```

The script pins the THM-2825 primary and path-operator scripts/outputs and
the THM-2859 globalization script/output.  It contains no executable Python
`assert`.  Normal and optimized modes produce byte-identical output.

