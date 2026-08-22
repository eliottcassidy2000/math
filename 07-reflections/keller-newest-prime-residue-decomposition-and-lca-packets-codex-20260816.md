# The strict-henselian factor atlas forgets its own symmetry

**Status: proved synthesis from THM-3530, THM-3533, THM-3535, THM-3538,
and THM-3539.  The actual newest-prime decomposition-group saturation and
all-level coordinate index zero remain OPEN.**

## 1. Inheritance pass

The closest proved mechanism is THM-3535's supported bottom transposition in
the full ternary wreath group.  The canonical hostile is THM-3537: at old
`L`, the predecessor is already ramified and two Newton packets collide, so
the newest-prime split-block model is unavailable.  The corrected near miss
is “full monodromy plus one inertia element determines local decomposition”;
it does not.  The least-used sidecar is the residue extension in the
decomposition exact sequence.

THM-3538 writes the prescribed-coordinate index as one internal factor for
each predecessor block plus one resultant for each unordered block pair.  At
depth `n`, that is

```text
m+binom(m,2)=m(m+1)/2,             m=3^(n-1),           (1)
```

apparently an exponential wall.  THM-3539 identifies both the maximal lawful
symmetry reduction and why it is not yet available.

## 2. The exact local sandwich

For the supported bottom transposition `t`, the full wreath calculation is

```text
N_W(<t>)=C_W(t)=Stab_W(ell_*),                         (2)
```

where `ell_*` is the third sibling in the ramified ternary block.  The Keller
decomposition group satisfies only

```text
<t>=I <= D <= C_W(t).                                  (3)
```

The quotient `D/I` is the residue Galois group.  Passing to the strict
henselization makes the predecessor cover split precisely by separably
closing that residue field.  Thus the construction that reveals every cubic
block simultaneously removes the descent datum saying which blocks remain
conjugate over the original newest divisor.

This is not a defect in THM-3538.  It is the exact price of its useful factor
atlas.  The missing object is now concrete: the residue decomposition cover,
not another global discriminant or primitivity calculation.

## 3. The quadratic LCA shadow

Let `r=n-1`, mark the ancestry block `b_*`, and let `H_r` be its stabilizer in
the predecessor wreath group.  On blocks, `H_r` remembers only the LCA depth
from `b_*`, giving `r+1=n` types.  On unordered pairs it remembers the rooted
three-leaf dendrogram:

- a pair containing `b_*` and one shell depth;
- two unequal shell depths; or
- one common shell depth plus the pair's mutual LCA depth.

There are exactly

```text
r(r+1)=n(n-1)                                          (4)
```

pair types.  Internal factors and resultants therefore total exactly

```text
n+n(n-1)=n^2                                           (5)
```

maximal-symmetry packets.

This is the useful compression.  It is also a floor: any actual decomposition
group is a subgroup of the stabilizer, so its orbits can only split these LCA
packets.  With no residue action beyond inertia, all factors in `(1)` remain
separate.

## 4. The right gate is weaker than `D=C`

The THM-3538 factors depend on predecessor blocks, not on permutations of the
three last-step roots inside each block.  Consequently full equality in `(2)`
is stronger than needed.  It is enough that the image of `D` on predecessor
blocks have the same point and unordered-pair orbits as `H_r`.

This “two-orbit saturation” is the next exact target.  Three possible attacks
are now sharply separated:

1. compute the residue splitting field of the predecessor cover along
   `P_(n-1)` and prove its group is `H_r`;
2. prove only point-and-pair orbit saturation, avoiding the invisible bottom
   base factors; or
3. find a residue factorization whose orbit splitting refutes saturation and
   identifies the extra arithmetic packet coordinate.

Any of the three is progress.  None may be replaced by the global identity
`G_n=W_n`.

## 5. What the packet quotient keeps and destroys

The LCA quotient keeps the marked ancestry, rooted distance shells, the
three-leaf dendrogram of a pair, orbit multiplicities, and--under the
decomposition gate--valuation/unit status.  It destroys the ternary address,
off-ray child label, strict-henselian section label, residue-field descent,
and the actual unit residue of each factor.

Therefore `(5)` does not prove the factors are units.  It converts a potential
exponential list into a quadratic list only after a local Galois theorem, and
then leaves `n^2` genuine arithmetic noncollision tests.  Producing one deep
resultant can still be expensive, so this is an orbit-type bound rather than
an algorithmic running-time theorem.

## 6. Updated concept board

| lane | new object | preserved | destroyed / still missing |
|---|---|---|---|
| Keller newest prime | residue decomposition group `D/I` | local conjugacy packets | strict henselization erases it |
| ternary tree | marked-leaf LCA dendrogram | rooted ancestry type | exact address |
| THM-3538 carrier | `n^2` conditional packets | valuation and unitness under the gate | representative arithmetic |
| old-`L` inertia | ramified predecessor hostile | need for cross-block resultants | newest split model |
| tournament/XOR lane | no intrinsic pair orientation appears | unordered pair relation only | orientation gauge and edge signs |

The last row is a useful negative connection.  The carrier indices form an
unordered-pair association scheme around a marked leaf, not a tournament:
there is no intrinsic binary orientation, and forcing one would discard the
LCA invariant rather than explain it.

## 7. Next decisive tests

The cheapest lawful next probe is not another level-five resultant expansion.
It is a residue factorization of the predecessor normal closure along one
newest divisor, recording decomposition-group orbits on blocks and pairs.
At finite levels `n=2,3,4`, THM-3538's existing witnesses can be revisited to
distinguish three statements that were previously conflated:

```text
carrier unitness;
point-and-pair orbit saturation;
full centralizer equality.                             (6)
```

Only the first is currently proved there.  Establishing or refuting the
second is the smallest test that changes the all-level workload.  Full
centralizer equality is optional extra structure.

No step here proves all-level index zero, a Keller-map classification, a
Jacobian-conjecture counterexample theorem, a physical current, or LRC(14).
