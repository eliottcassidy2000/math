# The allocation square becomes an allocation word along the address horn

**Status:** NON-CANONICAL SYNTHESIS around the `PROVED + VERIFIED-EXACT +
INDEPENDENTLY HOSTILE-AUDITED`
[THM-3285](../01-canon/theorems/THM-3285-same-ancestry-allocation-switching-horn.md).
The primary executable claim is frozen in
[`lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py`](../04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py)
and its matching
[`out`](../05-knowledge/results/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.out).
Its independent address-first audit is
[`lrc14_allocation_switching_same_ancestry_horn_independent_audit_20260803.py`](../04-computation/lrc14_allocation_switching_same_ancestry_horn_independent_audit_20260803.py).
No endpoint current, row exclusion, or LRC(14) conclusion is asserted here.

## Trigger and inheritance

Four nearby facts had looked mutually frustrating when read one address at a
time.

1. [THM-2791](../01-canon/theorems/THM-2791-full-arm-orbit-transfer-and-lower-central-chord.md)
   gives the two positive `tau=3` endpoint cylinders at
   `n_0=3454614` and `n_a=4143978`.  They lie on one literal Boolean rail
   ancestry sheet, but the theorem stops before endpoint allocation.
2. [THM-2807](../01-canon/theorems/THM-2807-positive-graded-address-two-simplex-and-allocation-lift-boundary.md)
   inserts `n_+=n_0+13` and obtains a positive address triangle from the
   exceptional `tau=12` target column.
3. The finite-exact addendum to
   [THM-2806](../01-canon/theorems/THM-2806-literal-fixed-sheet-central-allocation-scalar-law-and-endpoint-translation-no-go.md)
   proves that this `tau=12` triangle is target-only: its common carrier is
   already empty before the address cut.
4. The same addendum also leaves a distinct `9 x 9` label bank in which the
   common carrier `M` and right cofiber `R` are simultaneously nonempty but
   disjoint in every one of its `193` live cells.

The underused operation was therefore not another scalarization or Fourier
transform.  It was **restricting both disjoint allocation pieces along the
already known three-vertex address horn**.  The resulting object is not one
allocation square at one address.  It is an allocation-type word along an
address path.

## Exact reframe and finite census

Fix rail eight, physical clock one, both relative-present safeties, target
root one, delayed carry six, and the full THM-2806 common-label bank

```text
sigma in {0,1,2,3,8,9,10,11,12},
tau   in {3,4,5,6,7,8,9,10,11}.
```

Before address restriction let

```text
A = source one-sided carrier,
B = pulled target one-sided carrier,
M = A intersect B,
R = B minus A.
```

Push the `B`-side copies of `M,R` to the target coordinate and restrict them
at `(n_0,n_+,n_a)`.  Exhaustion of all `81` labels gives exactly three
patterns:

| labels | `M(n_0,n_+,n_a)` | `R(n_0,n_+,n_a)` |
|---|---:|---:|
| `sigma in {0,3,8,9,10,11,12}`, every `tau=3..11` | `010` | `101` |
| `sigma=1`, every `tau=3..11` | `000` | `001` |
| `sigma=2`, every `tau=3..11` | `000` | `101` |

Thus exactly

```text
1 clock x 7 sigma labels x 9 tau labels = 63 labels
```

carry the allocation-switching word

```text
R -> M -> R,
```

while the undecomposed target carrier has word `B=111`.  Every positive
vertex is one whole cylinder of weight

```text
27581135604,
```

with the identical exact mass and carry-six coefficient

```text
60781651775958960/371293,
790161473087466480.
```

The three exact target-coordinate translations are

```text
R(n_0) -> M(n_+) : 7/371293,
M(n_+) -> R(n_a) : 28553/28561,
R(n_0) -> R(n_a) : 371196/371293.
```

All three identities hold as translations of complete weighted intervals for
all `63` labels, not merely after integration.

## Why “same ancestry” is literal here

Numerically equal rail weights would not suffice after MISTAKE-281.  The scout
therefore rebuilds THM-2791's two raw contributor-wall sets.  All three open
cylinders lie strictly inside the same chamber

```text
[140890500190440,144190879112280).
```

The complete middle contributor set was enumerated independently.  Its two
factor counts are

```text
966606 and 28534,
```

their product is the rail weight `27581135604`, and their typed path-set digest
is the same THM-2791 endpoint digest

```text
15c804c7cea9f61feab3b641eccdc035d937142b446d1cc14e059210eb1534fd.
```

The explicit Boolean copy `(a,b,e')=(59162,26,56658)` is active at every
vertex, and the collision/root/deep labels remain `(5,6,12)`.  Because the
complete hull contains no contributor wall, the identity map—not an arbitrary
equinumerous bijection—identifies all literal contributor labels across the
three vertices.

This is the precise scope of “same ancestry”: the outer THM-2471 rail sheet is
fixed.  It does not say that every hidden semantic factor is covariant under
the translations.

## Connection contract

```text
source:
  the 81 clock-one label cells of the fixed rail-eight one-sided carrier bank,
  restricted to the three THM-2807 address cylinders;

target:
  a word in the disjoint allocation alphabet {M,R,empty} at those vertices;

map:
  B |-> (M=A intersect B, R=B minus A), followed by exact address restriction;

preserved:
  whole weighted cylinder, exact mass and delayed coefficient, rail ancestry,
  clock, sigma, tau, relative safeties, target root, carry, and open-boundary
  convention;

destroyed / not supplied:
  endpoint origin, endpoint current, the full bare/source/target/both K4,
  factorwise translation covariance, and a global packet action;

sidecar:
  the literal contributor-wall chamber and path-set digest;

decisive tests:
  two independently typed constructors agree on all 486 M/R address cells,
  and the three translation identities hold for every one of the 63 labels.
```

The mechanism is now transparent.  `B` is geometrically continuous along the
three selected address cylinders, but the source carrier `A` meets it only at
the middle.  The common/right allocation boundary is therefore crossed twice:
right cofiber, common overlap, right cofiber.  Scalarizing `B` erases exactly
the switching coordinate.

## Hostile controls and correction lineage

- The `sigma=1,2` rows are the nearest internal hostiles.  They show that the
  horn is not the whole `9 x 9` bank and locate the exact sigma boundary.
- `tau=12` remains outside the common target-label bank.  Its target-only
  simplex from THM-2806 is unchanged, so this result does not reverse the
  earlier no-go.
- MISTAKE-313 requires the physical clock to be part of the carrier.  The
  census is fixed at clock one and makes no seven-clock claim.
- MISTAKE-310 warns that a missing label need not mean empty physical support.
  Here emptiness and positivity are computed from literal half-open weighted
  intervals in the common target chart.
- MISTAKE-281 requires same-atom ancestry before composing controls.  The
  wall-chamber audit supplies it through `(u,v,a,b,e')`, but deliberately stops
  before endpoint atoms.
- MISTAKE-300 forbids identifying two selectors by shared vocabulary.  The
  direct constructor forms `A`, `B`, `M`, and `R` from their truth sets; the
  canonical constructor is only a second exact replay.
- The independent audit reverses the construction order and also tests three
  stronger interpretations.  Translation of the uncut full `R` carrier to
  the full `M` carrier fails, diagonal periodicity of the uncut full `R`
  carrier fails, and the complete fixed-clock factor section is not covariant
  under the first horn step.  Thus the theorem is exactly address-local and
  cannot be read as a global allocation or packet action.

## What this changes on the live board

- **THM-2791 rail sheet:** its remote endpoint chord is not isolated.  The
  same literal ancestry chamber contains a genuine common-carrier middle
  vertex.
- **THM-2806 allocation obstruction:** fixed-address allocation remains
  projectively flat, but allocation type is not constant along the address
  horn.  The useful invariant is a word/cocycle, not one scalar square.
- **THM-2807 address simplex:** the address triangle can be realized inside
  the genuine common-label bank only after allowing its allocation role to
  switch; it is not an all-common simplex.
- **Endpoint frontier:** the missing coordinate is narrower.  The next map no
  longer has to create a common physical cylinder from nothing.  It must attach
  endpoint-origin/allocation data to an existing `R-M-R` same-ancestry horn and
  test whether the two switches cancel or retain the `Z2^4079` transition.

## Cheapest next decisive test

Attach one explicit THM-2772 endpoint-origin fibre to a single horn label,
preferably `(clock,sigma,tau)=(1,0,3)`, and compute the four literal
bare/source/target/both states separately at all three vertices.  The test must
retain the middle `M` atom, both endpoint `R` atoms, the endpoint origin, and
the full-depth transition.  There are three informative outcomes:

1. the endpoint fibre is empty at the middle: an endpoint-origin obstruction;
2. the fibre exists but one translation fails: genuine allocation holonomy;
3. all three lift and the two allocation switches compose: a new physical
   endpoint cycle, still requiring a current/Fourier noncancellation audit.

The stopping certificate for this scout is the complete `81`-label pattern
census and its two sigma hostiles.  Another scan of the same three addresses
without adding endpoint-origin data would only reproduce the classified
allocation words.

## Reproduction and hashes

```bash
python3 04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py
python3 -O 04-computation/lrc14_allocation_switching_same_ancestry_horn_scout_20260803.py
```

Both modes byte-match the stored output.  LF-normalized hashes at freezing are

```text
c42d66498f460f2142ea375fe9d4047b82c62b872b35d5a1634d2bb4c80a68ee  script
e89dce3307e5d374e8583f92e1b2da1214e44929e52fdd42c6532d61adb3e246  output
```

The script has zero Python `assert` nodes and zero floating-point literals.
