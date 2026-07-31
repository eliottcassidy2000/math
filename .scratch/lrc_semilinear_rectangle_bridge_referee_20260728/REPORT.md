# Independent referee report: semilinear endpoint rectangles

**Verdict:** `PASS WITH CORRECTIONS APPLIED` for promotion as a finite exact
support theorem.  The six masks, all affine-group counts and searches,
nearest-map defects, Fourier covector transport, named-origin census, and
cyclotomic/Prony counts survive an independent replay.  The scratch report's
stale THM-2868 interpretation and overly broad PGL claim were repaired in
place.  No LRC(14) row is excluded.

## 1. Independent evidence

The referee companion pins the submitted script/report/output, then
reconstructs the six q0/q3/q11 step-2/step-68 masks with the separate
THM-2859 endpoint-coboundary evaluator.  The two evaluators agree on all
`6*169` address values.  Direct evaluation after the physical source/target
shift also confirms source-frame equals target-frame for all six states.
Every support is the claimed Cartesian rectangle.

The independent affine enumeration makes one pass over

```text
|GL_2(F_13)|  = 26,208,
|AGL_2(F_13)| = 4,429,152,
monomial GL_2 = 288.
```

Its `(maximum overlap, minimum Hamming defect, exact-map count)` table is

```text
q0 step2 -> q3 step2       (81,  9, 0)
q0 step2 -> q11 step2      (64, 34, 0)
q0 step2 -> q0 step68      (81,  0, 8)
q3 step2 -> q3 step68      (80, 20, 0)
q11 step2 -> q11 step68    (72, 18, 0)

one map on both q0 states -> both q3 states   (144,54,0)
one map on both q0 states -> both q11 states (128,68,0).
```

The typed subgroup counts also replay:

```text
exact character-3 linear stabilizer 156
axis-preserving linear maps         144
axis and character-3                 12
pointwise two-origin affine maps    156
both character-3 and two-origin      13.
```

Normal and optimized runs of the submitted companion byte-match its stored
output, SHA-256

```text
66c38cfd4c68f7c4ff3471d84ccad33c7329a142d5970c358d0d28ca07bc0ad3.
```

## 2. Mathematical checks

For each source rectangle, every matrix row using both input coordinates has
projection size thirteen.  This was replayed directly for all nonzero row
coefficients, and follows abstractly from Cauchy--Davenport:

```text
|rA+sB| >= min(13,|A|+|B|-1)=13.
```

Since each target coordinate factor has size nine or ten, every **exact**
affine rectangle map is monomial.  This proves the exact-map reduction; the
full AGL enumeration is still needed for the global nearest-map claims,
because a nonexact mixed map need not be monomial.

The selected nearest defects replay literally:

```text
q0 -> q3:
  target minus image = {9} x B3;

q3 edge under (a,b)->(a,2b):
  lost = A3 x {3},       gained = A3 x {12};

q11 edge under (a,b)->(a,6b+1):
  lost = A11 x {10},     gained = A11 x {11}.
```

For an affine state map `y=Mx+t`, pushforward sends a Fourier row covector
`k` to `kM^(-1)`; translation contributes only a scalar phase.  Thus the
q3 and q11 near maps send `(0,3)` to `(0,8)` and `(0,7)`, respectively.
Coefficient recharts `omega -> omega^2` and `omega -> omega^6` return those
characters to `3`.

The q0 edge has exactly eight affine maps.  Exactly two preserve the named
character `(0,3)`.  Their images of the ordered origins
`((0,0),(12,0))` are

```text
((0,5),(12,5)),       ((6,5),(7,5)).
```

No exact q0 edge map preserves the origins pointwise or setwise.

Finally, with `xi` of order `2366`, the Prony nodes `xi^13` and `xi^169`
have orders `182` and `14`.  For each rechart slope `s=2,6`:

```text
units u mod2366 with u mod13=s          78
such units fixing the order-182 node     0
such units fixing the order-14 node     13
such units preserving the unordered pair 0.
```

The unordered-pair conclusion also follows because a cyclotomic
automorphism preserves multiplicative order and hence cannot swap nodes of
orders `182` and `14`.

## 3. Required PGL scope correction

The exhaustive projective-line group has

```text
|PGL_2(F_13)|=2,184.
```

Requiring a Möbius denominator to be nonzero on the occupied source factor
is correct: otherwise the map leaves the affine endpoint plane.

No choice of a separate Möbius map on each coordinate gives a full
q0-to-q11 rectangle bridge at either step, or a full q3 or q11 direct edge,
in either direct or swapped orientation.  That is the promotable statement.

The original broader claim that coordinatewise PGL adds no exact factor
maps was false.  Non-affine individual factor maps exist:

```text
q3_A  -> q3_A:       8 maps, 6 non-affine;
q11_A -> q11_A:      6 maps, 4 non-affine;
q11_A -> q11_B2:     6 maps, 4 non-affine;
q11_B2 -> q11_A:     6 maps, 4 non-affine;
q11_B2 -> q11_B2:    6 maps, 4 non-affine.
```

None supplies the missing companion factor, so the full-rectangle no-go is
unchanged.  The source scratch report now states this narrower result.
The theorem must call the operation “coordinatewise
`PGL_2(F_13) x PGL_2(F_13)` Möbius maps”; it does not exhaust arbitrary
birational plane maps or `PGL_3`.

## 4. Correct relation to THM-2868 and THM-2874

THM-2868 is no longer reserved.  It is
`PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED`, and constructs a
signed two-origin coefficient projector, a 26-sample variable-offset Prony
atlas, and a projective Kummer coordinate.  It never claims that one
positive full endpoint mask is affinely transported while fixing both
origins.

Accordingly, the present affine no-go is a support-level boundary
**compatible with** THM-2868, not a warning against a proposed theorem.
It explains why THM-2868's successful route must live on signed coefficient
ratios rather than positive coordinatewise support bijections.

THM-2874 then proves the explicit `Q(zeta_91)`-rational clutch from that
projective frequency orbit to the centered endpoint Galois orbit, realizes
the sharp `C169` carry fibre over q3, and proves the Bockstein seam
transgression.  The hypothetical coefficient recharts `2` and `6` in this
audit address only the specific q3/q11 **near support maps**.  They are not
needed for THM-2868/2874, and the Prony-node calculation proves that the
existing normalized THM-2863 bank cannot simply be imported unchanged
under those recharts.

The proposed one-fibre correction is therefore a test for a new positive
common-support extension over q11/q7.  Failure would not retract or
reformulate the proved projective coefficient atlas.

## 5. Does the rectangle no-go explain the missing guard/q5 corner?

No lawful identification is currently available, and the natural quotient
gives an exact no-go.

The addendum's literal target-atom projection to `(guard,q5)` has image

```text
{(S,S),(D,S),(D,D)}
```

and misses `(S,D)`.  The referee independently reconstructs the complete
13-residue word at both step 2 and step 68; the words agree.  Moreover
endpoint-address representatives have coordinates zero in both the guard
and q5 slots, so this projection is constant over all 169 endpoint
addresses.

In particular, q0, q3, and q11 at both audited steps all map to `(S,S)`.
The natural map is

```text
pi(q,step,address)
  = (literal guard truth at q, literal q5 truth at q).
```

It destroys the endpoint address and sends both ends of each q3/q11
rectangle edge to the same chamber.  A minimal witness is the nonzero q3
boundary

```text
A3 tensor (delta_12-delta_3).
```

Under `pi` both fibres land at `(S,S)`, so this boundary maps to zero.  The
same holds for the q11 boundary
`A11 tensor (delta_11-delta_10)`.  Conversely, the absent `(S,D)` chamber
is not in the image at all.  Thus neither the one-row defects nor the full
rectangle no-go explains the missing chamber under the only intrinsic
factor-state quotient presently defined.

Any claimed connection must supply a different explicit map and show that
it preserves q-label, physical origin, Boolean positivity or augmentation,
and the guard/q5 predicates.  Factorwise PGL absence alone cannot provide
such a map.

## Promotion recommendation

Promote the finite affine theorem after carrying over the corrected report
wording and an independent audit companion.  State:

- exact mask bank and source/target equality;
- complete AGL counts, exact-map absences, and nearest defects;
- character/origin and cyclotomic-node boundaries;
- no **full coordinatewise Möbius rectangle bridge**;
- compatibility with the already proved THM-2868/2874 signed coefficient
  route; and
- no q7 transport, physical ancestry action, row exclusion, or LRC(14)
  decrement.

