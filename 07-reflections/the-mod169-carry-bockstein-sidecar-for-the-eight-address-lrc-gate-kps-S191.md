# The mod-169 carry/Bockstein sidecar behind the eight-address LRC gate

**kps-S191 synthesis, 2026-08-21.**  Literal statuses are kept separate:
THM-2334, THM-3534, THM-3593, THM-3654, and THM-3657 are **PROVED** at
their stated scopes; the chain maps proposed below are **OPEN**.

## 1. Inheritance pass

- **Closest proved mechanism.**  THM-3657 classifies the 169 two-current
  digit labels in THM-3593's rank-two correction quotient and proves that
  every source-defect cancellation must use one of eight exceptional labels.
- **Canonical hostile.**  THM-3585's equality of row spaces is static.  Its
  parent explicitly says that the two retained digits are not a THM-2334
  exact relation address or a physical current.  Equal dimensions and equal
  cardinalities do not repair that type mismatch.
- **Corrected near miss.**  THM-2334 already warns that a common abstract
  `F_13^2` label does not intertwine its 169 target twists with predecessor
  dynamics.  The new digit arithmetic identifies one concrete reason.
- **Least-used sidecar.**  The base-13 carry cocycle, equivalently the
  nonsplit extension `0 -> F_13 -> Z/169Z -> F_13 -> 0`, has not yet been
  retained in the current-to-address bridge.

## 2. The exact integer-coordinate atlas

The two-current parent assembles its retained inverse digits as

```text
a=r0+13 r1,                  0<=a<=168.                (1)
```

THM-3657's proved label classification becomes

```text
kernel Z = [0,11] union [78,90] union [157,168],       (2)

exception X = {12,25,48,71,97,120,143,156}
            = 84 plus_or_minus {13,36,59,72}.          (3)
```

Simultaneous branch reversal is

```text
(r0,r1)->(12-r0,12-r1),
a       ->168-a=-a-1 mod 169,                          (4)
```

with unique fixed residue `84`.  Thus `(2)` is a central block and two
boundary blocks, while `(3)` is four reflection shells.  The `37+124+8`
atlas is visibly order/carry geometry, not a translation-invariant partition
of a two-dimensional vector space.

The eight addresses also have a second proved provenance:

```text
six labels = THM-3534 endpoint residual,
two labels = (3,9),(9,3), the crossed Loc_(3,9) representatives. (5)
```

This is the first object in this route that simultaneously retains the
endpoint residual and the marked two-chamber location.  THM-3534 had to
quotient the endpoint line before its relative cospan closed; THM-3657 says
that a source correction cannot cancel unless the address data returns to
exactly that endpoint/chamber locus.

## 3. Why the two different 169s do not identify

THM-2334's target quotient is

```text
G=K_13/L_13 isomorphic to F_13^2,       |G^|=169.      (6)
```

The current-leg digits in `(1)`, when chronology composes them, form

```text
C_169=Z/169Z.                                          (7)
```

These groups have the same cardinality and are not isomorphic:

```text
exponent(C_169)=169,              exponent(F_13^2)=13. (8)
```

Consequently every additive map `C_169 -> F_13^2` kills `13 C_169` and has
image dimension at most one.  In particular, no additive bijection can turn
the two-digit ancestry object into the full two-dimensional target quotient.
This upgrades the old warning "the labels are merely both `F_13^2`" to an
explicit obstruction.

In digit representatives, addition in `(7)` is

```text
(r0,r1) plus (s0,s1)
 = (r0+s0 mod 13,
    r1+s1+kappa(r0,s0) mod 13),                       (9)

kappa(r0,s0)=1 if r0+s0>=13, and 0 otherwise.         (10)
```

The carry `kappa` is the missing two-cocycle.  Forgetting it replaces the
nonsplit cyclic extension by its split associated graded `F_13 direct_sum
F_13`.  A lawful bridge must therefore be one of the following:

1. a **filtered** map retaining the carry/Bockstein as a sidecar;
2. a deliberately **nonlinear** map whose defect from additivity is exactly
   `(10)`; or
3. a proof that the physical current is confined to a carry-trivial chart.

The third option is hostile: the kernel and exceptional sets in `(2)--(3)`
meet both boundary and central carry regions.

## 4. The proposed typed bridge

The desired composition is not an equality of 169-element sets.  It is the
following proof obligation:

```text
THM-2334 exact relation current
  -> split target character q in G^                         [proved]
  -> carry-sensitive lift (q,beta(q)) in a filtered C_169  [OPEN]
  -> two-current row a_(r0,r1) in THM-3585                 [OPEN]
  -> correction e(a) in E2                                 [proved]
  -> exceptional exposure X                                [THM-3657]
  -> endpoint/chamber response sidecar                      [THM-3534]
  -> fixed-branch rigidity test                             [THM-3654].
```

For any proposed arrow one must record:

| field | required answer |
|---|---|
| source | THM-2334 quotient character together with its exact-address orbit coefficient |
| target | the assembled mod-169 inverse-digit label, not a free pair of residues |
| map | a filtered lift whose additivity defect is the carry cocycle `(10)` |
| preserved | reversal `(4)`, nonzero source correction, and the endpoint/chamber mark `(5)` |
| destroyed | absolute relation address beyond mod 13, unless supplied as a sidecar |
| cheapest hostile | prove the lift cannot be additive by `(8)` and test whether it drops rank after reduction mod 13 |

If this chain exists and a closed physical word supplies a relation of the
form used by THM-3657, the 169-state noncancellation problem becomes an
eight-address exposure problem.  Nothing currently proves that hypothesis.

## 5. Concrete next computations

### 5.1 Generic amplitude atlas

THM-3657 proves that 124 nonzero rows share the `+1` projective line of the
quotient reversal.  Write, in its pinned RREF coordinates,

```text
e(a)=mu(a)(1,b)                  for generic a.         (11)
```

Exact reversal gives `mu(168-a)=mu(a)`.  The next probe should tabulate
`mu(a)` and test, with hostile controls, whether it is:

- constant on the carry-free intervals;
- a low-degree polynomial separately on those intervals;
- a discrete Green function with jumps exactly at `(2)--(3)`; or
- genuinely high-complexity.

A positive result would turn the projective gate into a scalar transport
law.  A negative result would identify the missing coordinate that the line
atlas suppresses.

### 5.2 Carry-defect rank

For every pair of labels define

```text
Delta(x,y)=e(x plus y)-e(x)-e(y),                     (12)
```

where `plus` is cyclic addition `(9)`, and compare it with the same formula
under split coordinatewise addition.  Compute the ranks and supports of the
two defect banks.  The decisive positive signal is that the cyclic defect
factors through the one-bit carry `kappa`; the decisive negative signal is a
rank larger than the carry sidecar can support.

### 5.3 Exceptional exposure against target twists

Do not identify THM-2334 twist indices with `(r0,r1)`.  Instead enumerate
all affine/linear identifications of the associated graded quotients that
respect reversal, then ask whether the nonconstant twist variance can be
supported away from the image of `X`.  This is a finite hostile test for any
candidate filtered lift, not evidence that such a lift is physical.

## 6. Cross-thread lesson from AMM 12592

The AMM dyadic frontier independently rejected every scalar recurrence for
the Rule-A horizon.  Its exact surviving state is a skew product: Sturmian
phase, normalized offset, clamp/front state, and a carry-like branch bit.
The analogy is precise at the level of research method:

```text
discard carry/phase -> attractive false scalar recurrence,
retain extension state -> exact but higher-dimensional transport. (13)
```

This does not transfer an AMM theorem to LRC.  It does identify the same
failure mode that repeatedly produced false simplifications in both lanes.
The next LRC move should preserve the mod-169 extension state from the start.

## 7. Status

**PROVED:** the atlas `(2)--(5)`, the quotient reversal, the source/two-row
projective disjointness, the eight-address necessary gate, and the elementary
group obstruction `(8)`.

**OPEN:** every current-to-digit chain map, every filtered/Bockstein lift,
physical exposure of `X`, scalar amplitude law `(11)`, and LRC(14).
