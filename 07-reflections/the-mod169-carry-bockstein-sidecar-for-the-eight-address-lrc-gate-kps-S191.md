# The mod-169 carry/Bockstein sidecar behind the eight-address LRC gate

**kps-S191 synthesis, 2026-08-21.**  Literal statuses are kept separate:
THM-2334, THM-3534, THM-3593, THM-3654, and THM-3657 are **PROVED** at
their stated scopes; the chain maps proposed below are **OPEN**.

**Exact update after THM-3658/3659.**  The proposed *one-bit* carry sidecar is
**REFUTED**.  The lawful linear split/cyclic character transform is now proved
in THM-3658, while THM-3659 computes the nonlinear defect exactly: on a carry
it is not one fixed vector but an address-dependent rank-two response with
55 values.  The strongest survivor is the full response table

```text
R(t0,t1)=e(t0,t1+1)-e(t0,t1),   0<=t0<=11, 0<=t1<=12. (0)
```

The discussion below is retained as the route that exposed this missing
coordinate; every occurrence of a bare carry-bit bridge is superseded by
the response-table formulation (0).

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

1. a **filtered** map retaining the carry and response table (0) as a
   sidecar;
2. a deliberately **nonlinear** map whose defect from additivity is exactly
   the address-dependent response in THM-3659; or
3. a proof that the physical current is confined to a carry-trivial chart.

The third option is hostile: the kernel and exceptional sets in `(2)--(3)`
meet both boundary and central carry regions.

## 4. The proposed typed bridge

The desired composition is not an equality of 169-element sets.  It is the
following proof obligation:

```text
THM-2334 exact relation current
  -> split target character q in G^                         [proved]
  -> carry/response-sensitive filtered lift               [OPEN]
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
| map | a filtered lift whose defect transports both the carry (10) and response (0) |
| preserved | reversal `(4)`, nonzero source correction, and the endpoint/chamber mark `(5)` |
| destroyed | absolute relation address beyond mod 13, unless supplied as a sidecar |
| cheapest hostile | prove the lift cannot be additive by `(8)` and test whether it drops rank after reduction mod 13 |

If this chain exists and a closed physical word supplies a relation of the
form used by THM-3657, the 169-state noncancellation problem becomes an
eight-address exposure problem.  Nothing currently proves that hypothesis.

## 5. Exact probes completed and the repaired frontier

### 5.1 Generic amplitude atlas

THM-3657 proves that 124 nonzero rows share the `+1` projective line of the
quotient reversal.  Write, in its pinned RREF coordinates,

```text
e(a)=mu(a)(1,b)                  for generic a.         (11)
```

Exact reversal gives `mu(168-a)=mu(a)`.  THM-3659 now tabulates this scalar.
It has exactly 16 values, but its interpolation degrees on the eight generic
intervals are

```text
(11,21,21,5,5,21,21,11),
```

and both full coordinate sequences have cyclic Fourier support `169/169`.
Thus the finite-state amplitude compression is real, while the proposed
low-degree and sparse-frequency compressions are **REFUTED**.  A useful law
must retain the 16-state partition or an equivalent address sidecar.

### 5.2 Carry-defect rank

For every pair of labels define

```text
Delta(x,y)=e(x plus y)-e(x)-e(y),                     (12)
```

where `plus` is cyclic addition (9), and compare it with the same formula
under split coordinatewise addition.  THM-3659 completes this computation.
The split and cyclic defect banks have exact records

```text
(rank, nonzero, distinct)=(2,26474,549), (2,26079,542).
```

Their difference vanishes without a carry and, with a carry, factors through
the split sum by (0).  The table has 156 cells, 139 nonzero cells, 55 values,
and rank two.  Consequently the one-bit factorization is **REFUTED**, but a
much smaller 156-cell factorization survives.

### 5.3 Exceptional exposure against target twists

Do not identify THM-2334 twist indices with `(r0,r1)`.  THM-3658 proves the
exact block Fourier change between the two character bases but also proves
that it is not a convolution or current map.  The repaired finite test is:

1. enumerate reversal-compatible associated-graded identifications;
2. transport the complete table (0), not just the carry bit;
3. ask whether a physically admissible coefficient law can place all
   response energy on generic reversal pairs and avoid `X`.

Arbitrary reversal-symmetric weights *can* avoid `X`; variance alone cannot
prove exposure.  The missing predicate must therefore come from chronology,
current conservation, positivity, or an exact-address sidecar.

### 5.4 Crisp next theorem target

Define the exceptional response energy for any candidate physical weight
`theta` by

```text
E_X(theta)=sum_(t in X or t+(0,1) in X)
             theta(t) ||R(t)||^2.                     (12a)
```

The norm here is only a placeholder until a characteristic-zero or
positive-semidefinite lift is typed.  The next useful theorem is not
`E_X>0` for arbitrary symmetric weights (false), but one of:

- an exact current identity forcing positive mass on response edges incident
  to `X`;
- a chronology restriction excluding every generic-only response cycle; or
- a dual certificate separating the physical coefficient cone from the
  generic-only kernel.

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

**PROVED:** the atlas (2)--(5), the quotient reversal, the source/two-row
projective disjointness, the eight-address necessary gate, the elementary
group obstruction (8), THM-3658's exact block Fourier transform, and
THM-3659's 55-value response-table factorization.

**REFUTED:** a constant one-bit carry response, low-degree generic interval
laws, cyclic Fourier sparsity, and exceptional exposure from arbitrary
reversal-symmetric variance alone.

**OPEN:** every current-to-digit chain map, every filtered response lift,
physical exposure of `X`, a useful law for the 16-state scalar automaton,
and LRC(14).
