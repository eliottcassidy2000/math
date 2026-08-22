# The `2 x 2` mixed-Haar kernel is a marked square cycle in `K4`

**Status: VERIFIED-EXACT representation sidecar; not independently audited.**
The companion is
`04-computation/tetrahedral_k4_haar_xor_cycle_bridge_20260816.py`.
It refines the tetrahedral representation atlas and the canonical U_full
minimum-joint-coordinate gate.  It does not construct the missing map from
actual endpoint cells to ancestry addresses.

## Four atoms and six pairwise relations

Label the four cells of a `2 x 2` coupling by

```text
00, 01, 10, 11 in F_2^2.                                  (1)
```

These are the four vertices of `K4`; its six edges are all pairwise
relations.  The three nontrivial Walsh characters are the three balanced
sign vectors

```text
chi_01=(+,-,+,-),
chi_10=(+,+,-,-),
chi_11=(+,-,-,+).                                         (2)
```

For each character, equal-sign pairs form one of the three perfect
matchings.  The remaining four edges form a square.  Thus the three objects
that had repeatedly appeared separately are different representations of
one tetrahedron:

```text
four joint atoms  -> three 2+2 Haar lines
                  -> three perfect matchings
                  -> three complementary four-cycles in six edge slots. (3)
```

## The opposite-face map turns Haar into holonomy

Let `Omega:C_0(K4)/1 -> H_1(K4)` be the oriented opposite-face map from the
tetrahedral atlas.  It satisfies

```text
Omega^T Omega = 4I-J.                                     (4)
```

Hence it is a factor-two isometry on the sum-zero vertex space.  With edges
ordered `(01,02,03,12,13,23)`, the actual U_full checkerboard direction is

```text
chi_11=(1,-1,-1,1),

Omega(chi_11)=2*(-1,1,0,0,-1,1).                         (5)
```

The primitive vector in (5) is exactly the oriented boundary

```text
00 -> 10 -> 11 -> 01 -> 00.                               (6)
```

Its two omitted diagonals are the equal-parity matching.  The other two
Walsh characters give the other perfect matching/complementary-square pairs.
The three primitive cycles are pairwise orthogonal and form a rational basis
of `H_1(K4)`.

This makes the U_full minimum theorem concrete.  Row/column marginalization
on four cells has rank three and kernel `span(chi_11)` over every odd
characteristic.  The one mixed-Haar number needed to repair that quotient is,
under `Omega`, one marked square holonomy.  This is why a four-atom problem
naturally asks for six pairwise edge slots, but only one of the three cycle
coordinates is selected by a fixed left/right factorization.

## `S4`, XOR, and tournaments

`S4` permutes the three Haar lines/perfect matchings through the familiar
quotient

```text
S4 -> S3,                 kernel V4.                         (7)
```

This is exactly the information loss seen in the Berggren variable-translation
language: the matching quotient forgets the `V4` origin and signed `H1`
orientation.  XOR chooses one Walsh line, not a canonical scalar invariant
of all four labels.

For a labelled tournament, the six edge signs are recovered over the
rationals from three independent score/cut coordinates plus the three Haar
cycle holonomies.  The corresponding integer determinant has absolute value
`32`; all `64/64` labelled tournaments are separated.  The non-unit index is
an integral lattice effect, not information loss over odd fields away from
two.

## Connection contract and boundary

| field | exact content |
|---|---|
| source | one formal `2 x 2` joint table on `F_2^2` |
| quotient | separate row/column margins |
| lost line | Walsh/XOR checkerboard `chi_11` |
| tetrahedral map | opposite-face `Omega` |
| target | one marked complementary-square class in `H_1(K4)` |
| preserved | nonzero mixed-Haar value / nonzero marked holonomy in odd characteristic |
| forgotten | absolute scale under projectivization, endpoint semantics, chronology and address |
| sidecar needed physically | one fixed map from actual endpoint records to the four joint atoms and their address |

The canonical U_full result explicitly proves that its present labelled
circle cells do not yet supply that last map.  The r=5 THM-2594 table supplies
a different common-base joint object, not an identification with U_full.
Accordingly (5)--(7) are a precise representation dictionary, not a physical
current, `D5` map, row exclusion, or LRC(14) closure.
