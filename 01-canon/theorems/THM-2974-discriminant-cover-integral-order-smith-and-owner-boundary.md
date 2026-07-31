---
id: THM-2974
title: "Discriminant-cover integral-order Smith and owner boundary"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.
  On THM-2769's divisor t=0, the THM-2971 discriminant-cover isomorphism
  preserves the common normalized rank-six algebra, double-transposition
  inertia, and all three V4 Kummer residues, whose common nonzero word is
  110.  It does not preserve the monogenic graph orders or a natural
  quartic-sheet owner.  The edge and orientation graph orders have
  discriminant exponents 6 and 18 versus maximal exponent 2, hence
  normalization colengths 2 and 8; their exact relative DVR Smith profile is
  (-1,0,0,1,2,4), so the orders are not nested.  Natural quartic inertia has
  no fixed sheet although both sextic actions have two.  The exact
  intertwiner determinant has boundary orders 6 at t=0, 1 above a simple
  residual-Delta root, and 2 above every simple J_or root.  This is an
  integral-order and affine-owner stopping theorem, not a Keller, SFC(4),
  JC(2), or degree-four exclusion.
source: codex-quartic-discriminant-integral-sidecar-2026-07-30
depends_on:
  - THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile
  - THM-2971-discriminant-cover-edge-orientation-sextic-algebra-intertwiner
related:
  - THM-2655-quartic-keller-resolvent-v4-quasietale-torsor-and-kummer-class-group-gate
  - THM-2685-equivariant-kummer-boundary-parity-completion-and-divisor-residue-gate
  - THM-2968-quartic-edge-and-oriented-cycle-s4-complements
script: 04-computation/quartic_discriminant_cover_integral_order_smith_thm2974.py
output: 05-knowledge/results/quartic_discriminant_cover_integral_order_smith_thm2974.out
script_sha256: 6a77c2f4408c0b3fcf554bb72508445fd5122069f65e35f95375d4cf5d66211a
output_sha256: 82fbbca57974a4cad62d3cf283c8345a0dddeb83443a984941a134a7fef1a20d
hash_basis: LF-normalized bytes
---

# THM-2974 -- discriminant-cover integral-order Smith and owner boundary

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Inheritance and statement

[THM-2971, the discriminant-cover sextic algebra
intertwiner](THM-2971-discriminant-cover-edge-orientation-sextic-algebra-intertwiner.md)
proves that the edge and oriented-cycle sextic algebras of a full-`S4`
depressed quartic become explicitly isomorphic after adjoining the square
root of the quartic discriminant.  It is a generic multiplication theorem.
It deliberately does not identify integral graph orders or reconstruct an
affine source-sheet owner.

This theorem computes those two missing coordinates on the sharp affine
hostile of [THM-2769](THM-2769-full-s4-pair-sum-affine-divisor-parity-hostile.md):

```text
f_t(X)=X^4-2X^2-8tX+1-4t.                              (1)
```

Its matching cubic and the two sextics are

```text
S_t(U)=U^3-4U^2+16tU-64t^2,
E_t(Y)=S_t(Y^2),
O_t(Z)=Z^6+c_4(t)Z^4+c_2(t)Z^2+c_0(t),                 (2)
```

with `O_t` in the exact THM-2971 orientation gauge.  The relevant invariants
factor as

```text
Delta=-4096t^2(27t^2-14t+3),
J_or=4096(t-1)(27t^3-25t^2+8t-1).                     (3)
```

Let `D=C(t)(v)`, `v^2=Delta`, and choose either completed prime of `D`
above `t=0`.  Since `v_t(Delta)=2` and the residue field is `C`, this
quadratic cover is etale and split at that divisor.  Write `R=C[[t]]` for
the chosen complete DVR and `M` for the common maximal rank-six `R`-order
inside the THM-2971 algebra.  Transport the orientation coordinate through
THM-2971 and define the two monogenic graph orders

```text
R_edge=R[Y]/(E_t),              R_or=R[T]/(O_t),
T=vA(Y),                        E_t'(Y)A(Y)=-8q mod E_t. (4)
```

Then:

1. the common maximal algebra has tame inertia `C2`, acting on either
   sextic as `2^2 1^2`;
2. the edge and orientation Kummer orders are respectively `(1,1,0)` and
   `(1,1,4)`, so both divisor-residue words are exactly `110`;
3. the graph-order discriminant exponents are `6` and `18`, while the
   maximal-order exponent is `2`, giving

   ```text
   length_R(M/R_edge)=2,             length_R(M/R_or)=8; (5)
   ```

4. in the two power bases, the relative DVR Smith profile of `(4)` is

   ```text
   (-1,0,0,1,2,4).                                      (6)
   ```

   Thus the graph orders are not nested in either direction;
5. natural quartic inertia is a double transposition and fixes no quartic
   sheet, although it fixes two sheets in each sextic action.  Therefore an
   inertia-fixed sextic sheet does not determine an inertia-fixed present
   quartic owner.

The normalized algebra, inertia, and Kummer residues survive.  The graph
order and natural affine owner do not.  In particular the preserved Kummer
word is the nonzero obstruction `110`, not a repaired Keller carrier.

## 2. Exact local inertia

Put `t=s^2`.  Expanding `(1)` about its two double roots at `t=0` gives

```text
X= 1+a s+O(s^2),                  a^2=3,
X=-1+b s+O(s^2),                  b^2=-1.               (7)
```

The loop `s -> -s` swaps both roots in each cluster.  Hence its natural
quartic inertia is a double transposition, with cycle type `2^2` and no
fixed quartic sheet.

The element is even, so THM-2968 and THM-2971 identify its two sextic
actions over the discriminant cover.  Direct relabelling gives cycle type

```text
2^2 1^2                                                     (8)
```

in each.  Because the residue characteristic is zero, the common normalized
rank-six algebra is tame.  Its maximal-order discriminant exponent is the
sum of the local different exponents:

```text
(2-1)+(2-1)=2.                                           (9)
```

Thus THM-2971 really does preserve maximal-order inertia.  The loss occurs
after choosing the two graph coordinates.

## 3. Newton polygons, Kummer residues, and graph-order indices

The coefficient Newton points and lower hulls at `t=0` are

```text
E_t: (0,2),(2,1),(4,0),(6,0),  hull (0,2),(4,0),(6,0);
O_t: (0,6),(2,2),(4,1),(6,0),  hull (0,6),(2,2),(6,0).  (10)
```

Therefore the edge and orientation root valuations are

```text
edge:        (1/2,1/2,1/2,1/2,0,0),
orientation: (1/2,1/2,1/2,1/2,2,2).                    (11)
```

Passing to the three binary radicands gives orders

```text
edge U=Y^2:        (1,1,0),
orientation W=Z^2:(1,1,4).                              (12)
```

Their reductions modulo two are both `110`.  This is also the divisor-wise
content of THM-2971's cubic Kummer identity

```text
F(U)=Delta U B(U)^2 mod S_t(U):                          (13)
```

the square contributes even order and `v_t(Delta)=2`.  The discriminant
cover preserves all three Kummer residues at the hostile divisor, hence
preserves the failure of THM-2685's zero-row condition.

The exact polynomial discriminants are

```text
Disc(E_t)=64q^2 Delta^2,
Disc(O_t)=2^66 q^12 J_or^4 Delta^3.                     (14)
```

At `t=0` their orders are `6` and `18`.  Both monogenic orders are integral
suborders of `M`; combining `(9)` and `(14)` with

```text
Disc(R_0)=Disc(M)[M:R_0]^2                              (15)
```

gives exactly the two colengths in `(5)`.

## 4. The relative Smith profile

The difference `8-2=6` agrees with the valuation of THM-2971's power-basis
determinant

```text
det(Phi)=2^30 q^5 J_or^2 v,
v_t(det(Phi))=5+0+1=6.                                  (16)
```

It would be incorrect to infer from `(16)` that one graph order contains
the other.  Write `v=tw`, where `w` is a local unit.  In the edge algebra,
the coordinate is `T=wT_0`, with

```text
T_0=
 Y((t-1)Y^4+(24t^2-12t+4)Y^2+144t^3-96t^2+16t)
 -------------------------------------------------------. (17)
             32t(27t^2-14t+3)
```

Column `k` of the power-basis matrix gains only the unit `w^k`, so its DVR
Smith exponents equal those computed from `1,T_0,...,T_0^5`.  Parity splits
that matrix into two `3 x 3` blocks.  The minimum orders of their minors are

| block | `delta_1` | `delta_2` | `delta_3` | Smith exponents |
|---|---:|---:|---:|---|
| even | `0` | `0` | `2` | `(0,0,2)` |
| odd | `-1` | `0` | `4` | `(-1,1,4)` |

Sorting their union gives `(6)`.  The negative exponent proves that neither
graph-order lattice contains the other.  The sum of the six exponents is
`6`, independently recovering `(16)`.

The exact companion verifies that `(17)` is not merely a fitted matrix: it
satisfies both

```text
E_t'(Y)T_0=64t^2 mod E_t,
(Delta/t^2)T_0^2=F(Y^2) mod E_t.                        (18)
```

Thus it is precisely the THM-2971 coordinate after removing the local unit
`w`, and the Smith calculation concerns the actual intertwiner.

## 5. Fixed sextic sheets are not quartic owners

All four roots in `(7)` remain finite over `t=0`; `(1)` is monic and defines
a proper finite root cover.  There is no omitted Jelonek sheet in this
hostile.  The double-transposition inertia fixes no quartic root.

It nevertheless fixes two edge sheets: the within-cluster pairs whose sums
tend to `+2` and `-2`.  These are fixed **pairs of swapped quartic sheets**,
not fixed quartic sheets.  The same inertia fixes two orientation sheets,
but by `(11)` both orientation coordinates collapse to zero with order two.

Hence the first failed implication is exact:

```text
inertia-fixed sextic sheet
    does not imply
inertia-fixed present quartic source sheet.             (19)
```

For an actual Keller Jelonek component, THM-2655 requires a finite affine
inverse quartic sheet fixed by inertia, and the semiregular `V4` kernel must
be invisible to that inertia.  Here the inertia is itself a nontrivial
`V4` element.  The example is not a Keller map and is not being used as a
Jelonek divisor; it is the minimal hostile to reconstructing the required
owner from sextic fixed points.

## 6. Every parameter-line boundary is non-unimodular for this coordinate

The factors in `(3)` are squarefree away from the displayed multiplicity
`t^2`, and the residual `Delta` and `J_or` factors are coprime.  Formula
`(16)` gives the exact orders of the coordinate determinant on the
discriminant cover:

```text
t=0:                           5 v_t(q)+v_t(v)=6;
simple root of 27t^2-14t+3:   v(v)=1;
simple root of J_or:           2v(J_or)=2.               (20)
```

Thus this exact edge/orientation coordinate change is a graph-order
isomorphism only on the generic `q Delta J_or !=0` open.  It fails to be an
integral graph-order isomorphism on every component where THM-2971 loses
separability or primitivity.  At a simple `J_or` root the map itself is
integral and `R_or` is a proper suborder of `R_edge`; the order-two
determinant says precisely that the change is non-unimodular.  One can
identify the two maximal orders after normalization, but normalization is
an extra sidecar and still cannot reconstruct the missing quartic owner.

## 7. Scope and next gate

The exact source-to-target contract is:

| item | outcome under the THM-2971 map at `t=0` |
|---|---|
| common maximal rank-six algebra | preserved |
| tame `C2` inertia and its sextic cycle type | preserved |
| three `V4` Kummer residues | preserved as the obstructing word `110` |
| monogenic graph-order lattice | not preserved; Smith profile `(6)` |
| present/omitted quartic-sheet owner | absent and unreconstructible |

The cheapest positive affine successor must therefore retain, before taking
the edge/orientation quotient, a common maximal order together with an
inertia-fixed natural quartic sheet and a zero divisor-residue row.  Another
discriminant identity, sextic fixed point, or field isomorphism cannot supply
those coordinates.

This theorem proves no Keller map, degree-four exclusion, SFC(4), `JC(2)`,
`DC(2)`, GMC, or LRC statement.  It does not say that a different affine
quartic chart with the required owner cannot exist.

## 8. Exact companion

Run

```text
python 04-computation/quartic_discriminant_cover_integral_order_smith_thm2974.py
python -O 04-computation/quartic_discriminant_cover_integral_order_smith_thm2974.py
```

Both modes byte-match after LF normalization the stored transcript

```text
05-knowledge/results/quartic_discriminant_cover_integral_order_smith_thm2974.out.
```

The companion checks the root-leading equations, both Newton polygons,
natural/edge/orientation inertia cycle types, fixed-sheet counts, the three
Kummer residues, maximal and graph-order discriminant exponents, both
normalization colengths, the exact inverse-derivative and squared-intertwiner
identities, all minors required for the relative Smith profile, the
determinant valuation, and every parameter-line boundary factor.  Every
truth-bearing check uses an explicit runtime requirement.  There is no
Python `assert` and no floating-point decision.

LF-normalized SHA256:

```text
script  6a77c2f4408c0b3fcf554bb72508445fd5122069f65e35f95375d4cf5d66211a
output  82fbbca57974a4cad62d3cf283c8345a0dddeb83443a984941a134a7fef1a20d
```

**QED (candidate pending independent hostile audit).**
