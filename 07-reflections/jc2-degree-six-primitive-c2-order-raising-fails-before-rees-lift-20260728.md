# The degree-six primitive quadratic-deck test fails before any Rees lift

**Status:** FINITE-EXACT REFUTATION OF THM-2206'S STATED CANONICAL
PRIMITIVE-ROWWISE TEST.  This note does not refute every possible integral
order-raising lattice, does not produce a split Keller survivor, and proves no
case of `JC(2)` or `DC(2)`.

Companion:

```text
04-computation/jc2_degree6_c2_primitive_order_raising_referee_20260728.py
05-knowledge/results/jc2_degree6_c2_primitive_order_raising_referee_20260728.out
```

Reproduce with

```bash
python3 04-computation/jc2_degree6_c2_primitive_order_raising_referee_20260728.py
python3 -O 04-computation/jc2_degree6_c2_primitive_order_raising_referee_20260728.py
```

The two modes are byte-identical.

## 1. Inheritance and the typed object

The closest proved mechanism is THM-2206: for

```text
A=Z_2[C_2],       I=(sigma-1),       I^j=2^(j-1)I,
```

an implication sending every anti-invariant defect from `I^j` to
`I^(j+1)` would iterate to zero.  THM-2206 names THM-2194's exact degree-six
Faber bank as the cheapest hostile test.

The canonical hostile is the exact square used as the positive algebraic
control in THM-2194's referee:

```text
T=(1+X)^2,        c=Lambda=Omega=0.                  (1)
```

The least-used relevant sidecar is the third Hamiltonian/Keller row `R_Q`.
It is not one of the boundary/`Phi`/`Psi` rows named in the test and is not
silently inserted below.

Keep the fields and modules separate, in the spirit of MISTAKE-308's typing
repair.  The Faber formulas are first rational polynomial formulas on one
chosen sheet.  The split integral deck is separately the algebra

```text
M=Z_2 x Z_2
```

with the swap action.  Its anti-invariant generator is

```text
epsilon=(1,-1),       IM=Z_2 epsilon,
I^jM=2^(j-1)Z_2 epsilon.                              (2)
```

An odd Faber seed changes sign between the sheets.  To keep the mate fixed,
its coefficient is also anti-invariant.  Thus an odd coefficient is
`k_j epsilon`; its boundary and second-flux products are invariant, while its
first-flux product is anti-invariant.  Treating `k_j` as deck-fixed would test
the wrong module.

## 2. The canonical primitive-row lattice

Use THM-2194's monic column normalization

```text
E6,E5,E3,E2,E1.                                      (3)
```

For degree `j`, its three chosen-sheet entries are exactly

```text
a^j B_j,       a^(j+1) F_j,       a^(j+2) G_j,       (4)
```

where `F` and `G` are first flux divided by four and second flux divided by
four.  Work formally over

```text
S=Z_(2)[a,c,Lambda,Omega].                            (5)
```

For a fixed rational row there is one canonical `Z_(2)`-line: multiply by a
power of two until all coefficients are 2-integral, then divide their common
2-content.  It is unique up to a `Z_(2)` unit.  With the columns (3) left in
their monic normalization, the exact primitive powers are

```text
boundary: 2^8,       F: 2^10,       G: 2^9.           (6)
```

The companion derives the rows from the Faber recurrence rather than copying
THM-2194's display.  Reducing (6) modulo two gives, in the full column order
(3),

```text
boundary  = [0,a^5,0,0,0],
F         = [0,a^6,0,0,0],
G         = [0,a^7 Lambda,0,0,0].                    (7)
```

Consequently the induced odd-column matrix, with columns `E5,E3,E1`, is

```text
[ a^5          0  0 ]
[ a^6          0  0 ]
[ a^7 Lambda   0  0 ].                               (8)
```

Over `F_2(a,c,Lambda,Omega)` it has rank one.  The `E3` and `E1` classes
survive.  If `a` is not a unit, even the `E5` conclusion can weaken; the
rank-one statement is already the most favorable generic reading.

Thus the canonical primitive-row test does **not** send every odd defect from
`I^j` into `I^(j+1)`.  Multiplying one row by another factor of two merely
kills its graded class; dividing it by two abandons the primitive integral
row.  Multiplication by a 2-adic unit does not change (7).

## 3. A primitive exact kernel at the square face

At (1), with `a=1`, the five exact rows are

```text
          B             F             G
E6        0             0             0
E5        3/256        -5/1024        0
E3       -1/16          3/128         0
E2        0             0             0
E1        1/2          -1/8           0.              (9)
```

Therefore

```text
128 E5 + 32 E3 + E1                                  (10)
```

lies in the exact characteristic-zero kernel of all three observables.  The
normalized degree-six bank

```text
E6 + 128 E5 + 32 E3 + E1                             (11)
```

has the same zero triple.  This is a statement about the three necessary
observables, not a claim that (11) is a Keller mate.

After the generic primitive row scalings (6), the nonzero specialized
odd-column matrix is

```text
[  3  -16   128 ]
[ -5   24  -128 ].                                   (12)
```

Its Smith factors are `1,8`; the row quotient has a `Z/8` torsion class.
Its primitive integer kernel is exactly

```text
(128,32,1).                                          (13)
```

In particular, for every `j>=1`,

```text
xi_j=2^(j-1)(128,32,1)epsilon
```

belongs to `I^j` but not `I^(j+1)`, because its `E1` coordinate is a unit on
the graded line.  All three observables of `xi_j` vanish exactly.  This is a
uniform all-grade failure of injectivity, not an artifact of reducing one
poorly scaled row.

Saturating the row lattice can remove the displayed `Z/8` cokernel torsion;
it cannot remove the exact kernel (13).  Excluding (13) by rescaling the `E1`
source column changes the monic Faber coefficient lattice and loses arbitrary
coefficient specialization.

## 4. What is refuted and what remains

The connection ledger is

```text
source:       odd split-deck Faber coefficient module;
target:       boundary/Phi/Psi observable rows;
map:          THM-2194's three rational Faber rows, made rowwise primitive;
preserved:    exact rows, monic seed normalization, and split parity;
destroyed:    the E3/E1 first augmentation classes;
hostile:      exact square and primitive kernel (128,32,1);
needed sidecar:
              an independently defined Rees/source lattice or another
              Hamiltonian observable, plausibly the R_Q/Keller row.
```

This closes THM-2206's explicitly named cheapest test negatively.  It does
not prove that no integral route can work.  A saturated Rees lattice carrying
Laurent pole order, coefficient weights, and ramification of the actual
quadratic root algebra would be **extra structure**, not inherited from
THM-2194's three rows.  It must be defined before reduction, shown
2-torsion-free, and shown stable under the required operations.

There is an additional specialization obstruction.  The chosen-sheet
coordinates include

```text
a=b/(2U),       c=D/b^2
```

and analogous rational expressions for `Lambda,Omega`.  THM-2194 works over
complex valued function fields and supplies no preferred 2-adic integral
model for these quantities.  Arbitrary complex specialization does not
preserve 2-divisibility.  A future all-degree argument therefore needs a
universal integral polynomial identity in original coefficients; pointwise
2-adic divisibility after choosing embeddings is not enough.

The cheapest repaired probe is to derive a denominator-cleared `R_Q`/Keller
row in the original quartic coefficients and test whether it detects (13)
without introducing 2-torsion.  If it does not, the three-row augmentation
route should be retired rather than repaired by ad hoc row or column weights.
