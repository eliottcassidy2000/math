---
id: THM-3144
title: "Mixed-depth selector persistence death barcode"
status: >
  PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM AUDIT.
  At support (1,3), bank I2, probability laws on the 41 physical depth-one
  and depth-two virtual-prefix currents exhibit an exact persistence split.
  The law 2/9 delta_(1)+2/9 delta_(2)+5/9 delta_(2,2) is strictly
  Hasse-positive for every degree 5 through 10 but fails at degree 11.  The
  degree-11 slice is itself nonempty, via a different exact four-state law.
  Nevertheless four upset rows at degrees 8, 10, and 11 have a positive
  integer combination strictly negative on all 41 states, so no single law
  persists through every degree 5 through 11.  This is a finite
  averaged-current theorem, not an original-response decomposition or a
  sequential stochastic pole flag.
source: root/multiscale-newton-flag/product-gamma-width3-2026-08-02
depends_on:
  - THM-3115-low-degree-monomial-fibre-newton-refinement-transport
  - THM-3120-row-pole-prefix-newton-flag-positivity
  - THM-3127-partition-refinement-strassen-upset-dual-and-filter-response
  - THM-3137-finite-stochastic-pole-selector-polytope-and-portability-wall
related:
  - THM-3136-one-sided-fixed-reference-elementary-tail-hasse-no-go
script: 04-computation/gmc_mixed_depth_selector_death_barcode_scout.py
output: 05-knowledge/results/gmc_mixed_depth_selector_death_barcode_scout.out
script_sha256: 085a8e833a8de818dbd29ef822618b63f37774e4b0d4fb608ca9279e7b226ac5
output_sha256: 8208f98fe059cfae54797da69310fef362c72254fb57bfb2d5084de900f835fe
hash_basis: LF-normalized bytes
---

# THM-3144 -- mixed-depth selector persistence death barcode

**PROVED CANDIDATE + VERIFIED-EXACT / AWAITING INDEPENDENT THEOREM
AUDIT.**

THM-3137 proves that probability laws on the eight physical depth-one pole
currents can be Hasse-positive through degree 9, but not at degree 10.  It
also proves that no probability law on the 33 physical depth-two currents is
Hasse-positive at degree 10.  The two fixed-depth faces nevertheless have a
positive join: a probability law mixing one- and two-pole prefixes is strictly
positive through degree 10.

That repair has a sharp persistence boundary.  The degree-11 selector slice
is nonempty, but no one law belongs simultaneously to all slices from degree
5 through degree 11.  Thus the object that dies is a common section across
degrees, not the last individual fibre.

## 1. The 41 physical stopping states

Use the THM-3120 reduced row at

```text
(a,b)=(1,3), bank I2.                                        (1)
```

Its physical reduced-pole multiset is

```text
(8,7,6,5,5,4,4,3,3,2,2,2,1,1,1,1).                         (2)
```

There are eight value-labelled depth-one states `(r)`, `1<=r<=8`.  There are
33 unordered depth-two submultisets `(r,s)` respecting the multiplicities in
`(2)`: all 28 distinct pairs and the five repeated pairs `(r,r)` for
`1<=r<=5`.  Let

```text
S=S_1 disjoint_union S_2,         |S_1|=8, |S_2|=33, |S|=41. (3)
```

For a state `sigma` in `S`, subtract its one- or two-letter virtual alphabet
from both the signed bank atoms and the distinguished residual alphabet `Q`,
exactly as in the fixed-`Q` current definition inherited by THM-3137.  Write

```text
G_N^sigma(mu)
 =Phi^sigma(h_N)m_mu[Q^sigma]
  -Phi^sigma(m_mu)h_N[Q^sigma].                              (4)
```

Every `G_N^sigma` has total mass zero.  For a probability law `lambda` on
`S`, put

```text
G_N(lambda)=sum_(sigma in S) lambda_sigma G_N^sigma.          (5)
```

By THM-3127, `(5)` is a nonnegative fine-to-coarse Hasse boundary exactly
when every coarsening upset has nonnegative mass.

For each degree, define the feasible selector slice

```text
P_N={lambda in Delta_40 : G_N(lambda) is Hasse-positive},     (6)
```

and the cumulative feasible set

```text
C_d=intersection_(N=5)^d P_N.                                (7)
```

The distinction between `P_N` and `C_d` is load-bearing below.

## 2. Exact finite row universe

For degrees `N=5,...,11`, the partition lattices have respectively

```text
10, 27, 47, 168, 573, 3588, 19542                           (8)
```

coarsening upsets, including the empty upset.  The numbers of nonzero
41-coordinate response rows are

```text
8, 25, 45, 166, 571, 3586, 19540,                            (9)
```

and the numbers having at least one negative coordinate are

```text
0, 1, 1, 7, 55, 871, 8642.                                 (10)
```

All rows in `(8)--(10)` are reconstructed from the THM-3115/3120 signed bank
and the physical multiplicities `(2)`.  No row is inserted as a numerical LP
fixture.

## 3. A strict common law through degree ten

Consider the three-state mixed-depth law

```text
lambda_(1)   =2/9,
lambda_(2)   =2/9,
lambda_(2,2) =5/9.                                           (11)
```

The repeated state `(2,2)` is physical because pole 2 has multiplicity three
in `(2)`.  Exact evaluation proves that `(11)` pairs strictly positively with
all

```text
8+25+45+166+571+3586=4401                                   (12)
```

nonzero response rows in degrees 5 through 10.  The sharp smallest positive
pairing occurs at degree 5 on the top upset `{(5)}`:

```text
9*minimum_upset_mass=802795680.                              (13)
```

Exact max flow equals exact demand in every degree 5 through 10.  Therefore

```text
lambda in C_10,                  so C_10 is nonempty.         (14)
```

The support of `(11)` meets both fixed-depth faces.  THM-3137 proves that
neither fixed-depth face meets `P_10`.  Hence mixing stopping depths, not only
mixing pole values, is essential.

## 4. The common law fails at degree eleven

At degree 11, law `(11)` has exact transport values

```text
flow   =37492336246227928606600,
demand =37492554756606147580488.                             (15)
```

Thus its deficit is

```text
218510378218973888,                                          (16)
```

with first unresolved negative type `(8,1,1,1)`.  In particular, `(11)` is
not in `P_11` and hence not in `C_11`.

Failure of this one law does not imply that `P_11` is empty.

## 5. An exact degree-eleven fibre law

Let

```text
D=366800496594452404486086550562442469686294.                (17)
```

Define a probability law on four states by the common-denominator numerators

```text
state    numerator
-----    ---------
(1)      109567364149254639329212553024928619727445
(1,2)     55718439327751924214134941704508925909799
(2,2)    200785902960670587548043927335254094230245
(7,8)       728790156775253394695128497750829818805.          (18)
```

The four numerators in `(18)` are positive and sum exactly to `D`.  Exact
evaluation of all 19,540 nonzero degree-11 upset rows gives nonnegative mass;
exactly three rows are equalities.  Exact max flow saturates demand.  Hence

```text
P_11 is nonempty.                                             (19)
```

The same law has the degree profile

```text
degree:  5     6     7     8      9      10     11
result: pass  pass  pass  fail   fail   fail   pass.          (20)
```

Thus the individual slices `P_N` are not nested.  A selector can leave the
feasible region and re-enter at a later degree.

## 6. Four-row cumulative Farkas certificate

Although `P_11` is nonempty, `C_11` is empty.  Consider the four coarsening
upsets

```text
F8       =P_8 \ {(1^8)},
F10      =P_10\ {(1^10)},
F11three =P_11\ {(1^11),(2,1^9),(2,2,1^7)},
F11one   =P_11\ {(1^11)}.                                   (21)
```

Let their exact 41-coordinate response rows be respectively

```text
R8, R10, R11three, R11one.                                  (22)
```

The positive integer combination

```text
R=461295 R8+3948 R10+R11three+22 R11one                     (23)
```

is strictly negative on every state in `S`.  The companion prints all 41
coordinates.  Their exact range is

```text
-213093302419995239808
 <= R_sigma <=
-340278874544980224 < 0.                                    (24)
```

Suppose `lambda` belonged to `C_11`.  Then its pairing with each of the four
necessary upset rows `(22)` would be nonnegative, and hence its pairing with
the positive combination `(23)` would be nonnegative.  But `(24)` and
`lambda in Delta_40` give

```text
<R,lambda><0,                                                (25)
```

a contradiction.  Therefore

```text
C_10 is nonempty,        C_11 is empty,       P_11 is nonempty. (26)
```

This is the exact persistence death barcode.

## 7. Exact verification

The companion uses only integers and exact rational arithmetic.  It:

1. enumerates all 41 multiplicity-valid physical states;
2. enumerates every coarsening upset in degrees 5 through 11;
3. reconstructs all 23,941 nonzero response rows and the sign census `(10)`;
4. checks strict positivity of `(11)` on all 4,401 degree-5-through-10 rows;
5. independently runs exact max flow in each positive degree;
6. verifies the exact degree-11 deficit `(15)--(16)`;
7. verifies every degree-11 upset on `(18)` and its three active equalities;
8. reconstructs the four rows `(21)--(22)` and all 41 coordinates of `(23)`.

Fresh normal and optimized runs byte-match the stored 20-line transcript.
The immutable LF-normalized hashes are recorded in the frontmatter.

## 8. Scope and first failed implication

The theorem is about probability averages of the derived virtual-prefix
currents `(4)` on one fixed support/bank fibre.  It does not prove:

- a positive decomposition of the original product-Gamma scalar response or
  operator;
- a sequential Markov process whose stopping distribution is `(11)`;
- a canonical transition map between the distinct selector slices `P_N`;
- compatibility across supports, banks, or later prefix depths;
- arbitrary-degree Hasse positivity, arbitrary-radial NC2, or the Gaussian
  Moment Conjecture.

The exact first failed implication is

```text
every relevant degree has (or nearly has) a positive selector
  does not imply
one common selector exists across those degrees.              (27)
```

The strongest survivor is not a longer fixed pole flag.  It is the finite
stopping-depth law `(11)` through degree 10, together with the separate live
degree-11 fibre `(18)`.  The missing sidecar is a lawful transport or gluing
map between selectors at adjacent degrees.

QED (candidate pending independent theorem audit).
