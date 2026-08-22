# The support-five degree-21 map is exactly maximal, but maximal is not full

Status: **PROVED structural rank ceiling and characteristic-zero rank** +
**FINITE-EXACT guarded sparse certificates**; full support five remains
**OPEN**.

This is canonized as
[THM-3323](../01-canon/theorems/THM-3323-cyclic-quartic-support-five-exact-degree-21-rank.md)
and continues
[THM-3321](../01-canon/theorems/THM-3321-hesse-moment-kernel-and-cyclic-quartic-support-four-exclusion.md)
without changing its coefficient chart.  Retain

```text
g=A*zbar+B*z^2+C*z*zbar^2+D*z^3*zbar+E*zbar^4
```

and `M_m=<g^m>` for `m=3,6,...,21`.  The inherited Hesse kernel

```text
mu(a,b)=2*a!*b!*C(a,b)/(a+b+2)!,
C(a,b)=[X^aY^b](1-X^3-Y^3-3XY)^(-1)
```

is the entry oracle throughout.

## 1. The first correction is structural, not computational

The proposed full-support degree-21 map has

```text
Phi_21 : direct_sum_(d=3,6,...,21) R_(21-d) -> R_21,
(h_d) |-> sum h_d M_d,
dim source = 13972,                 dim target = 12650.       (1)
```

It cannot have full column rank.  Here is a proof which does not assume a
generic-Hilbert-series conjecture.

First retain only the five forms of degrees `3,6,9,12,15`.  In the irreducible
parameter space of all five-tuples of forms of these degrees, matrix rank in a
fixed degree is generically maximal.  The locus of regular sequences is a
nonempty open set: the tuple

```text
x0^3, x1^6, x2^9, x3^12, x4^15
```

is one witness.  Hence the generic quotient has the complete-intersection
Hilbert series

```text
(1-t^3)(1-t^6)(1-t^9)(1-t^12)(1-t^15)/(1-t)^5,             (2)
```

whose coefficient at degree 21 is `1705`.  Every specialization has Macaulay
rank at most the generic rank, so every specialized quotient has degree-21
dimension at least `1705`.

Now add `M_18`.  Its degree-21 multiples have cubic multipliers, a
35-dimensional space.  But multiplication

```text
R_3 -> (R/(M_3,M_6,M_9,M_12,M_15))_21,    h |-> h*M_18     (3)
```

has the nonzero cubic `M_3` in its kernel.  Thus `M_18` adds rank at most 34.
Finally `M_21` has only its scalar multiple and adds rank at most one.  It
follows rigorously that

```text
dim (R/(M_3,...,M_21))_21 >= 1705-34-1 = 1670,
rank(Phi_21) <= 12650-1670 = 10980.                          (4)
```

This is the correct hostile control.  A degree-21 full-rank signal would be an
engine bug, not a proof of emptiness.

The formal product

```text
product_(d=3,6,...,21)(1-t^d)/(1-t)^5                       (5)
```

has coefficient `39` at degree 28 and `-354` at degree 29.  Therefore degree
29 is merely the first target **not ruled out by this formal count**.  No claim
about the generic Hilbert series of seven forms in five variables is being
made.

## 2. Memory-safe deterministic elimination

The companion

```text
04-computation/degree_four_support5_degree21_sparse_rank_scout_20260803.cpp
```

builds every entry from the Hesse coefficient sum, in exact arithmetic modulo
a guarded prime.  It orders fixed-degree monomials by descending lexicographic
order.  Each transient row uses one reusable dense byte accumulator; only
normalized pivot rows are retained, as sparse `(uint16 column, uint8 value)`
pairs.  This avoids the roughly 177 million entries of the dense rectangle.

Reproduction command:

```bash
g++ -O3 -DNDEBUG -std=c++20 \
  04-computation/degree_four_support5_degree21_sparse_rank_scout_20260803.cpp \
  -o /tmp/degree_four_support5_rank && \
/tmp/degree_four_support5_rank
```

At `p=101`, the rank additions in generator order are

```text
generator       M3    M6   M9   M12  M15  M18  M21
new rank      7315  2056  930   470  174   34    1
cumulative    7315  9371 10301 10771 10945 10979 10980.     (6)
```

The final store has `12,964,187` pivot nonzeros, its densest pivot has `5296`
nonzeros, and an instrumented optimized run used about 81.6 MB maximum resident
memory.  The algorithm is ordinary exact row reduction, with no random
projection, floating point, or probabilistic rank inference.

The same computation at `p=103` again gives rank `10980`; the per-generator
rank additions are identical.  Its final pivot store has `12,967,780`
nonzeros.  Different modular cancellations change the sparse fill, but not the
rank.

## 3. A maximal-minor certificate and the rational lift

Whenever a source row creates a new pivot, the program records that source row,
the pivot column, and the raw pivot before normalization.  In discovery order,
subtracting earlier selected rows makes the selected square matrix upper
triangular.  The product of the raw pivots is therefore the determinant in
pivot-discovery column order.  The program computes the permutation sign to
natural increasing column order and applies it explicitly.

At `p=101`, the output displays all `10980` selected rows and columns and gives

```text
raw pivot product                         = 82 mod 101,
pivot-to-natural column permutation sign = -1,
det(selected natural-column minor)        = 19 mod 101.       (7)
```

This proves rank at least `10980` over `F_101`; (4) proves the matching upper
bound.  Since `101>4*21+2=86`, every simplex-moment denominator is invertible.
The displayed modular minor is the reduction of a rational minor, so (7) also
proves rank at least `10980` over `Q`.  Combining with (4) yields

```text
rank_Q(Phi_21)=10980,
dim_Q (R/(M_3,...,M_21))_21=1670.                            (8)
```

Thus the degree-21 computation is not merely a scout: its exact rank lifts to
characteristic zero.

## 4. Controls

For both `p=101` and `p=103`, the same independent sparse engine rebuilds the
support-four delete-`E` system from the Hesse oracle and obtains

```text
2926 x 2024, rank 2024, nullity 0,                           (9)
```

matching THM-3321's dense FLINT certificate.  Omitting `M_15,M_18,M_21` gives
the hostile control

```text
2821 x 2024, rank 1978, nullity 46,                          (10)
```

at both primes.  The code also requests `p=83` and verifies that the
denominator guard rejects it before any rank computation.

Separate `-O2` and `-O3 -DNDEBUG` builds give byte-identical `p=101`
transcripts, including the selected minor and determinant, and both match the
`p=101` prefix of the frozen output.

## 5. What changed on the concept board

The support-five lane did not fail because the matrix was too large.  It failed
because degree 21 was the wrong proof object.  The correct output is a sharp
Hilbert slice:

```text
universal lower bound 1670
          +
explicit nonzero maximal minor
          =
exact degree-21 quotient dimension 1670.                    (11)
```

This is useful information about the open chart.  The first five moments use
all rank permitted by the complete-intersection comparison at this degree;
`M_18` uses all 34 ranks permitted after its forced `M_3*M_18` overlap; and
`M_21` contributes its one available rank.  There is no hidden early
degeneracy to exploit at degree 21.

The next move should not instantiate the vastly larger degree-29 rectangle
densely.  Two cheaper routes are now better motivated:

1. continue the sparse quotient degree by degree, reusing the degree-21 pivot
   basis and watching when the Hilbert function first departs from the formal
   count;
2. reinterpret the 1670-dimensional dual quotient through the Hesse surface
   `uv=r^3`, looking for a structural block or recurrence before forming the
   degree-29 map.

## 6. Boundary

**PROVED:** the universal degree-21 nullity bound `>=1670`, exact rational rank
`10980`, and exact rational degree-21 quotient dimension `1670`.

**FINITE-EXACT:** deterministic ranks and explicit maximal-minor certificate at
`p=101`, with an independent exact rank repetition at `p=103`, support-four and
deficient hostile controls, and the `p>86` guard.

**OPEN:** whether the seven moments have a common full-support projective zero;
whether a later Macaulay degree closes the chart; whether degree 29 attains full
rank; the rest of the cyclic degree-four eigenspace beyond this Hilbert slice;
and `FC(3)` in general.
