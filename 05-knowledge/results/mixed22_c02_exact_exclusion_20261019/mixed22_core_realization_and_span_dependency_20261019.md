# A realizable core-only obstruction and its short-span exclusion

## Scope

This note records a limitation of a boundary-generator method. The displayed
actual poset has balanced pairs. It is not a counterexample to the 1/3–2/3
conjecture, nor an unresolved seven-thin morphology: the candidate maximum-span
order is excluded by the published six-thin theorem stated below.

## 1. The six-point core and exact unrestricted-cone obstruction

Label the core `(x,y,b1,b2,a1,a2)=(0,1,2,3,4,5)`. Its predecessor masks are
`(0,12,0,4,1,21)`, or equivalently its Hasse relations are

- `b1 < b2 < y`,
- `x < a1 < a2`,
- `b1 < a2`.

There are 19 core linear extensions. The candidate majority order is
`b1,x,b2,a1,y,a2`. List its forward incomparable pairs in the order

`(x,y), (b1,x), (x,b2), (a1,y), (y,a2), (b1,a1), (b2,a1), (b2,a2)`.

Let a generator restrict a core extension to a prescribed prefix and suffix.
The three generators with boundary restrictions `(empty,empty)`,
`(b1,empty)`, and `(empty,a2)` have extension counts 19, 10, and 10.
Coefficients `1,5,5` give total mass 119 and forward comparison numerators

`108,90,80,80,90,111,80,111`.

Every resulting probability exceeds 2/3, with minimum `80/119`.
Equivalently, normalized generator-mixture weights are `19/119,50/119,50/119`.
This is an exact failure of the nonnegative generator cone screen.

## 2. Product structure alone does not fix the failure

Impose the rank-one fiber weight

`(1 + u*1[b1 is first]) * (1 + v*1[a2 is last])`, where `u,v >= 0`.

The denominator is `19+10u+10v+6uv`. For the three particularly restrictive
forward comparisons `x<b2`, `a1<y`, and `b2<a1`, the unnormalized `3N-2E`
polynomials are respectively

- `7 - 2u + v - 3uv`,
- `7 + u - 2v - 3uv`,
- `-8 + u + v + 3uv`.

At `u=v=15/11`, all eight forward deficits are positive. The normalized
minimum forward probability is `4635/6949`, and the minimum `3p-2` is
`7/6949`. The receipt contains all nine prefix/suffix generator columns,
including the longer `b1,b2` and `a1,a2` boundaries.

## 3. An actual 16-point realization

Take a five-point lower exterior `A0,...,A4` with predecessor masks
`(0,0,1,1,7)`. It has 11 extensions. Put every Ai below `x,b2,y,a1,a2`;
put `A0,A3` below b1 and leave the other three Ai incomparable with b1.
The six-point poset consisting of A and b1 has 26 extensions.

Add the order-dual upper exterior `B0,...,B4`, using the involution
`x<->y`, `b1<->a2`, `b2<->a1`, `Ai<->Bi`. Take transitive closure.
The resulting predecessor tuple, in the order of the six core vertices,
then A0,...,A4, then B0,...,B4, is

`(1984,1996,576,1988,1985,2005,0,0,64,64,448,59391,34783,34783,2047,2015)`.

For every fixed core extension w the number of ambient extensions is exactly

`(11+15*1[w starts b1])*(11+15*1[w ends a2])`.

Indeed, if b1 is not first, the first core vertex is x, and all A must already
occur. This gives 11 left interleavings. If b1 is first, A may interleave b1,
but must end before the next core vertex; the left interleaving count is 26.
The right factor follows by duality. Left and right choices are independent.

Thus the actual poset has 6,949 extensions and maximum incomparability degree
5. Its forward core comparison numerators are

`6273,5200,4635,4635,5200,6421,4645,6421`.

Every core comparison is biased by more than 2/3 in the candidate direction.
Nevertheless there are eight balanced pairs, all involving exterior vertices.
In particular

`p(A4<b1) = p(a2<B4) = 3949/6949`.

Consequently even genuine exterior realizability, factorial-moment conditions,
and rank-one boundary structure cannot imply an inequality that excludes
this core using only its eight forward comparisons. Such an inequality would
be false on the displayed actual poset. Additional hypotheses or exterior
comparison rows are essential.

## 4. Maximum-span anchors

Suppose hypothetically that this core comes from a poset with no balanced
pair and that x,y is a maximum-index-span incomparable pair in its strict
majority order. All exterior vertices are below both roots (A) or above both
roots (B). They cannot occur between the roots in the majority order.
Hence the full root span equals its projected core span, which here is 3.

The incomparable pair b1,a1 also has projected span 3. Any A between b1 and
x would increase its full span without increasing the root span, impossible.
Thus all A precede b1. Dually, all B follow a2. In particular the hypothetical
majority order requires `p(Ai<b1)>2/3` and `p(a2<Bi)>2/3` for every
incomparable contact, requirements violated by the explicit realization.

More precisely, only the last three A vertices can be incomparable with b1,
only the last A vertex with b2, only the first B vertex with a1, and only the
first three B vertices with a2. This follows immediately by comparing their
full index distances with the maximum span 3. These are structural deductions,
not an exhaustive new exterior-witness proof.

## 5. Existing theorem removes all short-span orders

For any root core `Q=I(x) union I(y)`, each vertex outside Q is comparable
with both roots. It lies below both or above both: either mixed direction
would imply the roots are comparable by transitivity. Therefore no exterior
vertex lies between x and y in any linear extension, including the strict
majority order. The projected root span is its full span.

If that span is at most 3 and the root has maximum span among all
incomparable pairs, each vertex can be incomparable with at most three
vertices before it and three after it. The entire ambient poset is six-thin.

Marcin Peczarski's *The Gold Partition Conjecture for 6-Thin Posets*,
Order 25 (2008), 91–103, proves the Gold partition conjecture for posets in
which each element is incomparable with at most six others; this also implies
the 1/3–2/3 conjecture. See the publisher's primary abstract:
https://link.springer.com/article/10.1007/s11083-008-9081-9
(DOI `10.1007/s11083-008-9081-9`, checked 19 October 2026).

Thus all projected root orders of span at most 3 may be discarded in the
seven-thin exclusion search, with this explicit published-theorem dependency.
This is not claimed as a new theorem or a new proof of the six-thin result.
The separate continuation searches only orders of span at least 4.

## Reproduction and exact controls

Run, with Python 3 and no external packages:

```sh
python -O verify_mixed22_core_realization_20261019.py
```

The script writes `mixed22_core_realization_receipt_20261019.json` and checks:
19 fixed-core fibers, all 9 generator columns, all 6,949 actual extensions,
and 120 pair counts independently by subset dynamic programming. It checks
all 8 balanced exterior pairs and the maximum incomparability degree 5.
All checks remain enabled under Python's optimization flag.
