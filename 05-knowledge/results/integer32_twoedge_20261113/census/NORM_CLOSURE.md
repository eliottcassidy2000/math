# Two operations reveal the norm quotient's missing ternary phase

Status: **PROVED IN NOTES after independent audit.** The exhaustive
64-state check is **FINITE-EXACT**. This is a finite-field observer model, not an
asserted realization by tournament edge moves.

## Inherited object and the new operation

[The existing sextic clock](../../sextic_subfield_clock_20261004.md) proves that
F=F2[t]/(t^6+t^3+1) is F64, b=1+t has order63, and

    pi(x)=(x^21,x^9) in F4* × F8*       (x nonzero)

is a surjection with kernel F4* of order3. Extend pi to0 by pi(0)=(0,0).
Its partition has22 cells: the singleton0 and21 scalar orbits of size3.
Multiplication M(x)=bx already descends through pi; its visible nonzero clock
has period21 and the actual state clock has period63. This inherited loss does
not make the multiplier quotient non-Markovian. The new ingredient here is
translation A(x)=x+1, used alongside M on the same64 states.

For a partition P, one simultaneous refinement replaces the label at x by

    (P(x), P(Mx), P(Ax)).

Repeated labels are identified. With only A or only M omit the other component.
The claim is exact:

| Allowed operations | Cell counts through strict refinements | Stable loss |
|---|---|---|
| M alone |22|all ternary scalar phases|
| A alone |22 →63|one pair in F4*|
| M and A |22 →63 →64|none|

## Proof candidate

Two nonzero elements have equal norm observations iff they differ by a scalar
in F4*. This follows because the common kernel in the cyclic group of order63
has order gcd(21,9,63)=3, precisely F4*. Let phi=b^21; then
F4={0,1,phi,phi^2}, phi^2=phi+1.

M preserves every equivalence of the initial norm partition: if y=cx then
My=cMx. Consequently M contributes no extra split in the first refinement.

Consider a nonzero initial cell F4* x. If x is outside F4, the vectors1,x are
linearly independent over F4. Suppose two translated representatives have equal
observations, so cx+1=lambda(dx+1) for c,d,lambda in F4*. Independence forces
lambda=1 and c=d. None of these translates is zero. Thus A splits each of the
20 such cells into three singletons.

In the remaining nonzero cell F4*, translation sends1 to0 and exchanges phi
with phi^2. Hence it splits this cell into {1} and {phi,phi^2}. Together with
{0}, this is63 cells. For A alone the last pair is invariant and all other
cells are singletons, so this partition is stable.

Because b is outside F4, b phi and b phi^2 are distinct elements outside F4.
They therefore lie in distinct singleton cells after the first refinement.
M splits the one surviving pair in the second refinement. The resulting
partition is discrete and stable. The first round really leaves two states
indistinguishable, so two strict rounds are necessary as well as sufficient.

Equivalently, the two states phi,phi^2 have equal initial observations and
equal observations after either one-letter operation, but the future word
A after M separates them. Refinement labels through depth k record all norm
observations of operation words of length at most k (induction on k). Thus
these are the coarsest refinements preserving the initial observation and
closed under the declared operations. No universal two-round bound follows.

## Typed connection and hostile controls

Source: the64 field states. Target: the initial22-cell norm observer, then its
operation-closed refinements. Map: pi and the displayed refinement. Preserved:
exact norm values, zero, and multiplication response. Initially lost: the
ternary scalar phase. Needed sidecar: enough translated future observations
to recover that phase. The cheapest decisive test is the pair phi,phi^2:
with multiplication alone it stays equivalent forever; with A after M it
separates. This gives an exact small comparison with merged-metagraph observer
closure, whose states and operations remain different.

The underlying arithmetic also separates32 from64: the Cartesian product
F4×F8 has32 elements, while F4 tensor_(F2) F8 is a six-dimensional binary
algebra with64 elements. Our22 is21 multiplicative norm fibres plus zero;
it is not a proposed meaning for every historical occurrence of22.

## Reproduction

From the repository root:

```sh
python3 05-knowledge/results/integer32_twoedge_20261113/census/check_norm_closure.py
python3 -O 05-knowledge/results/integer32_twoedge_20261113/census/check_norm_closure.py
```

The standard-library checker uses explicit exceptions. Normal and optimized
receipts are byte-identical, with12,434 gates. It tests every pair of states
for the norm/scalar-orbit equivalence, both separate operators, their joint
refinement, and the independent bounded-word description of each partition.
Polynomial-basis codes8 and9 are the last pair; their numbers are labels, not
invariant arithmetic distinctions. The inherited root t has order9 rather
than63, providing a hostile control for conflating field dimension with a
chosen element's multiplicative order.
