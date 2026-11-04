# Only the 59 output meets the exceptional two ready library

**PROVED arithmetic consequence of the inherited complete ready-root
classification; FINITE-EXACT numerical controls.** For integers m>=2 and
k>=1, the balanced-family value 15*2^m-1 equals one of

    3*2^k-1,  6*2^k+1,  5*2^k-1,
    13*2^k+3, 14*2^k+3

if and only if m=k=2 and the fifth expression is used. The common value
is 59.

These five expressions are the exact counts of the strong two-ready
roots with indegree two, from [READY sections 2 and 3](../../two_entry_funnels_20261103/READY.md).
The classification and its formulas are inherited. The present elementary
intersection calculation explains why the concrete 59 source exchange
cannot be repeated by merely changing the size of an exceptional template.

## Proof of the intersection statement

Equality with the first expression would give 2^k=5*2^m, and equality
with the third would give 2^k=3*2^m. Neither is possible: taking the
2-adic valuation forces k=m and then the odd factors disagree.

The second expression is 1 modulo 3, whereas 15*2^m-1 is 2 modulo 3,
so it cannot agree.

For the fourth expression the equation is

    15*2^m - 13*2^k = 4.

Both terms on the left are divisible by 2^min(m,k), so min(m,k)<=2.
If m=2, the equation requires 13*2^k=56, impossible. Otherwise m>=3,
so k is 1 or 2. These cases require 15*2^m to equal 30 or 56,
respectively. The first gives m=1, outside the stated range, and the
second is not divisible by 15.

For the fifth expression,

    15*2^m - 14*2^k = 4.

If m=2, this gives 2^k=4, hence k=2. If m>=3 and k>=2, both terms
on the left are divisible by 8, contradicting the right side. The only
remaining possibility is k=1, which requires 15*2^m=32 and is impossible.
This proves both directions, including the displayed positive control.

The condition m>=2 matters even for the numerical statement: m=k=1
gives 29 in the fourth expression. This is only an arithmetic boundary
control; the preceding session's actual C_m construction and count formula
were themselves stated for m>=2.

## Consequence for the missing source return

For every m>=3, any strong two-ready tournament of count 15*2^m-1 must
lie in the other branch of READY's complete classification: its ready
root has indegree one and is the newborn of a strong tournament with
exactly two Hamilton paths starting at its unique in-port. Thus an
H-preserving source exchange at the later balanced outputs cannot use
any of the five exceptional indegree-two families, regardless of how
their parameter k is chosen.

The same conclusion applies to a weak two-ready source after removing
its final transitive padding: [O1](../../ready_return_pentagon_20261103/OBSTRUCTION.md#o1-weak-sources-do-not-repair-a-positive-h-ready-obstruction)
shows that, for positive H>=3, the ready core has the same H and readiness.
That core must therefore be an indegree-one newborn of a two-start parent.

This is a necessary representation condition. It neither proves that
such parents exist at every m nor that their newborns can be reached
with the required named ports. The 59 construction is a valid base repair,
but uniform iteration needs a different source family or a method that
does not require a two-ready reset.

`exceptional.py` checks the five formulas at m=2,...,64 and k=1,...,66,
as well as the excluded numerical m=1 boundary. The finite check supports
the proof; no range of tested exponents substitutes for its valuation
argument. Normal and optimized output agree with `exceptional.json`.
