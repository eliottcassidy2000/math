# A cyclic wrap operator and the exact self-converse 7–21–63 ray

## The full wrap transfer

Let $W(T)=C_3[T,1,1]$, with quotient arcs $T\to p\to q\to T$.
It is strong; if $T$ is self-converse, the same anti-involution on
$T$ and the swap of $p,q$ make $W(T)$ self-converse. For the
ordered-slot polynomial $Q$ with empty slots allowed,
$$
 Q_{W(T)}(z)=z^2Q_T(z)+z\big(Q_T(z+1)-Q_T(z-1)\big).
$$
In particular $H(W(T))=H(T)+Q_T(2)$ for nonempty $T$.

Proof: in the [Grinberg–Stanley odd-cycle packing expansion](https://arxiv.org/abs/2307.05569),
an unused vertex contributes $z$ and an odd cycle contributes $2z$.
Without a mixed cycle, the two singleton factors give $z^2Q_T(z)$.
A mixed cycle has exactly the form $q\to P\to p\to q$, with $P$
an odd-order directed path on a nonempty subset $S$ of $T$. Its
contribution is $2z H(T[S])Q_{T-S}(z)$.
Finally $Q_T(z+y)=\sum_S Q_{T[S]}(y)Q_{T-S}(z)$ by grouping the
slots into two sets, and $Q_{T[S]}(-1)=(-1)^{|S|}H(T[S])$ by
the packing parity. Taking $y=1$ and $-1$ isolates the odd subsets.

## Why repeating the same wrap is not 3-adic refill

For odd-order input let $a=Q_T(1), b=Q_T(2), c=Q_T(3)$.
After three wraps the count is $4(a+b)+c\bmod9$. Every further
wrap doubles that residue. Thus from the third step on either
all counts are divisible by 9 or all have exact 3-adic valuation one.

For detail, one wrap acts on evaluations at $1,2,3$ by
$(a,b,c)\mapsto(a+b,4b+2c-2a,9c+3(d-b))$, with $d=Q_T(4)$.
Odd $Q_T$ gives $b=-a\bmod3,\ c=0\bmod3,\ d=a\bmod3$.
After the first wrap every evaluation is a multiple of 3,
and after the second the evaluation at 3 is 0 mod 9. At the
third step the values at 1 and 2 coincide mod 9 (both become
the sum of the previous values at 1 and 2), and then double.
Expanding that sum at step three gives
$-11a+13b+10c+6d=4(a+b)+c\bmod9$.

For a singleton seed, the residue at the third step is 6,
and the first and second counts are 3 and 15. So the entire
positive-depth tower has exact valuation one. For the order-nine
self-converse 567 seed below the initial evaluations are
$(567,13920,150405)$, giving residue 3 at the third step.
Simple repeated strong wrapping of this seed eventually loses
every divisibility layer above the first, rather than climbing the
$7\cdot3^j$ ray.

## Exact ray support with no strongness claim for the propagation

The [previous exact self-converse exclusion](../session315_dyadic_pairing_20261113/SPECTRUM_AND_SIEVE.md)
rules out 7, 21 and 63. The first two positive seeds on the ray are:
* $H=189$: Paley tournament on 7 vertices, arcs given by differences
  $1,2,4$, with anti-involution $i\mapsto-i\bmod7$.
* $H=567$: outgoing rows $(446,312,66,148,4,92,411,310,60)$,
  anti-involution $(0\ 4)(1\ 5)(2\ 6)(3\ 7)$ fixing 8.

Both seeds have exactly nine fixed Hamilton paths under their chosen
anti-involution composed with path reversal; the companion script
checks the literal masks, paths, and strongness independently.
Sandwiching a self-converse tournament between $C_3$ and its converse
in an ordinal sum multiplies $H$ by nine and $F$ by three
(the two outer paths are paired by reflection).
Consequently exactly the first three terms are missing, and for
every $j\ge3$ the ray has a paired witness with
$$
 (H,F)=(7\cdot3^j,\ 3^{\lceil j/2\rceil}).
$$
The two starting witnesses are strong. The ordinal sandwiches are
not, and do not solve strong self-converse ray coverage.
