# Every induced slot polynomial is a cyclic-triple invariant

For a fixed labelled ground set, the following data have exactly the
same equality fibres among tournaments:
1. the labelled set of cyclic triples;
2. the whole table H(T[S]), S a labelled subset;
3. the whole table Q_(T[S])(z), S a labelled subset;
4. the number of directed Hamilton cycles on each labelled odd subset
   of size at least three (cyclic order, not rooted).

This specializes the older [complete unsigned curl theorem](../session36_information_curl_contexts_20261023/context/CURL-AND-CLAN-COMPLETION.md).
That result already identifies the ordinary commutative walk series
and full substitution response with labelled cyclic triples.
Here the squarefree logarithm and the marked observer comparison are
the main point, not a new interval-inversion theorem.

## Exact proof and primitive coordinates

Work over the rationals with commuting variables x_i and x_i²=0.
Put x^S=product over S and H(empty)=1. Define
```
P_T = sum_S H(T[S]) x^S.
```
Choosing one disjoint subset per labelled slot, and then its path,
gives Q_(T[S])(z) = [x^S] exp(z log P_T), first for nonnegative
integers z and therefore as a polynomial. The exponential and
logarithm are finite because the positive-degree ideal is nilpotent.
Thus (2) and (3) are equivalent. The coefficients of log P_T are
the squarefree cumulants. The odd-cycle packing identity from
[Grinberg–Stanley](https://arxiv.org/abs/2307.05569) gives
```
log P_T = sum_i x_i + 2 sum_(S odd, |S|>=3) C_T(S) x^S
```
where C_T(S) is the Hamilton cycle count of the induced tournament.
Indeed exponentiation chooses disjoint primitive supports; the
factorials cancel the orders of the chosen factors, reproducing the
odd-cycle packing sum with weight two for each selected cycle.
This proves (2) equivalent to (4).

On three vertices H is 3 exactly on a cyclic triple, otherwise 1,
so (2) determines (1). Conversely, by [the Boussaïri–Ille–Lopez–Thomassé
inversion theorem](https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf),
two tournaments with the same labelled cyclic triples differ by a
sequence of module reversals. Every intersection of a module with an
induced subset is still a module. In any slot cover of an induced
subset, reverse each maximal run in the module. This is a bijection
to covers after the reversal, preserving the slot count. It proves
(1) implies (3) without identifying any named internal direction.

Consequently in the order-22 collision not even all induced odd-cycle
counts distinguish the pair. There is also a [strong order-24
collision](STRONG_SHADOW.md) satisfying the entire hereditary table
and all full-tournament first slot marginals. This does **not** say
that first marginals on every induced subset are also invariant.
For instance one named endpoint count on a two-vertex subset already
identifies its direction. Nor does the single univariate Q determine
the full labelled cyclic-triple hypergraph.

## A constructive interpretation and its limit

The primitive cycle coefficients are determined by the same local
triangle structure, but they are not independent choices. A collection
of triples declared cyclic may not come from any tournament. And even
among realizable triple fields, the high-order odd-cycle cumulants
enter the target count through compatibility of disjoint supports.
An additive “count one triangle, add two paths” strategy can therefore
fail when new triples affect higher primitives. Wedge observations
supply the missing signs, while the hereditary unmarked slot table
is the appropriate quotient when only scalar uniform continuations
are needed.
