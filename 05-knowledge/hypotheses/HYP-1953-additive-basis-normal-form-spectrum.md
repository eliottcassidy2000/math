---
id: HYP-1953
status: OPEN
source: codex-2026-06-01-S494
related:
  - HYP-1902
  - HYP-1920
  - HYP-1931
  - HYP-1941
  - HYP-1952
---

# HYP-1953: additive basis theorems form a coverage-density-normal-form spectrum

## Statement

Goldbach, Helfgott's ternary Goldbach theorem, the Hardy-Littlewood Goldbach
heuristic, Fermat's polygonal number theorem, and Zeckendorf's theorem are all
views of the same additive-basis machine.

For an atom set `A`, target `n`, and arity `r`, form the representation
hypergraph:

```text
R_A,r(n) = { {a_1,...,a_r} : a_i in A and a_1+...+a_r=n }.
```

The deep structure is not merely whether `R_A,r(n)` is nonempty.  It is the
triple:

```text
(representation count, local obstruction product, normal-form carry debt).
```

Goldbach lives at the sparse/probabilistic end: primes are thin atoms, the
binary edge count is conjecturally positive on even targets, and the
Hardy-Littlewood singular series predicts the local residue correction.

Helfgott's ternary theorem is the first smoothing lift: adding one prime turns
the thin binary edge problem into a three-dimensional surface with enough
averaging for a proof.

Fermat polygonal is the bounded-cover regime: polynomial atoms are structured
enough that allowing at most `k` `k`-gonal numbers beats all obstructions.

Zeckendorf is the unique-normal-form regime: Fibonacci atoms do not give many
representations after quotienting; instead, a local no-adjacent carry law gives
one canonical independent set in a path graph.

## Evidence

`04-computation/additive_basis_normal_forms_s494.py` computes four probes.

### Hardy-Littlewood layer

For even `n`, the script compares Goldbach pair counts with the rough unordered
Hardy-Littlewood heuristic:

```text
C_2 * n/log(n)^2 * prod_{p|n,p>2} (p-1)/(p-2).
```

For selected values:

```text
n=1000:  pairs=28,  HL_unordered=18.45, ratio=1.518
n=5000:  pairs=76,  HL_unordered=60.67, ratio=1.253
n=10000: pairs=127, HL_unordered=103.76, ratio=1.224
```

The main point is structural: Hardy-Littlewood is

```text
archimedean volume * p-adic local product.
```

This is the additive-prime analogue of the repo's repeated product/debt
ledgers.

### Helfgott layer

The ternary prime hypergraph has far more mass:

```text
n=99:   30 unordered triples
n=501:  283 unordered triples
n=999:  769 unordered triples
n=1999: 3105 unordered triples
```

The third summand is not cosmetic.  It increases the dimension of the
representation surface and gives the circle method enough averaging to prove
positive mass.

### Fermat polygonal layer

The script checks minimal `k`-gonal term counts up to `300`:

```text
k=3: max min-terms=3
k=4: max min-terms=4
k=5: max min-terms=5
k=6: max min-terms=6
k=7: max min-terms=7
k=8: max min-terms=8
```

This is the bounded-cover side of the spectrum: enough summands make the
structured atom set universal.

### Zeckendorf carry layer

For each Goldbach pair `p+q=n`, S494 overlays the Zeckendorf supports of `p`
and `q`, then compares the raw digit multiset with the canonical Zeckendorf
support of `n`.  The carry score is:

```text
L1(raw digits, target digits) + repeats + adjacent-index violations.
```

Selected results:

```text
n=50:   best pair (13,37), score=0
n=100:  best pair (11,89), score=0
n=500:  best pair (13,487), score=4
n=5000: best pair (673,4327), score=0
```

This suggests a new feature: among many Goldbach edges, some are compatible
with the target's Fibonacci normal form while balanced prime pairs may have
larger carry debt.  Abundance and normal-form compatibility are separate axes.

## Predictions

1. Goldbach-like searches should record not only representation count but also
   normal-form carry debt relative to a chosen recurrence basis.
2. Hardy-Littlewood singular series and Zeckendorf carry debt are dual local
   corrections: one is p-adic residue debt, the other is path-carry debt.
3. Helfgott's ternary lift is the additive analogue of adding one geometric
   dimension to turn a brittle edge problem into an averaged surface problem.
4. Fermat polygonal should be used as the bounded-cover sanity check for new
   atom systems: if a proposed atom set has no finite arity cover, it is
   Goldbach-like rather than polygonal-like.
5. Tournament Analysis can import this as a standard additive-basis profile:
   atom set, arity, representation entropy, local product, carry debt, and
   canonical normal form.

## Next Tests

- Scan even `n <= 10^5` for the minimum Zeckendorf carry score among Goldbach
  pairs and compare against Hardy-Littlewood singular-series bonuses.
- Replace Fibonacci normal form with polygonal and Ostrowski normal forms to
  see which basis makes Goldbach pairs lowest-debt.
- Build a "representation tournament" on all Goldbach pairs of a fixed `n`,
  orienting pair `A -> B` when `A` has lower carry debt, lower imbalance, or
  higher local singular contribution.
- Reuse this profile on LRC endpoint rows: representation count corresponds
  to possible repairs, singular product to endpoint debt, and carry debt to
  legal handoff/canonicalization cost.

## External Anchors

- Helfgott, "The ternary Goldbach conjecture is true", arXiv:1312.7748:
  https://arxiv.org/abs/1312.7748
- Hardy and Littlewood, "Some problems of 'Partitio numerorum'; III: On the
  expression of a number as a sum of primes", Acta Mathematica 44 (1923):
  https://link.springer.com/article/10.1007/BF02403921
- Fermat polygonal number theorem, MathWorld:
  https://mathworld.wolfram.com/FermatsPolygonalNumberTheorem.html
- Zeckendorf representation, MathWorld:
  https://mathworld.wolfram.com/ZeckendorfRepresentation.html

## See Also

- `04-computation/additive_basis_normal_forms_s494.py`
- `05-knowledge/results/additive_basis_normal_forms_s494.out`
- `07-reflections/goldbach-additive-basis-normal-forms-s494.md`
- HYP-1902
- HYP-1920
- HYP-1952
