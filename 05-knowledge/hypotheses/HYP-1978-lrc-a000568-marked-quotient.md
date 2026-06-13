---
id: HYP-1978
status: OPEN
source: codex-2026-06-01-S507
related:
  - HYP-1900
  - HYP-1901
  - HYP-1951
  - HYP-1967
  - HYP-1968
  - HYP-1970
  - HYP-1971
  - HYP-1972
  - HYP-1973
  - THM-357
  - THM-380
---

# HYP-1978: LRC is a marked A000568-style quotient problem

## Statement

The strongest form of the A000568 analogy is not:

```text
LRC = ordinary tournament isomorphism classes.
```

It is:

```text
LRC counterexamples are bad orbits in a rooted, endpoint-marked tournament
species lying over the ordinary A000568 tournament quotient.
```

The A000568 layer is the unrooted base: binary pair decisions modulo vertex
relabeling.  LRC adds data that ordinary A000568 intentionally forgets:

```text
rooted vertex:       the stationary observer,
safe/end labels:     the anchored 1/n endpoint clock,
pressure labels:     deletion-relief or blocker-debt arcs,
tie Hamiltonian path: increasing speed/order.
```

Thus the right proof object should be a quotient stack:

```text
marked LRC tournament class
  -> rooted tournament class
  -> ordinary unrooted A000568 class.
```

A proof of LRC would show that the bad marked classes are empty.

## Evidence

`lrc_a000568_isoclass_atlas_s507.py` compares four layers.

First, the ordinary quotient scale is exactly the A000568/Burnside one:

```text
n=6: 2^15 labelled tournaments, A000568=56, average orbit 585.14
n=7: 2^21 labelled tournaments, A000568=456, average orbit 4599.02
n=8: 2^28 labelled tournaments, A000568=6880, average orbit 39016.78
```

Second, rooting already changes the state space:

```text
n=3: A000568=2, rooted=4
n=4: A000568=4, rooted=12
n=5: A000568=12, rooted=48
n=6: A000568=56, rooted=296
```

This is the first warning that LRC cannot live on the unrooted quotient alone:
the observer is not fungible.

Third, `H` and score data are useful but coarse:

```text
n=6: A000568=56, H buckets=19, score buckets=22
```

So `H` is a quotient height or loneliness metric, not a complete class
coordinate.

Fourth, the initial LRC half-turn clock is a thin walk in the quotient:

```text
n=7: 24 clock cells visit 7 unrooted classes and 17 rooted classes,
     versus A000568(7)=456.
```

The LRC trajectory is not a random tournament census.  It is a structured path
through a tiny part of the quotient.

Finally, the hard `n=14` and `n=18` ladders show fiber behavior.  After the
first gate, the phase and safe-phase shadows stabilize while endpoint debt
continues to move:

```text
phase_half n14-s14 and n14-s28: same coarse rooted shadow,
  debt 168 -> 336, gap*debt = 5/11.

phase_half n18-s18 and n18-s36: same coarse rooted shadow,
  debt 352 -> 704, gap*debt = 1.
```

This looks exactly like motion inside a marked fiber over a stable tournament
base class.

## Predictions

1. There is a useful Burnside or cycle-index formula for the rooted/marked LRC
   species `(phase tournament, root, safe mask, endpoint labels, pressure
   labels)`.
2. The unrooted `phase_half` class controls global spread through `H`, while
   the rooted/marked fiber controls the actual `1/n` LRC obstruction.
3. Hard dyadic ladders are mostly fiber translations: after the first gate,
   ordinary phase shape stabilizes and endpoint labels carry the recursion.
4. A true counterexample class must satisfy a multi-layer bad predicate:

```text
no positive gap,
no unprotected endpoint,
safe_phase alarm,
pressure_k2 SCC carrying the endpoint core,
and non-harmless phase_H behavior.
```

5. The analogue of the A000568 odd-cycle Burnside filter is an LRC
   fixed-point-killing lemma: many would-be bad marked classes have no fixed
   representative because endpoint protection forces incompatible
   divisibility/gate labels.

## Proof Program

1. Define the marked species

```text
M_n = (T_phase, root 0, safe mask, endpoint-owner labels, pressure relation).
```

2. Quotient by permutations fixing the root.  Ordinary A000568 is recovered by
   forgetting all markings and then forgetting the root.
3. Define a bad predicate on `M_n` corresponding to THM-357 and THM-380:

```text
full open cover,
all endpoints protected,
terminal endpoint core nonempty,
pressure-realized SCC.
```

4. Derive a colored/rooted Burnside count or at least a fixed-point filter for
   bad classes.
5. Prove fixed-point killing or endpoint-private leaf peeling for every
   possible bad orbit.
6. Use the metric vector from HYP-1973 as the search fingerprint while the
   proof tracks exact endpoint labels.

## Sources

- `04-computation/lrc_a000568_isoclass_atlas_s507.py`
- `05-knowledge/results/lrc_a000568_isoclass_atlas_s507.out`
- `07-reflections/lrc-a000568-isoclass-analogy-s507.md`
- HYP-1900
- HYP-1901
- HYP-1970
- HYP-1971
- HYP-1972
- HYP-1973
- THM-357
- THM-380
