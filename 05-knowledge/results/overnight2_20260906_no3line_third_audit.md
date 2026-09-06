# Independent audit of the third collinear-triple factorial moment

**Status: INDEPENDENT AUDIT PASS.** The producer's conclusion remains
FINITE-EXACT at `n=8`; the contrast reduction is proved for the stated
short-edge model. No errors or source corrections were found.

Reviewed sources and result:

- [Native enumeration](../../04-computation/overnight2_20260906_no3line_third.cpp)
- [Rational weighting and controls](../../04-computation/overnight2_20260906_no3line_third.py)
- [Mathematical report](overnight2_20260906_no3line_third.md)
- [Frozen certificate](overnight2_20260906_no3line_third_certificates.json)

## Analytic and source audit

Each line event is a matching of three edges. Thus the union of two has
degree at most two, and adding a third matching preserves that bound
exactly when every new edge avoids saturated vertices. Previously present
edges must remain allowed. The native fast filter implements precisely this
predicate, while the literal mode independently rebuilds degrees. The
reported second pass is correctly described as a complete filter replay;
the two modes share a component-signature routine.

For connected degree-at-most-two bipartite graphs, edge count and shore
sizes distinguish odd paths, even paths with either shore larger, and
even cycles. Thus the packed component token is injective in the declared
universe. At most nine edges imply tokens below 32 and at most nine
components, so five-bit packing into a 64-bit integer is safe. The native
grid mask uses exactly the same row stride eight as the independent
determinant construction. Distinct rows and nonzero slope guarantee
distinct columns, excluding axis triples without losing any eligible event.

The skeleton counts are non-induced edge-subset counts, with no isolated
vertices. Dividing injective embeddings by shore-preserving automorphisms
therefore gives the correct containment probability. The factors are `L`
for an `L`-cycle, two for an even-edge path, one for an odd-edge path, and
the factorial of every repeated-component multiplicity. Extra skeleton
edges between incident vertices do not invalidate a selected copy.

Every unordered set of three distinct events gives exactly six ordered
tuples, even if events share edges. The factor six multiplies geometric
event multiplicity, rather than just the number of distinct union sets.

The inherited weighted short-cycle theorem allows only
`1,c4,c6,c8,c4^2` through nine edges. The two contrasts are consequently

```text
E=(M3(2C8)-M3(C16))/2,
F=(M3(4C4)-4M3(C4+C12)+3M3(C16))/12.
```

They annihilate every profile through seven edges. Restricting the native
census to eight- and nine-edge unions therefore recovers the complete
coefficients `E,F` of the full moment. It need not recover its remaining
three coefficients, and the report correctly says so. All seven partitions
of eight into parts at least two are included in the profile-law check.

## Independent exact reconstruction

I did not import any producer functions. For a single cycle, I enumerated
edge masks and decomposed each proper mask into cyclic runs. A run's length
and starting-edge parity determine its path type and larger shore. Full
masks give cycles. Convolution of these component banks gives the counts
for every skeleton. This replaces the producer's arbitrary-subset adjacency
traversal by a different exact mechanism.

All **1,050 counts** (150 retained geometric types in seven skeletons)
agreed. Recomputing the automorphism denominators and applying the raw
contrasts directly recovered:

```text
c8 coefficient = 172483/529200,
c4^2 coefficient = 11881/50400.

Positive / negative contributions:
c8:   25988161/2116800, -8432743/705600;
c4^2: 5120467/705600,   -4954133/705600.

Eight-edge contributions: 768463/2116800, 456371/2116800.
Nine-edge contributions: -26177/705600, 42631/2116800.

Raw contrasts: 172483/264600, 11881/4200.
```

In particular the aggregate signed cancellations are retained; a positive
individual profile is not used as a surrogate for the final coefficients.

I also independently enumerated all 576 row/column labelings of both
size-four skeletons, using literal integer collinearity determinants. The
histograms and first three factorial moments were:

```text
C8 histogram:   {0:16, 1:192, 2:272, 4:64, 5:32};
moments:       2, 61/18, 6.

2C4 histogram: {0:288, 2:128, 4:64, 6:64, 8:32};
moments:       2, 74/9, 104/3.
```

These independently confirm the ordered-event factor and label-copy
normalization. I reviewed the complete native loops and the matching
normal/optimized transcripts, but did not run a third full size-eight
geometry census. The geometric multiplicities in the frozen certificate
are the inputs to the independent reconstruction above.

## Reproducible cyclic-run audit

The following self-contained code can be copied into Python and run from
the repository root. It imports no repository module and writes no file.

```python
from collections import Counter
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import json

data = json.loads(Path(
    "05-knowledge/results/overnight2_20260906_no3line_third_certificates.json"
).read_text())

def need(ok):
    if not ok:
        raise RuntimeError("independent cyclic-run audit failed")

def cycle_bank(half):
    L, bank = 2*half, Counter()
    for mask in range(1 << L):
        m = mask.bit_count()
        if m > 9:
            continue
        if not mask:
            sig = ()
        elif m == L:
            sig = (3*L,)
        else:
            tokens = []
            for i in range(L):
                if mask >> i & 1 and not (mask >> ((i-1) % L) & 1):
                    q = 1
                    while mask >> ((i+q) % L) & 1:
                        q += 1
                    tokens.append(3*q + (0 if q % 2 else 2 if i % 2 == 0 else 1))
            sig = tuple(sorted(tokens))
        bank[sig] += 1
    return bank

cache = {h: cycle_bank(h) for h in range(2, 9)}

def bank(parts):
    result = Counter({(): 1})
    for h in parts:
        after = Counter()
        for a, x in result.items():
            for b, y in cache[h].items():
                sig = tuple(sorted(a+b))
                if sum(t//3 for t in sig) <= 9:
                    after[sig] += x*y
        result = after
    return result

parts = ((8,), (4,4), (3,5), (2,6), (2,3,3), (2,2,4), (2,2,2,2))
banks = {part: bank(part) for part in parts}
totals, checks = [Q(), Q()], 0
for row in data["profiles"]:
    signature = [tuple(q) for q in row["signature"]]
    sig = tuple(sorted(3*m + (1 if l > r else 2 if r > l else 0)
                       for m, l, r, cycle in signature))
    for part in parts:
        need(banks[part][sig] == row["skeleton_copies"][str(part)])
        checks += 1
    nr = sum(q[1] for q in signature)
    nc = sum(q[2] for q in signature)
    aut = 1
    for (m, l, r, cycle), multiplicity in Counter(signature).items():
        one = m if cycle else 2 if m % 2 == 0 else 1
        aut *= one**multiplicity * factorial(multiplicity)
    ordered = factorial(8)//factorial(8-nr) * (factorial(8)//factorial(8-nc))
    need(ordered % aut == 0)
    denominator = ordered//aut
    need(denominator == row["grid_copies"])
    E = Q(banks[(4,4)][sig] - banks[(8,)][sig], 2)
    F = Q(banks[(2,2,2,2)][sig] - 4*banks[(2,6)][sig] + 3*banks[(8,)][sig], 12)
    contributions = [6*Q(row["unordered_event_triples"])*v/denominator for v in (E,F)]
    need(list(map(str, contributions)) == row["ordered_weighted_contributions_E_F"])
    totals = [x+y for x,y in zip(totals, contributions)]
need(checks == 1050)
need(totals == [Q(172483,529200), Q(11881,50400)])
print("PASS", checks, "independent copy counts; totals", list(map(str, totals)))
```

The audited producer SHA-256 values, after CRLF-to-LF normalization, match
its report exactly:

```text
C++:          bc6b901088573ab2642cbca8b7cf4856110090869a18a5d6940d2110430ce83e
Python:       4a00e332070e7415cc6fe02999f226a1072dbcbbfbaac95120dceed706eb293e
certificate:  f2c566ac1b2bcedb530af72f3a290479841e7db87844b331558bed68c93ba727
both outputs: a8c1c4a7ecf28f19699f1b2700b73b72f9ccf0a893b53f0c0e3d1df73f708bc1
```

Acceptance is scoped to the fixed `n=8` result and the proved reduction.
It supplies no all-size sign law, zero-defect probability, or independence
claim. Only this audit document was written; producer files were unchanged.
