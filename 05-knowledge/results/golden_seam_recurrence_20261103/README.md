# A three coordinate Hamilton recurrence and its prime lanes

**PROVED all-size recurrences for the specified finite atom alphabets;
FINITE-EXACT bounded-value coverage and exclusion certificates. Same-agent
review with independent algorithms. Global spectrum completeness remains
OPEN.** This is a uniform construction family, not a new unrestricted
enumeration record or an extension of A000568's known coefficients.

The golden complement family now has a nonnegative integer3-by-3 recurrence
that retains each of its five Q4 presentations. If h(w) is the Hamilton-path
count after complementing the source word w, and a is a singleton source
component, its state is

    v(w) = (h(w), (h(aw)-h(w))/2, (h(aaw)-h(aw))/2).
    v(empty) = (1,0,0),       v(cw) = v(w) P_c.

All seven matrices P_c are displayed in [THEORY.md](THEORY.md). The recurrence
is linear in the number of source components in arithmetic operations.
Its smallest rational linear dimension is exactly3 for left attachments;
both-end attachments require exactly9, retaining the word's full matrix.
These minima are scoped to Hamilton observation and these operations, not
the full complement/re-marking metagraph or arbitrary joint ports.

## What the family reaches

The first coordinate never decreases under extension. An invariant cone
bounds the other two coordinates when h is bounded, so finite reachability
certificates cover **every word length and tournament order** at that value
bound. These are not searches cut off at a prescribed order.

| Source alphabet | Value bound | Odd values attained | Prime lane attained | 7p lane attained |
|---|---:|---:|---:|---:|
| Golden seven symbols |100000|30194 /50000|4907 /9590|851 /1674|
| Golden plus one marked five-vertex atom |100000|35661 /50000|6000 /9590|1042 /1674|
| All marked strong atoms through five vertices |10000|4990 /5000|1222 /1227|222 /223|
| Previous row plus eight marked six-vertex atoms |30000|14998 /15000|3243 /3243|586 /586|

Prime lane means odd prime p!=7. The7p lane means odd prime p!=3 with7p at
most the value bound. Every reported lane witness is checked to be strong.
The fourth row has exactly the two missing odd values7 and21. The inherited
global exclusion theorem agrees with these computations.

There are decisive limits and repairs:

- **Golden-only obstruction, every size:**37 and133=7*19 are absent. They
  are its first missing required prime and7p values. The small finite
  closed-state certificates are retained in [exclusions.json](exclusions.json).
- A single source atom r13, with labelled out-neighbor rows
  `(22,12,8,17,6)`, repairs both: `r13 a` has h=37, and `q0 r13` has h=133.
  The atom itself has H=11; its complemented presentation has H=13.
- All50 strong five-vertex charts give43 distinct counting matrices. This
  expanded family still misses47,83,103,491,709 in the prime lane and
  3227=7*461 in the7p lane, at **every size**.
- All903 strong six-vertex charts compile to544 distinct matrices. Eight
  specified charts suffice, together with the previous alphabet, to fill
  every remaining odd gap through10000. The same fixed65-symbol rule then
  passes a holdout extension through30000 without adding another chart.

[COVERAGE.md](COVERAGE.md) records the actual source words, guards and next
inverse problem. [THEORY.md](THEORY.md) also gives cubic recurrences for
every repeated word and an exact minimal-order12 recurrence for total
Hamilton mass over all golden words of each vertex size. That mass counts
marked presentations with Hamilton weights, not ordinary isomorphism classes.
It also proves that fixed left/right word wrappers cannot implement a
nonzero constant affine offset on the whole bank, excluding those wrappers
as uniform3H+2 or3H+4 gadgets for the incoming ternary induction target.

## Why the carrier needs its operation scope

`ba` and `q5` both have state `(5,4,5)`, so every left word context has the
same Hamilton response. Appending a singleton on the right gives9 and13.
The full matrix separates them. Its two-sided Hankel minor has determinant
of absolute value8192 and proves the minimum dimension9.

For weighted source banks modulo m, the smallest left-context observer is

    Z/m x Z/(m/gcd(2,m)) x Z/(m/gcd(2,m)).

The factors of2 are exact: future observations see H,2x,2y. Individual
physical words have additional realization constraints. A positive odd
coordinate vector such as(37,1,1) is not a realizable golden word merely
because it satisfies parity, primitivity and the invariant cone.

## Recovery and the requested numbers

The closest mechanism is the prior golden SCC codec and actual complement
construction, followed by the ordered-component expiry transfer. The hostile
is now the all-size golden gap37; the corrected near miss is assuming that
small recurrence dimension implies all odd outputs. The least-used sidecar
is which end of the source word is allowed to grow.

Anchor: prime and7p realization. Niche: the sharp3-versus9 context memory.
Wildcard: the earlier111=56+55 splittings. The live board is source words /
physical seams / endpoint pairing / matrix semigroup / value support /
operation-indexed memory.

[RECOVERY.md](RECOVERY.md) recovers the actual meanings of111,55,56:
six-vertex complement memory is56 visible plus55 hidden dimensions; the
square-sum graph at56 has111 edges, split into55 Hamilton-tree edges and56
chords; Fibonacci55 has a genuine terminal-flag codec. The two rank
decompositions are related by their kernel/image mechanism, without a
claimed natural isomorphism between their spaces. Golden words attaining
55 and111 are supplied;37 is absent although111=3*37 is attained.

Reproduce with Python's standard library:

```powershell
python -B -X utf8 05-knowledge/results/golden_seam_recurrence_20261103/verify.py
python -B -O -X utf8 05-knowledge/results/golden_seam_recurrence_20261103/verify.py
```

[AUDIT.md](AUDIT.md) specifies the independent paths, exact universes and
source hashes. The complete proof still needs an all-value reachability
theorem, or another uniform source construction, in the prime and7p lanes.
