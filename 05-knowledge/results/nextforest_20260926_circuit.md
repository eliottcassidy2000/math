# Six edges rule out every fixed additive binary-digit rank

**PROVED elementary theorem / FINITE-EXACT controls / Collatz OPEN.**
No claim of literature priority. The shortcut map is
`T(n)=(3n+1)/2` for odd n and `n/2` for even n.

## 1. A precise recovered connection

The previous [forest carrier](forest_20260926_excursions.md) excluded positive
combinations of digit-polynomial evaluations. That left arbitrary signed
weights and exceptions below a verification threshold. The older GMC result
[THM-3258, depth-two affine Farkas clutch](../../01-canon/theorems/THM-3258-depth-two-affine-farkas-clutch-and-complete-reset-distance-gauge-no-go.md)
suggests a stronger operation: search for a positive null combination of
the actual required inequalities. Its GMC theorem is not imported into
Collatz; its elementary certificate mechanism is reconstructed below.

Source: Collatz edges and their binary digit incidence vectors.
Map: `n->T(n)` becomes `d(n)=bits(T(n))-bits(n)` in the direct sum of
integer coordinate lines. Preserved: the change in every fixed additive
digit functional. Lost: which integer owns a bit, and which edges can follow
one another. Sidecar: actual source and target labels. Cheapest hostile:
check both the vector cancellation and the uncancelled integer-vertex boundary.

Unlike a superficial resemblance between small graphs, this map transports
an exact impossibility certificate. It does not transport an orbit cycle.

## 2. The six-edge certificate

For arbitrary real weights `w_0,w_1,...`, put

    V(n)=sum_s b_s(n) w_s.

Each sum is finite. No positivity, growth, computability, or lower-bound
assumption on the weights is needed. Write `Delta_n=V(T(n))-V(n)`.

| n -> T(n) | Delta_n | positive multiplicity |
|---|---|---:|
| 3 -> 5 | w_2-w_1 | 2 |
| 4 -> 2 | w_1-w_2 | 1 |
| 5 -> 8 | w_3-w_2-w_0 | 1 |
| 6 -> 3 | w_0-w_2 | 2 |
| 8 -> 4 | w_2-w_3 | 1 |
| 9 -> 14 | w_2+w_1-w_0 | 1 |

Direct addition gives the polynomial identity

    2 Delta_3+Delta_4+Delta_5+2 Delta_6+Delta_8+Delta_9=0.    (1)

If all six changes are nonpositive, every one must be zero. The six rows
have rank four; directly, the first and fifth rows give `w_1=w_2=w_3`,
the third gives `w_0=0`, and the fourth gives `w_2=0`.
Thus `w_0=w_1=w_2=w_3=0`.

For k>=2 the two ordinary edges

    2^k -> 2^(k-1),
    2^k-1 -> 3*2^(k-1)-1

have changes `w_(k-1)-w_k` and `w_k-w_(k-1)`. Monotonicity makes them
both zero. Consequently every weight vanishes.

**Theorem A.** If V(T(n))<=V(n) for every integer n>2, then V is identically
zero. Adding an arbitrary constant changes this to a constant potential.

The certificate does not contain the trivial 1<->2 cycle. Its cancellation
at actual integer vertices leaves

    {6,6,9} -> {2,5,14}.                                (2)

These multisets have identical aggregate binary digits. They are different
multisets and are not one legal orbit. In particular, their squared-value
sums differ by72. The obstruction belongs to additive digit observables,
not arbitrary functions of an integer or of its full digit word.

## 3. Finite exceptional sets do not repair the class

**Theorem B.** For any finite threshold N, if V(T(n))<=V(n) for every
integer n>=N, every weight still vanishes.

First apply the opposite-current pair above for all sufficiently large k.
It follows that `w_s=c` for every s>=K, for some finite K and real c.

Now fix any positive integer n, including a small n below N. Choose M
large enough that M-1>=K, `n,T(n)<2^(M-1)`, and the following padded source
is above N. If n is even, use

    x=2^M+n,           T(x)=2^(M-1)+T(n).                (3)

Its extra high contribution is c on both sides. If n is odd, use

    x=3*2^M+n,         T(x)=9*2^(M-1)+T(n).              (4)

The high binary blocks3 and9 each have two bits, all at positions>=K,
so their extra contributions are both2c. The size condition prevents any
overlap with the low blocks. The assumed inequality at x therefore gives
exactly `V(T(n))<=V(n)`. This transfers every small edge to the hypothesis.
Theorem A now applies. This is a proof for every N, not an extrapolation
from a finite census.

## 4. What the obstruction changes

The full digit word and its polynomial are lossless. Nevertheless, every
fixed linear functional of its coefficient vector fails as a nonconstant
eventual every-step rank. Thus adding more signed digit-position weights
cannot solve the problem. This covers arbitrary infinite weight sequences,
not merely a fixed-dimensional vector of statistics.

The theorem does **not** exclude:

- nonlinear functions of digit patterns or of exact carry streams;
- a bounded or unbounded non-additive correction to V;
- inequalities required only at an effectively selected subsequence;
- multi-step certificates that retain their original source and boundary.

This aligns with [THM-2163, radix relation-carry descent](../../01-canon/theorems/THM-2163-radix-relation-carry-descent.md):
a finite carry rule needs its actual terminal/source data. Here aggregate
digit conservation is exact, but discards ownership across distinct sources.
There is no transfer of that LRC theorem's target predicate to Collatz.

The constructive next move is to include interactions and an actual selected
return, and run the same positive-current test on that enlarged feature map.
A null certificate then diagnoses exactly which interaction remains absent;
failure to find one is only a search result, not a proof of a rank.

## 5. Reproduction

Run `python -B 04-computation/experiments/nextforest_20260926_circuit.py`
and the same command with `-O`. The assertion-independent script checks
the exact positive certificate and rational rank,255 scale pairs through
256 bits,12,285 padded ordinary edges, and456 template lifts. It also
checks the actual nonzero integer-graph boundary and nonlinear hostile.
The proof above supplies the all-height and all-threshold quantifiers.

[Script](../../04-computation/experiments/nextforest_20260926_circuit.py)
and [retained output](nextforest_20260926_circuit.out).
