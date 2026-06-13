---
id: THM-396
name: lrc-n14-small-pinch-universal-blocker
status: PROVED
date: 2026-06-02
depends_on:
  - HYP-2059
  - HYP-2060
---

# THM-396: A universal blocker of a small n=14 pinch is a sum-multiple shield

## Statement

Let `a,b` be positive speeds, set `D=a+b`, `g=gcd(a,b)`, and

```text
s = D/g.
```

Assume `s <= 14`.  Consider the pinch times

```text
t = m/D,  1 <= m < D,
```

at which the pair `(a,b)` itself is n=14-safe:

```text
||a m/D|| >= 1/14  and  ||b m/D|| >= 1/14.
```

If a third speed `c` is dangerous at every such pair-safe pinch time,

```text
||c m/D|| < 1/14
```

for all pair-safe `m`, then

```text
D | c.
```

Equivalently, the only single runner that can universally kill all safe
pinches of a reduced-sum-`<=14` pair is a multiple of the pair sum `a+b`.

## Proof

Write

```text
a = g alpha,    b = g beta,    alpha + beta = s,
```

with `gcd(alpha,beta)=1`.  Then `alpha` and `beta` are both units modulo `s`,
and `beta == -alpha (mod s)`.

At `t=m/D`,

```text
a t = alpha m / s,
b t = beta  m / s == -alpha m / s     (mod 1).
```

Hence the two pair distances are equal.  Since multiplication by `alpha` is a
permutation modulo `s`, and since `s <= 14`, the pair is n=14-safe exactly when

```text
m not == 0 (mod s).
```

Indeed, every nonzero residue modulo `s` has circular distance at least `1/s`,
and `1/s >= 1/14`.

Now suppose `c` is not divisible by `D`.  Put

```text
h = gcd(c,D),    q = D/h,    u = c/h.
```

Then `q >= 2` and `u` is a unit modulo `q`.  A residue `j mod q` is n=14-safe
from zero when

```text
min(j, q-j) / q >= 1/14.
```

Let `J` be the set of such safe nonzero residues modulo `q`.  For `q < 14`,
every nonzero residue lies in `J`.  For `q >= 14`, the unsafe nonzero residues
are only those with distance `d` satisfying `14d < q`, so

```text
|J| = (q-1) - 2 floor((q-1)/14) > floor((q-1)/2).
```

On the other hand, among residues `m=1,...,q-1`, at most `floor((q-1)/s)` are
multiples of `s`, and this is at most `floor((q-1)/2)` because `s >= 2`.

Multiplication by the unit `u` permutes residues modulo `q`, so the preimage
`u^{-1}J` has size `|J|`.  Since `|J|` is larger than the number of positive
multiples of `s` modulo `q`, there exists some `m` with

```text
m not == 0 (mod s)
```

and

```text
u m mod q in J.
```

Lifting this `m` to the original denominator gives a pair-safe pinch time
`m/D`, but

```text
||c m/D|| = ||u m/q|| >= 1/14,
```

contradicting universal danger.  Therefore `D | c`.

## Consequence for HYP-2060

The HYP-2059 target "find a reduced-sum-`<=14` clearing pinch" is too strong in
general.  However, THM-396 makes one obstruction exact: if a small pinch is
killed by one runner for all pair-safe residues, that runner must be a
sum-multiple shield.

The remaining proof residual is therefore not an unconstrained search over
times.  It is a finite blocker-cover problem: either a small pair has a
sum-multiple shield, or several non-shield runners must collectively cover the
pair-safe residue classes.

## Verification

`04-computation/lrc_n14_small_pinch_shield_s558.py` verifies the residue lemma
through pair sums `D <= 500`, records exact N3 counterexamples, and exhibits a
collective non-shield cover showing why the theorem cannot be strengthened to
"every failed small pinch has a universal shield."

## Independent verification (monad-reviewer-2026-06-02)

PROOF RE-DERIVED AND CONFIRMED. The argument is a clean elementary counting
proof; I checked each step from definitions:

1. `a t = αm/s`, `b t ≡ -αm/s (mod 1)` since `β ≡ -α (mod s)`, so the two pair
   distances are equal. Pair-safety `⟺ m ≢ 0 (mod s)` because every nonzero
   residue mod `s` has circular distance `≥ 1/s ≥ 1/14` (uses `s ≤ 14`). ✓
2. With `c` not divisible by `D`: `||c·m/D|| = ||u·m/q||`, `q = D/gcd(c,D) ≥ 2`,
   `u` a unit mod `q`. ✓
3. The two conditions are integer conditions on `m ∈ {1,…,q-1}` (where
   `m mod q = m`, and `m < q ≤ D` is a valid pinch index): c-safe `⟺ um mod q ∈ J`
   (count `|J|`), pair-safe `⟺ m ≢ 0 (mod s)` (at most `⌊(q-1)/s⌋ ≤ ⌊(q-1)/2⌋`
   excluded). ✓
4. KEY INEQUALITY `|J| > ⌊(q-1)/2⌋` for `q ≥ 14`, with `|J| = (q-1) - 2⌊(q-1)/14⌋`:
   independently verified by brute force for all `q = 2..2000` (0 failures), and
   provable since `2⌊(q-1)/14⌋ ≤ (q-1)/7 < ⌈(q-1)/2⌉`. The pigeonhole then forces
   an `m` that is both pair-safe and c-safe, contradicting universal danger. ✓

No pitfalls from MISTAKES.md apply (no µ-computation, no chained mod-arithmetic
overflow, no small-n-only coincidence — the `s ≤ 14` and `q ≥ 14` regimes are
both exercised). The theorem statement and its honest scope (single-blocker only;
collective covers remain open) are accurate. **Status confirmed PROVED.**
