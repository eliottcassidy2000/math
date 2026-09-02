---
id: THM-4332
title: "LRC(14) fixed-pool single-constraint implication rigidity"
status: >
  PROVED FINITE-EXACT + INDEPENDENT EXACT AUDIT. For the displayed
  thirty-label pool P, its complete 1/14-safe set is contained in the
  singleton safe set G_h exactly when h is already a label of P. Hence no
  subset of this unscaled pool can imply an outside singleton constraint.
  This rules out literal safe-set-implication entry only; it does not rule
  out adaptive, measure, scaled, or typed entry maps and does not prove
  LRC(14).
source: root + entry_scout + implication_auditor / LRC14 continuation session, 2026-09-01
depends_on: []
related:
  - THM-4156-divisor-complete-anchor-pool-haar-odd-tail-transfer
  - THM-4326-lrc14-rank-two-wall-graph-complete-typed-universe-closure
  - THM-4330-lrc14-affine-two-adic-root-types-and-anchored-pool-entry-sieve
script: 04-computation/lrc14_fixed_pool_implication_rigidity_thm4332.py
output: 05-knowledge/results/lrc14_fixed_pool_implication_rigidity_thm4332.out
independent_audit_script: 04-computation/lrc14_fixed_pool_implication_rigidity_thm4332_independent_audit.py
independent_audit_output: 05-knowledge/results/lrc14_fixed_pool_implication_rigidity_thm4332_independent_audit.out
verifier: 04-computation/verify_lrc14_fixed_pool_implication_rigidity_thm4332.py
script_sha256: 38f3969b80977ce534966073f9f203bbe7a8678abc801c16890f501a188b4dd8
output_sha256: 4de9a580755b6af5d518bcdb9f47817ce90a2b0e9ecfee0349e95e142e5ba031
independent_audit_script_sha256: 98edfdeaeddd1ee22b88bec2bda96d9b01d9d4954a2b0b7ccaff8daa9e734778
independent_audit_output_sha256: 9e9a8f6b839a56fc16490cc174001e3ee84158ddf7c4226cfb6115daa96d9241
hash_basis: raw LF bytes
audit: >
  PASS / ACCEPT. The primary combined arrangement and a separate fixed-pool
  interval implementation independently exhaust h=1,...,591 and return
  exactly the thirty pool labels. Normal and optimized independent runs
  agree. A third audit checked the analytic cutoff, every circle wall class,
  strict/closed boundary conventions, component endpoints, subset
  quantifiers, and the exact hostile phases.
---

# THM-4332 -- fixed-pool single-constraint implication rigidity

**PROVED FINITE-EXACT + INDEPENDENT EXACT AUDIT. LRC(14) REMAINS OPEN.**

## 1. Exact statement

Put

```text
P={8,10,15,16,20,30,40,42,60,63,80,84,85,88,95,
   120,126,132,143,145,168,170,176,190,193,240,252,
   264,286,290}.                                       (1)
```

For every finite positive integer set `A`, write

```text
G_A={x in R/Z:||ax||>=1/14 for every a in A},
G_h=G_{ {h} }.
                                                               (2)
```

> **Theorem.** For every positive integer `h`,
>
> ```text
> G_P subset G_h                 iff                 h in P.   (3)
> ```

Equivalently, with the strict danger comb

```text
D_h={x in R/Z:||hx||<1/14},                              (4)
```

one has

```text
D_h subset union_(p in P)D_p       iff       h in P.     (5)
```

Consequently, if `B subset P` and `G_B subset G_h`, then `h in P`. More
generally,

```text
B subset P and G_B subset G_C       imply       C subset P.  (6)
```

The conclusion in `(6)` is only necessary. For `B` properly contained in
`P`, the theorem neither says that every `h in P` is implied nor classifies
the implication relations among the labels of `P`.

## 2. Analytic cutoff

The forward implication in `(3)` for `h in P` is tautological. For the
converse, exact rational-wall reconstruction gives `G_P` exactly `150`
positive-length closed components. Its largest component is

```text
[393/1232,1301/4060],              length L=37/25520.   (7)
```

The singleton safe set `G_h` is the disjoint circle union of the `h` closed
intervals

```text
[(14k+1)/(14h),(14k+13)/(14h)],       0<=k<h,           (8)
```

each of length `6/(7h)`. If `G_P subset G_h`, the connected component in
`(7)` must fit inside one component in `(8)`. Therefore

```text
37/25520 <= 6/(7h)
                    iff 259h<=153120.                   (9)
```

Now

```text
259*591=153069,                  259*592=153328.         (10)
```

Thus every `h>=592` is excluded analytically. The boundary label `591`
belongs to the finite exact universe below.

## 3. Two exact finite paths

Both implementations exhaust every `1<=h<=591`, including all `561`
labels outside `P`, and return precisely the thirty labels in `(1)`.

### 3.1 Combined-arrangement path

For each `h`, the primary inserts the reduced rational walls of `P union
{h}`. It tests every circle wall class, using `0` for the duplicate class
`0=1`, and one midpoint of every open cell. On every atom inside `D_h` it
records the mask of pool labels whose danger combs contain that atom.

Every outside label has an atom with the zero mask. Hence even the union of
all thirty pool danger combs fails to cover `D_h`; the auxiliary nine-label
hitting-set branch is never reached for an outsider. The frozen summary is

```text
IMPLIED exactly P,
NEW_COUNT 0,
UNION_FAIL_COUNT 561,
OVER_CAP_COUNT 0.                                      (11)
```

### 3.2 Fixed-pool interval path

The independent implementation constructs `G_P` once, using only its
`7,134` distinct cut points including `0,1`. It obtains

```text
safe open cells       150,
safe wall points      300,
components            150,
largest length        37/25520.                        (12)
```

It separately verifies that the `300` safe wall points are exactly the
endpoints of the `150` safe cells. It then intersects those cells and points
directly with each strict `D_h`, without using the combined arrangement or
a hitting set. For every outside label it produces an exact point

```text
x_h in G_P intersect D_h.                              (13)
```

Selected hostiles are

```text
h=1,2,3,7: x=199/21280,
h=50:      x=1147/61600,
h=291:     x=15023/337560,
h=400:     x=14219/705600.                              (14)
```

The finite result together with `(9)--(10)` proves `(3)`. Taking complements
proves `(5)`. If `B subset P`, then `G_P subset G_B`; composing with a
hypothesized `G_B subset G_h` proves the first consequence. Apply that
argument to every `h in C` to prove `(6)`. **QED.**

## 4. Boundary and scope

- Safety is closed and danger is strict. Wall points and open cells are
  tested separately; equality is never counted as danger.
- Coincident walls are deduplicated as exact reduced rationals and
  membership is recomputed against every label.
- Every positive speed is dangerous near the circle cut, so no safe
  component wraps through zero.
- This theorem says that an outside singleton constraint is not redundant
  on all of `G_P`. It does not say that `P union {h}` is unsafe, that its Haar
  mass drops below `4/63`, or that `P` is optimal.
- THM-4326 uses retained Haar mass, not safe-set containment, so it is not
  contradicted. Literal implication from an unscaled pool subset cannot be
  its arbitrary-row entry map.
- THM-4330 is now proved and uses a different operation: a positive-rational
  common refactorization of the whole eleven-body, followed by Haar-measure
  invariance and the at-most-two-outsider pool chart. The present theorem
  forbids only unscaled pointwise containment `G_P subset G_h` for an outside
  singleton. Thus THM-4330's projective entry sieve and this literal
  implication no-go are compatible and neither subsumes the other.

## 5. Reproduction

From the repository root:

```text
python 04-computation/verify_lrc14_fixed_pool_implication_rigidity_thm4332.py
```

The verifier reruns the combined path, the independent path, and the
optimized independent path through `h=591`; checks both frozen outputs; and
prints the four raw artifact hashes. No arbitrary-row entry and no instance
of LRC(14) follows.
