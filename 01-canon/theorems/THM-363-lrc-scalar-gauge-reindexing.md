---
id: THM-363
name: lrc-scalar-gauge-reindexing
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S367
depends_on:
  - THM-358
---

# THM-363: LRC Scalar-Gauge Reindexing

## Statement

Fix integers `n>=2` and `k>=1`.  For an open cell `C` in the discontinuity
arrangement of

```text
alpha -> (floor(n {i alpha}))_{i=1..k},
```

write its constant floor vector as

```text
B_C(i) = floor(n {i alpha}),  alpha in C.
```

For `s,m mod n`, translation

```text
T_{s,m}(alpha) = alpha + s m / n  mod 1
```

permutes the open cells of the arrangement, and for every cell `C`,

```text
B_{T_{s,m}C}(i) == B_C(i) + s m i    mod n.
```

Consequently, if residue vectors are gauged by

```text
v'_i = v_i + m i  mod n,
```

then for each fixed `s`, the set of cells blocked by `v'` is the translate
under `T_{s,m}^{-1}` of the set of cells blocked by `v`.  In particular, scalar
addition preserves the number of blocked `(s, cell)` candidates in every full
micro-staircase cell system.

When `k=n-1`, THM-358 supplies that the zero vector blocks every open cell.
Therefore every scalar ramp `v_i=m i mod n` also blocks every open cell.

## Proof

The breakpoints of the `i`th coordinate are the rationals

```text
a/(n i),  a=0,...,n i.
```

If `alpha=a/(n i)`, then

```text
alpha + s m/n = (a + s m i)/(n i),
```

so translation by `s m/n` sends breakpoints for coordinate `i` to breakpoints
for coordinate `i`.  Hence it permutes the complement of the full breakpoint
set, i.e. the open cells of the common arrangement.

Let `C` be an open cell and choose `alpha in C`.  Put

```text
x = {i alpha},       c = s m i.
```

Then

```text
{i(alpha+s m/n)} = {x + c/n}.
```

Write `floor(n x)=q`.  Since adding `c/n` adds `c` to `n x` before reducing
modulo `n`, the integer floor changes by `c` modulo `n`:

```text
floor(n {x+c/n}) == q+c  mod n.
```

Thus

```text
B_{T_{s,m}C}(i) == B_C(i) + s m i  mod n
```

for every `i`.

Now compare the blocking condition for the gauged vector `v'` on `C`:

```text
s v'_i + B_C(i)
  == s v_i + s m i + B_C(i)
  == s v_i + B_{T_{s,m}C}(i)   mod n.
```

Therefore the forbidden-residue test

```text
s v'_i + B_C(i) in {0,n-1}
```

holds for some `i` exactly when the corresponding test for `v` holds on the
translated cell `T_{s,m}C`.  Translation is a permutation of cells, so blocked
candidate counts are preserved.

For `k=n-1`, the zero vector is the initial-segment Dirichlet equality
skeleton from THM-358: every open cell has some coordinate with
`B_C(i) in {0,n-1}`.  Applying the gauge reindexing to the zero vector gives
the full scalar-ramp family.

## Verification Record

`04-computation/lonely_runner_k13_scalar_gauge_s367.py` independently checks
the theorem in the `n=14,k=13` cell system by comparing raw coverage counts for
`200` random vectors and all `14` scalar shifts; it finds zero failures.  The
stored output is
`05-knowledge/results/lonely_runner_k13_scalar_gauge_s367.out`.

## Related

- THM-358: LRC initial-segment unit skeleton.
- HYP-1817: fourteen-runner micro-staircase lift.
- HYP-1818: scalar-ramp excision before composite micro-staircases.
- HYP-1823: fourteen-runner scalar-gauge quotient lemma.
- `07-reflections/lonely-runner-fourteen-runner-scalar-gauge-s367.md`.
