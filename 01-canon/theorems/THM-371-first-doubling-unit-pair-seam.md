---
id: THM-371
name: first-doubling-unit-pair-seam
status: PROVED
date: 2026-05-31
session: codex-2026-05-31-S453
depends_on: []
---

# THM-371: The first doubling seam for units and pairs

## Statement

Let `m` be odd.

1. The reduction map

```text
U(2m) -> U(m)
```

on residue units is a bijection.  Its inverse sends a unit `a mod m` to the
unique odd lift among `a` and `a+m`.

2. Consequently,

```text
phi(2m) = phi(m).
```

3. Later dyadic row steps are honest doublings:

```text
phi(2^(r+1)m) = 2 phi(2^r m),  r >= 1.
```

4. The maximum matching count has the complementary first-seam law

```text
floor(2m/2) = 2 floor(m/2) + 1.
```

For later row steps, the matching count doubles.

## Proof

For `m` odd, every residue class `a mod m` has two lifts modulo `2m`, namely
`a` and `a+m`.  Exactly one of these two lifts is odd.  Since `m` is odd,
coprimality to `m` is unchanged by adding `m`, and coprimality to `2m` is
equivalent to being odd and coprime to `m`.  Hence each unit modulo `m` has
exactly one unit lift modulo `2m`, and reduction is a bijection.

This proves `phi(2m)=phi(m)`.

For `r>=1`, the number `2^r m` is even.  Each unit modulo `2^r m` is odd and
coprime to `m`, and both lifts modulo `2^(r+1)m` remain odd and coprime to
`m`.  Thus every unit has exactly two unit lifts, giving
`phi(2^(r+1)m)=2 phi(2^r m)`.

For matching counts, when `m` is odd,

```text
floor(m/2) = (m-1)/2,
floor(2m/2) = m = 2((m-1)/2) + 1.
```

If the starting size is already even, say `n=2^r m` with `r>=1`, then
`floor(2n/2)=n=2 floor(n/2)`, so later row steps double exactly.

## Significance

The first row step `m -> 2m` is not a generic doubling.  LRC unit endpoints do
not double; they become the unique odd-lift copy of the odd root's unit
skeleton.  Tournament pair structure does the opposite-looking but equivalent
thing: the odd root's unmatched vertex gains its first twin, creating one extra
pair.

This is the exact arithmetic core behind the `r=0 -> r=1` seam in the repo's
two-mode tournament and LRC recursions.

## Related

- HYP-1905
- HYP-1881
- HYP-1891
- `07-reflections/adic-column-families.md`
- `04-computation/lrc_tournament_first_doubling_seam_s453.py`
