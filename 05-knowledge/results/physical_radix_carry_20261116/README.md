# Exact carry moves in physical tournament repairs

## Main result: an actual low digit and a high reset

[Radix four](RADIX4.md) gives a repeatable construction of the
entire low two-state cycle, then one reset, at arbitrarily large
counts. On one physical strong tournament spine the displayed
head-contact changes take

    (P,1) → (P,3) → (P+1,1)

in the literal readout H=4P+C, not just modulo four. The first
example is 7741→7743→7745. There are Θ(√X) distinct such windows
up to X on the proved fixed-cost subfamily. A second route,
marking just the two spacer vertices, lands at the same count
as the high reset while keeping both outer masks fixed.

This finishes one actual low cycle, not indefinitely many adjacent
high cycles. The high reset arrives at a square product. The note
gives a concrete same-cut barrier to reopening the low cycle at
some such squares, and leaves readers or cuts free to change.

## Two different ways to change the scale

- [Growing radix](VARIABLE_RADIX.md): on a related family the
  physical high product advances by one and the signed low
  remainder falls by 2^t−2, exactly. The full radix needs a
  specified phase (except at its first example). A subsequence
  gives arbitrarily large unshifted canonical radices, and
  symmetric marks inside the seam add a third adjacent count
  at every odd height. This does not supply the whole growing low bank.
- [Two-apex carry square](TWO_APEX.md): compose a fixed two-arc
  swap with a mutual arc flip, in either order. Two growing
  terms cancel to the exact signed residual 2−8ρ for an arbitrary
  binary word encoding ρ. With an empty word this closes to +2
  against a 3^t baseline. This selects a residual, not a target
  baseline; the address itself changes both ports' absolute counts.

[Low-digit budget](BUDGET.md) proves why merely increasing the
spacer, even allowing all its marks on the three exhibited outer
states, cannot fill the larger low bank.

Independent subset DP and direct grouped sums run with
`python3 check_carry.py`; the exact transfer-span certificate runs
with `python3 check_two_apex.py`. The new carry is physical,
but Chan–Pak still needs a way of choosing a viable window from
an arbitrary prescribed odd count.
