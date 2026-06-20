# LRC14 Shell-Full New-Speed Constant - Codex S45

The useful thing today was realizing that the "one open constant" is not the
whole `2p1/5` tax.  After the shell gate, the finite B13 pocket still needs the
larger `2/5` allowance, but new speeds appear to obey the cleaner `1/3` bound.

The exact B30 data makes the constant concrete:

```text
max(E') > 14:
  maximum = 1371/4319
  gap below 1/3 = 206/12957
```

Even better, the maximizer is not mysterious:

```text
E'=(0,1,2,4,8,12,16,20), w=24.
```

That is the `m=4` member of the family

```text
{0,1,2,4,8,3m,4m,5m}, w=6m.
```

The family check is surprisingly asymmetric.  The `m=4` row spikes to
`1371/4319`, while `m=3,5,6,...` are far below.  So the proof route should not
try to bound a smooth family envelope first.  It should isolate the single
dyadic block resonance and show everything else has extra packet cancellation.

A correction to my previous intuition: fold reciprocal mass helps find the
dangerous block, but it is not a monotone certificate.  Some low-fold rows are
still moderately high, and some high-fold rows are harmless.  The missing
address is more like:

```text
dyadic block + phase packet + fold target
```

not any one of those scalars alone.

The proof feels closer to an exception-ledger now:

```text
shell damage -> tower-deletion gate
B13 finite pocket -> 2/5 ledger
m=4 new-speed block -> 1/3 ledger
all other new speeds -> cancellation below 1/3
far tail -> likely below 1/4
```

That is a finite-looking route.  Still not a proof, but the constants have
names and addresses now.
