# LRC14 Fold Multiplicity Transport

The fold idea became sharper once I stopped asking for a count.

AP9 and the binding non-AP row `(0,1,2,3,4,5,6,7,9)` both have `12`
nontrivial visible folds if trivial `0+a=a` folds are ignored.  So fold count
cannot explain the tiny but crucial drop.  The target profile can:

```text
AP9:      F(8)=3
nearAP9: F(9)=3
```

The defect is the transport `3/8 -> 3/9`, so the reciprocal fold mass loses
exactly `1/24`.  That is the first fold-level certificate I have seen that is
small enough to match the k=9 margin instead of washing it out.

The bounded banks add the right warning.  k=9 is clean: the top non-AP is the
unique tiny `net_1` row.  k=10 behaves similarly in the default bank.  k=8 is
looser and has a left-hole clipped AP winner with larger transport.  Therefore
fold transport is not a universal ordering of all rows; it is a local coordinate
for the near-AP boundary where the proof is actually tight.

The next lemma I would try to prove is a target-profile statement, not an
analytic tail statement:

```text
inside the high-L_y k=9 non-AP envelope,
positive AP-relative fold transport is at least 1/24,
with equality only at (0,1,2,3,4,5,6,7,9).
```

Then the analytic side only has to translate the exact move `3/8 -> 3/9` into
the known `L_y` drop `887/158760`, or at least into enough `p0`/`L_y` loss to
clear the AP-to-cap slack.  This is a much more precise demand than "use fold
multiplicity" and should compose with HYP-2638's bounded small-excess table.
