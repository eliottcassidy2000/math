# LRC n=14 `Res_27` Fixed Bridge, S609

This session pushed on the n=14 proof by making the bridge between the 64
self-converse round classes and the `C=27` owner layer more explicit.

The first result is a sanity check with teeth: the exact Hamiltonian-path count
for every one of the 64 fixed classes is odd.  That verifies HYP-2160 on the
actual n=14 table.  But the scalar counts range from `1` to `3711175`, and the
regularity gauge disagrees with the HP-scarcity gauge on `2005/2016` pair
comparisons.  So the solved-face transfer is parity, not a scalar monotone.

The second result is the real n=14 improvement.  S578 scanned canonical
unit-spine slack through `42`; S609 extends the canonical slack fibre through
`81`, one full `3C` lift range for `C=27`.  Among `17550` slack rows there are
`9506` full D/U/N covers.  Of those, `9504` have an unblocked small pair, the
two no-cheap block-all rows are positive-measure controls, and there are zero
open residuals.  The only measure-zero full covers are still AP
`(3,6,9,12)` and `V*=(3,6,9,24)`, both discharged by `(1,13)` at `1/14`.

This sharpens the bridge lemma:

```text
64 fixed classes carry parity.
C=27 owner fibres carry certificates.
The proof must connect them without forgetting owner labels.
```

Incoming S606/HYP-2162 gives the sibling shell-orbit compression: the raw
`C=27` shell layer folds from `13` shells to the three gcd strata `{1,3,9}`.
That is why the S609 owner scan records both raw shell and gcd signatures; the
proof bridge should not carry thirteen independent shell obligations after the
twisted-involution normalization.

The newer upstream HYP-2163 efficient-pipeline result should be read as a
stronger formal route, and upstream HYP-2164 supplies a sharper pinch
certificate.  S609 is the compatibility layer showing that the old 64-class
parity scaffold and the `C=27` owner certificates line up with those routes
rather than competing with them.

The next proof move is not "read loneliness off the 64 class table."  It is:
prove every tight-boundary realization lifts into the `C=27` owner layer where
the cheap-pair/positive-measure dichotomy applies, or prove that failed lifting
already exposes one of those certificates.

Tournament Analysis stayed deliberately transitive in this pass.  The vertices
were fixed classes and slack rows, not runners; the switch was a certificate
ledger.  A nontrivial SCC should only appear after arbitrary unit
representatives or endpoint-owner fibres are reattached.

Artifacts:

- `04-computation/lrc_n14_res27_fixed_bridge_s609.py`
- `05-knowledge/results/lrc_n14_res27_fixed_bridge_s609.out`
- `05-knowledge/hypotheses/HYP-2165-lrc-n14-res27-fixed-class-bridge.md`
