# LRC14 Two-Far Live-Depth Kernel - codex S67

The useful overnight step was not another broad finite scan.  It was changing
the two-far deviation target from a seven-depth missed-sector problem to a
four-depth signed kernel problem.

For the HYP-2701 survival currency

```text
C = p1+p2+p3+p4-4p6
```

the deterministic depth coefficients are

```text
C(0)=0, C(1)=C(2)=C(3)=C(4)=1, C(5)=0, C(6)=-4.
```

With two far colors, a before-state of missed depth `t` can only move to
`t`, `t-1`, or `t-2`.  Hence depths `3` and `4` are exactly silent: their
possible after-depths stay inside `{1,2,3,4}`, where the currency is constant.
The actual hit law and the decorrelated iid death-chain law both contribute
`1`, regardless of the relation between the two far speeds.

This turns the HYP-2701 proof obligation into:

```text
sum_{t in {1,2,5,6}} signed_kernel_debt_t(B;u,v)
  <= boundary_margin(B,k).
```

That is a real reduction.  The bulk middle mass cannot spend the boundary
margin at all, so any proof that estimates all seven depths is paying for the
wrong object.  The remaining debt is concentrated in:

- shallow misses `t=1,2`, where relation can change hit/no-hit balance;
- high-tail misses `t=5,6`, where failing to hit enough sectors keeps the
  negative `p6` currency alive.

The compact S64 risk bank confirms the shape.  The k=9 and k=10 leaders both
spend margin through the same far pair `(15,16)`, but their `t=2` contribution
has opposite sign, while `t=5,6` repeat exactly.  That suggests the high-tail
part should be a finite far-pair phase atlas, and the shallow part should be a
signed live-state discrepancy bound rather than a scalar distance estimate.

Tournament Analysis: vertices are proof obligations/rows, not runners.  The
pairwise observable is smaller final slack, then more negative live-depth
deviation.  This quotient preserves the S64 risk order and is transitive in
the compact bank, so it is a valid proof-state compression for the current
true-wide branch.

The next sharp target is now narrow: prove the four live debts stay below the
already-positive two-far decorrelated boundary margins, with k=8 allowed to
use the cap dividend instead of the floor.
