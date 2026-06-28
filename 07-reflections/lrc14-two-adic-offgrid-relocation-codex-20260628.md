# LRC14 Two-Adic Off-Grid Relocation Reflection

Date: 2026-06-28

The prompt's resonant-transparency picture is valuable, but only after the
HYP-3418 correction.  The sentence "resonance puts danger on the grid while
the optimizer is off-grid" is not enough to prove the covering floor, because
the coprime-to-14 speeds alone usually choose `t=1/2`, and every even speed is
dead there.

The sharper reframe is to stop separating resonant from nonresonant speeds
globally and instead split the clock by the doubling map:

```text
S = O union 2E,     u = 2t.
```

Then even speeds become the smaller halved packet `E`, while odd speeds become
one of two branch filters.  The off-grid lift is `t=(u+1)/2`; it keeps odd
speeds near half-integers and lets the even half choose the real witness.  This
is exactly the "relocation" missing from the naive nonresonant witness.

The scout's result is small but clean: `24/24` audited covering rows have exact
two-adic relocation certificates, while the naive nonresonant witness fails
`24/24`.  The proof target is now concrete:

```text
E_safe(1/14) cap (odd_branch_0_good union odd_branch_1_good) is nonempty.
```

This target also integrates the concurrent owner-cut work.  HYP-3419 saw
`2:g2` as the even-cover label in the `10->20` charal tree, and HYP-3417 saw
the frontier current `{2:g2,11:g1,13:g1}`.  HYP-3422 says why that even label
matters analytically: it is the visible finite shadow of the two-adic descent
coordinate.

The next proof attempt should avoid two traps.  First, do not return to an
apex-7/Galois explanation for the covering floor; HYP-3418 puts that structure
off the critical path.  Second, do not scalarize the branch-overlap theorem too
early.  The branch-bad intervals are explicit rational intervals.  The proof
should try finite-ruler, Helly, or interval-piercing estimates on the halved
even safe set before invoking broad decorrelation language.
