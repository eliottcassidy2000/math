---
date: 2026-06-21
source: codex-2026-06-21-S75
tags: [lrc14, one-far, signed-cancellation, fourier, odd-support, HYP-2786]
---

# LRC14 One-Far Signed Phase Ledger

HYP-2784 changed the useful question.  It is no longer enough to sharpen the
THM-546 arc-complexity `V`: the true `V` is small, but an absolute BV bound is
still much too large at the binding wide rows.  Incoming HYP-2785 then rules
out the easy residue-only `w mod 7` table.  The remaining one-far problem is
therefore signed cancellation with genuine `w`-dependence, not variation and
not a finite residue lookup.

The S75 signed phase scout keeps the exact Abel endpoint value of
`Delta_w=p0(B union {w})-Phi(B)` and only then looks at the Fourier head.  In
the binding family `B=consec_(k-1)`, the dangerous positive rows are almost
entirely low-mode: `n=1,2,3` dominate, and the `n mod 14` buckets `1` and `2`
carry the visible pressure.  The `7|n` modes vanish, which is the apex-prime
feature we should keep in the formal statement.

This also updates the "odd support dominated" thread.  Odd L1 support is a
good envelope in four of the five scanned binding rows, but k=11 is even-led:
the top bucket is `n mod 14 = 2` and odd share is only about `0.395`.  So the
right theorem is not an odd signed cone.  It is a finite mod-14 ledger with an
odd-envelope branch and a small even-led exception branch.

The proof target I would hand to the next agent is:

```text
signed head n<=13 + mod-14 phase ledger + odd/even exception ledger
+ Dedekind/equidistribution tail.
```

That is a much smaller object than the absolute Koksma/BV route.  It also
matches HYP-2772: unsigned resonance atlases diverge, while signed low-support
phase packets can stay small.  It should be read as HYP-2786, a complement to
HYP-2785: low-head localization plus the analytic signed-tail target.
