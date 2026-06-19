# LRC14 Euler-copy squarefree profile

**Source:** codex-2026-06-19-S21, HYP-2629 / T877.

This is the Hill-row/crossing-gate complement to HYP-2628 / T876, which
establishes the exact-period `phi` packet interpretation of the same copy rule.

The user's copy rule is exactly the Euler totient:

```text
sum_{d|n} c(d)=n  =>  c(n)=phi(n).
```

The useful part is not the identity by itself.  The useful part is grouping the
totient copies by squarefree prime masks:

```text
copy_mass_N(M)=sum_{d|N, mask(d)=M} phi(d).
```

For squarefree `S`, adding a new prime `p` is a genuine recurrence:

```text
old masks stay,
new masks = (p-1) shifted copies of old masks.
```

That turns the LRC14 address chain

```text
{2,3} -> {2,3,5} -> {2,3,5,7}
```

into an Euler-copy recurrence rather than a loose modular slogan.

The sharp new readout is at `K_14`:

```text
P_14 = 5*6*6*7 = 1260
cr(K_14) = 315.
```

The raw product has full `{2,3,5,7}` copy mass `576`; the divided crossing value
has full copy mass `0`.  Dividing by four removes the whole dyadic gate.  This
is a stronger version of the S20 warning: the raw Hill product is the proof
ledger, and the crossing value is already too quotientized for the LRC14
mod-210/coimage transfer.

The detailed atom is also informative:

```text
full copy mass of 210  = phi(210)=48
full copy mass of 1260 = (2^2-1)(3^2-1)(5-1)(7-1)=576.
```

The repeated `6,6` blocks thicken the full mask by `3*4=12`, not merely by the
squareful quotient `6`.  Ordinary `phi(1260)=288` is a different projection:
it sees exponent thickening by `2*3`, while the full squarefree copy atom sees
all nonzero p-power divisor layers.

The Markov-Hurwitz equation remains an analogy, not a carrier.  In copy
language, `wxyz` counts four-block copy assignments, but for `(5,6,6,7)` the
diagonal energy is only `146`, far below `1260`.  The reframe clarifies the RHS
product, not the equality.

Assumption challenge: runners, raw divisors, divisor copies, squarefree masks,
Hill factors, crossing values, Markov-Hurwitz coordinates, mod-7 coimage
classes, and proof obligations were all considered.  The best quotient is the
Euler-copy squarefree profile: it preserves the mod-210 proof address and
exposes why `315` loses information.

Next useful test: apply this profile to the HYP-2626 repeated-residue k=10 tail
and see whether copy mass separates the quadratic-character cases better than
the raw `{7}` prime mask alone.
