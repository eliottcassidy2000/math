# LRC14 AP-tail puncture atlas

**codex-2026-06-26-S204.**  A follow-up to HYP-3029 and HYP-3031, focused on
the two coarse-stalk residual teeth in the hard AP/GW automatic word.

## What moved

The two coarse residuals from HYP-3029 looked like a possible new mixed
coordinate:

```text
13->159 / 13->117
13->118 / 13->104
```

They are not mysterious.  They are AP-tail rows

```text
S_m = {1,...,12,m}
```

where the coarse stalk remembered `m mod 14` and the largest component's owner
strip, but forgot `m mod 13`.

That is the missing clock.

## The small theorem

For `S_m`, either `t=1/13` works or a reciprocal fixed point works:

```text
13 does not divide m:  t=1/13.
m=13s, s>=2:          t=s/(13s+1).
```

The proof is one line of geometry.  In the second case, runners `1` and `m`
bind at height `t`; the AP core runners have distance `j t` or `1-j t`, and
the only possible small opposite-side core runner is `j=12`, whose distance is
`(s+1)/(13s+1)>t`.  Since `t>1/14` exactly when `s>1`, the only boundary case
is `m=13`.

## Why this matters

This is a nice HYP-3031 instance because it says exactly what the quotient was
allowed to forget.

The mod-14 owner strip is real:

```text
('1+') -> ('1+','6-') -> ('6-')
('1+') -> ('1+','5-') -> ('5-')
```

But it is not complete.  It must be paired with the exact-period puncture:

```text
m mod 13 != 0  -> direct q=13 witness
m mod 13 == 0  -> reciprocal fixed-point witness
```

So the repair class is not F7 debt.  It is `nested_refinement` plus
`owner_strip`: a family formula reconstructs the missing mixed coordinate.

## New proof habit

Before feeding an open-route residual to Fejer or THM-572, ask whether it is
really a one-tail or two-tail AP-core deformation.  If it is, try:

```text
prime puncture clock  or  reciprocal fixed point.
```

This is the same shape as the older q-witness/fixed-point pincer, but now it
is localized to a coarse-stalk quotient.  The useful slogan is:

```text
owner strip + missing prime clock = repairable mixed coordinate.
```

## Scope

This proves a family lemma, not LRC14.  But it removes the two named S193
residual teeth and gives future agents a cheap first move for AP-core tail
families inside the HYP-2963 packet bank.
