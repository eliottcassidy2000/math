---
source: codex-2026-06-03-S587
status: denominator-shield formalization plus Lemma A/B routing synthesis
tags: [LRC, Lemma-B, fold-sieve, denominator-gates, augmentation, HYP-2122]
---

# Fold-divisibility sieve: Lemma B after the visible fold

The session started with the user's instruction to work mostly on Lemma B but
keep Lemma A in earshot as noise and control.  The useful turn was to stop
treating `a+b=c` as the whole structural atom.

The real local fact is:

```text
D = a+b,  t = m/D  ==>  a*t + b*t = m.
```

So the pair `(a,b)` pinches at every denominator `D=a+b`, and any speed
divisible by `D` is at an integer at that same time.  The visible fold `c=a+b`
is the special shield `v=D`, but it is not the only shield.  This is exactly
why the `V*` floor row is less mysterious than it looked: it has no speed `12`,
but it has `24`, and `24` shields the `D=12` family.

That is the fold/sieve bridge in one line.

## What changed

Before S587, the story was:

```text
visible 3-term relation -> literal fold -> shield -> then sieve.
```

After S587, the sharper story is:

```text
pair denominator D=a+b -> any multiple of D is a shield -> then sieve.
```

The old visible fold is still real.  It is just the primitive case of a
divisibility family.  In THM-400 language, this is observer-coupled structure:
`v=q(a+b)` gives the unbalanced relation `q*a+q*b=v`, not a balanced
translation-invariant shadow.

The incoming HYP-2120/HYP-2121 perspective results say the same thing through
the tournament observer: LRC's useful recursion is rooted at a source
perspective, and full rooted classes still need incident threshold payload.
Denominator shields are rooted gates with exactly that payload.  They remember
which observer clocks have been killed; balanced energy is what remains after
forgetting the root.

## AP, V*, and shifted AP

The n=14 ledger has a clean calibration:

```text
AP:        all low D<n killed, D=n unshielded, M=1/14.
V*:        all low D<n killed, D=n unshielded, M=1/14.
unit AP:   low D killed, D=n shielded by 14, M-delta=+0.05357.
far AP:    no low denominator scaffold, M-delta=+0.28571.
```

So AP and `V*` agree at the level that matters: all lower pair-pinches are
stripped away, but the `D=n` clock remains.  Unit-shift AP is the warning
control.  It has lots of additive structure, but once `14` shields the
`D=14` clock, the floor witness is gone and the row loosens.

The doubled-apex row sits between the two worlds.  It kills the low
denominators and almost reaches the floor, but the remaining positive margin
points toward the endpoint/Phi residual rather than toward raw fold count.

## Why Lemma A stayed in the room

Lemma A is still the cleanliness test.  If there is no observer-coupled
denominator scaffold, then the row should behave like a discrepancy/gap row,
not a critical row.

The sample audit supports that reading.  Balanced energy was weak or negative
as a hardness predictor.  Unbalanced count was positive.  Circuit-free or
`u=0` controls kept substantial margins when they appeared.  The lesson is not
"fold-rich always tight"; it is "only observer-coupled denominator gates can
move a row into the terminal proof branch."

That matters because it prevents the fold proof from overfitting AP.  The
right split is:

```text
denominator gates present -> Lemma B / Phi / endpoint branch
denominator gates absent  -> Lemma A discrepancy branch
```

## Prime machinery position

The deletion-shadow audit explains where the prime-divisibility machinery
plugs in.  Removing a visible fold produces a shadow witness that is often
blocked when the fold is restored.  For AP and `V*`, every visible-fold shadow
tested is blocked.  But the residual denominators after blocking are not all
trivial: they include values such as `17,19,23,25,31,35`.

So the fold is not the whole proof.  It creates the residue ledger.  The
endpoint, `Phi`, CRT, and prime-fibre machinery should consume that ledger
after denominator shields have done the cheap local work.

## Tournament Analysis note

S587 deliberately did not use runners as tournament vertices.  It used proof
lenses:

```text
D_denominator_divisibility_shield
visible_fold_a_plus_b_eq_c
Phi_endpoint_prime_residual
D_eq_n_unshielded_delta_clock
augmentation_nonzero_count
circuit_free_A_margin
balanced_energy_background
```

The observable ranked delivery power, control status, divisibility power,
observer coupling, and maturity.  The tournament was transitive, with no
directed 3-cycles and one Hamiltonian path.

That quotient preserves the proof predicate we care about: which pair-pinches
are killed, and whether the `D=n` clock survives.  It destroys endpoint-owner
geometry and prime-fibre language on purpose; those are downstream gates.

## Handoff

The new proof target is a denominator-gate compression theorem.

For total `n`, show that if all low pair denominators `D<n` that occur are
shielded and `D=n` is unshielded, then the `1/n` clock survives as in AP and
`V*`.  If `D=n` is shielded, route to the multiple branch and prove positive
`Phi`.  If the denominator scaffold is absent, route to Lemma A and prove the
gap/discrepancy margin from lack of observer-coupled gates.

This does not yet prove the n=14 case, but it gives a sharper casing language:

```text
low D killed / D=n survives / D=n killed / no D-scaffold.
```

That casing is much closer to a finite certificate calculus than "many folds"
or "high energy."
