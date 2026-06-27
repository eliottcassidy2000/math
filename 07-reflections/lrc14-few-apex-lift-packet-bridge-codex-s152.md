# LRC14 Few-Apex Lift Packets: Boundary-Moment Bridge, Not Another Scalar

The useful correction in this pass is small but structural: the covering
residual is not "many multiples of 14" anymore.  THM-571 already closed
`|14Z cap S|>=7` by the apex-majority gamma descent, modulo the below-frontier
LRC input.  So the live branch should be read as few-apex:

```text
S = R union 14Q,  1 <= |Q| <= 6.
```

HYP-2961 still had the older `|Q|>=7` language, which would send us back into a
closed branch.  HYP-2968 updates the proof route: the remaining work is the
scale-separated / finite-core side where only a few apex multiples are present
and the residual speeds have enough mass to threaten every lifted phase.

The exact packet is `u=14t`.  A multiple `14m` becomes `m u`, independent of
the lift.  A residual speed `r` forbids rational intervals in lift `k`:

```text
u in ((14n-1)/r-k, (14n+1)/r-k) cap [0,1].
```

That is the boundary-moment bridge in a usable form.  It is not merely a Haar
measure number.  It has owners: the Q comb, fourteen lift labels, residual
interval deletions, and a possible K33/state-lift route if the deletion pattern
becomes nonunit and zero-open.

S152 audited `8190` structured AP-replacement rows with `qdiv>14` and
`1<=|14Z cap S|<=6`.  Every row had positive exact lift-safe mass.  The
smallest was `563/105105`, from the familiar `12->84` hard covering row.
The exact-M fallback on the four tightest rows gave `7/89`, `2/23`, `2/23`,
and `14/173`, all safely above `1/14`.

The pattern is encouraging because increasing the number of 14-multiples does
not move the bank toward zero-open behavior.  The minimum number of open lifts
increases from `2` at `k14=1` to `8` at `k14=5,6` in this bank.  That suggests
the obstruction is not "Q has too many comb teeth"; the real danger would have
to be a coordinated residual-lift deletion pattern, exactly the kind of
labelled non-scalar packet the HYP-2956 theorem is supposed to retain.

Assumption challenge:

```text
Considered vertices:
  runners, Q multipliers, residual speeds, fourteen lifts, q-divisor
  obligations, interval components, boundary events, Fourier modes,
  and proof obligations.

Chosen vertices:
  lift packets and proof carriers.

Preserved:
  qdiv>14, Q/R split, |14Z cap S|, exact lift-safe mass,
  and positive-open witness existence.

Destroyed:
  exact Farey maximizer and fine residual owner identities.
```

The next theorem should be a fixed-margin few-apex lift theorem:

```text
primitive qdiv>14, 1<=|14Z cap S|<=6
  -> positive regular-open lift packet
     or K33/HYP-2908/THM-572 state lift.
```

That theorem would plug directly into the labelled packet theorem: F5 is
positive covering boundary-moment strictness, F6 is the zero-open
non-migration kernel, and S152 found no F6 packet in this structured bank.
