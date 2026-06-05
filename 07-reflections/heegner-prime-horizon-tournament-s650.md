# Heegner Prime Horizon Tournament S650

This is a refinement of the incoming S649/HYP-2225 Euler-Heegner boundary
packet.  The extra point here is narrower: THM-410 gives an exact tournament
model for the `p-2` interior witnesses, and the pair-witness split for the
reversed long edge is `(sigma, lambda, delta)=(0,p-2,0)`.

The user's phrasing "first `n-2` inputs" is the important signal.  The prime
generating set

```text
{2,3,5,11,17,41}
```

is the Euler/Rabinowitsch set for

```text
Q_p(x) = x^2 + x + p.
```

The polynomial is prime at `x=0..p-2`, and the first forced composite is

```text
Q_p(p-1) = p^2.
```

If we treat `x=0` as a boundary input, the positive interior run is exactly
`x=1..p-2`, hence `p-2` witnesses.  That is the tournament-shaped part.

The Heegner connection is also typed.  The six `p` values are not themselves
the full Heegner list.  They map by

```text
d = 4p - 1
```

to

```text
{7,11,19,43,67,163},
```

the Heegner tail.  The missing base values `{1,2,3}` are separate base fields,
not failed data.  This matters because the repo has repeatedly been burned by
letting a scalar migrate between roles.

THM-410 supplies the clean tournament side.  In a transitive tournament, reverse
one long edge `0 -> p-1`.  The only directed triangles are

```text
0 -> v -> p-1 -> 0
```

for the `p-2` interior vertices.  So the same boundary/interior count appears:

```text
#C3 = p-2.
```

This is not a statement about global Hamiltonian-path counts.  It is a local
pair-witness horizon.  For the reversed boundary edge, all interior witnesses
fall into the cyclic slot:

```text
sigma=0, lambda=p-2, delta=0.
```

That pins the analogy to the already-recorded identity
`sigma+lambda+delta=n-2`.

The wild but plausible proof use is this: class number one is a no-hidden-ideal
condition, and the Euler prime horizon is a no-early-interior-collapse
condition.  The tournament model says a marked boundary reversal can force all
interior witnesses to be accounted for with no extra leakage.  If this analogy
is real, then the useful object is not a tournament on primes.  It is a
tournament on proof obligations or witness slots:

```text
boundary prime
interior input
residue obstruction
discriminant/class group
forced square boundary
```

The speculative next theorem would not say "Heegner numbers are tournaments."
It would say something like:

> A class-number-one quadratic has a single-boundary-collapse witness ledger:
> every interior slot is discharged before the boundary square, and failure of
> class number one is witnessed by an interior composite slot.

The non-Heegner controls in S650 fit that mood: `p=7,13,19,23,...` often fail
immediately at `x=1`, long before the boundary.  That is exactly what a hidden
side channel should do when unique factorization is gone: an interior slot
collapses early.

For LRC and unit distance, the lesson is familiar but useful.  The hard object
is again a retained side channel, not the visible scalar.  LRC has shell
modulus, gcd strata, carry, and owner labels; unit distance has spine, bulk,
direction support, and cap endpoints; the prime-polynomial row has boundary,
interior slots, discriminant, and class group.  All three are "safe until a
boundary" problems only after the right side channel is kept.
