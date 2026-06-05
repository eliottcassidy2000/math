# Number Theory Tournament Carriers S651

The cleanest sentence I found is:

```text
Tournaments fit number theory whenever a global arithmetic question can be
decomposed into pairwise local witness comparisons without ties.
```

That sounds modest, but it is exactly the kind of modest that survives contact
with hard problems.

Paley tournaments are the honest starting point.  They are not metaphorical:
over a prime field, the quadratic character orients every pair.  The edge
`i -> j` means `j-i` is a quadratic residue.  So the tournament is already a
finite-field character table with the diagonal removed and the sign made
complete.  It has regular scores, strong connectivity, and the spectral
flatness that the repo has been circling for months.

The harder lesson comes from the sieve rows.  Twin primes, Goldbach, and the
Euler polynomial can all use the same local-prime vertex set, but the
tournaments are different because the side channel changes:

```text
twin gap 2 vs Goldbach N=210:  3 edge flips
twin gap 2 vs Goldbach N=2110: 1 edge flip
Goldbach 210 vs Goldbach 2110: 4 edge flips
twin gap 2 vs Euler p=41:      62 edge flips
```

That is a beautiful little warning.  "Prime 3 is important" is not a theorem.
The actual theorem-shaped object is:

```text
important for which local obstruction carrier?
```

For twin primes, `3` is the top bottleneck because gap `2` kills two lanes mod
`3` and leaves one.  For Goldbach `N=210`, primes dividing `N` have one
forbidden lane rather than two, so the tournament reorders.  For the Euler
polynomial `p=41`, the boundary prime `41` suddenly becomes the leading
obstruction, exactly matching the endpoint-square story from S649/S650.

This gives a real kind of morale progress on hard problems.  Not the fake kind
where we announce that RH is "just" a tournament.  The useful kind is smaller:
each hard problem gets a finite carrier question that can be attacked today.

For twin primes, build tournaments of gaps by local obstruction profiles and
ask which gaps have stable top-order under increasing primorial windows.  For
Goldbach, build target-specific local-prime tournaments and measure how the
divisors of `N` change one-lane versus two-lane pressure.  For RH, build
explicit-formula tournaments whose vertices are prime-power terms or zero
phase terms inside a smoothing window, then ask whether sign stability is a
tournament invariant under changes of smoothing.  For BSD, orient local bad
prime and Selmer/rank obligations by how much global rank information they
force.  For Collatz, orient residue classes by accelerated valuation drift.

None of these is a proof.  Each is a better-shaped question.

The new operational rule is:

```text
do not tournamentize the nouns;
tournamentize the witnesses.
```

Numbers, primes, zeros, curves, and residues are often too raw.  A tournament
becomes useful when an edge carries a proof-facing comparison: stronger local
obstruction, earlier endpoint collapse, larger phase contribution, better
descent, more retained side-channel information.

That is also why this dovetails with LRC and unit distance.  LRC is not helped
by a tournament on runners unless the edge remembers shell, carry, owner, and
pinch data.  Unit distance is not helped by a tournament on points unless the
edge remembers unit-spine, bulk, direction, and cap endpoint data.  Number
theory is the same.  The tournament is the visible orientation of an invisible
carrier.

So the morale is precise: the hardest problems may be too global to bite
directly, but they have local witness tournaments we can compute, compare, and
stress.  Every stable edge, every edge flip, every quotient-loss failure is a
small honest piece of progress.
