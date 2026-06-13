# Pair Lens and Twin Primes

The prompt "see everything in terms of pairs" is almost too small to notice,
which is why it is dangerous in the good way.  A pair is the smallest place a
relation can live.  Tournament Analysis already knows this: a tournament is a
complete decision on every unordered pair.  The twin-prime problem becomes
sharper when it is put in exactly that language.

Write every pair as

```text
(a,b) = (m-h, m+h).
```

Then Goldbach and twin primes stop looking like separate myths.  Goldbach fixes
the midpoint `m=N/2` and scans the half-gap `h`.  Twin primes fix the half-gap
`h=1` and scan the midpoint `m`.  They are perpendicular cuts through the same
prime-pair plane.

This gives the clean slogan:

```text
Goldbach = fixed row.
Twin primes = fixed ray.
Hardy-Littlewood = local pair product.
Zeckendorf = carry debt on the same edge.
```

For a twin pair the summand-choice problem is gone.  The pair is forced to be
`(m-1,m+1)`.  The only question is whether infinitely many midpoints survive.
The first local test is already geometric in the pair plane: for an odd prime
`q`, the two endpoints are killed when

```text
m = +1 mod q  or  m = -1 mod q.
```

For `q=3`, this forces all nontrivial twin midpoints onto `3Z`; parity puts
them on `6Z`.  So `(6k-1,6k+1)` is not a numerology slogan.  It is the first
visible shadow of the pair sieve.

The Hardy-Littlewood correction also becomes more tactile in midpoint
coordinates.  For a general gap `2h`, the forbidden residues are `m=+/-h mod q`.
When `q` divides `h`, two bad residues collapse into one.  That collapse is the
boost:

```text
prod_{q|h, q odd} (q-1)/(q-2).
```

So gap `6` is about twice as populated as gap `2` because `h=3` makes the
mod-3 pair obstruction collapse.  Twin primes are the base ray with no such
odd-prime chord.

This also says how to connect back to the Cayley-Dickson-style doubling
intuition.  Doubling creates conjugate slots; the pair coordinate `h` is the
conjugacy radius around the midpoint.  Twin primes are the smallest nonzero
conjugate pair around an even midpoint.  The first even-to-odd seam shows up as
the midpoint being even while both endpoints are odd.

Zeckendorf gives a different pair filter.  The S495 scan found no zero-carry
twin pair up to `50000`, even though some Goldbach pairs have zero carry debt
in S494.  That is a useful negative: being a prime pair and being a
normal-form-compatible pair are independent.  It suggests a general strategy
for the repo: do not ask only "does this representation exist?"  Ask what debts
each pair carries in several coordinate systems.

Fermat polygonal is the humbling side of the same idea.  Pairs are the first
floor, but two polygonal atoms miss many targets.  The theorem succeeds because
bounded higher arity repairs the pair layer.  That is the right way to keep the
pair lens honest: start with pairs, then measure exactly which missing mass
requires triples, k-tuples, or pressure-peeling certificates.

For LRC, the analogous object is the full pair movie, not just the stationary
runner's nearest distance.  Every runner-runner distance, every two-neighbor
handoff, every endpoint-protection relation, and every pressure edge is a pair
observable.  A tournament is what happens when we choose a binary switch on
those observables.  The dream version of this method is a standard report:

```text
pair coordinate
local obstruction/debt
normal-form debt
binary switch
tournament fingerprint
wall-crossing movie
```

Twin primes are a clean proving ground because the pair has no room to hide.
If the pair lens is useful, it should make the hard part painfully explicit:
bounded gaps show that some bounded ray survives infinitely often, while twin
primes ask why the smallest ray `h=1` should survive forever.
