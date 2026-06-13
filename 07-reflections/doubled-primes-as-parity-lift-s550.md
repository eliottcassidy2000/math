---
source: codex-2026-06-01-S550
status: synthesis + hypothesis
tags:
  - primes
  - goldbach
  - lemoine
  - parity
  - odd-cycles
  - doubled-primes
  - residue-wheel
---

# Doubled primes as the parity lift of prime cycles

The clean way to line up the pieces is:

```text
odd cycle        = parity obstruction / not 2-colorable
even cycle       = parity-compatible / two-sheeted
prime p          = irreducible odd cycle core
doubled prime 2p = even two-sheet lift of a prime core
```

That makes the Goldbach/Lemoine split feel much less arbitrary.

## Odd and even cycles

Odd cycles are the first real obstruction to a two-color or transitive picture.
They cannot alternate colors consistently.  In this repo's tournament language,
directed odd cycles are the objects in the conflict graph `Omega(T)`; they are
what make `H` larger than the transitive baseline.

Even cycles are different.  They can alternate between two sheets.  They may
still carry structure, but they are parity-compatible.  A cycle of length `2p`
is especially sharp: it is the bipartite double cover of the odd prime cycle
of length `p`.

So `2p` has exactly two pieces of information:

```text
the factor 2  = the parity sheet / bipartite lift,
the prime p   = the irreducible odd cycle core.
```

That is the whole reason doubled primes matter.  They are even without losing
their prime identity.

## Even Goldbach and odd Lemoine

For large primes, every prime is odd.  Therefore:

```text
odd + odd = even.
```

Goldbach for even `N` is the symmetric prime convolution:

```text
N = p + q.
```

The parity channel is already solved:

```text
1 + 1 = 0 mod 2.
```

For an odd target, two odd primes cannot work.  The Lemoine/Levy form is:

```text
N = p + 2q.
```

Now the parity channel is:

```text
1 + 0 = 1 mod 2.
```

The doubled prime is exactly what lets a binary prime formula survive on the
odd side.  Without `2q`, the natural prime-only route for odd numbers is
ternary Goldbach.  With `2q`, the ternary form is compressed to a diagonal:

```text
N = p + q + q.
```

So the doubled prime is the repeated-prime pair `(q,q)` carried as one even
term.  It is a diagonal Goldbach atom.

## The mod-6 wheel

On the `6k+-1` wheel, large primes occupy:

```text
P = {1,5} mod 6.
```

Doubling sends that prime set to:

```text
2P = {2,4} mod 6.
```

Then the two additive laws are mirror images:

```text
P + P   = {0,2,4} mod 6   all even residue channels,
P + 2P  = {1,3,5} mod 6   all odd residue channels.
```

Goldbach is the same-sheet convolution of the prime necklace.  Lemoine is the
cross-sheet convolution of the prime necklace with its doubled image.  The
small-prime wheel does not prove either conjecture, but it explains why the
local residue obstruction disappears in both cases: the relevant channels are
complete.

## Why doubled primes are not ordinary composites

A generic even composite has too much structure:

```text
12 = 2^2 * 3,
18 = 2 * 3^2,
30 = 2 * 3 * 5.
```

A doubled prime has exactly one dyadic factor and exactly one odd prime core:

```text
2p = 2^1 * p.
```

It is the first-even row of the add/multiply grid from HYP-1984.  Horizontal
motion `x+2` moves among odd cores; vertical motion `x*2` moves between parity
sheets.  A doubled prime is the point where those two motions meet cleanly:

```text
p --x*2--> 2p.
```

This is why `2p` should be read as a parity lift, not as a loss of primality.
It is the minimal even carrier of prime information.

## What this says back to cycles

The doubled prime lets an odd prime cycle participate in an even/bipartite
operation.  It is the arithmetic analogue of taking the double cover of an odd
cycle:

```text
C_p      = odd, irreducible, parity-obstructing;
C_{2p}   = even, bipartite, still remembers p.
```

That makes the Lemoine term `2q` conceptually precise.  It is not "prime plus
some even number."  It is:

```text
prime odd cycle + bipartite lift of a prime odd cycle.
```

The odd target is reached by combining one obstruction with one lifted
obstruction.  The lift pays the parity debt.

## Working slogan

```text
Goldbach:  even numbers are covered by P + P.
Lemoine:   odd numbers are covered by P + 2P.
2P:        the prime set on the parity-lifted sheet.
```

The doubled primes are important because they preserve the prime cycle while
making it parity-compatible.  They are the seam between odd-cycle arithmetic
and even-cycle transport.

## LRC boundary addendum

The n=4 LRC character formula adds a useful warning.  A doubled-prime leg is
not just an even speed; it is an even speed whose odd core is still prime.
In the mod-4 pair ledger, even reduced cofactors are `chi4`-silent, so the
bridge can be pairwise quiet.  But the same dyadic move can destroy the
all-odd `t=1/4` boundary witness.

That makes doubled primes boundary-active but pairwise quiet.  The right next
tournament vertices are bridge cells `(N,p,q)` for `N=p+2q`, compared against
Goldbach cells `N=p+q`, with endpoint-owner and residue labels retained rather
than projected away.
