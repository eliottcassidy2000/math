# Doubled primes as first-even parity gates

This pass started from a simple parity puzzle:

```text
even number = prime + prime          (Goldbach)
odd number  = prime + doubled prime  (Lemoine/Levy)
```

The doubled prime looks like a hack if we only read it additively.  But inside
the LRC add/multiply stack it has a cleaner meaning.

All odd primes live in the odd-core layer.  The number `2` is not another odd
atom; it is the dyadic transport operator.  So `2q` is not just "a composite
summand."  It is:

```text
prime core q, transported exactly one dyadic step.
```

That makes `p+2q=N` the first-even bridge version of `p+q=N`.  Goldbach uses
two odd prime atoms to hit an even target.  Lemoine uses one odd prime atom and
one first-even prime-core atom to hit an odd target.  The doubled prime is
where the parity bill gets paid.

## The cycle analogy

The tournament side already teaches the same lesson.  In the OCF, the
surviving pieces are odd directed cycles:

```text
H(T) = I(Omega(T), 2) = 1 + 2 alpha_1 + 4 alpha_2 + ...
```

Even-cycle structure tends to be cancellation, pairing, or quotient scaffolding.
Odd cycles are where obstruction survives.  For primes, the odd primes are the
surviving atoms; `2` is the parity sheet-change.  A doubled prime is therefore
an odd atom carried across the parity sheet.

That is the important distinction.  Arbitrary even numbers can have many
branches and many prime factors.  Doubled primes are the smallest even
numbers that still remember a single prime core.

## Why LRC should care

The current LRC stack already separates:

```text
x + 2  : odd-core horizontal motion
x * 2  : dyadic vertical motion
```

The doubled prime sits exactly at their interface.  It is additive because it
appears in `p+2q`.  It is multiplicative because the even leg has address
`2^1 * q`.  It is prime-sensitive because the odd core `q` has not been
blurred into a generic composite.

This matters most in the boundary/tight cases.  HYP-2043 says that for
`n=4`, every all-odd triple has the closed witness `t=1/4`, so any tight
triple must contain an even speed.  THM-393 says the pairwise correction uses
the mod-4 odd character, so even reduced cofactors are pairwise silent.  A
doubled-prime leg therefore has a very specific role:

```text
it can break the all-odd witness,
but it does not necessarily create pairwise chi4 debt.
```

So the doubled prime is a boundary-active but pairwise-quiet gate.  That is
probably the answer to "why are doubled primes important?" in this LRC frame:
they move difficulty out of the easy pairwise layer and into the full-support
or endpoint-boundary layer where LRC actually lives.

## Formal object to try next

Do not make runners the default vertices.  The more natural vertex set here is
the representation bridge:

```text
b = (N, p, q) with N = p + 2q.
```

For comparison, Goldbach bridges are:

```text
g = (N, p, q) with N = p + q.
```

Each bridge should carry:

```text
target parity,
v2 of the even leg,
odd core of the even leg,
representation abundance in the target fiber,
endpoint-debt proxy phi(N)/N,
pairwise chi4-neutrality,
and any endpoint-owner labels available from the LRC row.
```

A Tournament Analysis switch can orient one bridge over another by lower
boundary debt or higher repair potential, with ties broken by increasing
`(N,q,p)`.  The fingerprints to report are the usual repo ones: score
histogram, directed cycles, SCCs, Hamiltonian-path count, and edge flips
against the existing S513 add/multiply denominator tournament.

## Working synthesis

Goldbach is additive smoothing inside one parity class.  Lemoine is additive
smoothing plus a dyadic gate.  The doubled prime is the minimal object that
does both at once.

For LRC, that makes doubled primes candidate bridge events: first-even,
prime-core-preserving, pairwise-quiet, and boundary-active.  That is exactly
the kind of object the current endpoint/pressure/marked-fiber program wants
to keep instead of projecting away.
