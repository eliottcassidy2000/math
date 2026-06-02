# LRC n=17 prime-gate attempt: the skip-8 packet

codex-2026-06-02 S559.

The request was to make an attempt at `n=17`.  In repo convention this means
16 moving speeds and threshold `1/17`.

## What the prime row gives immediately

For prime `17`, the top denominator-sieve branch is clean.  If no speed is
divisible by `17`, then for every speed `v` and every unit `a mod 17`,
`va mod 17` is nonzero.  Hence at `t=a/17`,

```text
||v t|| >= 1/17.
```

So all 16 unit walls `a/17` are closed witnesses.  Any open-cover
counterexample must contain a `17`-multiple.

This is not new technology; it is THM-369 in prime-row dress.  The useful part
is what comes after the forced gate.

## What the first repair does

The one-gate swaps

```text
{1,...,16} - {d} + {17q}, q <= 32
```

all stayed positive-gap.  The closest non-sieve-complete row was
`drop16_add17x1` with `gap/th=0.025641`; the closest sieve-complete one-gate
rows bottomed at `gap/th=0.027574`.  Representative one-gate rows also had
empty endpoint cores after peeling.

That says a single prime gate repairs the unit wall but exports the obstruction
to finer denominators rather than building an open cover.

## The new pressure family

The sharper family was

```text
{r} union {17*q : 1 <= q <= 16, q != s}.
```

Here a primitive breaker `r` keeps gcd one, while the `17q` speeds carry the
gate/sieve burden.  The closest exact rows all skipped `s=8`:

```text
gap/th = 0.003676...
forbidden length = 81463/82654
boundary count = 450
```

The witnesses land at denominators `9248 = 17*544`; examples include
`609/9248`, `479/9248`, and `65/9248`.  This is the first real n=17 object from
the session.  It looks like a half-gate packet: skipping `8` leaves the
middle divisor of the `1..16` gate payload, and the remaining `17q` ladder
nearly covers the circle but still leaks a small positive corridor.

## Assumption challenge

I did not use runners as the Tournament Analysis vertices.  The candidate
vertices considered were:

- runners;
- unit walls `a/17`;
- skipped gate labels `s`;
- primitive breaker choices `r`;
- gate multiples `17q`;
- endpoint leaves;
- whole repair rows;
- proof obligations such as "has positive gap" or "has private endpoint leaf".

I used repair rows because the n=17 predicate under attack is row-level:
does this repair of the killed unit wall produce an open cover?  This preserves
gap size, sieve status, and endpoint debt for candidate configurations.  It
destroys runner-level incidence and should not be treated as a final proof
certificate.

The endpoint-audited row tournament was transitive with one Hamiltonian path.
That is useful mostly as a sanity check: the chosen scalar repair gauge is a
ranker, not a cyclic obstruction detector.

## Honest status

No n=17 proof.  No open-cover candidate either.

The best next theorem is narrower:

```text
For n=17, in the primitive-breaker plus 15-gate family,
the skipped gate s=8 is extremal and still leaves an explicit positive gap.
```

The next computation should derive a formula for the `skip 8` gap instead of
merely printing it.  If the formula is clean, compare it with n=16 dyadic
half-gate packets and n=18 square-core gate ladders.
