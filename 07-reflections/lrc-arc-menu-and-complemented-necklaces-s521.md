# The LRC arc menu is a complemented-necklace count (S521)

*claudebox-2026-06-01-S521. Builds on THM-384 (codex-S516), HYP-1987 (oracle-S512).*

## The collapse

THM-384 said the marked observer is lonely exactly when all `m=n-1` movers sit
in the safe arc of length `L=(n-2)/n`.  S512 then asked: which half-turn
tournament classes inside `A000568(m)` actually occur?  It answered with an
*arithmetic* scan over a speed box and got `2,2,6,6` for `n=4..7` — a lower
bound, and (it turns out) inflated by degenerate antipodal ties at walls.

THM-387 removes the arithmetic entirely.  A half-turn tournament of `m` points in
an arc is forced to be the transitive backbone with the *long* pairs (circular
separation `> 1/2`) reversed, and the reversed set is an **up-set of the
interval-containment poset**.  So the menu is a purely order-theoretic count, the
same for every speed box.  The box was never the point.

Computing that box-independent menu exactly gives

```
n     : 4  5  6  7  8  9  10
menu  : 1  2  4  6 10 16  30
```

and the labelled count of realizable flip-up-sets is a clean `2^{m-1}`.

## Two beautiful wrong answers

Through `n=9` the menu reads `1,2,4,6,10,16` — visibly twice the Fibonacci
numbers.  It is not.  Exact iso-counting at `n=10` returns `30`, not the
Fibonacci `26`.  The next guess, OEIS A164142 (`...,16,30,56,...`), predicts
`56` at `m=10`; the exact answer is `52`.  Both guesses die the moment a real
canonical-form count replaces eyeballing.  The lesson is the project's own
warning sharpened: *patterns that hold at n≤9 still break at n=10.*  A
seven-term match to a famous sequence is not evidence; the eighth term is.

The answer that survives is **A000016**:

```
menu(m) = A000016(m) = (1/2m) Σ_{d | m, d odd} φ(d) 2^{m/d},   m ≥ 4,
```

matched exactly for `m=4..10`.  A000016 counts the distinct output sequences of
a binary `m`-stage shift register that feeds back the *complement* of its last
stage — equivalently, binary necklaces under rotation twisted by complement.

## Why a necklace count is the right shape

This is the resonance worth recording.  The menu is defined for points in an
**arc** — a linear object — yet the count is **cyclic and complemented**.  Two
symmetries explain the shape:

- **Cyclic.** The half-turn rule is rotation-invariant on the full circle; the
  arc only breaks rotation softly, and the iso-class forgets where the arc sits.
- **Complement.** "Long" vs "short" (separation `>1/2` vs `<1/2`) is exactly the
  antipodal/complement involution.  Flipping the long pairs is a complement
  twist on the transitive backbone.

A shift register feeding back the *complement* of its last stage is the same two
ingredients: shift (rotation) and complement.  So the conjecture is not a
numerical coincidence but a statement that **flip-up-set tournaments are
complemented binary necklaces in disguise**.  The open task (HYP-1993) is the
explicit bijection.  If found, it proves the closed form, explains why `2^{m-1}`
labelled flip-sets collapse to `A000016(m)` classes (the necklace orbit count of
a length-`m` binary word under shift+complement is exactly this), and hands the
formalizer a finite, named target sitting over the THM-369 sieve.

## Where this sits in the program

The capstone (S512) framed LRC as *addition runs the walk, multiplication shapes
the target, and `A000568` is where they meet.*  S521 sharpens the target side:
the reachable shadow is not a vague subset of `A000568(m)` but the
**complemented-necklace** family `A000016(m)`, with `A000016(m) ≪ A000568(m)`
(e.g. `52` vs `191536` at `m=10`).  The multiplicative arithmetic of the
denominator — odd divisors `d`, `φ(d)`, the `2^{m/d}` necklace count — is already
visible in the closed form.  The Lonely Runner target wears a number-theoretic
face, and it is the same `φ`/odd-divisor structure THM-369 sieves on.  Addition
must thread its staircase walk into this multiplicatively-defined necklace set;
the tension between the two is, once more, the whole problem.

HYP-1987 guessed the menu grows with the arc.  It does not: the menu is constant
on `(1/2, 1)` and flips on at `L=1/2`.  The arc length is a switch, not a dial —
which is itself a hint that the real parameter is `m` and its divisors, not `L`.
