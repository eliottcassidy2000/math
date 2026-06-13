# Irreducible Prime Carriers and Certificate Tournaments

Source: codex-2026-06-12.  Companions: arXiv:2411.18366,
arXiv:2410.15880, HYP-2447, HYP-2448, OPEN-Q-070,
`04-computation/irreducible_prime_carrier_tournament_codex.py`.

The user's sentence is unusually clean:

```text
irreducible polynomials can produce primes;
primes can produce irreducible polynomials.
```

That is close to a heartbeat.  A prime is an atom of `Z`.  An irreducible
integer polynomial is an atom of `Z[x]`.  Evaluation `f -> f(n)` is the bridge,
but it is a violent quotient: it turns a structured object into one integer.
Bunyakovsky says that if the polynomial atom is primitive enough, then this
quotient should hit integer atoms infinitely often.

Singh's paper is valuable because it runs the arrow backward.  A large value
with few prime factors limits how many irreducible factors the polynomial can
have.  Cohn's digit criterion is the same move in base notation: if the digit
value is prime, the digit polynomial cannot split.  This makes "prime values"
not merely output data but certificates.

Iravanian's real-factor recombination paper supplies the missing algorithmic
middle.  Factoring over `R` produces linear and quadratic blocks.  A rational
integer factor must choose a subset whose traces and later coefficients land
back in `Z`.  The first trace test is a subset-sum filter, not a proof by
itself.  The atlas's little warning is perfect: `x^4-10x^2+1` is irreducible
but already has two false integer-trace subset candidates.  So recombination is
a tournament of surviving candidate subsets, with false edges that must be
killed by richer coefficient/support data.

This is where the repo connection becomes more than metaphor.  We have kept
learning the same lesson in other dialects:

```text
scalar quotient passes; retained support channel decides.
```

In LRC14, "q is blocked" is too scalar; the proof-bearing object is the
resource vector over shell-27 classes, 13-clock data, divisor fibers, and
owner-deletion/Bprime exceptions.

In the `[72,36,16]` problem, the Type II weight enumerator is healthy.  That is
not the code.  The code is a support/matroid/design realization problem.  A
scalar enumerator is like one value of a polynomial: necessary evidence, but
not the factorization itself.

So the tournament should not live on the nouns.  Do not tournamentize primes,
polynomials, runners, or codewords first.  Tournamentize the witnesses:

```text
forward prime-hit witnesses;
reverse Singh/Murty/Cohn certificates;
local fixed-divisor obstructions;
real-factor recombination survivor sets;
Newton polygon support data;
code support/matroid construction moves;
LRC resource vectors.
```

The stored proof-carrier tournament is nontransitive.  That is the right
answer at this stage.  A transitive ranking would mean one scalar proxy had
quietly taken over.  The three directed cycles say the carriers trade off:
Bunyakovsky is strong forward and weak reverse; Singh/Cohn are strong reverse
and weak asymptotic; recombination is algorithmically strong and proof-noisy;
support gates are the best bridge to LRC14 and length 72.

The possible infinite tournament is then clear enough to build.  A vertex is a
certificate state `(f, retained data, cutoff X)`, not just `f`.  Orient one
state toward another when it has smaller unresolved factorization ambiguity
after normalizing fixed divisor, degree, and certificate depth.  As `X` grows,
edge flips become data.  Bunyakovsky predicts that admissible irreducible
vertices keep receiving forward prime-hit witnesses forever; Singh, Cohn, and
Iravanian provide the reverse edges that keep the state from dissolving into
wishful numerology.

This does not prove a famous conjecture.  It gives the project a better-shaped
object to compute.
