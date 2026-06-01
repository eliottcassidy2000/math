# The clearance cascade: transitivity's two facts, and where it deadlocks (S521)

*claudebox-2026-06-01-S521. Framing LRC as a cascade = a product of conditional
clearances, structured by the TWO facts of tournament transitivity: forward
chaining (Fact 1) and the hidden sideways-pruning (Fact 2). The cascade clears the
transitive backbone freely; it can only deadlock at the wrapping 3-cycles — which
are exactly the apex/seam = the observer — connecting to the project's odd-cycle
(OCF) core.*

## The two facts of transitivity

- **Fact 1 (forward chaining).** `X->Y and Y->Z  =>  X->Z`. In a cascade this
  propagates a clearance FORWARD along an order: once `X` precedes `Y` and `Y`
  precedes `Z`, `X` precedes `Z` for free.
- **Fact 2 (hidden, sideways pruning).** `X->Y  =>  not (Z->X and Y->Z)` for EVERY
  third node `Z`. A single arc excludes a JOINT configuration of every other node —
  it is the acyclicity dual of Fact 1 (it forbids the 3-cycle `X->Y->Z->X`). In a
  cascade this is constraint-propagation: one established clearance prunes the
  options of all remaining runners.

Both are equivalent to acyclicity; the point is the cascade USES them differently —
Fact 1 to chain, Fact 2 to prune.

## The cascade = product of conditional clearances

Loneliness is "clear every runner into the safe band `[1/n,1-1/n]`." Decompose it as
a chain (product of conditional clearances): order the runners along the band and
clear them one at a time,
`mu = P(clear r_1) * P(clear r_2 | r_1) * ... * P(clear r_{n-1} | rest) * (seam closes)`.
Within the band the half-turn tournament is the **transitive backbone** (a linear
order), so Fact 1 chains the order and Fact 2 prunes consistently — the backbone
runners clear "for free" in sequence. The cascade is a CSP that propagates.

## Where it deadlocks: the wrapping 3-cycles = the seam = the observer

The only place transitivity FAILS in the round tournament is the 3-cycles, and
(proved, verified 200/200):
> **a triple is a 3-cycle  <=>  it WRAPS the circle (no empty semicircle, every gap
> `< 1/2`)  <=>  it involves a LONG pair (separation `> 1/2`) = the apex/seam.**
The within-semicircle triples are transitive (the backbone); the wrapping triples
(through the long/apex pairs) are the 3-cycles. The OBSERVER sits at the seam (the
apex arc, where the band closes through `0`). So:

> the clearance cascade runs freely along the transitive backbone (the band order);
> it can only **deadlock at the wrapping 3-cycles at the seam — exactly the
> observer's clearance step.** LRC = the seam closes = the observer is cleared = no
> frustrating wrapping cycle traps it.

Fact 2 is precisely the statement that no single wrapping arc creates a 3-cycle; its
FAILURE (a wrapping triple IS a 3-cycle) is localized at the long/apex pairs, i.e.
at the observer's seam. The cascade's success is the absence of a frustrating cycle
through the observer.

## Connection to the OCF (the project's odd-cycle core)

A cascade deadlock = a **frustrated cycle** — a directed cycle that cannot be
linearized, i.e. an ODD cycle (3-cycle, 5-cycle, ...). These are exactly the
project's **Odd-Cycle Collection** (`H(T) = I(Omega(T), 2)`, Redei/OCF). The
observer is clearable (a source, lonely) iff it lies on no frustrating odd cycle.
So:

> **the clearance cascade completes  <=>  the observer is on no frustrating odd
> cycle through the seam  <=>  the observer can be made a source (lonely).**

This ties the cascade directly to the OCF: the conditional clearances are the
linear-extension steps; the obstruction is the odd-cycle structure `Omega` restricted
to cycles through the observer. The hidden fact (no 3-cycle from one arc) is the
base case of "no frustrating odd cycle"; loneliness is its global version.

## What the framing buys

- **Localization.** Of the `n-1` conditional clearances, the backbone ones are free
  (Fact 1/2 hold); the entire difficulty is the single **seam closure** (the
  observer's clearance), governed by the wrapping odd cycles. This is the
  fiber-only / apex localization from another angle.
- **A CSP / cascade structure.** LRC = a constraint-propagation cascade completes;
  transitivity gives the propagation (chain + prune); the failure mode is a
  frustrated wrapping odd cycle. This invites constraint-satisfaction / cycle-space
  tools (and connects to the even/odd-subgraph = cycle-space framework `E_n`).
- **The product form.** `mu = prod (conditional clearances)`; the backbone factors
  are bounded below (room along the band); the seam factor is the observer's gap
  (`>= 1/n`?), the one term that can vanish (the tight/frustrated case).

## Honest assessment and lead

This is a structural/conceptual synthesis, not a proof: it recasts loneliness as a
clearance cascade whose only deadlock is a frustrating wrapping odd cycle at the
observer's seam, rigorously localizing the difficulty (3-cycle <=> wrapping, proved)
and tying it to the OCF odd-cycle core. The lonely measure factors as a product of
conditional clearances; the backbone factors are controlled; the seam factor is the
LRC content.

Lead: make the cascade a genuine inductive bound — order runners along the band,
lower-bound each conditional clearance by the available gap (Fact 1/2 guaranteeing
the order and pruning), and reduce LRC to the seam-closure factor being positive,
i.e. to the observer lying on no frustrating odd cycle. Combined with the OCF
machinery (`H = I(Omega,2)`, the odd-cycle collection through the observer), this
is the cascade route to the apex/fiber-only LRC, now in the project's native
odd-cycle language.
