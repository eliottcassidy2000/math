# Investigating the speed bound B(n): the fast-runner engine and why it is the whole problem (S521)

*claudebox-2026-06-01-S521. The speed bound is the external ingredient the
multiplicative-walk methodology needs (lrc-methodology-needs-a-speed-bound-s521.md).
This investigates where it comes from, how it couples to the methodology, and why
it is the genuine hard core. Honest: it does not produce a usable B(14).*

## Standard reductions (free)

LRC(n) for reals reduces to rationals (closed condition) reduces to integers
(clear denominators) reduces to primitive (gcd 1, scale `t`).  So WLOG: distinct
positive integers with gcd 1.

## Coupling to the methodology (proved + verified)

- **Loneliness at `q` depends only on speeds mod `q`.**  Hence "lonely at some
  `q ∈ Q`" depends only on speeds mod `L = lcm(Q)`.  So that statement is a finite
  check over speed sets with entries in `[0, L)` — about `L^{n-1}` of them.
- **Denominator grows like `log(speed)`** (verified): the extremal family is
  `(1, lcm(2..Q))`, speed `~ e^Q`, first lonely at `q ~ Q+1`.  So with a speed
  bound `B`, the methodology's modulus is `Q ~ ln B` (small), and `L = lcm(2..Q)
  ~ B`.  The finite check is over `~ B^{n-1}` systems.

Conclusion: the methodology contributes a fast per-check (small denominators), but
the **number** of checks is `~ B(n)^{n-1}`.  **`B(n)` is the entire bottleneck.**

## Where B(n) comes from: the fast-runner bound

The engine of every LRC speed bound (and the `n ≤ 7` proofs):

> If the other `n-2` movers have a common-safe time window `I` (all simultaneously
> in `[1/n, 1-1/n]`) and the fastest mover `v` satisfies `v·|I| ≥ 1`, then `v`
> sweeps the whole circle inside `I` and lands in the safe band there — so the
> system is lonely.  Hence a counterexample needs `v_max < 1/|I|`.

Demonstrated for `n=4`: `v3 < 1/|I(v1,v2)|`, and empirically `|I| ~ 1/(2 v_2nd)`,
so `v_max ~ 2 v_2nd` **when the window exists**:

| (v1,v2) | \|I\| | v3 must be < |
|---|---|---|
| (1,3) | 1/6 | 6 |
| (7,9) | 1/18 | 18 |
| (11,13) | 1/26 | 26 |

The min window over coprime slow pairs with speeds `≤ S` is `~ 1/(2S)`, giving
`v3 < ~2S` — a clean ratio bound *conditional on a nonempty window*.

## Why it is hard: the window can be empty

For `n ≥ 4` the `n-2` slower movers can have an **empty** common-safe window —
their bad sets (each of measure `2/n`) can tile the circle (measure bound
`1 - (n-2)·2/n ≤ 0`).  Then `|I| = 0` and the fast-runner argument gives **no
bound** on `v_max`.  Those configurations require separate arguments, and the
recursion does not close cleanly because deleting the fastest runner does not
yield a standard sub-instance (the threshold `1/n` does not change).  This
casework is exactly why:

- the proven cases stop at `n ≤ 7` (the casework was carried through there);
- `B(n)` is **super-exponential** in general (each level of the inductive bound
  multiplies); the best general finiteness (Tao, 2018) gives velocities of size
  `~ n^{O(n)}`, enough to make LRC(n) a finite computation but astronomically
  beyond reach for `n = 14` (`B(14)^{13}` is hopeless).

## The deep tool: covering systems

A genuine counterexample is an **interval-cover of the circle** by the bad-set
clusters `B_i` (runner `i` contributes `v_i` evenly-spaced intervals of width
`2/(n v_i)`).  This is the continuous cousin of an Erdős **covering system** of
congruences.  Bounds on covering systems — Hough's minimum-modulus theorem (2015)
and the Balister–Bollobás–Morris–Sahasrabudhe density results — are the natural
deep machinery for bounding the moduli `v_i` of such a cover, and this is the
project's existing covering-systems thread.  A covering-systems-style bound on the
interval-cover would be a principled route to `B(n)`; adapting the discrete
theorems to the interval setting is the open technical step.

## Verdict

- The methodology is the correct **reduction** and gives a `log`-size modulus, but
  the finite check it reduces to is `~ B(n)^{n-1}` systems.
- `B(n)` is produced by the fast-runner / inductive-interval engine, with the
  empty-window casework as the obstruction; it is explicit but super-exponential,
  done through `n ≤ 7`, hopeless at `n = 14`.
- So "investigate the speed bound" lands on the honest truth: **`B(n)` essentially
  *is* the Lonely Runner problem for fixed `n`.**  The methodology does not shrink
  it; it reframes the per-instance test.  The only principled hope for a better
  `B(n)` is the covering-systems connection.

## Seed

Make the interval-cover ↔ covering-system correspondence precise and import a
Hough-type minimum-modulus bound to bound `min_i v_i`, then bootstrap the rest via
the fast-runner ratio bound `v_{i+1} ≲ 2 v_i` on the (nonempty-window) branches.
The empty-window branch is where a new idea is required — and it is the same
"fully covered / rational coincidence" core that every earlier S521 reframing
landed on.
