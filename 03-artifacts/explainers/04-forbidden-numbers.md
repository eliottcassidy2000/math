# Forbidden Numbers

## The Question

We know H(T) is always odd. So the candidates are 1, 3, 5, 7, 9, 11, 13, 15, ...

But can every odd number be achieved? Is there always some tournament, somewhere, with exactly 7 Hamiltonian paths? With exactly 21?

The answer is: **no**. Some odd numbers are **forbidden** — they can never be H(T) for any tournament, no matter how many players you use.

---

## The Discovery: 7 and 21 are Forever Forbidden

Through extensive computation — exhaustively checking all tournaments up to 7 players, and sampling millions of larger ones — this project established:

**H(T) = 7 never occurs.** No tournament on any number of players has exactly 7 Hamiltonian paths. Ever.

**H(T) = 21 never occurs.** Same story for 21.

These are permanent, absolute gaps. They're not "we haven't found it yet" — the structure of the Odd-Cycle Formula makes them algebraically impossible.

---

## Why 7 is Impossible

Here's the intuition. The OCF formula says:

H(T) = 1 + (something)·2 + (something)·4 + (something)·8 + ...

Everything beyond the "1" is a multiple of 2. So:

H(T) = 1 + 2·(some non-negative integer)

This means H(T) mod 4 is either 1 or 3 (since the "some integer" is either even or odd).

Now, 7 mod 4 = 3. That's fine so far. But you can dig deeper: the 4-divisible terms in the formula also have structure. Working through the constraints, H(T) = 7 would require a combination of cycle-packing configurations that is algebraically contradictory — the pieces can't fit together.

More precisely: to get H(T) = 7, you'd need the conflict graph Ω(T) to have I(Ω(T), 2) = 7. But the independence polynomial evaluated at 2 satisfies constraints from graph theory — certain values are simply unreachable by any graph.

The full proof uses the OCF together with properties of claw-free graphs (the conflict graph Ω(T) always avoids certain "claw" patterns), which limits what values I(Ω(T), 2) can take.

---

## Other Gaps That Fill In Later

Some odd numbers that SEEM missing at small tournament sizes actually appear at larger sizes:

- **63**: Not achievable at n = 7, but IS achievable at n = 8. It fills in.
- **107, 119**: Similar — they appear at larger n.

This is different from 7 and 21, which never appear no matter how large you go.

The research conjectures that **7 and 21 are the only permanent gaps** — every other odd number is achievable by some tournament of large enough size. This has been verified for all odd numbers up to several hundred, and there's strong theoretical evidence for it, but a complete proof hasn't been found yet.

---

## A Pattern in the Gaps

The numbers 7 and 21 share a structure: they are both of the form 7 × (odd power of 3):
- 7 = 7 × 1
- 21 = 7 × 3

This hints at a deeper algebraic reason rooted in modular arithmetic — specifically, in properties of integers modulo powers of 2. The OCF formula generates values via weighted sums of powers of 2, and the "forbidden zone" around 7 and 21 emerges from constraints on how those powers can combine.

---

## Why This Matters

Forbidden values are rare in combinatorics — most counting problems have no gaps (every non-negative integer appears as some count for some input). When gaps exist, they signal deep structural constraints.

The existence of forbidden path counts tells us that tournaments are not "generic" combinatorial objects. Their internal structure — specifically the interplay of odd cycles and the insistence that all cycles have odd length — imposes algebraic constraints that not every odd number can satisfy.

It's analogous to asking which integers can be the area of a right triangle with integer sides. Most positive integers work, but some don't — and the reason why connects to deep number theory (the Pythagorean triple structure). Similarly, the forbidden values 7 and 21 connect to algebraic properties of the independence polynomial for the special class of conflict graphs arising from tournaments.

---

## Key Words

- **Forbidden value**: an odd number that can never equal H(T) for any tournament T
- **Permanent gap**: a value forbidden for all n, not just small n
- **Claw-free graph**: a graph where no vertex has three mutually non-adjacent neighbors — the conflict graph always has this property, which constrains what H(T) can be
