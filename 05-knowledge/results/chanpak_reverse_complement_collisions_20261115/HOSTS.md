# Prescribed absolute counts: certified carrier hosts for 12- to 16-digit primes

**Status: DEMONSTRATION (each instance PROVED by the closed form). This is not a coverage theorem; Chan–Pak remains
OPEN.**

## The host and the certificate

The host is the two-tape carrier of THM-5192 (A):
- a marked guard;
- an X-tape;
- a strong marked core C of order 3–5;
- a Y-tape;
- an unmarked guard;
- one apex.

Its Hamiltonian-path count is the PROVED closed form

    F_C(X, Y) = qXY + tX + sY + a(p(X) + p(Y)) + n + f·[X odd] + e·[Y odd].

So a solution (C, X, Y) of F_C(X, Y) = N is a complete certificate that the explicit tournament 𝒯_C(X, Y) has
exactly N Hamiltonian paths.

`hostsolve.c` re-derives each hit a second time. It runs the independent cut-transfer recurrence through every
tape vertex (guard, X-digits, block, Y-digits, guard). It aborts if either count differs from N.

## Results

| target N | digits | core (order; mark mask; out-masks) | X | Y | vertices | time |
|---|---:|---|---:|---:|---:|---:|
| 999,999,999,989 | 12 | 3; 1; 2 4 1 (directed triangle, one mark) | 550,071 | 151,495 | 42 | < 1 s |
| 9,999,999,999,971 | 13 | same | 101,405 | 8,217,845 | 44 | 0.9 s |
| 99,999,999,999,973 | 14 | same | 2,162,403 | 3,853,736 | 48 | 3.8 s |
| 999,999,999,999,989 | 15 | same | 1,206,503 | 69,070,121 | 52 | 10.6 s |
| 9,999,999,999,999,937 | 16 | 4; 6; 2 12 9 1 | 13,568,255 | 23,031,700 | 54 | 243 s |

- The targets are the largest primes below 10^12, ..., 10^16.
- The previous exhaustive range of the 168-form bank is [333, 10^9] (THM-5192 (D)).

**End-to-end controls** (`bf_check.py`): for the targets 10007, 20011, 65537 and 100003 the literal tournaments
(15–18 vertices) were built and their Hamiltonian paths counted by subset DP (`hp2.c`). All four agree. For
example, the Fermat prime 65537 is H of an 18-vertex carrier with the one-mark triangle, X = 65 and Y = 83.

## The solver and its cost

Fix C and the column X. The row c Y + a p(Y) + e[Y odd] = M, with c = qX + s and M the remaining target, confines Y
to the interval [(M − e)/(c + a), M/c], because 0 ≤ p ≤ Y. The solver scans the columns where this interval has
O(1) length, from X ≈ √(aN)/q up to 4√(N/q), and tests each candidate Y exactly.

The first moment (THM-5192 (D)) predicts about (2/q)·ln(X_hi/X_lo) hits per form in this window. That is about one
hit for the triangle forms, in line with the order-3 forms sufficing up to 15 digits.

Finding one hit costs about √(aN). This is the same order as the concurrent deterministic inverse (O(√T log²T),
[TWISTED_LATTICES.md](../chanpak_cone_collision_placement_20261115/TWISTED_LATTICES.md) §4, credited), which also
certifies the absence of a solution. Neither column scanning nor the 2-adic row recursion beats √N per expected
hit, for the following reasons:
- the row recursion is (c, f, M) ↦ (2c + a(1−b), a(1−b), M − (c+f)b), the carrier forms' closure under peeling one
  low digit;
- its live branches peak at the square root of the window;
- the hit probability per column scales the same way.

## What this lane shows, and what it does not

- *Shown.* A prescribed absolute count, chosen far beyond every census, is hit by an explicit host found by
  arithmetic alone. The host's two tapes are fixed by the target through the proved closed form, not by search
  over tournaments.
- *Not shown.* That every odd N has a hit; that is HYP-9722. The reverse-complement family in RC_COLLISIONS.md
  controls the *second* moment of exactly these rows. A pointwise first-moment statement is what Chan–Pak needs,
  and it remains open.
