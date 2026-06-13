---
source: claude-2026-06-03-S607
status: three verified efficient finite targets + an assembled finite reduction (n=14 not closed)
tags: [lonely-runner, n14, finite-target, efficient, pinch, nerve, clock-symmetry, sieve, involution]
---

# Three efficient finite targets for LRC@14

S606 concluded the certifying involution must be *combinatorial* (every clock
symmetry is measure-preserving, hence sign-preserving). This session turns that
into objects that are actually **finite and efficient** to compute on, and shows
how they assemble into a finite reduction.

## Target 1 — the clock-symmetry obstruction is itself finite [verified]

The clock maps that respect the speed structure are `t ↦ u·t` for units `u` with
`u² ≡ 1 (mod 2n)`. For `n=14` (`2n=28`) this is the **finite** set
`{1, 13, 15, 27}` (`27 ≡ −1` = negation). All are measure-preserving, hence
sign-preserving on `p₀`. So the search for a *clock* involution is a finite, dead
search — confirming the certifier must act on the nerve. (The point of stating it
as a finite set: there is nothing more to try on the clock side.)

## Target 2 — the realized nerve is poly-sized, not `2ⁿ` [verified]

The `(★)` sum nominally runs over `2ⁿ` subsets, but only the **realized
close-sets** `{i : t ∈ F_i}` occur — one per breakpoint interval, so `O(Σv)` of
them. For the `n=14` AP: **63 faces** vs `2¹³ = 8192`, with size profile (13
singletons, 32 pairs, …, one full `{1..13}` at `t→0`) and a single maximal facet
(the full set). `p₀` = total length of the depth-0 strata. The pivot/Morse
involution (S606) and its critical cells (the two-block Helly witnesses,
codex-S599) live on *this* poly complex — so "reduce `p₀` to its Helly core" is a
poly-time operation, not a `2ⁿ` enumeration.

## Target 3 — the pinch oracle: loneliness at `O(n²)` times [verified, THM-369]

THM-369 says the optimal time is always a **pair-sum pinch** `t = m/(v_a+v_b)`.
Verified: `M(V) = max over pinches of min_i ‖v_i t‖` matches the exact `M` on
every test, and `M ≥ 1/14 ⇔ LRC@14`. The optimal pinch sits at **bounded
`t·n = O(1)`** — up to the negation mirror `t ↔ 1−t` (the AP's optimum is
`13/14 ≡ 1/14`; `(1,3,4,5,9)` at `1/6`; `(1,5,11,24,25)` at `4/13`, i.e.
`t·n ≈ 1.85`). So a **bounded window** `(0,c/n) ∪ (1−c/n,1)` holds the witness,
giving `O(n²)` candidate pinches — a finite, poly oracle per config. (The exact
`c` is the open S556o first-window conjecture; any fixed `c` keeps it `O(n²)`.)

## The assembled finite reduction

By the sieve (THM-401, modulus `2n−1 = 27`) a config's pinch pattern depends only
on `V mod 27`, so there are **finitely many pinch types**. Then:

> **LRC@14 ⇔ every pinch type is cleared by a bounded-window pinch** (each check
> `O(n²)` on a poly nerve), **except** the single residual type that is the
> measure-zero **AP wall** (witnesses `t = j/(2n)`, S551).

Each ingredient is finite/efficient; what is *not* finite is the AP wall — the
one type whose loneliness is witnessed only at measure-zero points and whose
clearance is the observer-coupling residual (S580o). That is exactly where the
combinatorial sign-reversing involution (S606) would have to act.

## Honest status

Not a proof of LRC@14. What is new and verified: the three targets are finite and
efficient (clock side closed as a finite dead-end; the nerve and the pinch oracle
poly), and they assemble into "finitely many pinch types + the AP wall." This
makes the program concretely computable and isolates the one infinite residual.

## Next (finite, efficient)

1. **Enumerate the pinch types `V mod 27`** for `n=14` (finite) and run the
   `O(n²)` bounded-window oracle on each; tabulate which clear with room and which
   approach `1/14` (the near-AP types). A complete clear of all non-AP types
   would reduce LRC@14 to the AP wall alone.
2. **Pivot-Morse on the realized nerve** of each near-AP type: compute the
   critical cells and check they match codex-S599's singleton/pair walls — a poly
   certificate per type.
3. **The AP-wall residual**: bound the observer-coupling defect (S580o) on the
   single residual type via the nerve's critical cells — the finite form of the
   S556o first-window LP.
