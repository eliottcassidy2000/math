# Odd-prime Paley objects and a corrected small-prime comparison for LRC(14)

> **CORRECTION (MISTAKE-228).** The first S215 synthesis promoted an analogy
> into a law.  Only the odd-prime Paley dichotomy below is exact.  The prime
> `2` is outside that theorem; the Paley adjacency spectrum is a shifted,
> half-scaled Gauss sum and also has a principal eigenvalue; `Φ₃` and `Φ₆`
> differ by a sign twist; and no small-prime label in the LRC table implies an
> LRC statement.  The original exploration is retained here as a disciplined
> source of prompts, not as canon.

*boxeph-2026-07-21-S215. Replayed by
`04-computation/understanding_primes_via_paley_tournaments_boxeph_S215.py`.*

## Exact kernel: the odd-prime Paley dichotomy

For an odd prime `p`, put `i → j` when `j-i` is a nonzero quadratic
residue modulo `p`.

- If `p ≡ 3 (mod 4)`, this is a Paley tournament.  Its adjacency
  eigenvalues are the principal value `(p-1)/2` and
  `(-1 ± i√p)/2`, the latter each with multiplicity `(p-1)/2`.
- If `p ≡ 1 (mod 4)`, the relation is symmetric and gives a Paley graph.
  Its adjacency eigenvalues are the principal value `(p-1)/2` and
  `(-1 ± √p)/2`, again with the usual equal nonprincipal multiplicities.
- If `g_p` is the quadratic Gauss sum, then `g_p²=(-1)^((p-1)/2)p`
  (up to the conventional choice of sign for `g_p`).  The nonprincipal
  eigenvalues are `(-1 ± g_p)/2`: they are not `g_p` itself.
- `p=2` is outside this construction.  Calling reversal an involution can be
  useful elsewhere, but it does not produce a Paley law or an `i√2` datum.

The replay checks the tournament characteristic polynomials for `p=3,7,11`,
the graph/tournament dichotomy for `p=3,5,7,11,13`, and the numerical Gauss
sums.  In particular,

`char_A=(x-(p-1)/2)(x²+x+(p+1)/4)^((p-1)/2)`

for the checked tournament primes.  This exact kernel is the survivor of the
original reflection.

## The `Φ₃` / `Φ₆` boundary

Paley-3 is the directed 3-cycle, whose nonprincipal eigenvalues satisfy
`Φ₃(x)=x²+x+1`.  The LRC parameter recalled by S215 uses
`Φ₆(x)=x²-x+1`, and

`Φ₆(x)=Φ₃(-x)`.

That sign twist is essential.  A shared quadratic field does not identify the
two cyclotomic objects or transfer an extremal statement between them.
Likewise, the Paley-7 quadratic field `Q(√-7)` is not the real cubic field
`Q(cos(2π/7))`; placing both under the label “7” loses type and degree.

## Historical small-prime table, now typed as heuristic

| prime | exact datum used in the comparison | historical prompt | status for LRC(14) |
|---|---|---|---|
| `2` | `Φ₂=x+1`; reversal is an involution | chirality | analogy only; not a Paley object |
| `3` | Paley-3 is a 3-cycle; `Φ₆(14)=183` | atom / argmax | separate facts joined only by `x↦-x` |
| `5` | Paley graph; `g_5=√5`; golden field | Fibonacci foil | heuristic; no tight/slack theorem follows |
| `7` | Paley-7; `x¹⁴-1≡(x-1)⁷(x+1)⁷ mod 7` | apex | the congruence is exact; the label adds no proof |
| `11` | Paley-11; `rank ker_Z(1,…,12)=11` | rank / scarcity | numerical coincidence; no map or forcing result |

The script also checks that `2cos(2π/7)` satisfies
`x³+x²-2x-1` and that `11` is its only positive multiple at most `14`.
Those are true observations, but neither one bridges a Paley object to the
lonely-runner predicate.  In particular, “scarce” does not force speed `11`,
and equality of ranks is not a map between relation lattices.

## What would turn a prompt into mathematics?

For any proposed Paley/LRC connection, supply all of the following before
promotion:

1. the LRC-side vertex set and the Paley-side vertex set;
2. an explicit map between them, not a shared number or field label;
3. the LRC predicate preserved by that map;
4. the information destroyed (phase, coefficient signs, saturation, or
   multiplicity) and a sidecar carrying it;
5. a hostile example distinguishing correlation from implication.

The historical “periodic table” remains useful as a question generator:
quadratic-residue orientations may expose symmetry, and the mod-`4` split may
suggest comparisons.  Its strongest warranted conclusion is narrower:

> Odd primes have a clean Paley tournament/graph dichotomy, and Gauss sums
> control the nonprincipal spectrum after a shift and half-scaling.  Whether
> any such object preserves a useful LRC(14) predicate is open in this thread.

Links: HYP-8860, THM-1830, THM-2043,
[[the-rank11-ap-core-is-the-achiral-vertex-where-the-rank-or-euler-frontier-meets-boxeph-S214]],
[[what-an-lrc14-disproof-must-be-and-why-fibonacci-is-the-foil-boxeph-S206]].
