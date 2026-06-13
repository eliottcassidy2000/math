---
source: codex-2026-06-03-S592
status: computation + synthesis (closes the non-circulant small-census gap in HYP-2136 through m=9)
tags: [LRC, regular-tournament, dichromatic-number, chi, AP, Paley, vertex-transitive, tournament-analysis, HYP-2136, THM-401]
---

# Full regular census: chi=2 is the AP/interval orbit through nine vertices

The S581o handoff left one concrete gap: the scripts separated the AP/interval
regular tournament from Paley by dichromatic number, but had only computed the
circulant part of the regular census at `m=9` and had not assigned chi to the
non-circulant regular class at `m=7`.

S592 closes that small-census gap.

## Certified census

`lrc_regular_full_census_chi_s592.py` counts labelled regular tournaments exactly
and discovers unlabeled classes until the orbit-size sum

```
sum_T m! / |Aut(T)|
```

equals the exact labelled count.  That makes the non-circulant part a completed
census, not a sample.

| m | labelled regular tournaments | regular classes | unique chi=2 class | non-AP chi=2 classes |
|---|---:|---:|---|---:|
| 5 | 24 | 1 | AP/interval | 0 |
| 7 | 2640 | 3 | AP/interval | 0 |
| 9 | 3230080 | 15 | AP/interval | 0 |

The key small rows:

| m | class | aut | VT | chi | H | c3 |
|---|---|---:|:--:|---:|---:|---:|
| 7 | AP/interval | 7 | yes | 2 | 175 | 14 |
| 7 | Paley/QR | 21 | yes | 3 | 189 | 14 |
| 7 | non-circulant middle | 3 | no | 3 | 171 | 14 |
| 9 | AP/interval | 9 | yes | 2 | 3267 | 30 |
| 9 | all 14 non-AP regular classes | 1,3,9,81 | mixed | 3 | 3159..3357 | 30 |

So the first regular class beyond AP/Paley exists already at `m=7`, but it has
`chi=3`, not `chi=2`.  At `m=9`, every one of the 14 non-AP regular classes also
has `chi=3`.

## What this answers

1. **Does chi add beyond vertex-transitivity?** Yes.  At `m=7`, AP/interval and
   Paley are both vertex-transitive, regular, self-converse, and have the same
   number of 3-cycles (`14`), but chi separates them: `2` versus `3`.
2. **Is there a tight regular config that is regular but not AP/Paley?** In the
   completed small regular census through `m=9`, there are many regular non-AP
   classes, but none are `chi=2`.  Combined with S581b's tight-speed audit through
   `n=8`, the tight regular orbit remains exactly AP/interval.
3. **Does its chi differ?** The non-AP regular classes differ upward: `chi=3`.
   The LRC-tight regular orbit is the minimally cyclic regular class, not the
   H-maximal or Paley-like class.

This is a useful correction to the phrase "maximally cyclic."  Regularity makes
the 3-cycle count constant, but chi still sees how those cycles pack.  AP is
regular yet dichromatic-2, hence closest to transitive.  Paley and the other
regular classes are more robustly cyclic.

## Link back to the pair-sum modulus

THM-401 pins the additive face arithmetically: the pair-sum sieve's floor
modulus is `C=2n-1`.  S591 translates the same point into tournament language:
the LRC comparator is the interval/round beat, not the multiplicative QR beat.
S592 reinforces that translation on the full small regular census.  The
additive interval/AP class is not merely one regular option; through `m=9` it is
the only regular option with `chi=2`.

So the slogan sharpens to:

> the additive face is the `2n-1` pair-sum modulus, and its regular tournament
> shadow is the interval/AP orbit; Paley is a multiplicative symmetry class, not
> the LRC beat.

## Tournament Analysis note

The TA vertices in S592 are regular isomorphism classes, not runners, arcs, or
gap sections.  The quotient preserves class-level invariants (`chi`, `H`, `Aut`,
regularity) and destroys labelled observer/source geometry, wall ownership, and
pinch denominators.  The challenged assumption is exactly the tempting one:
regular or vertex-transitive does not imply Paley-like, and high Hamiltonian
path count is not the LRC tightness signal.  The class-level gauge "lower chi,
then higher H" makes AP the unique source through `m=9`; relative to pure
`H`-ordering it flips one edge at `m=7` and eight edges at `m=9`.

## Honest limit

This is not an all-order theorem.  It proves the completed census through
`m=9`; the open theorem is still:

```
regular tournament has chi=2  iff  it is the interval/AP regular orbit.
```

The LRC-specific extension is also still open: every tight speed configuration
should realize a `chi=2` near-AP/round tournament, while regular classes with
`chi>=3` should be strictly lonely.

Namespace note: the S581o log used `HYP-2135` for this chi question; the live
hypothesis is renumbered to `HYP-2136`.  The same-day S591 sumset-support
calculus is now live as `HYP-2141`.

Artifacts: `04-computation/lrc_regular_full_census_chi_s592.py`,
`05-knowledge/results/lrc_regular_full_census_chi_s592.out`.
