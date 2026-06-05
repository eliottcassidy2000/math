# HYP-2227: Number-Theory Tournaments Are Local-Witness Carrier Quotients

**Status:** OPEN method hypothesis with finite atlas evidence.

## Claim

Tournaments fit number theory best when their vertices are not raw integers,
but local witnesses, residues, phases, valuations, or proof obligations.  A
tournament is useful when each edge records a tie-free local comparison:

```text
which side channel preserves more obstruction information?
which local prime is the stronger bottleneck?
which proof obligation explains the other without scalar collapse?
```

This is the common carrier behind several repo threads:

```text
Paley/character tournaments       -> finite-field characters as orientations
local sieve tournaments           -> primes/residue lanes ranked by obstruction pressure
endpoint horizon tournaments      -> boundary/interior witness slots, as in S650
valuation drift tournaments       -> residue classes ranked by descent pressure
explicit formula tournaments      -> prime-power/zero terms ranked by phase contribution
local-global rank tournaments     -> BSD-style local obligations ranked by retained data
```

The morale-progress use for hard problems is to replace a global, seemingly
opaque conjecture with finite side-channel ledgers:

```text
name the vertex set;
name the pairwise observable;
name the switch/gauge;
measure the quotient loss;
only then ask for a global theorem.
```

## Evidence From S651

`04-computation/number_theory_tournament_atlas_s651.py` builds three concrete
families.

### Paley Character Tournaments

For primes `p == 3 mod 4`, the Paley orientation

```text
i -> j iff j-i is a quadratic residue mod p
```

is a literal finite-field tournament:

```text
p=3:  score_hist={1:3},  c3=1,   H=3
p=7:  score_hist={3:7},  c3=14,  H=189
p=11: score_hist={5:11}, c3=55,  H=95095
p=19: score_hist={9:19}, c3=285, H skipped
```

This is not analogy.  It is number theory making a tournament from a
multiplicative character.

### Local Sieve Obstruction Tournaments

S651 orients local primes by obstruction weight `-log(survivors/prime)`.
The same vertex set of local primes reorders when the side channel changes:

```text
twin gap 2 ranking:       [3,2,5,7,11,13,17,19,...]
Goldbach N=210 ranking:   [2,3,5,11,13,7,17,19,...]
Goldbach N=2110 ranking:  [3,2,7,5,11,13,17,19,...]
Euler p=41 ranking:       [41,43,47,53,61,2,3,5,...]
```

The edge-flip counts between these local-prime tournaments are:

```text
twin gap 2 vs Goldbach N=210:  3
twin gap 2 vs Goldbach N=2110: 1
Goldbach 210 vs Goldbach 2110: 4
twin gap 2 vs Euler p=41:      62
```

The large twin-vs-Euler flip count is the point: the primes are not the whole
object.  The side channel changes the tournament.

### Hard-Problem Transfer Tournament

S651 ranks proof-carrier methods by local witness retention, side-channel
retention, finite probe quality, hard-problem reach, and overclaim risk.  The
resulting Tournament Analysis is transitive with one Hamiltonian path:

```text
local_sieve_obstruction_tournaments
> character_phase_paley_tournaments
> endpoint_horizon_witness_tournaments
> valuation_drift_tournaments
> explicit_formula_balance_tournaments
> local_global_rank_tournaments
> raw_conjecture_numerology
```

Fingerprints:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}
directed_3cycles=0
scc_sizes=[1,1,1,1,1,1,1]
hamiltonian_paths=1
```

## Hard-Problem Morale Moves

This is not a claimed solution to RH, BSD, twin primes, Goldbach, Collatz, or
Heegner uniqueness.  It is a way to make progress that is small, checkable,
and resistant to numerology.

For twin primes:

```text
vertices = local primes or residue lanes;
edge = stronger obstruction to a fixed gap;
next finite task = compare gap tournaments and singular-series side channels.
```

For Goldbach:

```text
vertices = local primes for the target N;
edge = one-lane/two-lane obstruction pressure;
next finite task = track how divisors of N reorder the local-prime tournament.
```

For RH:

```text
vertices = prime powers, zero phases, or explicit-formula terms;
edge = larger contribution to a smoothed error window;
next finite task = find stable cancellation tournaments under smoothing changes.
```

For BSD:

```text
vertices = bad primes, root-number factors, Selmer/rank proof obligations;
edge = which local datum forces more global rank/parity information;
next finite task = local-global rank tournament for small elliptic-curve families.
```

For Collatz:

```text
vertices = residue classes or valuation jumps;
edge = stronger descent/drift after one accelerated step;
next finite task = compare residue drift tournaments across moduli.
```

For Heegner prime polynomials:

```text
vertices = boundary/interior proof slots;
edge = no-early-collapse witness priority;
next finite task = extend S650's endpoint horizon to non-Heegner failure profiles.
```

For LRC and unit distance:

```text
vertices = shell/spine/bulk/carry/endpoint side channels;
edge = retained proof information;
next finite task = prove no-leak lemmas after a side-channel jackknife.
```

## Assumption Challenge

For number theory, the tempting but weak vertex set is "numbers" or "primes."
S651 explicitly challenges that.  Useful vertex sets include:

- residues;
- local primes;
- residue lanes;
- valuations;
- characters;
- prime powers;
- zero terms;
- class groups;
- discriminants;
- boundary/interior witness slots;
- local-global proof obligations;
- quotient-loss side channels.

The quotient preserved by the S651 atlas is local obstruction order.  It can
destroy actual arithmetic dependence across primes, analytic cancellation, and
global existence data.  Therefore the tournament should be treated as a
side-channel audit, not as a replacement for the theorem.

## See

`04-computation/number_theory_tournament_atlas_s651.py`;
`05-knowledge/results/number_theory_tournament_atlas_s651.out`;
`07-reflections/number-theory-tournament-carriers-s651.md`;
HYP-2226, HYP-2225, HYP-2224, HYP-2223, HYP-2217, HYP-2216, HYP-2215,
THM-410, THM-115.
