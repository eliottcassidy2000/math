---
id: HYP-2254
title: Easier Frontier Targets For LRC And Unit Distance
status: OPEN
source: codex-2026-06-06-S678
type: comparative proof-target atlas
---

# HYP-2254: Easier Frontier Targets For LRC And Unit Distance

The question "is there an easier value than LRC `n=14` or unit-distance
`n=21`?" has two different answers depending on the predicate.

For a full theorem, the smaller values are the honest easier values:

```text
LRC: total n <= 13 is already proved in the repo/literature window.
Unit distance: small exact-witness rows through the one-slab Moser rows are
much easier than the two-slab frontier.
```

For local carrier lemmas, several later values look easier than the headline
frontiers because their retained side channels are cleaner.

## LRC Split

LRC `n=14` is hard not merely because it is the first open total runner count.
It is hard in the repo's current quotient because

```text
C = 2n - 1 = 27 = 3^3,
```

so the THM-401 shell has `gcd` strata `1,3,9`, the HYP-2169/HYP-2177
doubling-sporadic wall is active, and the carry/owner lift has an apex-debt
side channel.

S678's scorecard finds later values with cleaner *local* `C=2n-1` carriers:

```text
n=15: C=29 prime
n=16: C=31 prime
n=18: C=35 squarefree and no triadic Vstar wall
n=19: n prime and C=37 prime
n=22: C=43 prime
n=24: C=47 prime
```

The prettiest carrier lab is `n=19`: the total denominator is an odd prime and
`C=37` is prime, so both the odd-prime shortcut direction and the repo's
pair-sum shell direction are clean.  This does not imply that proving LRC(19)
is globally easier than finishing LRC(14); it says LRC(19) may be the best
place to prove and stress-test carrier lemmas that avoid the `Res_27` triadic
wall.

## Unit-Distance Split

For the unit-distance side, the first split is between lower-bound/spine
certificates and exact maximum upper bounds.

THM-408 gives explicit unit-spine Moser families:

```text
P_m^+ : n = 8m+6, E = 27m+6 for m>=1
P_m^- : n = 8m+5, E = 27m+3 for m>=1
```

Thus `n=13=P_1^-` and `n=14=P_1^+` are easier one-slab targets than
`n=21=P_2^-`.  They retain the same unit-spine section but avoid the second
full slab and its ear/core side channel.

The centered Eisenstein shell `n=19` is a different kind of easier target:
it is not a THM-408 row, but it has full triangular-lattice symmetry and
therefore looks like the cleanest place to test "unit spine from symmetry"
statements before mixing in Moser ears.

For exact upper-bound work, `n=22=P_2^+` is not easier than `n=21`, despite
having an explicit spine.  In the repo's `n=22` thread, the hard question is
precisely whether the endpoint/ear side channel can add a `61`st edge beyond
the `60`-edge Moser lane.

## S678 Evidence

S678 adds `04-computation/lrc_ud_easier_frontier_targets_s678.py` and stores
`05-knowledge/results/lrc_ud_easier_frontier_targets_s678.out`.

The script ranks LRC values `3 <= n <= 24` by local `C=2n-1` arithmetic:
factorization, folded `<2,-1>` orbits, `gcd` strata, odd-prime `n`, prime `C`,
the triadic shell, and the Vstar doubling-sporadic wall.

It ranks unit-distance values `3 <= n <= 30` by Moser family membership,
slab count, unit-spine lower-bound difficulty, exact-upper difficulty, and
centered Eisenstein shell status.

The route tournament is transitive:

```text
split_local_carrier_from_full_theorem
> LRC_clean_C_later_values
> UD_one_slab_spine_values
> UD_centered_Eisenstein_shells
> raw_smaller_n_is_easier
> literal_14_21_numerology
```

with `directed_3cycles=0` and one Hamiltonian path.

## Next Moves

Use these easier targets as proof laboratories:

1. For LRC, try to prove a "prime `C` no-new-wall" carrier lemma first on
   `n=15`, `16`, `19`, `22`, or `24`, then identify exactly which clause fails
   at `C=27`.
2. For unit distance, prove the one-slab upper/lower side-channel ledger on
   `n=13` and `n=14`, then ask which endpoint-compatible ear invariant first
   appears at `n=21`.
3. For the bridge, use `n=19` on both sides as a clean symmetry lab: LRC has
   prime `n` and prime `C`, while unit distance has the centered Eisenstein
   hex shell.

The practical motto is:

```text
Do not ask only "what is the next n?"
Ask "which side channel is the proof trying to learn?"
```
