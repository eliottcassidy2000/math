# Easier Frontier Targets For LRC And Unit Distance (S678)

Prompt: see if there are values of the LRC or unit distance problem that would
be easier to prove than `14` or `21` respectively.

The useful answer is not a single number.  It is a split between full theorem
difficulty and carrier-lemma difficulty.

## LRC

For the full Lonely Runner theorem, `n<=13` is the genuinely easier/proved
range, and `n=14` remains the immediate frontier.  But the `n=14` carrier is
not arithmetically generic.  It sits at

```text
C = 2n - 1 = 27 = 3^3.
```

That means:

- the folded THM-401 shell has `gcd` strata `1,3,9`;
- the HYP-2169/HYP-2177 Vstar doubling wall is active;
- the carry lift has the `14|v iff k==r mod 14` apex-debt side channel.

So if the target is a clean carrier lemma, not the whole theorem, several
later values are easier laboratories:

```text
n=15: C=29 prime
n=16: C=31 prime
n=18: C=35 squarefree, no Vstar triadic wall
n=19: n prime and C=37 prime
n=22: C=43 prime
n=24: C=47 prime
```

The standout is `n=19`.  It is bigger than `n=14`, but it is cleaner in both
directions: the total denominator is an odd prime, and the repo's pair-sum
modulus `C=37` is also prime.  That makes it a promising "no side-channel"
control problem.  If a proposed LRC14 proof mechanism cannot first say
something sharp for `n=19`, the mechanism may be overfitted to `Res_27`.

## Unit Distance

For unit distance, "easier" splits even more sharply.

For explicit lower-bound or unit-spine certificates, `n=13` and `n=14` are the
clean values.  THM-408 gives

```text
P_1^- : n=13, E=30
P_1^+ : n=14, E=33
P_2^- : n=21, E=57
P_2^+ : n=22, E=60
```

The first two are one-slab rows.  The `n=21` row is the first minus-family
two-slab row, so it has the same `27` quantum but a second slab's worth of
bulk/ear side channel.

For exact upper-bound work, `n=22` is not easier than `n=21`, even though its
spine is explicit.  The repo's `n=22` problem is precisely the question of
whether an endpoint-compatible ear can push the `60`-edge Moser lane to `61`.

The third useful unit-distance value is `n=19`, not because it is a Moser row,
but because it is a centered Eisenstein shell.  It is the clean symmetry lab
for traceability and unit-spine statements before Moser ears enter.

## Shared Nineteen

The nicest surprise is that `19` is a good side target in both problems.

On the LRC side:

```text
n=19, C=37 prime.
```

On the unit-distance side:

```text
n=19 = 1 + 3*2*3,
```

the centered Eisenstein hex shell of radius `2`.

So `19` is not a replacement for the `14/21` frontiers.  It is a symmetry
control.  It can test which parts of a proof are genuinely structural before
the `27` side channels start talking over everything.

## Tournament Analysis

Vertices are proof obligations, not runners or points.  I considered runners,
residues mod `2n-1`, `gcd` shells, fixed clocks, Moser points, slabs, unit
directions, ears, centered Eisenstein shells, and exact-upper-bound
obstructions.

The target-selection route tournament is transitive:

```text
split_local_carrier_from_full_theorem
> LRC_clean_C_later_values
> UD_one_slab_spine_values
> UD_centered_Eisenstein_shells
> raw_smaller_n_is_easier
> literal_14_21_numerology
```

That ordering feels right.  The scalar coincidence `14/21` is the least useful
part.  The useful part is knowing what payload the proof is allowed to forget.

## Next

The next experiment should not try to "prove LRC19" wholesale.  It should prove
one named clause in a prime-`C` lane:

```text
prime C no-new-wall lemma
```

Then run the same clause at `C=27` and watch exactly where it breaks.  If it
breaks only at triadic carry/apex debt, it plugs directly into HYP-2253.

For unit distance, the analogous move is:

```text
one-slab upper/lower side-channel ledger at n=13,14
```

Then add the second slab and isolate the first invariant that is no longer
local.
