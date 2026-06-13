# HYP-2209: Hard Counting Sequences Should Be Extended by Shadow Recursions, Not Only by Next Terms

**Status:** OPEN synthesis, supported by S633 exact shadow-sequence scout.

## Claim

When a difficult counting sequence is hard to extend directly, useful progress
often comes from computing nearby sequences that retain one structural side
channel:

- fixed layer under an involution;
- merged layer under the same involution;
- nonfixed-pair layer;
- q-deformed Burnside pressure;
- parity or fixed-point bisection;
- skinny quotient preserving a geometric/proof predicate;
- transporter quotient `Anti(X)=Iso(X,JX)`.

The point is not that these shadows replace the original hard sequence.  They
give recursive coordinates around it, so a hard value can be approached through
which side channel is missing.

## Evidence

S633 adds `04-computation/sequence_shadow_lab_s633.py` and stores
`05-knowledge/results/sequence_shadow_lab_s633.out`.

The motivating fixed sequence is self-converse tournament classes:

```text
SC(1..10) = 1, 1, 2, 2, 8, 12, 88, 176, 2752, 8784.
```

The script records the mirror table:

```text
T(n)        = A000568(n)
SC(n)       = fixed layer under converse
M(n)        = (T(n)+SC(n))/2
N(n)        = (T(n)-SC(n))/2
```

It also reuses the HYP-2074 bisection:

```text
SC(2m)   = A(m,4)
SC(2m+1) = sum_{odd lambda |- m} 2^ell(lambda) * 4^c(lambda) / z(lambda).
```

So the odd bisection is not mysterious once the fixed vertex is named: it is
the even base-4 Burnside shadow with a `2^(#parts)` tax.

## Recursive Feel

The dominant singleton-partition heuristic gives:

```text
SC(2m+2)/SC(2m)       ~ 4^m/(m+1)
SC(2m+3)/SC(2m+1)     ~ 2*4^m/(m+1)
SC(2m+1)/SC(2m)       ~ 2^m
```

This is not a proof-quality asymptotic claim for small `m`; it is a useful
growth currency.  It says the fixed-point tax is multiplicative and partition
local, not a random correction.

## LRC And Unit-Distance Shadows

The same scout records thinner shadows connected to the repo's LRC and
unit-distance threads:

- round/self-converse LRC shadow: `2^floor((m-1)/2)`;
- folded shell orbit counts under `<2,-1>` on `Z/(2n-1)`;
- centered-hex markers `1,7,19,...` as geometry-side shell events.

At `n=14`, the transporter shell shadow again exposes the familiar
`C=27` gcd strata:

```text
gcd reps: [1, 3, 9].
```

This connects S633 to S632/HYP-2208: difficult counts become more usable when
the transporter and quotient-loss side channel is named.

## Assumption Challenge

Alternate vertices considered: raw values, partitions, automorphism groups,
anti-cosets, LRC shells, q-colors, unit-distance shells, and proof obligations.
S633 chooses extension methods as Tournament Analysis vertices.

Preserved predicate: each hard value is surrounded by fixed, merged,
bisection, q-shadow, skinny, or transporter companions.

Destroyed data: individual representatives, exact embeddings, and full
Hamiltonian-path spectra.

The challenged assumption is that progress on a sequence means computing only
the next term of that sequence.  Sometimes the better next term is a nearby
sequence whose recursion exposes the missing side channel.

## Next Tests

1. Apply the fixed/merged/nonfixed/shadow recipe to strong-tournament
   `H`-spectra rather than class counts.
2. Add a q-shadow table for round/LRC necklaces and compare it to A000568's
   odd-cycle Burnside q-shadow.
3. For unit-distance search beams, record `total/fixed/merged` under available
   geometric symmetry groups before pruning.
4. Turn S632's transporter records into a reusable helper: given `X` and `J`,
   report fixed count, merged count, transporter count, and quotient-loss notes.
