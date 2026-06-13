# HYP-2218: Goldbach/Lemoine Same-Pair Projection Has a Diagonal 14/21 Shadow

**Status:** OPEN synthesis, supported by S642 and parent doubled-prime work
HYP-2044, HYP-2049, HYP-2051, HYP-2215, HYP-2216, and HYP-2217.

## Claim

Even Goldbach and odd Lemoine/Levy should be compared in prime-pair space, not
only target-number space.  For odd primes `p,q`, define two projections:

```text
E = p + q
O = p + 2q
```

The even projection `E` is swap-blind.  The odd projection `O` is ordered
because it remembers which coordinate was doubled.  Together they reconstruct
the ordered pair:

```text
q = O - E
p = 2E - O
```

So an even `E` and odd `O` use the same ordered prime pair exactly when
`O-E` and `2E-O` are odd primes.  Necessarily `E < O < 2E`.

This makes the user's duplicate pairs precise.  The diagonal `p=q` maps

```text
p -> (2p, 3p).
```

The prime `7` gives the visible bridge:

```text
7+7 = 14
7+2*7 = 21.
```

Thus the pair `(14,21)` is a Goldbach/Lemoine diagonal shadow of the prime `7`.
This sharpens, but does not replace, HYP-2217's `27`-quantum section/bulk
bridge for LRC `n=14` and unit-distance `n=21`.

## S642 Evidence

`04-computation/goldbach_lemoine_pair_projection_s642.py` stores
`05-knowledge/results/goldbach_lemoine_pair_projection_s642.out`.

Through `LIMIT=300`, S642 verifies:

```text
Goldbach misses among even 4..300: []
Lemoine misses among odd 7..300: []
Odd-prime-compatible Lemoine misses: [7]
q=2-only odd exceptions: [7]
```

The last line matters.  Lemoine representations using `q=2` are an exceptional
prime-2 channel: they can prove an odd row, but they do not correspond to an
even same-pair Goldbach row with `E=p+q`, because one coordinate is even.

For `E=14`, the same-pair shadows are:

```text
E=14: O=17 via (11,3)
      O=21 via (7,7) diagonal
      O=25 via (3,11)
```

The off-diagonal odd shadows are paired by reflection; the diagonal is fixed.

## Reflection Law

Swapping the prime pair fixes `E` and reflects `O`:

```text
(p,q) <-> (q,p)
O      <-> 3E - O.
```

Equivalently:

```text
O + O' = 3E.
```

The duplicate diagonal is the fixed locus:

```text
O = 3E/2.
```

So the unordered/ordered distinction is not cosmetic.  It is a small
anti-coset: the even projection quotients by the swap, while the odd projection
chooses one side of the swap.  Adding both projections recovers the pair.

## Transfer Reading

This is the additive-number-theory analogue of the repo's recurring
two-shadow carrier:

```text
sum/product for pi,e         -> reconstruct unordered algebraic pair
E/O for Goldbach/Lemoine     -> reconstruct ordered prime pair
source/sink or carry/owner   -> reconstruct LRC side-channel data
spine/bulk                   -> reconstruct unit-distance carrier data
```

The diagonal `p=q` is the fixed-point locus under swap.  It is where the
Goldbach even row, the Lemoine odd row, and the multiplicative doubled-prime
row agree:

```text
p+p = 2p
p+2p = 3p
```

For `p=7`, this puts `14` and `21` in the same diagonal packet.  Combined with
HYP-2217, the current picture is two-layered:

1. `14/21` is a prime-7 diagonal shadow in the Goldbach/Lemoine pair map.
2. The LRC/unit-distance proof difficulty still lives in retained side
   channels: `C=27` shell/carry on the LRC side and `57=20+37` spine/bulk on
   the unit-distance side.

HYP-2220 adds a third compatible shadow: in the triangular pair-count carrier,
`C(14,2)=91=7*13` has aliquot shadow `s(91)=21`, while the same row has Vieta
shell root `sqrt(8*91+1)=27`.

## Tournament Analysis

S642 uses proof lenses as vertices.  The tournament is transitive with one
Hamiltonian path:

```text
invertible_EO_pair_projection
swap_reflection_O_to_3E_minus_O
diagonal_duplicate_fixed_locus
exceptional_prime_2_channel
degree_and_shadow_price_ledger
raw_goldbach_lemoine_counts
```

The ranking says the best object is not the target number alone.  It is the
invertible pair projection plus the swap/reflection side channel.

## Next Tests

1. Build a larger same-pair graph and study even-degree versus odd-degree
   scarcity after separating the prime-2 exception channel.
2. Compare the diagonal packet `p -> (2p,3p)` against LRC rank-one rows
   `n=2p` and tournament forbidden packets `H=3p`.
3. Add a local residue ledger: for each modulus `m`, record which residues of
   `(E,O)` preserve the invertible pair condition.
4. Use the reflection law as a toy anti-coset model for HYP-2208: `E` is the
   quotient, the two odd shadows are the anti-coset orbit, and the diagonal is
   the fixed locus.

## See Also

`04-computation/goldbach_lemoine_pair_projection_s642.py`;
`05-knowledge/results/goldbach_lemoine_pair_projection_s642.out`;
`07-reflections/goldbach-lemoine-same-pair-projection-s642.md`;
HYP-2220; HYP-2217; HYP-2216; HYP-2215; HYP-2211; HYP-2208; HYP-2051;
HYP-2049; HYP-2044.
