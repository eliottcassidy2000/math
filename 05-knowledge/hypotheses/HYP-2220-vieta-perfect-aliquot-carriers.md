---
id: HYP-2220
name: vieta-perfect-aliquot-carriers
status: OPEN
date: 2026-06-04
session: codex-2026-06-04-S644
depends_on: [HYP-2219, HYP-2218, HYP-2217, HYP-2216, HYP-2211, THM-401, THM-361]
---

# HYP-2220: Vieta/Perfect/Aliquot Carriers Expose a 14/21 Shadow

## Claim

The triangular pair-count carrier

```text
A = C(n,2)
```

has an exact Vieta/discriminant shadow:

```text
8A + 1 = (2n - 1)^2.
```

Thus the LRC shell clock `C=2n-1` is not only a modulus chosen by the
antipodal-witness method.  It is the odd square root of the pair-count
discriminant, equivalently the Vieta carrier for

```text
x^2 + x - 2A = 0.
```

Perfect numbers are fixed controls inside the same carrier: even perfect
numbers are exactly the triangular rows

```text
N = C(2^p,2) = 2^(p-1)(2^p - 1)
```

when `2^p-1` is prime, and they satisfy `s(N)=N` for the proper-divisor-sum
map `s(N)=sigma(N)-N`.

The LRC `n=14` row is not perfect, but it has the exact aliquot shadow:

```text
C(14,2) = 91 = 7*13
s(91) = 1 + 7 + 13 = 21
sqrt(8*91 + 1) = 27 = 2*14 - 1.
```

So `21` is a divisor-sum side shadow of the same triangular carrier whose
Vieta square root is the LRC `n=14` shell clock `27`.

## S644 Evidence

`04-computation/vieta_perfect_aliquot_carriers_s644.py` stores
`05-knowledge/results/vieta_perfect_aliquot_carriers_s644.out`.

It verifies:

1. The only integer `<=20000` that is both a doubled prime and a tripled prime
   is `6`.  The proof is exact: if `2p=3q`, then `3|p`, so `p=3`, then `q=2`.
2. The same `6` is the unique distinct positive product-sum resonance from
   THM-361:

```text
6 = 2*3 = 1+2+3 = 1*2*3.
```

3. Even perfect numbers are triangular pair-count fixed points:

```text
6=C(4,2), 28=C(8,2), 496=C(32,2), 8128=C(128,2), ...
```

4. The `n=14` row is:

```text
A=C(14,2)=91, factor(A)=7*13, shell=27, s(A)=21.
```

5. More generally, if `p` and `q=2p-1` are prime, then

```text
C(2p,2) = p(2p-1)
s(C(2p,2)) = 1 + p + (2p-1) = 3p.
```

The first rows are:

```text
p=2  -> n=4  -> C(n,2)=6    -> s=6
p=3  -> n=6  -> C(n,2)=15   -> s=9
p=7  -> n=14 -> C(n,2)=91   -> s=21
p=19 -> n=38 -> C(n,2)=703  -> s=57
p=31 -> n=62 -> C(n,2)=1891 -> s=93
```

This turns the `p=7` diagonal packet from S642 and HYP-2219/S643's duplicate
branch locus into an aliquot-family member, not an isolated coincidence.

## Interpretation

The user asked for "Vieta carriers."  Here the carrier is the hidden root data
of a quadratic.  The pair count `A` alone is a scalar, but the discriminant
square root recovers `2n-1`, the same odd shell that THM-401 uses to organize
LRC witness clocks.

Perfect numbers now have a disciplined role.  They are not claimed to explain
forbidden tournament `H` values.  They are fixed points of the divisor-sum side
channel inside the triangular carrier.  That gives a control family against
which the `n=14` row can be compared:

```text
perfect control:       s(C(n,2)) = C(n,2)
n=14 carrier shadow:   s(C(14,2)) = 21
```

The small number `6` is the base crossing of several channels at once:

```text
6 = doubled prime = tripled prime
6 = first perfect number
6 = C(4,2), shell 7
6 = unique distinct product-sum resonance {1,2,3}
```

This is why `6` should be treated as a carrier crossing/fixed point, not just
as a small-number curiosity.

HYP-2219 supplies the companion-graph side of the same story: an unordered
odd-prime pair `{p,q}` casts one or two odd companions `E+p,E+q`, the duplicate
branch folds at `(E,O)=(2p,3p)`, and the prime-2 Lemoine channel is an apex
boundary.  HYP-2220 adds the divisor-sum observer on the triangular row:

```text
duplicate branch:      p -> (2p,3p)
aliquot semiprime row: C(2p,2) -> 3p, when 2p-1 is prime.
```

## Proof-Route Use

For LRC, this suggests an auxiliary audit over triangular pair-count shadows:

1. Keep `C=2n-1` as the primary shell clock.
2. Add the divisor-sum side channel `s(C(n,2))` as a cheap arithmetic shadow.
3. Flag semiprime rows `n=2p`, `2p-1` prime, where `s(C(n,2))=3p`.
4. Compare these rows with rank-one LRC rows `n=2p`, S642 diagonal packets
   `(2p,3p)`, and unit-distance/tournament `21` guardrails.

For `n=14`, this packages three facts in one carrier:

```text
n=14 = 2*7
C=2n-1 = 27 = 3^3
C(n,2)=91 -> s=21 = 3*7
```

The proof burden remains the LRC `27` carry/lift/owner quotient from HYP-2217
and HYP-2167.  The new contribution is a compact arithmetic observer that may
help decide which clocks deserve proof attention.

## Tournament Analysis

S644 uses proof lenses as vertices.  The proof-lens tournament is transitive:

```text
aliquot_shadow_C14_pair_count_to_21
triangular_vieta_square_root_C_2n_minus_1
even_perfect_numbers_as_triangular_fixed_controls
semiprime_n_2p_family_s_Cn2_equals_3p
doubled_tripled_prime_six_seam
goldbach_lemoine_diagonal_packet
raw_perfect_number_numerology
```

The pairwise observable is: which lens better preserves the hidden carrier
needed for LRC-style proof transfer: pair-count square root, aliquot
fixed/shadow data, prime `2/3/7` seams, and side-channel risk.

## Next Tests

1. Add a `triangular_shadow_ledger(n)` helper to LRC scripts, returning
   `C(n,2)`, `2n-1`, factorization, aliquot class, and the semiprime
   `n=2p,2p-1 prime` flag.
2. Compare known LRC easy/hard rows by aliquot class: perfect, deficient,
   abundant, semiprime-shadow, and prime-shell.
3. Study whether `s(C(n,2))` predicts any channel damage or quotient size in
   the `Res_(2n-1)` ledgers.  The null result would still be useful: it would
   certify the aliquot shadow as a clean metaphor rather than a proof engine.
4. Link the `p=7` packet across S642 and HYP-2217 as:

```text
HYP-2219/S643: (2p,3p)=(14,21) is the duplicate branch
S644:          s(C(2p,2))=3p=21 is the aliquot side shadow
HYP-2217:      LRC n=14 proof carrier remains C=27
```

## See Also

`04-computation/vieta_perfect_aliquot_carriers_s644.py`;
`05-knowledge/results/vieta_perfect_aliquot_carriers_s644.out`;
`07-reflections/vieta-perfect-aliquot-carriers-s644.md`;
HYP-2219; HYP-2218; HYP-2217; HYP-2216; HYP-2211; THM-401; THM-361.
