# HYP-2216: Local Obstruction Product Ledgers Refine Side-Channel Jackknife

**Status:** OPEN synthesis, extending HYP-2215 with a concrete helper target.

## Claim

HYP-2215 names the shared transfer move: side-channel jackknife over a raw
scalar quotient, hidden witness, and retained side channels.  HYP-2216 refines
that move into a reusable ledger: price each local obstruction channel as a
factor, record what scalar would be unsafe if the channel were forgotten, and
then isolate the transverse branch that the product ledger cannot decide.

The common structure behind the current cluster

```text
H=21 impossibility,
LRC n=14,
unit-distance n=21,
A000568/self-converse counts,
Schanuel/pi-e algebraic shadows,
twin primes,
Goldbach
```

is not a single theorem.  It is a proof-design schema:

```text
scalar quotient + retained side channel + residual obstruction.
```

The mistake pattern is to identify the scalar with the object.  The progress
pattern is to make the side channel first-class.

## Atlas

| domain | tempting scalar | retained side channel |
|--------|-----------------|-----------------------|
| `H=21` | Hamiltonian path count `21` | strong-component norm, `c3`/odd-cycle budget, OCF packet |
| LRC `n=14` | `C=2n-1=27` | gcd strata `1/3/9`, unit/nonunit shells, lift/carry/owner data |
| unit distance `n=21` | `21` vertices, `57` edges | `20`-edge spine plus `37=C_hex(3)` bulk, traceability, ears |
| A000568 | unrooted tournament count | Burnside partitions, fixed/merged/nonfixed, q-deformation, transporters |
| Schanuel/pi-e | rational/algebraic shadow | algebraic independence, trace/norm/discriminant, transverse shadows |
| twin primes | prime gap `2` | admissible tuples, singular series, sieve weights, parity barrier |
| Goldbach | `N=p+q` or `N=p+p+p` | major/minor arcs, local congruences, singular series, density |

This table should be read operationally: each row asks which side channel must
survive before a quotient is allowed.

## S640 Evidence

S640 adds `04-computation/resonance_product_ledger_s640.py` and stores
`05-knowledge/results/resonance_product_ledger_s640.out`.

The script factors recurring integers:

```text
21 = 3*7
27 = 3^3
57 = 3*19
37 = 37
246 = 2*3*41
456 = 2^3*3*19
```

The point is not numerology.  It is prime-local side-channel pressure:

- prime `2` is the parity/traceability seam: Redei oddness, twin gap `2`,
  Goldbach parity, LRC even-row first doubling, dyadic carry data;
- prime `3` is the first fragmentation prime: LRC `C=27=3^3`, ternary
  Goldbach, and the `Phi_3` carrier;
- prime `7` is the forbidden-scalar warning: `H=7` and `H=21` are impossible
  as tournament path counts, while unit-distance `n=21` is legal only with
  spine/bulk side channels retained.

S640's Tournament Analysis uses next-proof routes as vertices, not problem
domains.  The majority tournament is transitive:

```text
local_obstruction_product_ledger
scalar_twins_side_channel_atlas
transverse_shadow_fallout_helper
discriminant_branch_labels_for_LRC
A000568_prime_sieve_shadow
H_spectrum_same_scalar_fibers
unit_distance_ear_singular_series
Schanuel_conditional_transfer
raw_scalar_numerology
```

The top route is the important one.

## Concrete Next Build

Build a reusable `local_obstruction_product_ledger` helper with a common schema:

```text
domain
local modulus / prime
forbidden residue classes
surviving residue classes
local weight
lost side channel
```

First instances:

1. LRC `C=27` shell strata under `<2,-1>` and gcd classes `1/3/9`.
2. Twin-prime admissible tuples: for each prime `p`, forbidden residues of
   `n` such that one translate is `0 mod p`.
3. Goldbach local obstructions: parity and odd-prime residue survival for
   `N=p+q` or `N=p+p+p`.
4. Unit-distance direction masks: local direction/support admissibility before
   embedding search.
5. A000568/Burnside cycle types: local cycle-part constraints as a tournament
   analogue of admissible residue data.

The goal is to make "singular series" a repo-native object, not only an
analytic-number-theory phrase.

## Applications

### Same-Scalar Twin Atlas

Search for pairs agreeing in a scalar and disagreeing in a side channel:

- same tournament `H`, different OCF packet;
- same unit-distance edge count, different unit spines/direction masks;
- same LRC `C=2n-1` shell, different carry owner;
- same A000568 raw count layer, different fixed/merged/nonfixed companion;
- same algebraic shadow status, different transverse fallout.

This gives small counterexamples to unsafe quotienting before the quotient is
used in a proof.

### Transverse-Shadow Helper

Generalize S638: once one scalar descends, identify which transverse shadows
are forced to remain non-tame unless the hidden source also descends.  The
pi/e power-sum recurrence is the prototype.

### A000568 Sieve Shadow

Treat Burnside cycle types as local obstruction data.  Fixed, merged,
nonfixed, q-deformed, and transporter companions are the tournament analogue of
admissible residue side channels in twin-prime and Goldbach sieves.

### Unit-Distance Ear Singular Series

For candidate ears from the exact `n=21` cores, record local direction/support
admissibility before full embedding search.  This mirrors admissible k-tuples:
first kill impossible local masks, then spend geometry on survivors.

## Assumption Challenge

Alternate vertices considered: problem domains, integers, primes, tournaments,
unit-distance points, LRC shells, Burnside partitions, prime tuples,
circle-method arcs, algebraic shadows, and proof routes.

S640 chooses proof routes as Tournament Analysis vertices.

Preserved predicate: whether a proposed quotient still certifies the domain's
hard yes/no property.

Destroyed data: representatives, embeddings, lift choices, branch labels,
local residue weights, and which transverse shadow carries the obstruction.

## See Also

`04-computation/resonance_product_ledger_s640.py`;
`05-knowledge/results/resonance_product_ledger_s640.out`;
`07-reflections/resonance-product-ledger-s640.md`;
HYP-2215; HYP-2214; HYP-2213; HYP-2212; HYP-2211; HYP-2210; HYP-2209;
HYP-2206; HYP-2200; HYP-2164; THM-407.
