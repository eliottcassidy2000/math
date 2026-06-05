# LRC14 Owner-Derivative No-Leak (S666)

Prompt: spend a session making progress on the LRC `n=14` proof.

S666 follows the S665/HYP-2240 rule:

```text
scalar/product collision -> address coordinate -> derivative sum
```

For LRC `n=14`, the scalar/product collision is the visible `Res_27` shadow:
AP, Vstar, and `2AP` share the same `C=27` gcd/product shell with their local
`+27` carry perturbations.  The address coordinate tested here is not the full
carry vector.  It is a small owner-deletion derivative paired to each residue.

## New Evidence

S663 had checked local carries of Hamming weight `1` and `2`.  S666 extends the
exact maximin audit to Hamming weight `3`.

For the three floor atoms AP, Vstar, and `2AP`:

- probes: `1134`
- floor atoms: `3`
- strict local carries: `1131`
- new strict rows at weight `3`: `858/858`

Minimum scores:

| Family | Weight 1 | Weight 2 | Weight 3 |
|---|---:|---:|---:|
| AP | `1/13` | `1/12` | `1/11` |
| Vstar | `2/25` | `1/12` | `1/11` |
| `2AP` | `1/13` | `1/12` | `1/11` |

So the local-carry tax now has a thicker finite base: isolated wraps do not
stay at the `1/14` floor through weight `3`.

## Owner-Deletion Repair

The checked owner derivative is deliberately tiny.  For each speed `v`, pair
its `Res_27` residue with:

```text
(number of D/U/N obligations covered by v,
 private-owner flag)
```

The private-owner flag is the deletion derivative: it is true exactly when
deleting `v` would uncover at least one D/U/N obligation.

Projection audit:

- visible shadow alone: `3` mixed buckets, each `1` floor + `377` strict;
- visible + cheap pair: still `3` mixed buckets;
- visible + paired owner cover count: only `2` singleton leaks remain;
- visible + paired private-owner flag: `0` mixed buckets.

The two leaks before the private-owner bit are exactly:

- `AP:w1:carry(11)`;
- `Vstar:w1:carry(11)`.

That is the best kind of failure: small, named, and repaired by a single
interpretable derivative bit.

## Interpretation

Cheap-pair existence is too coarse.  Every bounded row in this local atlas has
an unblocked small pair.  What matters is the paired owner ledger: which
residue owns which part of the D/U/N certificate, and whether deleting that
residue breaks the cover.

The result does not prove LRC `n=14`.  It does give a cleaner local lemma
target:

> In the `Res_27` fiber of AP/Vstar/2AP, any nonzero local carry that preserves
> the visible shell must change the paired owner-deletion ledger, unless it is
> part of a globally coherent scalar floor lift.

That wording is important.  S611 showed scalar unit lifts can have nonzero
global carries and still remain floor by scaling invariance.  The emerging
distinction is not "zero carry versus nonzero carry"; it is local incoherent
carry versus globally coherent carry.

## Tournament Analysis

Vertices are repair channels.  The observable is
`(mixed fibers, compression, owner alignment, pairedness, carry independence)`.
The repair-channel tournament is transitive:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1}`;
- `directed_3cycles=0`;
- `hamiltonian_paths=1`;
- top channel: `visible+owner_private_flag`.

This is a good sign for proof work.  Unlike the broad S665 atlas, this local
question has a clear best finite repair channel.

## Next Move

Extend from local Hamming-weight carries to coherent carry spaces:

1. scalar unit lifts, which are globally coherent floor rows;
2. two-block carry patterns, which should expose the first nonlocal obstruction
   family;
3. HYP-2165 owner-route lifts of the 64 fixed classes.

The proof should aim to show that any row invisible to the `Res_27` shell but
not globally coherent must pay either an owner-deletion change or a positive
loneliness tax.
