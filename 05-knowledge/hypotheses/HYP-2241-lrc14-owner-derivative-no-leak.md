---
id: HYP-2241
status: OPEN finite evidence for LRC14 owner-derivative no-leak
source: codex-2026-06-05-S666
related:
  - HYP-2240
  - HYP-2239
  - HYP-2237
  - HYP-2231
  - HYP-2230
  - HYP-2222
  - HYP-2167
  - HYP-2165
  - HYP-2164
---

# HYP-2241: LRC14 Owner-Derivative No-Leak

## Claim

HYP-2240's address-coordinate derivative repair grammar becomes concrete over
the LRC `n=14` `Res_27` carry seam:

```text
fixed odd-wall packets
+ fixed C=27 gcd/product shell
+ paired carry/owner/deletion derivatives
=> AP, Vstar, 2AP, or strict looseness.
```

The anticipated obstruction is not a new scalar invariant.  It should be a
no-leak statement: if a perturbation of the AP/Vstar/2AP floor atoms keeps the
same visible `Res_27` and product/gcd shell data, then either the paired
owner-derivative ledger reconstructs the original carry cocycle or the row pays
positive loneliness tax.

## S666 Evidence

S666 adds `04-computation/lrc14_owner_derivative_repair_s666.py` with stored
output in
`05-knowledge/results/lrc14_owner_derivative_repair_s666.out`;
and reflection
`07-reflections/lrc14-owner-derivative-no-leak-s666.md`.

The script builds the local carry atlas around the three known floor atoms
`AP`, `Vstar`, and `2AP`.  For each base row, it exhaustively adds `+27` to
one, two, or three coordinates.  These perturbations keep the same visible
`Res_27` residue shadow and `C=27` gcd/product shell as their base row, so the
visible quotient has three mixed buckets.

Exact maximin evidence:

- `1134` probes total: `3` floor atoms and `1131` strict local carries.
- New S666 radius: every one of the `858` weight-3 local carries is strict.
- AP minima: weight `1,2,3` give `1/13`, `1/12`, `1/11`, achieved by moving
  `(13)`, `(12,13)`, `(11,12,13)`.
- Vstar minima: weight `1,2,3` give `2/25`, `1/12`, `1/11`, achieved by moving
  `(13)`, `(13,24)`, `(11,13,24)`.
- `2AP` minima: weight `1,2,3` give `1/13`, `1/12`, `1/11`, achieved by moving
  `(26)`, `(24,26)`, `(22,24,26)`.

The owner-deletion side channel is intentionally small.  For each speed `v`,
pair its `Res_27` residue with:

```text
(number of D/U/N obligations covered by v,
 whether v is the unique owner of at least one D/U/N obligation).
```

The second coordinate is a deletion derivative: deleting `v` uncovers an
obligation exactly when this private-owner bit is true.

Projection audit:

| Repair key | Groups | Mixed floor/strict groups | Max bucket |
|---|---:|---:|---:|
| visible shadow + gcd shell | `3` | `3` | `378` |
| visible + cheap pair | `112` | `3` | `106` |
| visible + paired owner cover count | `976` | `2` | `2` |
| visible + paired private-owner flag | `1067` | `0` | `2` |
| visible + paired private-owner count | `1134` | `0` | `1` |
| visible + carry L1/support | `1134` | `0` | `1` |

The key new leak diagnosis is precise: paired owner cover counts alone fail
only on `AP:w1:carry(11)` and `Vstar:w1:carry(11)`.  Adding just the
private-owner deletion bit repairs those two remaining leaks and stays
carry-free.

## Tournament Analysis

Vertices are repair channels, not runners.  The pairwise observable is
`(mixed fibers, compression, owner alignment, pairedness, carry independence)`;
the switch is majority, with channel-list order as the tie Hamiltonian path.

Fingerprints:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1}`
- `directed_3cycles=0`
- `scc_sizes=[1,1,1,1,1,1]`
- `hamiltonian_paths=1`

The top channel is `visible+owner_private_flag`: it is the smallest checked
carry-free owner-deletion repair that separates every local carry through
weight `3`.

## Assumption Challenge

S666 did not assume the right vertices are runners.  Candidate vertices included
runners, residues, gaps, fixed circle sections, section boundaries,
wall-crossing events, pair-sum denominators, owner obligations, deleted-speed
cards, carry coordinates, and proof obligations.  The selected quotient
preserves the finite predicate "does this side channel separate known floor
atoms from bounded local carry perturbations?" and deliberately destroys phase
order and raw speed identity.

## Next Lemma Target

Prove a local no-leak lemma:

> In the `Res_27` fiber of the AP/Vstar/2AP floor atoms, any nonzero local
> carry that preserves the visible shell must change the paired
> owner-deletion ledger, unless it belongs to a globally coherent scalar floor
> lift.

The immediate next computation is to extend the same audit from local
Hamming-weight carries to coherent carry subspaces: scalar unit lifts,
two-block carry patterns, and HYP-2165 owner-route lifts of the 64 fixed
classes.
