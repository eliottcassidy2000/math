---
id: HYP-2625
title: LRC(14) modular recurrence address - mod 6 is the universal skeleton, mod 30 is the first non-universal address, and mod 210 is the divisor-profile interface to the support-6 mod-7 tail
status: OPEN
source: codex-2026-06-19-S18
depends_on:
  - HYP-2598
  - THM-538
  - HYP-2617
  - HYP-2619
related:
  - HYP-2627
  - HYP-2621
  - HYP-2622
  - HYP-2624
  - OPEN-Q-108
---

# HYP-2625 - LRC(14) Modular Recurrence Address

## Claim

The useful modular arithmetic for the LRC(14) proof should be organized as a
squarefree divisor-profile recurrence, not as a one-dimensional runner
recurrence.

The hierarchy is:

```text
mod 6   = denominator-2/3 universal-center skeleton (proved in HYP-2598)
mod 30  = add denominator 5 as the first non-universal recurrence address
mod 210 = add the support-six mod-7 tail address from THM-538/HYP-2617
```

Only the first row is a standalone universal-center proof tool.  The `5` and
`7` layers are exact modular addresses for recurring intervals and coimage
tails, but their center gaps are too small to certify the large-spread cluster
by themselves.  The proof target is therefore a coupled object:

```text
prime-mask inclusion-exclusion over {2,3,5,7}
  + signed projective mod-7 coimage tail
  + finite low-height wall ledger.
```

## Computation

Script:

- `04-computation/lrc14_modular_recurrence_probe_codex_s18.py`
- output: `05-knowledge/results/lrc14_modular_recurrence_probe_codex_s18.out`

The script enumerates `P subset {1,...,13}` exactly and checks rational centers
`a/b` by the strict condition

```text
||v a/b|| > 1/14  for every v in P.
```

It also computes the prime-mask inclusion-exclusion recurrence for center
primes `{2,3,5,7}`.

## Exact Sequence Rows

The known HYP-2598 survivor sequence is recovered as the `{2,3}` row:

```text
1, 11, 47, 109, 156, 146, 91, 37, 9, 1, 0, 0, 0, 0
```

with formula

```text
C(7,s) + C(9,s) - C(5,s).
```

Adding denominator `5` gives the mod-30 address row:

```text
1, 13, 72, 223, 437, 581, 545, 366, 174, 56, 11, 1, 0, 0
```

with inclusion-exclusion terms:

```text
C(7,s)+C(9,s)+C(11,s)-C(5,s)-C(6,s)-C(7,s)+C(4,s),
```

where the `C(11,s)` term is "avoid multiples of 5".

Adding denominator `7` gives the mod-210 address row:

```text
1, 13, 78, 280, 658, 1066, 1231, 1030, 623, 266, 76, 13, 1, 0.
```

The complementary "hits every prime in `{2,3,5,7}`" residual is:

```text
0, 0, 0, 6, 57, 221, 485, 686, 664, 449, 210, 65, 12, 1.
```

These are finite exact sequences, not asymptotic guesses.

## Denominator Guardrail

The exact center-denominator table through `b<=30` shows:

- `b=2,3` are the only universal cluster centers, because `1/b > 2/7`;
- `b=5,7,10,14,15,21,30` supply many safe center addresses for small parts, but
  none have a wide enough center gap to solve the large-spread cluster alone;
- treating those larger denominators as universal witnesses is the wrong
  overclaim.

This explains why mod `30` looked structurally real in HYP-2621/HYP-2622 but
was not enough.  It is a recurrence address, not a proof by itself.

HYP-2624 adds a useful downstream constraint after this recurrence address was
formed: height-2 one-large coimage walls already address every nonzero class for
`k=8,9`, and they reduce the `k=10` tail to repeated-residue packets.  Thus the
mod-210 row should be read as the shared address space for HYP-2624's wall
ledger and the remaining signed reciprocal tail, not as a separate proof.

HYP-2627 adds an external four-factor check on the same address: the raw
Harary-Hill product for `K_14` is `7*6*6*5=1260`, with squarefree core
`210={2,3,5,7}`.  The divided crossing value `315` loses the dyadic coordinate.
Moreover THM-523's half-clash satisfies
`15/36-2/5-1/70-1/504=1/(2*1260)`, so the known `1/1260`
lonely-measure scale is exactly the reciprocal raw product after symmetry
doubling.  This supports the rule used here: keep the raw denominator/profile
ledger until the squarefree recurrence and coimage projection have both been
recorded.

HYP-2628 refines this address by replacing the bare squarefree mask with
exact-period `phi` packets.  HYP-2629 adds the Hill-row scan: adding a
squarefree prime `p` appends a shifted copy layer multiplied by `p-1`, so
`{2,3}->{2,3,5}->{2,3,5,7}` is an honest copy recurrence.  At `K_14`,
`P_14=1260` has full `{2,3,5,7}` copy mass `576`, while the divided crossing
value `315` has full copy mass `0`; quotienting by four loses the dyadic gate
before the mod-210/coimage address has been recorded.

## Proof Route

Use the prime-mask recurrence as the finite transfer state for small parts:

```text
mask(v) = {p in {2,3,5,7}: p divides v},  v in {1,...,13}.
```

The mask histogram is:

```text
0000: 3, 0001: 3, 0010: 2, 0011: 2, 0100: 1, 0101: 1, 1000: 1.
```

The proposed proof split is:

1. Keep the HYP-2598 mod-6 cases as the universal-center skeleton.
2. Use the mod-30 row as the first recurrence address for non-universal
   recurring intervals.
3. Couple the mod-210 row to the THM-538 support-6 floor and HYP-2617 projective
   coimage fibers.
4. Delete finite low-height wall rows as in HYP-2616/HYP-2619, then prove the
   signed mod-7 reciprocal tail on the remaining non-null coimage classes.

## Tournament Analysis

The computation uses proof obligations as vertices, not runners.  The route
tournament is transitive, with Hamiltonian path:

```text
support6_coimage_tail
> mod6_universal
> mod210_address
> mod30_address
> primorial_height_escape
> all_small_denominators
> raw_runner_residues
```

Pairwise observable: preservation of the LRC predicate, proof legitimacy,
residual compression, and recurrence/address value.  This quotient preserves
the modular proof predicate while discarding individual witness-time geometry.

## Status

Open.  HYP-2625 does not prove LRC(14).  It identifies the missing modular
recurrence layer behind the existing mod-6 and mod-30 clues, and it gives a
specific next proof object: a squarefree divisor-profile recurrence coupled to
the signed support-6 mod-7 coimage tail.
