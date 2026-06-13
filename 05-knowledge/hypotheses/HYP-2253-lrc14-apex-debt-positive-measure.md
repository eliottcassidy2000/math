---
id: HYP-2253
status: OPEN sharpened LRC14 proof target
source: codex-2026-06-06-S677
related:
  - HYP-2252
  - HYP-2241
  - HYP-2168
  - HYP-2167
  - HYP-2164
  - THM-406
---

# HYP-2253: LRC14 Apex-Debt Positive Measure

## Claim

HYP-2252's target

```text
no new normalized p_0=0 wall in the Res_27 carry/owner fiber
```

should be split into three cases:

1. **No-multiple rows:** harmless, because the six `n`-clock witnesses
   `j/14`, `j in {1,3,5,9,11,13}`, survive.
2. **GCD-scaled walls:** nonprimitive scalar reductions such as `2AP`, handled
   by scaling and gcd-scaled endpoint witnesses.
3. **Primitive apex debt:** rows with `gcd(V)=1` and some speed divisible by
   `14`.  This is the real residual, and it should have `p_0(V,1/14)>0`.

In carry coordinates `v=r+27k`, the apex-debt congruence is exact:

```text
14 | v  <=>  k == r (mod 14)
```

because `27 == -1 (mod 14)`.  Thus the improved theorem target is:

> Every normalized primitive `Res_27` carry/owner row with apex debt has
> positive Lebesgue safe measure.

Equivalently: the only `p_0=0` walls are no-multiple endpoint walls or
nonprimitive scalar endpoint walls.

## S677 Evidence

S677 adds `04-computation/lrc14_apex_debt_lebesgue_s677.py` with stored output
in `05-knowledge/results/lrc14_apex_debt_lebesgue_s677.out`.

It computes exact rational `p_0` for `720` coherent carry probes over AP and
Vstar:

- `single_apex`: carry one coordinate to its first multiple of `14`;
- `apex_plus_one`: one apex-debt coordinate plus one ordinary `+27` side carry;
- `interval_block`: contiguous carry blocks of height `1` or `2`;
- `affine_mod14`: small affine laws `k=a*r+b (mod 14)`.

Primitive multiple branch:

- primitive multiple probes: `414`;
- primitive multiple walls: `0`;
- primitive multiple positive-measure rows: `414`;
- minimum primitive-multiple safe measure:

```text
p_0 = 181/28028
```

achieved by AP with speed `6` carried to `168`:

```text
(1,2,3,4,5,7,8,9,10,11,12,13,168).
```

This row has `gcd=1`, no `n`-clock witnesses, no gcd-scaled endpoint witnesses,
and `4` positive safe components.

## Endpoint Split

The endpoint buckets in the S677 audit make the proof split visible:

- `(multiple=False, gcd=1, n_clock=6, route=positive)`: no-multiple rows keep
  the ordinary six `n`-clock witnesses.
- `(multiple=True, gcd=1, n_clock=0, route=positive)`: primitive apex-debt
  rows lose the ordinary clock but open positive measure.
- `(multiple=True, gcd=28, gcd_clock=168, route=wall)`: a nonprimitive scalar
  wall, not a primitive counterexample lane.

So HYP-2252's "no new wall" should not spend energy on AP/Vstar no-multiple
walls.  It should attack primitive apex debt.

## Tournament Analysis

Vertices are proof filters, not runners.  Candidate vertex sets considered
included runners, `n`-clock endpoints, gcd-scaled endpoints, apex congruence
sites, carry vectors, safe intervals, owner obligations, and proof
obligations.

The filter tournament is transitive:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`;
- `directed_3cycles=0`;
- `scc_sizes=[1,1,1,1,1,1,1]`;
- `hamiltonian_paths=1`;
- order:
  `apex_congruence_debt > primitive_multiple_positive_measure >
  owner_private_derivative > gcd_scaled_endpoint_wall >
  no_multiple_n_clock > raw_res27_shadow > raw_first_moment`.

## Next Lemma Target

Prove an apex-debt measure lemma:

> Let `V` be a normalized primitive LRC14 row in the `Res_27` carry/owner
> fiber.  If some carry site satisfies `k_i == r_i (mod 14)`, then
> `p_0(V,1/14)>0`, unless the row is a nonprimitive scalar wall after gcd
> reduction.

The likely proof channel is HYP-2241's owner-private deletion ledger, now
localized at apex-debt sites.  The owner route should show that creating a
multiple of `14` removes the old endpoint witness set but necessarily opens a
positive safe interval elsewhere.
