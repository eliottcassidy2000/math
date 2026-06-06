---
id: HYP-2252
status: OPEN proof-partition evidence for LRC14
source: codex-2026-06-05-S676
related:
  - HYP-2241
  - HYP-2164
  - HYP-2167
  - HYP-2155
  - THM-406
---

# HYP-2252: LRC14 Lebesgue Wall Certificate

## Claim

For total runner count `n=14`, Lebesgue safe measure should be used as a
proof filter:

```text
p_0(V, 1/14) > 0  => interval certificate, easy
p_0(V, 1/14) = 0  => endpoint wall, needs a set-level certificate
```

After quotienting by global scalar symmetry and the `Res_27` floor ledger,
the only checked endpoint walls are AP, Vstar, and nonprimitive `2AP`.  Thus
the LRC14 proof target sharpens to:

> Prove there are no new lifted `p_0=0` walls in the `Res_27` carry/owner
> fiber beyond scalar floor orbits of AP/Vstar/`2AP`; then verify the named
> wall atoms by their explicit endpoint witnesses.

This is a strengthening of the HYP-2241 no-leak route: positive maximin tax is
upgraded to positive Lebesgue safe measure.

## Exact Evidence

S676 adds `04-computation/lrc14_lebesgue_wall_s676.py` with stored output in
`05-knowledge/results/lrc14_lebesgue_wall_s676.out`.

The script sweeps the exact rational danger-arc subdivision.  It does not use
grid sampling.

At `delta=1/14`:

- AP has `p_0=0`, no positive safe components, and `6` endpoint witnesses
  `1/14,3/14,5/14,9/14,11/14,13/14`.
- Vstar has the same six endpoint witnesses.
- `2AP` has `p_0=0` and `12` endpoint witnesses, the two preimages of the AP
  witness set under doubling.

The three wall atoms also have the same one-sided opening:

```text
p_0(1/14 - eps) / eps = 23324/6435
```

for `eps = 1/700, 1/1400, 1/2800, 1/5600` in the exact audit.

## Carry And Swap Evidence

Local `+27` carry perturbations through Hamming weight `3` over AP, Vstar, and
`2AP`:

- probes: `1134`;
- `p_0=0` walls: exactly the three unperturbed floor atoms;
- positive safe measure: `1131`.

Minimum positive `p_0` values in the local carry atlas:

| Family | Weight 1 | Weight 2 | Weight 3 |
|---|---:|---:|---:|
| AP | `7/858` | `4193/190476` | `55875/2074072` |
| Vstar | `426/35035` | `4045/168168` | `553771/18666648` |
| `2AP` | `7/858` | `530221/23711688` | `29590369/735062328` |

One-swap AP rows with new speed `<=60`:

- rows checked: `611`;
- zero-measure rows: exactly one, `12->24`, namely Vstar;
- all other one-swaps have positive safe measure.

The smallest positive one-swap safe measure found is `1/1260` at `12->36`.

Global scalar floor orbits explain why the finite wall statement must be
quotiented: all `36` AP/Vstar unit scalar rows modulo `27` still have `p_0=0`
by measure-preserving scaling.

## Interpretation

Lebesgue measure alone cannot prove LRC14, because the hard wall has
`p_0=0`.  But it can cleanly classify cases:

1. Positive measure rows have an open interval of lonely times and require no
   endpoint finesse.
2. Endpoint-wall rows require exact boundary witnesses.
3. Counterexamples would have to be point-saturated, not merely
   measure-saturated.

So the remaining theorem should not be phrased as "find positive measure for
everything."  It should be:

```text
every lifted row is either positive-measure strict
or a scalar/quotient copy of AP, Vstar, or 2AP with explicit endpoint witnesses.
```

## Tournament Analysis

Vertices are proof routes, not runners.  Candidate vertex sets considered
included runners, danger arcs, breakpoints, safe components, endpoint
witnesses, `Res_27` rows, carry vectors, owner obligations, and proof
obligations.

The route tournament uses the observable
`(exactness, endpoint handling, lift relevance, actionability, compression, safety)`.
It is transitive:

- `score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1}`;
- `directed_3cycles=0`;
- `scc_sizes=[1,1,1,1,1,1,1]`;
- `hamiltonian_paths=1`;
- order:
  `endpoint_wall_certificate > local_carry_measure_tax >
  res27_floor_quotient > global_scalar_floor_orbit > one_swap_wall_sieve >
  positive_measure_interval > raw_first_moment`.

## Next Lemma Target

Turn HYP-2241 into a measure-wall theorem:

> In each normalized `Res_27` carry/owner fiber, any row not belonging to a
> scalar AP/Vstar/`2AP` floor orbit has `p_0(V,1/14)>0`.

Then the endpoint certificate on AP/Vstar/`2AP` completes that fiber.
