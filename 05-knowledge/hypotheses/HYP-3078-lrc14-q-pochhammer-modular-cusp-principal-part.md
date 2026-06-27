# HYP-3078: LRC14 q-Pochhammer / Modular-Cusp Principal-Part Gate

Status: SYNTHESIS / exact q-series proof-interface scout; not a proof.
Source: codex-2026-06-27-S246.

## Claim

Every usable q-series or modular analogy in the LRC14 packet bank should be
represented as a q-cusp ledger `F_P(q)` with:

- a bounded finite principal part, meaning only finitely many negative
  q-powers;
- every polar term assigned to a named residual exit such as AP/GW boundary,
  q-witness, K33/F7, THM-572, or a covering/Haar/Ramanujan certificate;
- a product/tail representation whose nonpolar coefficients are generated,
  positive/certified, or discharged by existing packet sidecars.

Equivalently, the modular-function fact "meromorphic at the cusp has finitely
many negative q-powers" becomes an LRC quotient guardrail: a quotient may forget
an infinite positive product/partition tail only after all polar debt is finite
and named.  An infinite negative q-tail is uncontrolled residual debt, not a
legal proof object.

## Exact Scout

`04-computation/lrc14_q_pochhammer_modular_cusp_s246.py` verifies the local
series controls through finite exact arithmetic:

- `(q;q)_infty = prod_{n>=1}(1-q^n)` has Euler pentagonal support through
  `q^30`.
- `1/(q;q)_infty` produces the partition tail `p(n)` through `n=20`.
- `q d/dq log((q;q)_infty) = -sum sigma_1(n) q^n`, joining the divisor-function
  lane used in the repo's sigma/mu/phi work.
- `Delta(q)=q(q;q)_infty^24` is a cusp zero.
- `j(q)=E4(q)^3/Delta(q)` has finite principal part `q^-1`.

The deliberately illegal comparison object is `q^-1 + q^-2 + ...`: an infinite
polar tail.  That is the modular analogue of a quotient that hides infinitely
many unpaid residual obstructions.

## Packet Fields To Add

Future HYP-2963 packet rows should be able to carry:

- `q_cusp_ledger_id`
- `q_pochhammer_product_tail`
- `principal_part_order`
- `polar_exit_word`
- `polar_debt_coefficients`
- `partition_tail_certificate`
- `log_derivative_divisor_channel`
- `modular_transform_status`
- `finite_principal_part_status`
- `tail_nonnegativity_certificate`
- `illegal_infinite_polar_tail_flag`

## Proof Work

Build q-cusp ledgers for AP/GW, q=23, K33/F7, C27, covering, and route-state
closure rows.  The target theorem is not "modularity proves LRC14" but a
controlled-forgetting theorem:

> A q-series quotient is theorem-safe only if its principal part is finite and
> every polar coefficient is explained by a named packet exit.

If every HYP-2963 row admits such a ledger and all polar exits discharge through
the existing sidecar stack, the modular-cusp rule becomes a finite obstruction
gate rather than another analogy.

Incoming HYP-3075 (now the Hurwitz-Markov-Pell modular-cusp carrier) is
adjacent rather than conflicting mathematically: it treats q-series and rare
scalar coincidences as legal only after finite cusp, continued-fraction,
quadratic-unit, endpoint-shell, and carry addresses are retained.  HYP-3078 is
the narrower q-expansion audit: rare or polar q-terms are legal only after
their principal-part address and packet exit are named.

## Tournament Analysis

Vertices are proof carriers and q-expansion sidecars, not runners.  The scout
orders them by retained finite polar debt, product/tail address, divisor-channel
address, and LRC exit:

`labelled_packet_sheaf > modular_cusp_principal_part > j_single_pole_guardrail
> q_pochhammer_product_tail > delta_cusp_zero_boundary
> log_derivative_divisor_channel > partition_recursive_tail
> ramanujan_exact_period_projector > route_state_closure_median
> raw_q_series_numerology`.

The tournament is transitive in the scout with singleton SCCs and no directed
3-cycles.  Raw q-series numerology is last because it can hide illegal polar
tails.
