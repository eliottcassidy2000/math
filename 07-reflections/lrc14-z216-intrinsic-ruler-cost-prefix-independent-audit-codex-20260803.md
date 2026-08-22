# Independent audit of the two `z1=216` intrinsic-ruler cost prefixes

**Status: FINITE-EXACT AUDIT, PROJECTED-ATLAS SCOPE ONLY.**  This audit
independently confirms the two bounded cost-prefix scouts and their combined
ledger change

```text
373184 -> 373179 -> 373161,
z1=216 wall rows: 380 -> 375 -> 357.
```

It imports neither new prefix scout and trusts none of their expected row,
rank, screen, taxonomy, or ledger constants.  It reconstructs the universe
through the older pinned THM-3281 and THM-3270 atlas routes, then redoes the
selection and exact screens.  This remains only the maintained projected
`k=3` necessary atlas: it has no endpoint, owner, phase, current, arbitrary
`k<=1`, rung, physical-cover, or LRC(14) conclusion.

Companions:

- [`lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py`](../04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py);
- its [frozen output](../05-knowledge/results/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.out).

## What was rebuilt

Two separately routed old atlas parsers agree on the same `480` rows, split as
`447` wall and `33` order rows, with row digest
`53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649`.
The audit derives, rather than imports, the prior partition:

```text
19 gcd-eight = 17 low-cost + 2 costly,
 1 gcd-eighteen,
33 order,
47 rows in the three earlier natural ruler families.
```

These four sets are pairwise disjoint.  The pinned closure transcripts give
the exact lineage

```text
373284 - 1 - 17 - 33 - 47 - 2 = 373184,
447 wall - 19 gcd8 - 1 gcd18 - 47 natural = 380 wall.
```

The two costly addresses are independently recovered as atlas indices
`238,370`, rather than copied from either new prefix script.

## Family-prefix audit

On all `380` remaining wall rows, the audit groups by the complete intrinsic
fibre `(gcd(216,L),L)` and ranks fibres by `sum L*r`.  The first four ranks are

```text
(cost, rows, gcd, L, indices)
(243936, 1, 72, 11088, (210,))
(443520, 1, 24, 18480, (50,))
(517440, 3, 24,  5880, (18,140,299))
(609840, 1, 72, 27720, (134,)).
```

Thus the complete prefix through the first nonsingleton fibre has exactly the
claimed five rows.  Removing them and reranking all `375` survivors gives

```text
(609840,  1, 72, 27720, (134,))
(655200,  1, 72, 32760, (121,))
(1361808,16,72,  3528, (17,27,46,56,79,89,126,138,
                         205,248,258,286,298,347,376,425))
(1729728, 1,72, 72072, (157,)).
```

Hence the second complete prefix contains exactly eighteen rows.  Both
prefix boundaries are strict in cost, every selected fibre contains every
live atlas row with its key, and the two selected sets are disjoint.

## Screen and status-cut audit

All twenty-three old exact upper screens were rerun.  Their closure totals are

| prefix | states | crude | exact status | residual |
|---|---:|---:|---:|---:|
| first five rows | 42 | 38 | 4 | 0 |
| next eighteen rows | 750 | 479 | 271 | 0 |

For every one of the `275` status instances, the audit reconstructs
`(q,cofactor,marginals,capacities,histogram)` directly from the divisor tuple
and checks the returned rational Farkas certificate.  It then independently
tests the three solver-free event templates:

- coordinate union: the tail event is `union_(j in S){j=1}`, so its mass is
  at most `sum_(j in S)m_j`;
- zero-reduced union: a zero marginal first forces every incident binary
  pattern to mass zero, after which the same union bound applies;
- two-fan: `B union (A intersect (C union D))` has mass at most
  `m_B + min(m_A,m_C+m_D)`.

The first prefix has four coordinate-union statuses (five alternative cuts).
The second has the independently recovered taxonomy

```text
248 coordinate union,
  1 zero-reduced union,
 11 two-fan,
 11 weighted exact-Farkas only.
```

Thus `260/271` second-prefix statuses have a short solver-free monotone-event
proof, while the remaining eleven are exactly checked but are not claimed to
have minimal proof complexity.  The per-row screen and taxonomy census is
frozen in the output, so the aggregate cannot conceal a row swap.

## Boundary and reproduction

The only direction used is

```text
empty exact necessary upper screen
    => selected projected wall row is empty.
```

Family cost is a deterministic work order, not a safety invariant.  After the
audit, `357` wall rows in `33` complete families remain.  The next rank head is
the `gcd72/L72072` singleton followed by the three-row `gcd24/L25872` fibre.

Normal and optimized Python runs have identical transcripts, empty standard
error, no `assert` truth gates, and no floating-point literals:

```text
python3 04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py --processes 3
python3 -O 04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py --processes 2

semantic sha256  4bb4eca9949eb4f07b618f9e49bbb5de63492364a29a4a5abfc84a6acf105e11
script sha256    d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3
output sha256    e035531d0af15998a8a07d083042ff2486c80c7fdd729d902922179409d1725a
```
