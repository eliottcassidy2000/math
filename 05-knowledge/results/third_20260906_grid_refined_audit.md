# Independent audit of the refined full-word grid consumer

**ACCEPTED consumer reduction and independent exact census, relative to the
pinned full-word maximum table.** The table's upper-bound proof and complete
independent optimization audit are separate dependencies. This audit accepts
[the refined consumer](third_20260906_grid_refined.md): its three complete
sets are reproduced independently, leaving **8,202 scales**, maximum
**16,704**, in the final relaxation. No mathematical repair to that
consumer was required. LRC(14) remains OPEN.

## 1. Dependencies and independence

The starting domain is the complete 8,301-scale set from the
[independently audited baseline](third_20260906_grid_bootstrap_audit.md).
The audit reads its own independently generated array, with compact-JSON
SHA256

```text
a25d83f0eeb630bb82e84cdfac4e3cf7312f892f6c426d6affd5239a064e4b58.
```

The full-word maximum supplier is
[the separate optimization theorem](third_20260906_grid_full_words.md),
with [its independent optimization audit](third_20260906_grid_full_words_audit.md).
Its exact ordered `[t,M(t)]` table is pinned by

```text
ca6b6f562db1fc3632f8b7570b89a16020a981ae8aa130be200dc1bdcb4264ca.
```

Every table clock must match the independently reconstructed baseline,
with no omission, duplicate, or reordered clock. The maximum semantics
exclude diagnostic owner order, prefix-visit counts, and the supplier's
additional global rows. The supplier's successful replay marker is also
required. The audit checks all supplied owners' divisibility, primitivity,
and literal ceiling cost; these checks do not substitute for the
separate proof that their scores are upper bounds.

Neither producer is imported. The refined auditor explicitly reuses
its own pinned
[baseline audit module](../../04-computation/third_20260906_grid_bootstrap_audit.py),
whose source SHA256 is

```text
ec841b2a3e9898dcf675d4c0212b78e3d16e21b8d13a05aa54b535a231226cc0.
```

The reused helpers regenerate the strict atlas by a prime sieve,
recover all 5,855 component-length multisets by literal integer interval
intersection, and supply the bounded-knapsack routine. A new C++ main
performs the refined enumeration. The old audit's main is not executed.
The producer output is consulted only after the independent arrays have
been constructed, for exact comparison.

## 2. The analytic consumer is sound

For each baseline clock t, the full-word supplier maximizes the exact
seven-state excess over all allowed divisors of t satisfying every
projected complement word. An actual profile-surviving row has its
seven sheet states in this finite set. Therefore its true excess obeys
`E_actual<=M(t)`.

The inherited pair theorem supplies some actual edge with `e<=6`.
For its primitive ratio `(p,q)`, the two sheet states are exactly

```text
d_p=e*gcd(t/e,p), d_q=e*gcd(t/e,q).
```

The baseline's conditioned capacity maximum `E_pair` is another upper
bound for the same actual excess, provided these forced values are
present and their two uses fit the capacities. Consequently

```text
E_actual <= min(M(t),E_pair).
```

This minimum is valid without a common maximizing word. It is generally
a relaxation of the joint profile optimization conditioned on the pair;
the consumer does not assert that it equals that joint optimum.

The exact individual-component credit C is a lower bound on actual
pair overlap on every translated grid. Therefore `C>min(M,E_pair)`
forces a weak-safe lift by the inherited grid inequality. The strict
inequality is necessary; equality belongs to the survivor predicate.
Every divisor `e|t`, `1<=e<=6`, is retained, and all actual atlas ratios
remain available. A discarded clock has no possible candidate for the
actual edge supplied by the theorem. A surviving clock only supplies
an abstract compatible candidate for the declared inequalities.

The earlier weak-safe six-component phase input, clock bound, and
actual-edge theorem retain their original hypotheses and statuses.
No grid endpoint is promoted from weak safety to strict clearance.

## 3. Independent finite enumeration

For each of the 8,301 input clocks, the independent engine:

1. Reconstructs the exact allowed divisor capacities and ceiling weights.
2. Solves a seven-slot bounded knapsack to verify `0<=M(t)<=E_bag(t)`.
3. Checks every admissible e. Its aggregate prefilter minimizes the
   exact best measure separately at each component count J, using the
   literal atlas; it imports no two-line envelope implementation.
4. Retains all ratios needed after that valid prefilter and computes
   the sum of the individual component ceilings.
5. Reserves both forced pair uses and solves a five-slot bounded
   knapsack, then compares C with `min(M(t),E_pair)`.

The conditioned capacity bound is compared with the original bag
maximum, not incorrectly required to be at most M(t). Only their
minimum bounds the refined test. All arithmetic is integral, including
every ceiling and exact division of seven weight sums.

A positive existential witness permits an early exit. At a rejected
clock all candidate divisors and ratios are excluded explicitly, or by
an aggregate bound which the individual-component credit dominates.
Thus the computation reconstructs the entire declared finite sets:

| Stage | Count | Maximum |
|---|---:|---:|
| Full-word cost with aggregate credit | 8,288 | 21,600 |
| Full-word cost with individual components | 8,202 | 16,704 |
| Individual components with both upper cost bounds | 8,202 | 16,704 |

All three complete arrays agree exactly with the refined producer.
The final two arrays are identical. Exactly 99 baseline clocks are
removed, and the only survivor greater than 14,904 is 16,704.
The final compact-JSON array hash is

```text
4f29481d984ead40d0144556ce1c45dce210e30b964bb65835a7904ca6353e59.
```

## 4. The retained upper boundary survives genuine pair conditioning

The independently checked boundary is

```text
t=16704, e=4, (p,q)=(3,308),
d=(12,16,72,58,64,9,9).
```

The ratio belongs to the exact strict atlas. Each of the seven states
divides t, their gcd is one, and all 126 nonempty proper-subset
complement words pass the pinned full profile table. Their literal
ceiling excess is 188. The pair formula forces precisely `(12,16)`,
which appears in the word. The literal interval engine gives component
credit **172**.

The supplied global maximum is `M(16704)=188`. This word attains it
while containing the forced pair states. The full-profile maximum
conditioned on this pair is therefore also exactly 188: it is at least
the owner's score and at most the global maximum. No separate search
or unattained upper estimate is needed for that equality.

Thus the pair remains in the declared conditioned relaxation with
`172<188`. Rejecting an earlier pair with different forced states
cannot justify removing the entire clock. This is a stopping witness
for the proposed finite consumer improvement, not a physical decoder
realization, an unsafe row, or a proof that the actual scale ceiling
cannot be improved by additional information.

## 5. Reproduction

The [independent refined source](../../04-computation/third_20260906_grid_refined_audit.py)
and [full retained output](third_20260906_grid_refined_audit.out) contain
all three arrays, every removed clock, and the boundary checks. The
script compiles its exact C++17 engine in a temporary directory.

```bash
python3 -B 04-computation/third_20260906_grid_refined_audit.py
python3 -B -O 04-computation/third_20260906_grid_refined_audit.py
```

Normal and optimized outputs match byte for byte. The run passes
3,642,410 C++ checks, 16,748 refined Python checks, and 17,567 checks
in the reused independent geometry auditor. Raw LF SHA256 values:

```text
source 023e1ec58595b268e28f0049ad748e9b9f392965fe9534ca7dc9b547d8c9fb7f
output 1a3a537a8e31ccc2a6e0741646d50d6a2b49c0240b73c4d71e9df8af169d97d1
semantic 449bc42ba2e701ebd6ea8fbf6fd3fd90a0981000c928c957e516993d1b01755d
```

This audit verifies the consumer of the pinned maxima. The mathematical
upper-bound status of that table comes from its separate complete
optimization proof and audit, not from agreement of downstream counts
or the existence of maximizing owners.

## Additional scope audit: every qualifying connected complement

A separate analytic review of the global
[joint-shadow supplier](lrc14_joint_shadow_empty_core_next_sep06.md) and its
[recursive source](lrc14_recursive_gcd_empty_core_next_sep06.md) accepts the
explicit connected-complement corollary in the final producer note. Their
necessary profiles apply to every primitive thirteen-speed strict failure,
including arbitrary proper subsets and unbounded speed height. No decoder
rank, equality relation, or physical box enters that supplier.

For any chosen six/seven partition, divide each part by its own gcd t and g.
Global primitivity gives gcd(t,g)=1 and both normalized shapes are primitive.
The global profiles give g<=90 and every projected word used by the sheet
state compiler. Connectedness of the seven-label strict atlas supplies the
actual edge whose sheet gcd is at most6. The six-label safe phase needs no
inert connectivity on that side. The remaining grid and finite-cost arguments
are unchanged, and use no equality or height assumption. Hence every
qualifying partition of a strict failure has t in the same8,202-element set.
This audit does not extend the bound to partitions without a connected
seven-label strict atlas. No source, output or finite certificate changes.
