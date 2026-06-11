# HYP-2442 - Erdos-Moser tower two-step recurrence is a support gate

**Status:** OPEN, with exact small support evidence.
**Source:** codex-2026-06-11-P5.
**Companions:** THM-483, THM-455, HYP-2413, HYP-2430, THM-487.
**Computation:** `04-computation/erdos_moser_support_gate_codex_p5.py`;
stored output `05-knowledge/results/erdos_moser_support_gate_codex_p5.out`.

## Statement

Let `T_{2^k-1}` be the skew-Sylvester Mersenne tower tournaments and let
`a_k = trans(T_{2^k-1})`. The known exact data now read

```text
k:   2  3  4  5   6   7   8
a_k: 2  3  5  7  11  15  23
```

matching the two-step recurrence

```text
a_{k+2} = 2 a_k + 1.
```

The next case is not "find the next datum" but a support-realization problem.
For `T511`, the recurrence predicts `a_9 = 31`, while current stored work gives
only

```text
30 <= trans(T511) <= 47.
```

The sharpened target is:

```text
Find X subset T127 with trans(X)=15 and trans(D(D(X))) >= 31,
```

or prove no such local support packet exists.

Equivalently classify packets by

```text
q(X) = trans(D(D(X))) - 2 trans(X).
```

Pure transitive chains have `q=0`; the desired `T511` recurrence witness needs
some marked support with `q>=1`.

## Evidence

The P5 support-gate atlas verifies:

```text
m   trans(T_m)  trans(D(D(T_m)))  full_bonus
3       2              5              1
7       3              7              1
15      5             11              1
31      7             15              1
```

But pure chains lose the bonus:

```text
t   trans(D(D(TT_t)))  defect
2          4              0
3          6              0
5         10              0
7         14              0
11        22              0
15        30              0
```

This exactly explains the recent `T511` witness-lift obstruction: a maximum
`TT15` chain in `T127` gives `30`, and the stored one-extra/two-extra local
packets also stay at `30`.

## Interpretation

The recurrence has two layers:

- scalar shadow: `2a+1` is the attractive two-step formula;
- support gate: the extra `+1` is not contained in a naked maximum chain.

This is the Erdos-Moser analogue of the length-72 Type II code lessons from
HYP-2430 and the order-5 fixed-projection gate HYP-2439..2441: scalar
feasibility is not support realization. The real theorem is the hidden support
object that realizes the scalar shadow.

## Next Proof Target

Do not brute-force `T511` first. Instead mine `T127` for marked packets:

- maximum chains plus incidence-word classes of outside vertices;
- border-twin packets inherited from the `T63 -> T127` step;
- orbit packets under the skew-Walsh/tower automorphism group;
- minimal packets whose double-double lift has `q=1`;
- proof of the chain-only lemma `trans(D(D(TT_t))) = 2t`.

If a packet with `q>=1` exists, it proves `trans(T511)>=31`. If all local
packets have `q=0`, the recurrence mechanism is nonlocal and the next search
should move to global zigzag automata rather than larger chain neighborhoods.

## Tournament Analysis

For this session the tournament vertices are proof routes/support packets, not
the vertices of `T_m`. The pairwise observable ranks routes by support locality,
proof leverage, computational tractability, and reuse. The resulting proof-route
tournament is transitive:

```text
marked_support_gate
  > two_step_recurrence
  > border_twin_accounting
  > chain_only_lift
  > direct_T511_sat
  > random_baseline
```

Preserved predicate: existence of a transitive witness after two skew doublings.
Destroyed data: the full tower arc algebra outside the selected support packet.
The T511 failure says that destroyed data may contain exactly the missing `+1`.
