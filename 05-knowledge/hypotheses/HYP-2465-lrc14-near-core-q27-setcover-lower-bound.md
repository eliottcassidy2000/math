# HYP-2465 - LRC14 near-core rows cannot block Q27 in the carry window

**Status:** OPEN synthesis with exact bounded set-cover certificate.

**Source:** codex-2026-06-13.  Strengthens HYP-2463 and sharpens OPEN-Q-082 into a nine-core compression target.

## Claim

Inside the HYP-2444 carry window `1 <= v <= 13*84 = 1092`, a primitive LRC14 replacement row that retains at least nine of the twelve speeds in

```text
CORE = 7*{1,...,12}
```

always has a Q27 witness.

More explicitly, let `D subset CORE` with `|D|=e <= 3`, and let `A` be any set of `e+1` added speeds from `1..1092` excluding original core speeds.  If

```text
S = (CORE \ D) union A
```

has 13 distinct speeds and is primitive, then some `q in Q27={d*m: d in {1,2,7,14}, 1<=m<=27}\{1}` has a strict witness.

Equivalently: a would-be primitive Q27 blocker in this carry window must either leave the near-core regime by deleting at least four core speeds, or escape the bounded carry window / known compression frame.

## Method

For a partial row `B=CORE\D`, call a twist `(q,a)` an obligation if it is safe for `B`.  An added speed `v` covers the obligation if `v` is dangerous at `(q,a)`.  Blocking Q27 is exactly a set-cover problem over these obligations.

The script adds two constraints:

- cardinality: choose at most `|D|+1` added speeds;
- primitivity: choose at least one added speed not divisible by `7`, since every core speed is divisible by `7`.

The resulting binary MILP asks whether a primitive Q27 cover exists.  Infeasibility is the desired witness that all such rows have a Q27 opening.

## Exact Evidence

Script: `04-computation/lrc14_near_core_q27_setcover_codex.py`

Stored output: `05-knowledge/results/lrc14_near_core_q27_setcover_codex.out`

The candidate universe has `1080` speeds: all integers `1..1092` except the original 12 core speeds.

Q27 set-cover infeasibility table:

```text
delete_count e=0: cases=1,   add_budget=1, infeasible=1,   feasible=0, unknown=0
delete_count e=1: cases=12,  add_budget=2, infeasible=12,  feasible=0, unknown=0
delete_count e=2: cases=66,  add_budget=3, infeasible=66,  feasible=0, unknown=0
delete_count e=3: cases=220, add_budget=4, infeasible=220, feasible=0, unknown=0
```

The obligation counts grow rather than collapse:

```text
e=0: 126..126 Q27 obligations
e=1: 182..462
e=2: 252..728
e=3: 364..1008
```

This is why the near-core blocker cannot exist: deleting core speeds creates more Q27 obligations than `e+1` primitive additions can jointly cover.

## Plain Shell Noise

The plain shell `q<=27` is much noisier.  A direct scan over one deletion and two arbitrary additions finds `877` primitive rows that still block all plain `q<=27` witnesses.

But every one of those `877` rows has a Q27 witness.

So the proof should not aim to classify all plain-shell blocker pairs.  The correct object is the fibered Q27 set-cover ledger.  Plain shell blockers are numerous scalar shadows; Q27 infeasibility is the stable support statement.

## Proof Target

HYP-2463 said the eight named hard residues do not stack.  HYP-2465 says something stronger and less residue-specific:

> In the bounded carry window, retaining nine core speeds is already enough to force Q27, regardless of the added residues.

The remaining LRC14 compression theorem should therefore try to prove:

1. A primitive Q27 blocker can be normalized into the HYP-2444 carry window, or else it opens a divisor/carry/Bprime channel.
2. If it stays in the carry window and retains at least nine core speeds, HYP-2465 gives a Q27 witness.
3. If it deletes four or more core speeds, then the row is no longer a small perturbation of the 7-core and must pay a different descent tax: low clocks, AP/Vstar/2AP descent, owner-private deletion, or a support-load contradiction.

## Assumption Challenge / Tournament Analysis

Candidate tournament vertices considered:

- runners,
- added speeds,
- deleted core speeds,
- individual twists `(q,a)`,
- Q27 obligations,
- candidate-cover signatures,
- proof obligations.

The selected tournament uses proof obligations/set-cover gauges as vertices.  This preserves which proof route is load-bearing and discards raw runner geometry.  Stored Tournament Analysis is transitive with leader `primitive_Q27_set_cover`, then `near_core_compression`, then `plain_shell_cover_noise`.

The challenged assumption is that the eight hard residues are the only important obstruction atoms.  They are not: plain-shell blockers can be built from many non-hard residues.  The real invariant is not the residue name; it is the inability of a small primitive set-cover to cover all Q27 obligations.

