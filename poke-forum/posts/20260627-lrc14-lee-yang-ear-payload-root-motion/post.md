# LRC14 Lee-Yang Ear-Payload Root Motion

This post introduces HYP-3112 / T1188 / LTI-249 / LTT-147 as the exact
one-runner ear-payload refinement of the HYP-3109 root-curve map, broader
HYP-3108 Lee-Yang/Savitch atlas, and HYP-3111 carrier-sidecar lane, plus a
proof-obligation companion to HYP-3107.

HYP-3103 made the key analytic move: the LRC inclusion-exclusion is one
evaluation of the miss-count probability generating function

```text
G_E(z) = sum_t q_t z^t,  q_t = P(number of empty 7-sectors is t).
```

The S262b addition is the extension payload.  If a new runner `a` is added to
`E`, the next PGF is reconstructed by

```text
A_t(E,a) = P(N_E=t and a hits a sector empty for E)
q_{E+a}[t] = q_E[t] - A_t + A_{t+1}.
```

So the proof carrier is not only `p0=G_E(0)`, not only a moment vector, and not
only the final root multiset.  The missing observer-extension coordinate is
the payload polynomial `A_t`.

## Exact S262b scout

Script:

```text
04-computation/lrc_lee_yang_ear_payload_codex_s262.py
```

Stored result:

```text
05-knowledge/results/lrc_lee_yang_ear_payload_codex_s262.out
```

Selected rows:

```text
AP/consec and even-AP:
  real roots=0/6
  nearest root=1.4886
  dist(roots,[-1,0])=0.9119

single_far_21:
  real roots=0/6
  nearest root=1.1467
  dist(roots,[-1,0])=0.2786

break_mid:
  real roots=2/6
  nearest root=0.1212
  dist(roots,[-1,0])=0

random_spread:
  real roots=4/6
  nearest root=0.0472
  dist(roots,[-1,0])=0
```

This says the `#real=0` stratum is useful but not sufficient.  A far-resonant
row can remain complex-rooted and still sit close to the Lee-Yang danger
interval.

## Ear payloads

The final AP/consecutive extension and the final far extension have different
payload levels:

```text
nested AP/consec final +7:
  A_mean=1.965291
  dist(roots,[-1,0])=0.912
  reconstruction=ok

single-far final +21:
  A_mean=2.993492
  dist(roots,[-1,0])=0.279
  reconstruction=ok
```

Working heuristic:

```text
low-level nested ears stabilize roots away from [-1,0];
high-level far ears keep roots near [-1,0];
nonnested/broken ears create real-root contacts or named proof debt.
```

## Translation of the ear-decomposition cue

The graph-theory ear facts are best used as grammar, not as imported theorems:

```text
directed ear -> retained one-runner payload A_t
odd ear      -> parity split of A_t
nested ear   -> AP/consecutive legal refinement
nonnested ear -> root collision or named sidecar debt
```

This is the LRC counterpart to HYP-2879's strong tournament one-vertex-ear
calculus.  The tournament vertices for the next pass should be payloads,
root-motion events, and proof obligations, not runners.

## Next ledger

Build `lrc14_lee_yang_ear_payload_ledger` over HYP-2963 and the THM-573
`<=6` multiples-of-7 residual.  Each row should carry:

```text
miss_count_pgf_coefficients
miss_count_pgf_root_multiset
lee_yang_negative_interval_distance
root_axis_gap_deg
fugacity_winner_profile
last_legal_ear
ear_payload_A_vector
ear_payload_mean_level
ear_payload_parity_bias
root_motion_reconstruction_status
nested_ear_status
negative_interval_contact_status
destroyed_coordinate
terminal_exit_or_named_debt
```

The proof-facing conjecture is a separation statement:

```text
If a root approaches or meets [-1,0], then the packet has high-mean far-ear
payload, nonnested ear debt, component-bound debt, first-obstruction debt,
K33/THM-572 state-lift debt, or AP/GW boundary status.
```
