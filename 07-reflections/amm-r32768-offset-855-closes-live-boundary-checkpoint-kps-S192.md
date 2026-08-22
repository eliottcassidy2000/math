# AMM R=32768: offsets 855 and 856 close -- live boundary checkpoint

**kps-S192/S193/S195/S196 checkpoint, 2026-08-21.  Three FINITE-EXACT state
runs pin the adjacent fixed-Rule-A split, and a disjoint state-only engine has
now independently hostile-audited the lower endpoint.**

The pinned THM-3644 Rule-A engine

```text
sha256=8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad
```

was run without trace at the next dyadic epoch.  The exact stable records are

```text
R=32768,

D0=855:
  (855,'CLOSED',20238,20449,32551,40044,16383,8132,None),
  elapsed=4701.792 seconds;

D0=856:
  (856,'CLOSED',20233,20450,32549,40045,16383,8134,None),
  elapsed=4223.111 seconds.
```

The tuple fields are inherited literally from THM-3644's `stable_result`
contract; this checkpoint does not reinterpret them or infer monotonicity in
the offset.  In particular it proves only that the fixed Rule A closes at
the two pairs `(R,D0)=(32768,855),(32768,856)`.

The predicted split between 855 and 856 is therefore **REFUTED**.  The first
full-diagnostic run at `D0=854` ended inconclusively with a Python
`MemoryError` at

```text
junkL1 += abs(v-u).
```

This variable is an observer: it is reset and accumulated for trace reporting
but is never used in a branch, clamp, state update or terminal decision.  The
failure therefore supplies **no** Rule-A status and is not evidence of death,
closure or mathematical infeasibility.

An exact state replay with only the two `junkL1` observer updates removed by
an AST transformation has now completed.  The pinned dynamics are otherwise
byte-identical, and the transformed engine reproduces the original
`stable_result` on the hostile/positive controls

```text
(R,D0)=(64,0),(64,1),(512,4),(512,5).
```

The exact lower record is

```text
D0=854:
  (854,'DIE',8246,20448,25379,40043,4394,8130,19875),
  elapsed=4076.338 seconds.
```

Thus `854:DIE` and `855:CLOSED` pin a one-step local fixed-Rule-A transition.
No monotonicity assumption is being used: the result says nothing about
offsets below 854 or above 856.  In particular THM-3026's admissible-profile
monotonicity is not being conflated with this fixed heuristic engine's
status.

The completed run used a scratch copy whose SHA-256 is

```text
78bed4b5304a64dfaea3942d7d658ee676c9867853d254cb9b64dfb78dcfc7cb.
```

It differs from the durable runner pinned below by exactly one terminal blank
line.

The independent audit did not reuse the AST transformation.  It wrote a new
state-only engine which streams the coefficients of `2G_R`, retains only the
degree, sparse junk state, feed flag and `c_0`, and contains no trace or
`junkL1` observer.  Direct source inspection confirms that the two production
`junkL1` updates feed only `trace[*]["junkL1_bits"]` and cannot influence a
clamp, sparse-state update, branch, `minus2`, `T6b`, death row or overflow
bit count.  The independent ledger contains

```text
45 Python-integer production/state-only controls,
5 FLINT integer controls spanning DIE and CLOSED outcomes,
65,536 Arb-certified degree floors with zero ambiguous intervals,

FLINT target:
  (854,'DIE',8246,20448,25379,40043,4394,8130,19875),
  elapsed=1835.3803317 seconds.
```

Thus the lower endpoint is independently hostile-audited.  The artifact
remains called a checkpoint only because no separate theorem ID has been
reserved; its fixed-instance finite-exact content no longer has a pending
execution gate.

Reproduction command:

```bash
python3 -B 04-computation/amm12592_r32768_offset_boundary_probe_kps_s192.py 855 856
python3 -B 04-computation/amm12592_r32768_offset_boundary_lowmemory_probe_kps_s193.py 854
```

The original durable runner imports
`04-computation/amm12592_R16384_offset_transition_thm3644.py` and calls its
source-pinned `stable_result` at the requested offsets.  The low-memory runner
source-pins the same engine, AST-selects the same three functions, removes
exactly the two diagnostic `junkL1` additions, and cross-checks four terminal
controls before the long replay.  The disjoint audit instead reconstructs the
state transition directly and supplies the independent implementation gate.

Checkpoint artifacts:

```text
runner sha256:
65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0

stored output sha256:
a3c1c856979fd7b0b21ca12f93d08b3340345fd175187b9c724a14bc91b97710

low-memory state runner sha256:
962577096266bdc6df7b004eb38412b85a8a2b47eb9dd000f4666e9829043207
```
