# AMM R=32768: offsets 855 and 856 close -- live boundary checkpoint

**kps-S192/S193 checkpoint, 2026-08-21.  Two FINITE-EXACT independent
closure runs; the lower adjacent state replay and independent audit are
pending.**

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

An exact state replay is now running with only the two `junkL1` observer
updates removed by an AST transformation.  The pinned dynamics are otherwise
byte-identical, and the transformed engine reproduces the original
`stable_result` on the hostile/positive controls

```text
(R,D0)=(64,0),(64,1),(512,4),(512,5).
```

If 854 dies, the pair 854/855 will pin a one-step local Rule-A transition.  If
854 also closes, the search must continue downward; no monotonicity assumption
is being used.  In particular THM-3026's admissible-profile monotonicity is not
being conflated with this fixed heuristic engine's status.

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
controls before the long replay.  A durable theorem package will replace this
checkpoint only after the adjacent state result and an independent replay are
available.

Checkpoint artifacts:

```text
runner sha256:
65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0

stored output sha256:
24703605154f9f55d5a9e0dc7600c2873006b774ac3a77ea71b2b4bc9ac8dc25

low-memory state runner sha256:
962577096266bdc6df7b004eb38412b85a8a2b47eb9dd000f4666e9829043207
```
