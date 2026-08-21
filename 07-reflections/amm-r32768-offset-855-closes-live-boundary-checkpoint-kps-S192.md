# AMM R=32768: offsets 855 and 856 close -- live boundary checkpoint

**kps-S192 checkpoint, 2026-08-21.  Two FINITE-EXACT independent runs;
lower adjacent boundary and independent audit pending.**

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

The predicted split between 855 and 856 is therefore **REFUTED**.  The newly
compelled hostile run at `D0=854` remains in progress.  If 854 dies, the pair
854/855 will pin a one-step local transition.  If 854 also closes, the search
must continue downward; no monotonicity assumption is being used.

Reproduction command:

```bash
python3 -B 04-computation/amm12592_r32768_offset_boundary_probe_kps_s192.py 855 856
```

The durable runner imports
`04-computation/amm12592_R16384_offset_transition_thm3644.py` and calls its
source-pinned `stable_result` at the requested offsets.  A durable theorem package
will replace this checkpoint only after the adjacent hostile result and an
independent replay are available.

Checkpoint artifacts:

```text
runner sha256:
65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0

stored output sha256:
193250d6f2d9de4c960ea02b03c25c86b3d93bc835f61ea16b738a8f1e65ec3b
```
