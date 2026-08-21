# AMM R=32768: offset 855 closes -- live boundary checkpoint

**kps-S192 checkpoint, 2026-08-21.  FINITE-EXACT single-run result; adjacent
boundary and independent audit pending.**

The pinned THM-3644 Rule-A engine

```text
sha256=8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad
```

was run without trace at the next dyadic epoch.  The exact stable record is

```text
R=32768,
D0=855,
stable_result=
  (855,'CLOSED',20238,20449,32551,40044,16383,8132,None),
elapsed=4701.792 seconds.
```

The tuple fields are inherited literally from THM-3644's `stable_result`
contract; this checkpoint does not reinterpret them or infer monotonicity in
the offset.  In particular it proves only that the fixed Rule A closes at
this one pair `(R,D0)=(32768,855)`.

The adjacent run at `D0=856` remains in progress, and the newly compelled
hostile run at `D0=854` has been started.  If 854 dies, the pair 854/855 will
pin a one-step local transition.  If 854 also closes, the search must continue
downward; no monotonicity assumption is being used.

Reproduction command:

```bash
python3 -B 04-computation/amm12592_r32768_offset_boundary_probe_kps_s192.py 855
```

The durable runner imports
`04-computation/amm12592_R16384_offset_transition_thm3644.py` and calls its
source-pinned `stable_result(run_fast,32768,855)`.  A durable theorem package
will replace this checkpoint only after the adjacent hostile result and an
independent replay are available.

Checkpoint artifacts:

```text
runner sha256:
65b6fab453a0adb156bc79e1f9168fd46d95037fd345e2322861ec1f76f15ba0

stored output sha256:
3d76ff740c034d13200dd9c77b1952a931eef3795599e7bda9a2db0d87cbb2a5
```
