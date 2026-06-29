# LRC14 Gate-Escape Transversal Router

HYP-3453 is a small but important join: it connects HYP-3438's exact survivor
gate words to HYP-3450/HYP-3451's component-cover obstruction graph.  Rebased
over HYP-3452 and HYP-3439, it should be read as the bank-level companion to
the AP-with-`84m` phase audit and rescue-core bridge: HYP-3452 gives the
canonical tail phase, HYP-3439 joins one-branch rescue cores to survivor
routes, and HYP-3453 separates clean exits from gate exits across the broader
covering bank.

The new information-theoretic point is that "low endpoint rank" is still a
compressed statement.  It forgets whether the escape is a clean component
where the dead-cover obstruction is absent, or a gate sitting against bad-core
blocks where the dead-cover/Menger graph can actually attach.  That is a
failure of compression akin to the earlier overlap-tax loss: the scalar or
component summary is true but not composable.

Exact readout:

```text
rows_with_low_rank_component=135/135
rows_with_low_rank_gate=133/135
rows_with_dead_components=130/135
rows_with_dead_components_and_low_rank_gate=130/130
rows_dead_without_low_rank_gate=[]
dead_zero_clean_only_rows=[random_covering_044, random_covering_053]
```

So the finite lemma can split cleanly:

```text
dead_components(row) > 0  =>  rank <= 2 survivor gate exists.
no rank <= 2 survivor gate => dead_components(row) = 0.
```

This is stronger than merely reprinting HYP-3450.  It says that any attempted
full saturation must expose a graph-composable gate, not just a survivor
component somewhere.  The clean-only exceptions are already discharge cases,
because they have no dead components at all.

The AP-tail row keeps the expected corridor-fence shape: four rank-`2` gates
with labels `B1:7|E:84`, `E:84|B1:5`, `B0:5|E:84`, and `E:84|B0:7`, all with
adjacent cover delta `(1,1)`.  At `84x05`, two gate exits remain and two clean
`E:420|E:420` exits appear.  That is a useful warning for the all-`m` proof:
do not require every low-rank escape to be a gate; require the gate only when
the dead-cover obstruction is nonempty.

Next proof move: prove the gate-transversal implication without the sampled
bank.  The best route seems to be a finite interval lemma:

1. If `E_safe` has a dead-both component, the boundary between dead and
   nondead components must cross a mixed component.
2. Any such mixed component has an exact HYP-3438 gate word.
3. Endpoint-spine or corridor-fence parity bounds force one gate endpoint rank
   to be at most `2`, unless the obstruction graph disconnects and becomes a
   clean exit.

That would turn HYP-3451's Green/Menger/conductance graph from a heuristic
router into a formal local-to-global proof skeleton.
