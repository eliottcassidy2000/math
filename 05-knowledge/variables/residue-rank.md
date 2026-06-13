# Variable: residue rank

**Symbol:** `rank_res`
**Type:** integer statistic on a projected support family
**Defined in:** HYP-1780 / codex-2026-05-30 S355

This registry entry uses `rank_res` for support-family independence rank, not
for older local variables named `rank_res` that denote matrix restriction rank
inside beta-complex scripts.

## Definition

Let `F(T)` be a support family attached to a tournament, such as the directed
odd cycles used to build `Omega(T)`.  Let `pi` be a projection or forgetful map,
and let `R_pi(T)` be the family of supports that survive the projection.

For the deletion projection at vertex `v`:

```text
R_v(T) = {C in F(T) : v notin C}.
alpha_k(R_v) = number of k-tuples of pairwise disjoint survivors.
rank_res(v) = max {k : alpha_k(R_v) > 0}.
```

The max-loss deletion residue uses a vertex

```text
v* in argmax_v #{C in F(T) : v in C}
```

and records `rank_res(v*)` together with the small alpha profile of `R_v*(T)`.

## Values at Small Examples

| Tournament | Max-loss deletion residue |
|---|---|
| H=63 single-core `1001100` | `alpha=(0,0)`, `rank_res=0`, exact kill |
| H=63 single-core `1100110` | `alpha=(0,0)`, `rank_res=0`, exact kill |
| THM-025 n=9 | `alpha=(2,1)`, `rank_res=2`, `I(R,2)=9` |
| Paley T7 | `alpha=(20,1)`, broad rank-2 residue |
| Interval T7 | `alpha=(16,2)`, broad rank-2 residue |

## Equations Involving This Variable

- `rank_res(v)=0` iff the deletion projection kills all supports.
- `rank_res(v)>=2` iff at least two surviving supports are disjoint.
- For deletion residues on `n<=9` tournaments, the observed OCF evaluation is
  `I(R_v,2)=1+2*alpha_1(R_v)+4*alpha_2(R_v)` because at most two disjoint odd
  cycles can fit after deleting one vertex.

## Relationships

- Refines [projection defect](projection-defect.md): `loss_v` and `keep_v`
  measure size; `rank_res(v)` measures the independence shape of the survivor.
- Related to [Omega(T)](omega-graph.md): the rank is the independence number
  layer visible inside the induced residue subgraph.
- Related to [alpha_k](alpha-k.md): residue rank is the highest nonzero
  alpha layer after projection.

## Sources

- `04-computation/tournament_tda.py`
- `04-computation/residue_rank_probe_s355.py`
- `05-knowledge/results/residue_rank_probe_s355.out`
- HYP-1780 and HYP-1785

## Tags

#residue-rank #projection-defect #omega #ocf #tournament-tda
