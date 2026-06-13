# Variable: projection defect

**Symbol:** `delta_proj`
**Type:** family of integer / normalized statistics
**Defined in:** opus-2026-05-29-S12

## Definition

A projection defect measures how much tournament structure is lost under a
forgetful projection. S12 used three concrete versions:

1. **Support multiplicity defect**

```text
delta_support(T) = |V(Omega(T))| - |{ vertex supports of odd cycles }|.
```

2. **Old-projection deletion loss at a vertex**

```text
loss_v(T) = #{ odd directed cycles containing v }
keep_v(T) = #{ odd directed cycles not containing v } = |V(Omega(T-v))|
rho_v(T) = loss_v(T) / |V(Omega(T))|.
```

An old-projection kill vertex has `keep_v(T)=0`.

The S355 refinement also records the independence rank of the surviving
residue:

```text
rank_res(v) = max {k : alpha_k({C : v notin C}) > 0}.
```

Thus `keep_v` measures residue size, while `rank_res(v)` measures whether the
survivor has disjoint-cycle structure.

3. **Odd-n cycle-space projection**

For odd `n`, the projection `T -> T_cycle` is the GF(2) cycle-space component
`(I + L(K_n))*T`. At even `n`, the projection is ambiguous because
`Cut(K_n) ∩ Cycle(K_n)` is nonzero.

## Values at Small Examples

| Tournament | Key defect |
|---|---|
| H=63 cid 2519 | core vertex 3 has `loss=31`, `keep=0`, `rho=1` |
| H=63 cid 3285 | core vertex 0 has `loss=31`, `keep=0`, `rho=1` |
| THM-025 n=9 | vertex 3 has `loss=92`, `keep=2`, `rho=0.979`; the residue has `alpha=[1,2,1]` |
| Paley T7 | all vertices have `loss=60`, `keep=20`, `rho=0.75`; support defect = 44 |
| Interval T7 | all vertices have `loss=43`, `keep=16`, `rho=0.729`; support defect = 23 |
| S355 residue rank probe | H=63 has max-loss `rank_res=0`; THM-025 has max-loss `rank_res=2` with `I(R,2)=9` |

## Equations Involving This Variable

- `keep_v(T) = |V(Omega(T-v))|`.
- `loss_v(T) + keep_v(T) = |V(Omega(T))|`.
- `rank_res(v)=0` iff the deletion residue is empty.
- `rank_res(v)>=2` iff two surviving odd cycles are disjoint.
- If `keep_v(T)=0` and every odd cycle contains `v`, then `Omega(T)` is complete
  whenever all odd cycles pairwise meet at `v`.

## Relationships

- Related to [Omega(T)](omega-graph.md): measures what remains of the odd-cycle
  conflict graph after forgetting a vertex or directed-cycle multiplicity.
- Related to [r_core(s)](single-core-cycle-count.md): single-core complete-Ω
  tournaments are exactly old-projection kills with transitive deletion.
- Related to [residue rank](residue-rank.md): refines deletion loss by the
  independence shape of the surviving cycles.
- Related to `HYP-408` / ghost cycles: the path-homology version asks whether
  through-v-only cycles vanish into `im(d_4)` under old projection.

## Sources

- `04-computation/projection_defect_bridge_s12.py`
- `04-computation/residue_rank_probe_s355.py`
- `05-knowledge/results/projection_defect_bridge_s12.out`
- `05-knowledge/results/residue_rank_probe_s355.out`
- `07-reflections/projection-defect-as-common-residue.md`

## Tags

#projection-defect #omega #old-projection #even-graphs #path-homology
