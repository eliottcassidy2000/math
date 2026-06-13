---
id: HYP-1961
status: OPEN
source: codex-2026-06-01-S491
related:
  - HYP-1895
  - HYP-1903
  - HYP-1921
  - HYP-1930
  - HYP-1940
  - HYP-1942
  - HYP-1950
  - HYP-1960
---

# HYP-1961: LRC pressure DAG peel layers are endpoint-certificate data

## Statement

When an LRC two-neighbor pressure search returns a DAG, the result should be
read as a positive certificate object.

The mobile pressure relation at a time `t` compares deletion reliefs:

```text
relief_i(j) = nearest_distance_i(after deleting j) - nearest_distance_i.
```

Orient `j -> i` when deleting `j` gives more relief to `i` than deleting `i`
gives to `j`.  A strict pressure SCC is the disproof-like signal: it says
blocking dependence closes into a labelled loop.  A strict pressure DAG says
the opposite.  The pressure dependencies can be topologically sorted, giving
source and sink peel layers.

The proposed proof use is:

```text
pressure DAG + endpoint-private rows -> ordered peel certificate;
no compatible peel row -> labelled pressure SCC obstruction.
```

S503 makes the owner-cycle obstruction formal.  By THM-379/THM-380, if a
nonempty terminal endpoint core is pressure-realized, then a pressure DAG
excludes it.  Thus the peel-layer program should now record not only pressure
layers, but also which endpoint-core protection incidences each pressure edge
realizes.

## Evidence

`lrc_pressure_dag_s491.py` audits the pressure relation on the current hard
rows:

```text
n14 initial
n14 d=7
n14 d=14
n18 initial
n18 d=3
n18 d=9
n18 d=18
```

Across the bounded pressure search windows:

```text
case          times  cyclic  max_scc  max_tri  max_chain
n14 initial      42       0        1        0          1
n14 d=7         113       0        1        0          3
n14 d=14        183       0        1        0          4
n18 initial      40       0        1        0          1
n18 d=3         100       0        1        0          3
n18 d=9         233       0        1        0          3
n18 d=18        425       0        1        0          4
```

So every sampled strict pressure graph was a DAG.  Representative gap-midpoint
peels are explicit:

```text
n14 d=7 source peel:
  {1,14,49} -> {7,35,56,77,91} -> {0,84}

n18 d=18 source peel:
  {1,36,90,162,270} -> {18,54,126,180,198,288,306} -> {0}
```

These are not empty observations.  They are ordered dependency certificates.
The endpoint-debt rows remain positive:

```text
n14 d=7:  gap/th=5/924,  unprotected=84,  product=5/11
n14 d=14: gap/th=5/1848, unprotected=168, product=5/11
n18 d=9:  gap/th=1/176,  unprotected=176, product=1
n18 d=18: gap/th=1/352,  unprotected=352, product=1
```

## Predictions

1. Known near-miss LRC rows should continue to produce pressure DAGs unless
   they are genuinely counterexample-like.
2. The first serious disproof signal is not a smaller scalar gap; it is a
   labelled pressure SCC after endpoint-private leaves are removed.
3. A branch-and-bound proof should store source/sink pressure layers and try
   to match them to endpoint-private rows.
4. If all bounded perturbations of an `n=18` ladder remain pressure-DAG, then
   `n=18` is better treated as a mixed-torsion proof lab than as a disproof
   search target.

## Sources

- `04-computation/lrc_pressure_dag_s491.py`
- `05-knowledge/results/lrc_pressure_dag_s491.out`
- `07-reflections/lrc-pressure-dags-s491.md`
- `04-computation/lrc_endpoint_pressure_formal_s503.py`
- `05-knowledge/results/lrc_endpoint_pressure_formal_s503.out`
- `07-reflections/lrc-endpoint-pressure-formalization-s503.md`
- HYP-1930
- HYP-1942
- HYP-1950
- THM-379
- THM-380
