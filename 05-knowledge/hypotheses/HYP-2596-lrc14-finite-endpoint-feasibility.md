# HYP-2596: LRC(14) Finite Endpoint Feasibility

**Status:** negative feasibility result for brute force; positive route map for
structural finite cores.

HYP-2595 would certify all `V >= 711` if the colored resonance discrepancy bound
and the HYP-2593 `Sigma` floor are promoted to theorems.  The remaining
question is whether the finite endpoint `V < 711` can simply be computed.

Script: `04-computation/lrc14_finite_endpoint_feasibility_codex.py`  
Output: `05-knowledge/results/lrc14_finite_endpoint_feasibility_codex.out`

## Findings

The raw S3 coordinate count below `V=711`, writing
`S=P union {V-e:e in E}` with `0 in E`, is

`1675268836111283636384624399`.

Even the `k=3` layer alone has `16070952040` raw shapes before covering,
primitivity, and exact-M pruning.  Thus a literal all-shapes THM-524 exact
endpoint is not feasible by brute force.

The tempting replacement "just check the constructive modulus `q=14V`" also
fails at low `V`.  A binary MILP over forbidden residue classes searched for
primitive q-covering 13-sets containing `V` whose classes cover all residues
modulo `14V`.  It found:

`V=15`,
`S=(1,2,3,4,5,8,9,10,11,12,13,14,15)`.

This row has no `q=210` CRT witness, but exact THM-524 gives

`M(S)=2/23` at `tau=4/23`,

so LRC is comfortably true.  Therefore the large-V colored modulus is not the
whole finite endpoint method.

The same MILP proved infeasibility for selected `V=14,16,20,28,42,56`, and hit
time limits at `V=70,84,98,120`.  This is useful as an obstruction classifier,
but not yet an endpoint proof.

A small calibrated exact-M core scan for `k=3`, all large speeds `<=42`, scanned
`61161` covering primitive S3 rows, exact-checked `7` near-threshold rows, and
found `0` below `1/14`.  The worst exact checked row was

`S=(1,2,3,5,7,8,9,11,12,13,23,30,42)`,
`M=3/37` at `tau=11/37`.

This supports the existing lesson from the proved `k=2` slice: exact-M checking
is tractable only after a structural theorem has reduced the hard core.

## Tournament Analysis

Vertices were proof strategies rather than runners:

`naive_shape_enumeration`, `q14V_milp`, `k2_finite_core`,
`arc_width_dropmax`, `colored_resonance_largeV`, `endpoint_exactM_core`.

The pairwise observable was `(proved scope, endpoint cost, remaining proof gap)`,
with larger scope, smaller cost, and smaller gap winning.  The resulting
tournament was transitive:

- score histogram `{0:1,1:1,2:1,3:1,4:1,5:1}`;
- directed 3-cycles `0`;
- SCCs all singleton;
- Hamiltonian paths `1`;
- leader `k2_finite_core`.

Readout: the proven `k=2` finite-core proof is the model to imitate; the
colored-resonance route is the best large-V lever; direct enumeration is the
sink.

## Proof Target

Do not attempt the endpoint as `V<711` brute force.

Instead:

1. Prove HYP-2595's `Delta <= 8*(k+c_GP)+1`.
2. Prove the HYP-2593 uniform `Sigma` floor.
3. Classify low-`V` `q=14V` obstruction rows by exact-M denominator families.
4. Generalize the k=2 finite-core pattern: drop one large speed, prove a
   lower floor for the residual safe-arc width except in a bounded hard core,
   then exact-check that core.

This does not prove LRC(14), but it prevents a false endpoint plan and gives a
more realistic final architecture.

See also `HYP-2595`, `HYP-2594`, `HYP-2593`, `HYP-2589`, `HYP-2581`, and
`OPEN-Q-108`.
