# Mixed-(2,2), c=4: exact analytic census and sparse residual probes

## Scope

This packet proves the stated analytic exclusions and verifies 37 selected
residual cores. It does not exclude the entire mixed-(2,2), c=4 morphology.
The remaining residual cores require further certificates. A relaxed-cone
obstruction, if found in later work, would be a method limitation rather than
a poset counterexample without an additional realizability argument.

The ambient assumptions are a finite seven-thin poset with no balanced pair,
a strict-majority total order, and a maximum-index-span incomparable root x,y.
The root has two lower slaves, two upper slaves, and four common neighbors.
The core labels are `x,y,b1,b2,z1,z2,z3,z4,a1,a2 = 0,...,9`.

## Degree-valid catalogue

Every possible core is represented by a naturally labeled poset on its eight
nonroot vertices, in lower/common/upper group order, together with `bi<y` and
`x<aj`. Reverse group comparabilities would contradict the defining root
incomparabilities. Label within each group by a linear extension.

There are 2,800,472 naturally labeled eight-vertex nonroot posets. For a common
vertex with C nonroot comparabilities, its core incomparability degree is
`9-C`, so it requires C>=2. For a lower or upper slave, the degree is `8-C`,
so it requires C>=1. Both roots have degree exactly seven. Thus the all-vertex
seven-thin core filter is precisely:

- each of the four common vertices has at least two nonroot comparabilities;
- each of the four slaves has at least one nonroot comparability.

This leaves 1,165,395 degree-valid cores and rejects 1,635,077. The discovery
catalogue was streamed using predecessor ideals, with an exact resumable DFS
cursor. Independent replay constructs successor-closed subsets in reverse
vertex order and checks equality with the entire stored core set using exact
28-bit relation codes. It checks uniqueness and transitivity as well as the
counts. No probabilistic set-equality test is used.

All four common neighbors lie between x and y in majority order by the
maximum-span condition. Hence the root span is at least five; the published
six-thin short-span theorem does not remove any candidate here.

## Analytic exclusions

For incomparable u,v use the directed defect

`delta(u,v)=|D(v)\D(u)|+|U(u)\U(v)|`.

The included injection lemma gives `p(u<v)<=2/3` when delta<=1. Seed the roots
and every core-saturated vertex, propagate the forced lower-exterior direction
upward and the forced upper-exterior direction downward, and use the bound

`delta_P(u,v)<=delta_Q(u,v)+alpha(u)+beta(v)`.

Here alpha(u)=0 when every lower exterior vertex is forced below u and is
otherwise bounded by `7-degree_Q(u)`; beta is the dual upper-exterior bound.
The included side-specific proof justifies these bounds for arbitrary root
cores. Two directed bounds at most one give a balanced pair. A single such
bound forces the opposite strict-majority direction.

The exact analytic classification is:

| class | cores |
|---|---:|
| directed-defect balanced-pair witness | 733,492 |
| cycle of forced majority arcs | 216,636 |
| remaining residual cores | 215,267 |
| all degree-valid cores | 1,165,395 |

The cycle arcs include poset comparabilities, x<y, x<zi<y, and directions
forced by the defect bound. The compact audit stores an actual directed cycle,
not merely a numerical cycle-detection status. Every witness edge is replayed
independently. Residual records are checked to have no cycle of all these
forced arcs, which also rules out an omitted two-sided defect witness.

## Sparse residual probes

Probe the first 16 residual cores and every 10,000th residual core thereafter,
where the residual order is the discovery catalogue's deterministic order.
In zero-based residual indices the latter are 9,999, 19,999, ..., 209,999.
There are 37 probes in total. Each probe covers all admissible projected
majority orders of its core, using analytic directed bounds and exact positive
boundary-generator certificates.

All 37 pass independent replay, with 366 generator-certified orders and no
LP-derived certificate. This is selected evidence only; it does not certify
the other 215,230 residual cores.

The generator reduction is the same exact nonnegative prefix/suffix increment
argument used in the included saturation note: lower exterior vertices can
interleave only before the first root, and upper exterior vertices only after
the second root. Their counts telescope into nonnegative generator weights.
A saturated core vertex makes every increment after it on the left, or before
it on the right, vanish. For a certificate with positive integer weights w_i,
comparison counts N_i, and generator size E, the verifier checks

`3*sum(w_i*N_i)<=2*E*sum(w_i)`

on every retained generator, then checks that every admissible candidate
majority order is covered by forward comparisons in one such certificate.
All arithmetic is exact and all checks remain enabled under Python `-O`.

## Reproduction and preservation

From this packet directory:

```sh
python -O verify_c4_analytic_and_probes_20261019.py
```

Only Python's standard library is needed. The verifier checks the complete
catalogue equality, every analytic witness, the exact residual-inventory
correspondence, the specified probe schedule, and every selected positive
certificate. It writes `independent_census_replay_receipt_20261019.json`.

The original 13,152,413-byte gzip audit and original probe stream are preserved
in the discovery stage. The compact classification retains the entire core
set and analytic witnesses; the repack receipt gives the source and compact
hashes. A later exhaustive residual search must cover this exact 215,267-core
inventory before claiming the full c=4 exclusion.
