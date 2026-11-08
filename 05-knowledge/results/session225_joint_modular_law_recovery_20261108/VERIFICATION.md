# Fresh verification of the reconstructed packet

Status: FINITE-EXACT for the finite universes below. These are new checks, not recovered old logs.

## Author's portable exact graph run

Executed the exact bytes of check.py with Python3 from stdin and no arguments. Exit code 0. Python standard library only; no file/network access. Output:

q=3 states=45 SCCs=6 PASS first-density, second-law, every-SCC output and line; L(0)=34/27
q=5 states=225 SCCs=10 PASS first-density, second-law, every-SCC output and line; L(0)=151/125
q=7 states=637 SCCs=15 PASS first-density, second-law, every-SCC output and line; L(0)=400/343
q=9 states=1377 SCCs=27 PASS first-density, second-law, every-SCC output and line; L(0)=11/8
q=15 states=6525 SCCs=61 PASS first-density, second-law, every-SCC output and line; L(0)=5179/3375
q=21 states=18081 SCCs=90 PASS first-density, second-law, every-SCC output and line; L(0)=27389/18522

Reproduce from this directory: python3 check.py
Optional larger controls: python3 check.py 1 25 35 49 63

The optional extended command was attempted but its process could not be created because the execution transport disconnected. It is NOT counted as a fresh completed run. The default run above had already completed successfully.

Universe: q^2*(2q-1) normalized states (D,E,U,d), E in {0,d}, deduplicated at d=0. The code constructs actual O/Z edges, verifies indegree2, computes directed SCCs independently, and checks closure. It uses exact Fraction arithmetic for projector masses and integer accumulation for the finite-sum prediction. It checks every first-J state, every output residue, and output/transversal counts in every SCC, not merely a fitted zero factor. Controls include nonprimitive p=7 and composite q=15,21.

## Independent fresh controls

The independent reviewer reran exact BigInt rational branch computations in memory:
- all four kernels, every residue at q=3,5,7,9,25;
- each branch sums over r to q;
- all branches at q=15 equal the CRT products;
- L_15(0)=5179/3375, product of marginals=5134/3375, difference=1/75.

The reviewer also freshly built the restricted q=63,d=21,g=1 graph. It has 4536 vertices and six SCCs of 756 vertices. Every vertex has indegree2 and every edge stays in its SCC. In every SCC all 63 F-residue counts equal12. On each E=0 slice, for every slope a modulo63, all 63 counts of D-aU equal6. This is a scoped q63 check, not a fresh full-q63 graph run.

Regenerate the latter: enumerate D,U in0..62 with gcd(21,U,63)=1 and E in {0,21}; build O=(D,0,2U), Z=(D-E+U,21,2U-21) modulo63; run two-pass directed DFS/Kosaraju, then the stated exact counts.

## Limits

The old broader graph runs are not presented as fresh. No full repository test suite, Lean proof, formal axiom audit, uniform-in-q estimate, height-window census, scalar-coverage verification, or CI pass is claimed. The analytic all-q theorem is independently reviewed; these finite controls supplement, not replace, its proof.
