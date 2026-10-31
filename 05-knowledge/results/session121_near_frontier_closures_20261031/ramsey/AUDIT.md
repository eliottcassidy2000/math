# Audit and reproducibility

2026-10-31. Scope: all-parameter proof plus finite exact controls; no
formal proof assistant certification or unrestricted Ramsey record.

## Independent mathematical review

A separate automated mathematical review checked the all-parameter proof
without relying on numerical interpolation. It found no flaw and
confirmed:

* deficit stabilization follows from the inherited sharp union block;
* children at join/full-C5 roots have strictly smaller clique budgets;
* flattening top-level unions gives the additive-star spectrum exactly;
* the active-component semilinear-star formula retains base activation;
* the fixed-A Presburger compilation and threshold search terminate;
* the least eventual period obeys q|d|L;
* finite jumps do not bound the location of the final jump.

The proof explicitly cites the classical effective Presburger/semilinear
equivalence. It does not describe that equivalence as new or implemented.

## Checker

Run from this directory:

    python check_translation.py > normal.json
    python -O check_translation.py > optimized.json
    cmp normal.json optimized.json

Both executions completed and the receipts are byte-identical. All
checks use explicit exceptions, never optimizable Python assertions.

The primary exact dynamic program computes F for A<=6,B<=16 by all
three root operations. It uses only coordinatewise maximal feasible
C5 budget vectors: every feasible vector can be increased to one, and
capacity is monotone. A separate full-vector DP through A<=4,B<=6
compares the resulting exact values without this pruning. The first
four all-B closed formulas are checked through B=16. Literal adjacency
matrices for all24 optimal plans with A<=4,B<=6 are independently
checked by subset-recursive maximum clique and independent-set search.
This validates realized graphs, not just scalar budget arithmetic.

Four explicit semilinear families are checked in [0,18]^2 by comparing
literal additive closure of their actual points with the activated
component formula. All generators are nonnegative, so bounding both
coordinates cannot omit an intermediate path to an in-box point.
One test specifically excludes a period vector occurring without its
base component. These are finite controls of the algebraic formula,
not its proof.

Three delayed-plateau hostile functions are checked for superadditivity
on x,y=1,...,119. They show why no observed numerical plateau is promoted
to an infinite-tail certificate. Row5 and row6 numerical entries are
exact only on the tested rectangle; their apparent translation laws
outside that rectangle remain unclaimed.

## Provenance

Read from isolated source HEAD
`0578f16726e4e535189fb6074cc2ddb024c95d82`:

`05-knowledge/results/ramsey_neighbourhood_lifts_20261030/GRAPH_CAPACITY.md`

Its G2 supplies the capacity upper bound and sharp actual L-attainer;
its G4 explicitly names eventual exact translation as OPEN. The new
argument closes that target for every positive A. No changes to its
proof of the G2 input were needed.
