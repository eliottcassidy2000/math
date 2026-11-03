# Boundary carry: matching covers and uniform ternary returns

**PROVED scoped / FINITE-EXACT; full {7,21} classification OPEN.**
Same-agent proof review with independently implemented exact algorithms;
no subagent audit is claimed.

From the named H=31 source with its old pins retained, an infinite family

    h_m = 8*3^m + 6*2^m - 5, m>=1,

now supports **all four two-step ternary branches**, with explicit native
programs at every m:

    h -> 3h+2 -> 9h+8 or 9h+10,
    h -> 3h+4 -> 9h+14 or 9h+16.

The repaired steps keep every original vertex and clone address. Temporary
births, native complements, deletions and maximum order have separate integer
budgets. Source labels and named pins survive; their response profiles are
not claimed invariant. These are raw H counts, not H/Aut.

The structural advance is [BC1](PROOFS.md): the largest exponential base in
a transitive-clone family is exactly one plus a bipartite matching number
across the root's two neighbourhoods. The 3^m clock therefore has an exact
two-vertex boundary cover. [BC2](PROOFS.md) then compiles changes outside that
root star with a constant helper budget, uniformly over the clone size.
[BC7](PROOFS.md), prompted by incoming guard/profile work, proves the exact
d+K-dimensional function space for a crossing-matching cap d and queries
using at most K path components. This is a linear moment budget, not a
finite-state native realization theorem.
[BC8](PROOFS.md) integrates the incoming arbitrary-child law: the same
programs work on any selected-path child S but implement 3H+cH(S).
The exact additional obligation for a constant offset is c(H(S)-1).

The failed short programs had a real size obstruction. The four second-step
functions have minimum marked base orders **9,9,8,10**, certified by complete
rooted-neighbourhood censuses through nine and explicit upper witnesses.
This says nothing about minimum orders for just one scalar value.

## Inheritance and concept board

Closest proved mechanisms: the common-frame clone lift in
[adaptive LIFTING](../adaptive_four_returns_20261103/LIFTING.md), the ordered-run
kernel [THM-3121](../../../01-canon/theorems/THM-3121-path-cover-walk-content-substitution-kernel.md),
and the literal one-helper [native compiler](../adaptive_metagraph_reformulation_20261103/METAGRAPH.md).
The named H31 switch is inherited from [BRANCH](../source_selection_closure_20261103/BRANCH.md).
The canonical hostile is [H701's missing 1-ready interface](../return701_algebra_20261103/RETURN701.md).
The corrected near miss is confusing an attained scalar, a retained name,
or a formal phase matrix with a live legal return interface. The least-used
sidecar recovered is the root's crossing matching, including its actual
endpoints rather than only their number.

| Lane | Live object | Result / next decisive test |
|---|---|---|
| Anchor | Actual 3H+2/+4 returns | All four two-level branches lifted for all m; arbitrary depth remains open. |
| Boundary | Matching and two-cover | Exact growth theorem; test a bounded active-cover transfer under native repairs. |
| Niche | Older H31/H33 score examples | Explicit source pair has a changing two-clock defect, not a constant sign bit. |
| Niche | Nice square-sum pairs | Interface closure is the relevant analogy; scalar Hamiltonicity is weaker. |
| Wildcard | 32+-1, 13+-1 and golden arithmetic | Exact split/nonsplit clocks with scalar fibres, but no free native action. |
| Hostile | Short one-birth continuation | Size lower bounds force a larger source representation; record its operation debt. |

The board was revisited after: the first fourteen uniform programs; failure
of their one-birth reentry; the complete order-eight obstruction; order-nine
partial recovery; and the two-cover-filtered order-ten repair. Each stage
changed the proposed interface rather than merely extending the scalar scan.

See [PROOFS](PROOFS.md), [connections](CONNECTIONS.md), [audit](AUDIT.md), and
the exact programs in [returns.json](returns.json).

## Reproduce

From the repository root, with Python 3 and a C++20 compiler:

```sh
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/neighbourhood_bank.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/scout.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/scout.py --input 05-knowledge/results/boundary_carry_returns_20261103/next_input_0.json --output next_0.json
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/scout.py --input 05-knowledge/results/boundary_carry_returns_20261103/next_input_1.json --output next_1.json
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/next_level.py
g++ -O3 -std=c++20 05-knowledge/results/boundary_carry_returns_20261103/catalog.cpp -o /tmp/boundary-carry-catalog
g++ -O3 -std=c++20 05-knowledge/results/boundary_carry_returns_20261103/audit_catalog.cpp -o /tmp/boundary-carry-audit-catalog
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/reproduce_censuses.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/returns.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/arithmetic.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/boundary.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/moment_budget.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/audit.py
python3 -B 05-knowledge/results/boundary_carry_returns_20261103/child_audit.py
```

The final six commands check the saved certificates without
repeating the searches. `next_level.py` uses three local processes, not agents.
The small orbit banks contain 1,1,2,4,12,56 classes at orders 1 through 6;
their completeness is independently checked by literal relabelling orbits.

## Handoff: the next proof target

Develop a carrier retaining the two-cover vertices, ordered path endpoints,
the distinct 3^m and 2^m coefficients, and the full pinned source relation
needed by the repair word. Prove an invariant source bank on which **both**
ternary offsets return to that bank with charged helpers. First hostile test:
the four new target families, especially the 295 family whose minimum order
already exceeds a two-birth plan from the 97 source. Do not substitute a
formal observer matrix or an isomorphic unpinned source for a legal return.

The most useful small question is now: can two-cover sources be normalized
by a bounded set of active-boundary rules while leaving the exterior only
through its exact ordered-cover responses? BC2 supplies reachability between
supplied same-star targets; it does not supply the target's existence. A
positive answer needs both halves. The spectral classification remains open.

[BC6](PROOFS.md) records the immediate boundary test: the four outputs have
constant exterior ready defects 0,0,none,17, respectively. None is 1-ready
or 2-ready at its base, even after choosing another vertex. The all-m
defects are certified by the four-run bound and an independent word count.
