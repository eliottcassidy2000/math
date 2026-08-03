# Independent audit: full S364 blockers are exactly scalar ramps at n=14

Status: **FINITE-EXACT / INDEPENDENTLY VERIFIED.**  No flaw found.  This is a
classification inside the S364 residue/open-cell model, not a proof of
LRC(14).  No claim about `n=13` is made here.

## Claim audited

For `n=14`, `k=13`, let `b` range over the open-cell floor patterns

```text
b_i = floor(14*{i alpha}),  i=1,...,13,
```

and let `s` range over `Z/14`.  A residue vector `v` is a full blocker when
every `(s,b)` has some coordinate satisfying

```text
s*v_i + b_i in {0,13} mod 14.
```

The audited classification says that the full blockers are exactly

```text
v_i = m*i mod 14,  m in Z/14.
```

The forward direction is already PROVED for all `n` by
`THM-364-lrc-scalar-ramp-cell-blocking.md`.  The new content is the finite
converse at `n=14`.

## Logical audit

### Open-cell universe

For coordinate `i`, all discontinuities in `[0,1]` are exactly
`a/(14*i)`, `0<=a<=14*i`.  Their exact rational union has `813` points and
therefore `812` consecutive atomic open intervals.  Independent one-third and
two-thirds samples agree on every interval.  The `812` resulting patterns are
all distinct.  Thus the historical S364 pattern deduplication loses no cell at
`n=14`; `14*812=11368` is the exhaustive candidate count.

The same audit at `n=15` gives `961` breakpoints, `960` atomic intervals, and
`960` distinct patterns.

### Gauge sign and iff

`THM-363-lrc-scalar-gauge-reindexing.md` uses the action

```text
v_i -> v_i + m*i.
```

Choosing `m=-v_1` gives the normalized vector

```text
w_i = v_i-v_1*i,  so w_1=0.
```

There is no lost direction: `w=0` iff `v_i=v_1*i` for every `i`.  Therefore a
normalized nonzero full blocker is exactly the same thing as an ungauged
nonscalar full blocker.  The sign in the classifier is consequently correct.

### One-hot classifier under audit

The original one-hot encoding has `13*14=182` variables.  For each coordinate
it uses one width-14 at-least-one clause and all `C(14,2)=91` pairwise
at-most-one clauses.  Hence the exact-one layer has `13*(1+91)=1196` clauses,
including exactly `1183` binary clauses.  Adding the gauge unit, the normalized
nonzero clause, and `11368` candidate cover clauses gives the reported `12566`
clauses.

Given `v_1=0`, its nonzero clause ranges over coordinates `2,...,13` and
residues `1,...,13`; it excludes exactly the zero normalized vector.  A cover
clause lists precisely the coordinate/residue choices for which
`s*v_i+b_i` is `0` or `13`.  Both directions are correct once exact-one is
enforced.

### Independent encoding

The audit script does not import S364, S371, or the one-hot classifier.  Each
residue is represented by four bits; binary values `14,15` are explicitly
forbidden.  Equality atoms are constrained by a five-clause biconditional with
the bits.  A complete 16-value truth-table audit checks both directions of
every biconditional and the domain clauses.  There are no one-hot clauses.

Candidate clauses are ORs of these proved equality atoms.  MiniSat22 and
Lingeling, neither used in the original run, both return UNSAT for normalized
nonzero blockers.  Direct ungauged enumeration then decodes exactly `14`
models, verifies each against all `11368` unreduced candidates, and finds the
scalar multipliers `0,...,13` with no nonscalar model.

## Smaller shift certificate

Shift `0` is universally blocked, independently checked on all `812` patterns.
There is also an exact involution

```text
b -> 13-b,   s -> -s.
```

Indeed, if `y=s*v_i+b_i`, the reflected residue is `-y-1`, which swaps `0`
and `13`.  The cell involution is induced by `alpha -> 1-alpha`.  Consequently
the complete shift-`s` and shift-`-s` clause batches are identical after cell
reindexing.  Only shifts

```text
1,2,3,4,5,6,7
```

are needed.  This reduces the `n=14` cover layer from `11368` to `5684`
raw representative candidates.  Removing identical supports created by the
nonunit shifts leaves `4142` distinct cover clauses, and the independent
normalized CNF has `234` variables and `5083` clauses.  This is a symmetry and
syntactic-deduplication reduction, not a minimal unsatisfiable core.  No
smaller inclusion-minimal shift subset was established; a bounded drop-one
search was solver-hard and was stopped without a claim.

For `n=15`, the same reflection leaves the seven representatives `1,...,7`.
Its `6720` raw representative candidates give `5606` unique cover clauses.
The independent normalized CNF has `266` variables and `6675` clauses;
MiniSat22 and Lingeling both return UNSAT.  This agrees with the targeted
`n=15` result but is not an ungauged all-model enumeration in this audit.

## New generalization lane: a rank-one affine-cover rigidity problem

The exact result suggests a clean all-`n` tangent:

> **HYPOTHESIS.** For every `n>=2`, scalar ramps are the only full blockers in
> the `k=n-1` micro-staircase residue-cell system.

A dependent one-hot scout proves the normalized system UNSAT for every
`3<=n<=12`; `n=2` is trivial.  Together with the independently audited
`n=14,15` cases, the only hole through `15` is `n=13`, whose direct SAT run was
solver-hard and is deliberately left open.  Reproduce the small cases with

```text
python 04-computation/lrc_microstaircase_scalar_rigidity_small_n_scout_opus_20260803.py
```

and compare
`05-knowledge/results/lrc_microstaircase_scalar_rigidity_small_n_scout_opus_20260803.out`.

There is a useful two-dimensional shadow.  Sample the cell immediately to the
right of `alpha=a/n`.  Its pattern is `b_i=a*i mod n`, so a full blocker must
cover every `(a,s) in (Z/n)^2` by the affine tests

```text
< (a,s),(i,v_i) > in {0,-1} mod n                    (A)
```

for some `1<=i<=n-1`.  Thus the columns `(i,v_i)` define a two-level affine
hyperplane cover of the finite rank-two module.  For a scalar ramp
`v_i=m*i`, all columns lie on one rank-one line and `(A)` is automatic: put
`t=a+ms`; if `t=0` use residue `0`, if `t` is a nonunit choose
`i=n/gcd(t,n)` to obtain `0`, and if `t` is a unit choose `i=-t^(-1)` to
obtain `-1`.

The boundary shadow is necessary but not sufficient.  At `n=14`, the
nonscalar normalized vector

```text
(0,1,0,0,0,0,0,0,0,0,0,0,0)
```

passes all `14*14` right-adjacent tests `(A)` yet misses `336` candidates in
the full `812`-cell arrangement.  This hostile identifies the missing
ingredient: a human proof must combine finite-module cover rigidity with the
interior Farey-cell order, not classify projective directions alone.  A
promising formulation is therefore: classify equality cases of `(A)`, then
show that every nonscalar equality case acquires curvature on some interior
atomic interval.  This is a finite-geometry/circular-cover problem that may
generalize without reproducing the SAT search.

## Reproduction

From the session worktree root:

```powershell
python 04-computation/lrc14_microstaircase_binary_audit_microstaircase_audit_agent.py
python -O 04-computation/lrc14_microstaircase_binary_audit_microstaircase_audit_agent.py
```

Both are full runs, including the `n=14` ungauged enumeration.  Their semantic
stdout transcripts are byte-identical after LF normalization:

```text
stdout sha256 = 35f5bfd2f3f4a15ab5e12e2c0017696fe295c3e0a55df1a8dc5f7698fbb54569
script sha256 = 8f342b69ed1b9e7abf93ec5779b87ed6f5add2b9a9a1802a0634a870dc1b02aa
```

One normal run, one full `-O` run, and a second normal run produced this same
LF-normalized hash.  Native conflict/decision counters were observed to vary
between processes and are intentionally excluded from semantic stdout.

The run used CPython `3.13.14`, `python-sat 1.9.dev5`, MiniSat22, and
Lingeling on Windows 11.  The frozen transcript is
`05-knowledge/results/lrc14_microstaircase_binary_audit_microstaircase_audit_agent.out`.

## Boundary

This closes HYP-1818 test-plan item 3 only for the stated finite `n=14`
residue-cell system.  It does not supply the lifted prime-grid descent, a
tight-lift lemma, or a lonely time for arbitrary thirteen-speed sets.  The
separate all-`n` rigidity suggestion remains a hypothesis; in particular this
audit did not complete `n=13`.
