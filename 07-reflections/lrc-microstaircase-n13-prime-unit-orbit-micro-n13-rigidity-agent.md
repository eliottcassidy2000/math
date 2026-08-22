# The missing n=13 micro-staircase case closes by prime-unit orbit splitting

Status: **FINITE-EXACT.**  In the `n=13`, `k=12` residue/open-cell
micro-staircase system, the full blockers are exactly the thirteen scalar
ramps

```text
v_i = m i mod 13,  m in F_13.
```

Equivalently, after the THM-363 gauge `v_1=0`, the zero vector is the only
full blocker.  This is a finite classification, not a speed-to-residue lift,
not LRC(13), and not LRC(14).

## Inheritance pass

- **Closest proved mechanism:** THM-363 makes scalar addition a genuine gauge
  action, and THM-364 proves for every `n` that every scalar ramp is a full
  blocker.  Only the converse at `n=13` was missing.
- **Canonical hostile:** at composite `n=14`, the normalized vector
  `(0,1,0,...,0)` passes all `14^2` right-adjacent affine-cover tests but
  misses `336` full-cell candidates.  Thus boundary data alone is not a
  generally valid replacement for the Farey interior.
- **Corrected near miss:** the original S363 statement forgot the scalar-ramp
  equality family; S364/THM-364 repaired it.  More recently, a direct `n=13`
  normalized-nonzero SAT call was deliberately left open because it was
  solver-hard.  The present result changes the quotient before solving; it
  does not reinterpret that stalled call as evidence.
- **Least-used sidecar:** the binary audit's rank-two boundary shadow
  `a i+s v_i in {0,-1}` becomes especially rigid over the field `F_13`.
  Its low-support consequence is recorded below, although the full
  boundary-only classification remains unresolved.
- **Evidence guardrails:** following MISTAKE-267, every truth gate is an
  explicit `require`, so `python -O` retains it.  Following MISTAKE-331, no
  solver-selected proof basis, conflict count, or timing enters the semantic
  transcript; the deterministic branch CNFs are what receive hashes.

The live concept board was: scalar gauge; field-unit scaling; first-support
flag; cell reflection; rank-two affine cover; composite CRT/torsion strata.
The field-unit action is the connection that changed the search.

## Exact universe

For each `1<=i<=12`, the discontinuities are exactly

```text
a/(13 i),  0<=a<=13 i.
```

Their rational union has `599` points and therefore `598` atomic open
intervals.  Samples at one-third and two-thirds of each interval give the
same floor pattern, and all `598` patterns are distinct.  Hence the
unreduced universe has exactly

```text
13 * 598 = 7,774
```

shift/cell candidates.  Shift zero is universally blocked.  The exact
reflection

```text
(s,b) -> (-s,12-b)
```

preserves every blocker support, so shifts `1,...,6` represent all nonzero
shifts.

## The new prime-unit reduction

Start with an arbitrary full blocker `v`.  THM-363 gives the normalized
vector

```text
w_i = v_i-v_1 i,  so w_1=0.
```

For each nonzero `c in F_13`, multiplication `w -> c w` preserves full
blocking.  Indeed, for every fixed cell pattern `b`, the test for `c w` at
shift `s` is exactly the test for `w` at shift `cs`; multiplication by `c`
permutes the complete shift set.

If `w` is nonzero, let `j in {2,...,12}` be its first nonzero coordinate and
put `c=w_j^(-1)`.  Then `c w` has

```text
w_1=...=w_(j-1)=0,  w_j=1.
```

Thus every normalized nonzero orbit meets one of only eleven disjoint
first-support branches.  The map is:

```text
source       normalized nonzero blocker w
target       branch (j, c w) with first nonzero value 1
map          c = w_j^(-1)
preserved    full blocker predicate
forgotten    absolute nonzero residue scale
sidecar      first nonzero coordinate j
```

This reduction is unavailable verbatim at composite modulus: a first nonzero
residue can be a nonunit.  There the natural replacement is a gcd/CRT orbit
split, connecting this result directly to HYP-1832's torsion-leak grammar.

## Independent encoding and verdict

The standalone computation does not import the previous one-hot classifier,
the `n=14` binary audit, S364, or S371.  Each of the twelve residues is encoded
by four bits and thirteen equality atoms.  Exhaustive sixteen-value truth
tables verify both directions of every equality biconditional and the clauses
excluding binary values `13,14,15`.

The common formula has:

```text
variables                  204
base clauses               816
representative cover clauses 3,588
```

Every cover clause has exactly two permitted residues at each of twelve
coordinates.  There are no duplicate representative clauses.  Each of the
eleven branch formulas receives its own DIMACS SHA-256 digest.  Under a
`250,000`-conflict bound per branch, both CaDiCaL 1.9.5 and MapleSAT return
UNSAT for every first-nonzero coordinate `2,...,12`.

The positive direction is checked independently against the unreduced
universe: all thirteen scalar ramps block all `7,774` candidates.  The
nonscalar normalized hostile

```text
(0,1,0,0,0,0,0,0,0,0,0,0)
```

misses `242` candidates; its first miss is shift `1` on cell `311`,

```text
alpha in [27/52,61/117),
b=(6,0,7,1,7,1,8,2,8,2,9,3).
```

Combining the UNSAT converse with THM-364 proves the stated finite
classification.

## Boundary-shadow tangent

Immediately to the right of `alpha=a/13`, the floor pattern is

```text
b_i = a i mod 13.
```

Therefore a full blocker must cover every `(a,s) in F_13^2` by

```text
a i+s v_i in {0,-1}
```

for some `i`.  An exact scan of all normalized low-support vectors finds:

```text
support 1: 132 scanned,   0 boundary survivors
support 2: 7,920 scanned, 0 boundary survivors
```

So any nonscalar boundary-cover equality case at `n=13`, if one exists, has
support at least three.  This contrasts sharply with the composite `n=14`
support-one boundary hostile.  A bounded million-conflict attempt at the full
boundary-only normalized classification did not return a verdict and was
stopped; no sufficiency claim is made.

The next useful mathematical target is consequently smaller than the
598-cell SAT problem:

> Classify two-level affine covers of `F_p^2` by the labelled directions
> `(i,v_i)`, first for `p=13`; if nonscalar equality cases survive, prove that
> an interior Farey cell exposes each one.

For primes the field-unit orbit split should be applied before that
classification.  For composites, replace it with valuation/CRT strata.  This
prime-versus-composite bifurcation is a plausible organizing principle for an
all-`n` proof.

## Consequence through n=15

The earlier dependent scout closed `3<=n<=12`, `n=2` is trivial, and the two
independent classifiers close `n=14,15`.  With the formerly missing `n=13`
case now closed, scalar ramps are FINITE-EXACTLY the only full blockers for
every

```text
2 <= n <= 15.
```

This is evidence for the all-`n` rigidity hypothesis, not a proof of it.

## Reproduction

```powershell
python 04-computation/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.py
python -O 04-computation/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.py
```

Both modes produced the stored semantic transcript exactly on CPython
`3.13.14` with `python-sat 1.9.dev5`.

Artifacts:

- `04-computation/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.py`
- `05-knowledge/results/lrc_microstaircase_n13_prime_unit_orbit_micro_n13_rigidity_agent.out`
- THM-363: `01-canon/theorems/THM-363-lrc-scalar-gauge-reindexing.md`
- THM-364: `01-canon/theorems/THM-364-lrc-scalar-ramp-cell-blocking.md`
- HYP-1823: `05-knowledge/hypotheses/HYP-1823-lrc-scalar-gauge-quotient.md`
- `07-reflections/lrc14-microstaircase-binary-audit-microstaircase-audit-agent.md`
