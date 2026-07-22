# Canon: admission, scope, and correction

`01-canon/` contains the repository's adjudicated mathematical claims. A file
being in canon does not erase its status: proved, cited, finite-exact, verified,
conditional, refuted, retired, and superseded records can all be valuable when
their scope is explicit.

Current corrections in [`ACTIVE-GUARDRAILS.md`](ACTIVE-GUARDRAILS.md) and the
full [`MISTAKES.md`](MISTAKES.md) override older prose. Bare theorem numbers are
not stable addresses because legacy collisions exist; cite ID plus slug/path.

## Status vocabulary

| Status | Meaning |
|---|---|
| `PROVED` | Complete in-repo proof, with dependencies stated |
| `CITED` | Exact external theorem imported with hypotheses/version |
| `FINITE-EXACT` | Exhaustive proof over a stated finite universe |
| `VERIFIED` | Reproducible evidence; no general proof claimed |
| `CONDITIONAL` | Correct implication with a named open input |
| `OPEN` | Precise unresolved statement |
| `RESERVED` | Unproved namespace or provisional candidate under audit; never a result/dependency |
| `REFUTED` | False; witness and surviving repair retained |
| `RETIRED` / `SUPERSEDED` | Historical pointer; use the linked replacement |

Do not use `PROVED` for a finite census outside its universe, `CITED` for an
abstract-level paraphrase, or `VERIFIED` for a theorem whose consequence was
not actually computed.

Current high-signal status: THM-2058--2083 are proved in their stated scopes;
THM-2061/2065/2073 are reductions, THM-2069 is a deletion-code filter,
THM-2074 is density-one only, THM-2082 is a method boundary, and none proves LRC(14).
THM-2084--2086 remain reserved outside the proof graph. MISTAKE-235 repairs the
LRC/GMC weighted-fiber claim, MISTAKE-236/237 narrow the Jacobian claims, and
MISTAKE-238/239 repair dyadic transfers; MISTAKE-240 quarantines a vacuous Lean class predicate.

## Theorem record template

```markdown
---
id: THM-NNNN
title: Exact descriptive slug
status: PROVED | CITED | FINITE-EXACT | VERIFIED | CONDITIONAL | OPEN | RESERVED | REFUTED | RETIRED
source: instance/date or full external citation
depends_on: []
related: []
scripts: []
outputs: []
formalization: []
---

# THM-NNNN — title

## Statement
[Exact domain, hypotheses, conclusion, and quantifiers.]

## Scope and quantifiers
[Uniform/fixed-n; primitive/all; labeled/isomorphism; open/closed boundary.]

## Logical role
[Necessary, sufficient, iff, reduction, certificate, or obstruction.]

## Mechanism
[Why the statement holds; the coordinate in which it becomes natural.]

## Proof
[Complete proof, exact citation import, or explicitly delimited finite proof.]

## Equality and failure boundary
[Extremals, first failure, sharpness, and canonical hostile examples.]

## Dependencies and consumers
[Minimal inputs; what this enables; what it does not prove.]

## Quotients and sidecars
[Preserved predicate, destroyed information, and restoration data.]

## Verification / formalization record
[Universe, filters, controls, command, hashes, build, axioms, imports.]

## Correction lineage
[Earlier claim, MISTAKE entry, repaired statement, supersession history.]
```

Not every section needs equal length; every section relevant to the claim must
be answered. A short theorem can remain short.

A pure namespace reservation must say `RESERVED / UNPROVED EMPTY STUB`, leave
`depends_on: []`, and contain no claim. A reserved proof candidate must name
its missing audit/proof obligations and candidate prerequisites explicitly;
it still cannot be cited as a result or as a dependency of proved canon.

## Admission gate

Before promotion:

1. search the exact statement, constants, quantifiers, synonyms, and IDs;
2. audit types and necessary/sufficient/iff direction;
3. audit symmetries, orbit representatives, quotient loss, and non-vacuity;
4. attack equality, boundary, degenerate, and structured adversarial cases;
5. verify the conclusion—not only an intermediate statistic—with controls;
6. state external source version and exact imported theorem;
7. for Lean, build, audit axioms/hypotheses, and confirm root-import reach; and
8. run the namespace collision check before claiming an ID.

See [`../00-navigation/RESEARCH-PROTOCOL.md`](../00-navigation/RESEARCH-PROTOCOL.md)
for the full validity gate.

## Corrections

When a claim is wrong, repair current canon promptly and preserve history:

- add a `MISTAKE-*` record with the minimal witness or first invalid step;
- retain the strongest surviving statement and its proof;
- link the repaired theorem and correction in both directions;
- update rolling frontier/guardrails if startup truth changed; and
- open a court case only when the mathematical disagreement remains genuine.

Refutation is progress when it leaves a reusable failure genus, repaired
statement, and new question.
