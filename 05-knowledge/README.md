# Persistent knowledge and reproducible evidence

`05-knowledge/` stores the research graph below canon: hypotheses, raw results,
variable dictionaries, and external reference maps. The goal is cumulative
understanding, not indiscriminate accumulation.

## Structure

- `hypotheses/` — every precise tested claim, positive or negative. The bounded
  `hypotheses/INDEX.md` is current routing; its linked historical ledger is
  provenance. A detail file carries statement, evidence, status, and repair.
- `results/` — frozen outputs paired with their generating scripts. The index
  records which conclusion each output supports.
- `variables/` — definitions, domains, equations, known values, transformations,
  and cross-links for recurring quantities.
- `reference/` — exact literature imports. Start with
  [`reference/CORE-PAPERS.md`](reference/CORE-PAPERS.md).

The indexes are discovery surfaces, not automatically current truth. Resolve a
claim through current canon and corrections before relying on it.

## Hypothesis record

A hypothesis should state:

```text
claim and exact quantifiers
status and last test
why it was plausible
universe / assumptions / inherited filters
positive and hostile controls
result and reproduction command
if true: mechanism, boundary, dependencies, generalization
if false: minimal witness, first failed implication, surviving core,
          repaired statement, failure genus, new next question
consumers and related concepts
```

“No examples found” is evidence only over the named universe. “No OEIS match”
is a search result, not a novelty theorem.

## Computational result header

Every load-bearing output/index entry should make these fields recoverable:

```text
claim tested
exact universe and representation
filters and why each is sound
positive controls / rediscovery gates
adversarial controls
consequence printed
source script and matching output
reproduction command and environment
source/output hashes
independent path or optimized/unoptimized replay
```

An input model must not print its own assumption and call it validation. A
negative result inherits every filter. If pruning is used, include a future-
witness/multiple audit and a less-pruned cross-check for the load-bearing path.

## Preserve distinctions

Do not collapse:

- indexed sequences into support sets without a collision tax;
- scalar equations into coefficientwise identities;
- labeled outputs into isomorphism invariants;
- signed coordinates into absolute magnitudes;
- samples into exhaustive universes; or
- a standalone Lean build into root-imported formal coverage.

Results are never discarded, but stale or unsound results receive a prominent
correction pointer so future searches cannot mistake them for live evidence.
