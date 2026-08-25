# arXiv:2608.13306v1 matching-logic source audit

**Status: PREPRINT CLAIM / UNDER AUDIT. No computational or machine-checked
certificate. THM-4090 is a separate proved two-sort sharpening of one
displayed counterexample.**

## Identity and freshness

- Xiaohong Chen and Grigore Roşu,
  [*Completeness and incompleteness of basic matching logic*](https://arxiv.org/abs/2608.13306v1),
  arXiv:2608.13306v1 `[cs.LO]`, submitted 2026-08-13, 26 pages;
  [PDF](https://arxiv.org/pdf/2608.13306v1).
- Frozen PDF SHA-256:
  `c364f86abc5dba910a506b5aec62ef9541b2022a658aa2efa2178f7637ede132`.
- The acknowledgements disclose use of Claude Opus 5 for constructions,
  proofs, and drafting, and ChatGPT Sol for referee reports; the authors say
  they verified/corrected the work and accept responsibility. The paper has
  no accompanying code or machine proof and proposes mechanization as future
  work.

This source is unrelated to p-adic zeta values despite the adjacent user
prompt. Its object is matching-logic proof theory.

## Positive theorem ledger

The positive scope is: one sort; arbitrary finitary signature, possibly
infinite; nullary symbols allowed; set-valued symbol interpretations; element
variables denote singletons; no definedness, fixpoints, or free set variables;
closed `Gamma,phi`; arbitrary possibly infinite `Gamma`; the displayed
standard Hilbert system.

Let `E={(sigma,i)}` be the coordinate set, `[p]` the box associated to a word
`p in E*`, and

```text
Delta_Gamma={[p]gamma: gamma in Gamma, p in E*}.           (1)
```

The paper claims:

- Lemma 8: `[[Delta_Gamma]]` is the largest backward-closed subset of
  `[[Gamma]]`;
- Theorem 13: `Gamma models phi` iff `Delta_Gamma models_loc phi`;
- Theorem 14: `Gamma proves phi` iff `Delta_Gamma models_loc phi`;
- Corollary 15: `Gamma models phi` iff `Gamma proves phi` (global
  completeness);
- Corollary 16: the definedness extension is conservative on
  definedness-free `Gamma,phi`.

The model-theoretic mechanism is a double cover. From a local countermodel and
backward-closed core `C`, Definition 10 forms

```text
N=C x {0,1},                                               (2)
```

forbids mixed-sheet positive-arity tuples, and duplicates constants. Lemmas 9
and 11 show that from inside `C`, all outside values carry only one in/out bit
per element variable. Corollary 12 converts the local countermodel to a global
one.

The proof imports two load-bearing results rather than proving them:

- `(L)` strong local completeness, cited to Chen's thesis/report;
- `(S)` soundness, cited to prior Chen--Roşu work.

The arbitrary-signature extension additionally uses ordinary first-order
compactness. A bounded structural read found no immediate gap in Lemmas 4--14
conditional on these inputs; this is not a specialist or formal audit.

## Negative theorem ledger and exact boundaries

- Theorem 19 claims that validity for one-sorted definedness-free matching
  `mu`-logic is not recursively enumerable already over one unary and two
  binary symbols with no constants. The mechanism is an MRDP reduction using
  a conjunctive least fixpoint. It excludes sound weakly complete calculi with
  r.e. proof relation; it does not exclude arbitrary undecidable/set-theoretic
  proof relations. The aconjunctive, unary-only, and applicative-minimal
  boundaries remain open.
- Proposition 30 gives a satisfiable three-sort `Gamma models phi` with
  `Gamma` not deriving `phi`, refuting global completeness of that displayed
  calculus, not
  axiomatizability of semantic consequence.
- Theorem 33 says every well-founded calculus built from hypotheses/valid
  leaves and localization-respecting rules can derive only localized
  consequences. Corollary 36 applies this to `H(forall)` with a nominal. The
  paper explicitly does **not** claim that `H(forall)` lacks every complete
  calculus.

[THM-4090](../../01-canon/theorems/THM-4090-two-sort-matching-logic-global-completeness-obstruction.md)
sharpens Proposition 30 relative to the paper's cited soundness input: two
sorts `a,b` and one unary `f:b->a` suffice. The
theory

```text
Gamma={forall x:b forall y:b. f(x and y)}                 (3)
```

forces `M_b` to be a singleton, hence globally entails `forall x:b.x`; sort
flow prevents an `a`-sorted hypothesis from entering a `b`-sorted derivation.
The proof was independently audited rule by rule.

## No current mathematical interface to the headline repo problems

| frontier | missing source object/map | exact firewall |
|---|---|---|
| planar JC | no polynomial map, Keller determinant, function field, curve, or degree owner | completeness of an encoding would not prove its semantic premise or validate the encoding |
| LRC(14) | no torus maximum, integer-speed reduction, lonely time, width, or finite termination | global completeness constructs no quantitative witness; recursive encodings may require `mu` |
| p-adic zeta | no p-adic measure, period, valuation estimate, linear form, or irrationality exponent | the natural-number arithmetic in Theorem 19 serves only MRDP hardness |
| tournaments | coordinate reachability may have loops, missing pairs, and both directions | neither reachability nor sort flow is an intrinsic tournament orientation |

A finite tournament can be encoded as a relational model, but Corollary 15
would concern consequence in the encoding; it supplies no Hamiltonian count,
Pfaffian, discriminant, SCC inequality, or Paley theorem.

The lawful lesson is proof-system hygiene: local consequence can lose an
unreachable global carrier, and adding fixpoints can destroy effective proof
search. Any future use needs a semantics-preserving one-sorted,
fixpoint-free, nominal-free encoding plus translation soundness. No such
encoding currently exists for the repo's mathematical frontiers.
