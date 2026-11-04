# Scope and adversarial audit

## Proved

- Edge identity for every finite live binary DAG and every certified
  output-avoiding constant region satisfying all-constant-parent closure.
- Full equality topology for the certified region plus its immediate
  boundary, including merges and arbitrary external sharing.
- Literal-substitution deletion certificate in the permissive-output model,
  and its sharp one-gate correction in strict B2.
- Safe telescoping only when regions are measured in the current DAG and
  every counted object is actually removed in that step.

## Quantifier and typing checks

- Sources are induced-region sources, not the chosen set of known constant
  gates. The comparable-chain counterexample proves why this matters.
- Edges count distinct child gates, not pins. Repeated pins are allowed;
  their one-parent constant nodes contribute f, not a fictitious merge.
- C can contain semantic constants, but the application must certify them.
  The theorem does not decide constancy or imply polynomial-time detection.
- Liveness and exclusion of the output gate from C are essential to the
  nonnegative fanout excess. Semantic nonconstancy is a sufficient reason
  for that exclusion, but does not follow merely from choosing output outside
  C. The deletion corollaries explicitly assume a nonconstant restricted
  output, so their one possible strict terminal gate is a NOT, not a newly
  charged constant. A dead C gate or an output gate in C would have no child.
- Identity (1) remains an incidence theorem when a B gate happens also to
  be semantically constant but was not put in C. It still has one C parent
  and is deletable. The strongest semantic application simply takes all
  known constant gates and closes them first.
- Equality of the certificate size does not prohibit more deletions from
  later alias propagation or resynthesis.
- A free affine input change is not used. The two-essential-input criterion
  applies to strict B2 in the stated ordinary primary-input coordinates.

## Hostile and sharp cases

- Comparable chain: r=1,f=1,E=0 gives three deletions, rather than treating
  two comparable constants as two independent sources.
- Merge tree: r constant leaves, r-1 constant binary merges, one boundary
  give exactly 2r certified gates, for all r.
- Fork: one constant with two distinct boundary children has E=1 and three
  certified gates. A common later child is not charged twice.
- Repeated parent: (G,G) creates f=1; counting pins as two independent
  incoming edges would falsify the identity.
- Terminal NOR: a=0 in NOR(a AND b,x) gives NOT x. Both old gates can be
  removed in the signed-output model, but strict B2 needs one terminal gate.

## What the computation does and does not certify

The topology enumeration permits every closed region, even when arbitrary
simultaneous seed values might not be attainable by a Boolean restriction.
That is a stronger incidence universe, not a claim of semantic realization.
The separate truth-table census uses actual restriction-induced constants.
It evaluates the normalized function independently on the restricted inputs
and counts the reachable normalized DAG, retaining asymmetric truth tables
and the terminal sign. No source candidate program is imported.

Normal and optimized Python runs match exactly. The program uses explicit
exceptions, so Python optimization cannot remove its mathematical checks.
The all-size proof is elementary counting and substitution, not an
extrapolation of these finite experiments or a formal proof-assistant result.

## Not claimed

No proof of the provisional elimination algorithm's global coverage; no
improved linear circuit coefficient; no forced linear frequency of f or E;
no minimum-circuit theorem for long GP/P2 sequences; no literature priority.
