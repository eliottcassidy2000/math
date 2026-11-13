# Independent audit: the one-birth defect budget on the 405 ray

**Audit conclusion: no mathematical correction required.** The all-order
claim in [BIRTH_BUDGET.md](../returns/BIRTH_BUDGET.md) follows from its inherited
identities, the correctly bounded finite head, and the displayed exclusions.
The independent literal recount below agrees with every supplied profile.

## Quantifiers and exact inequalities

The statement allows any actual tournament at the original full H, any order,
and any chosen root after arbitrary same-H preparation. The next operation is
specifically one newborn beaten only by that root. Its indegree is one and
its readiness defect F_new-H_new equals the old start count s. This does not
restrict multi-birth words or arbitrary endpoint replacements.

For the requested next count2h-c, the target equation is F-h=s-c. With
Delta=F-h-q+s, it becomes q+Delta=2s-c. Since q>=s, this implies
s>=c+Delta>=c+d(d-1). Both inequality directions are correct. At c=s=71,
q>=s and Delta>=0 force q=71 and Delta=0, not merely F=h. The inherited
weighted-cut identity then forces d<=1. A source root would have s=h>71,
so d=1. Deleting that root produces exactly the parent conditions
H(U)=71, S_U(a)=0, F_U(a)=405*2^k. The newborn formulas have the correct
root direction and empty-second-path convention.

## The finite reduction is genuinely all-order

Primality of71 and multiplicativity across ordered SCCs give one nontrivial
strong H71 component and otherwise singleton components. This uses H=1 iff
the tournament is transitive, so a strong H1 component is a singleton.
The possible cover counts at a root before, inside, or after that component
are respectively powers of2 times71, a marked core F, or the core P(2).
These are exactly the inherited ordered-two-bin SCC factors. Ignoring the
additional zero-start condition enlarges the candidate family and is safe
for exclusion.

The published bound was checked directly in
[Busch, Corollary1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/):
the minimum strong count at order9 is75, and its three residue-class formulas
increase thereafter. Thus every strong H71 core has order at most8.

The input's completeness is inherited, not re-proved by this audit.
[Session212's census proof](../../../../session212/cycle-spectrum/PROOF.md),
Section6, describes generation of all6880 order-eight classes by all one-vertex
extensions and exact canonicalization. Padding any smaller strong H71 core
by transitive sinks gives an order-eight class, so extraction captures an
isomorphic copy of every possible core. Calling the eleven outputs
**presentations**, rather than assuming a new canonical class list, is sound;
duplicates would not invalidate the exclusion.

## All eleven literal presentations checked independently

`birth_budget_audit.py` imports none of the original counting code. For each
of the eleven supplied outgoing-mask arrays it checks tournament orientation
and strong connectivity, then literally enumerates every permutation of every
vertex subset. These counts reconstruct all Hamilton totals, all81 marked
(s,q,F,indegree,Delta) records, and every ordered two-bin total P(2).
It also checks sum_root F_root=P(2)-H, including empty-bin conventions.

All534,293 literal permutation tests reproduce the supplied records exactly.
The eleven F lists are the lists in BIRTH_BUDGET.md, including both404 values,
both71 equality entries, and the maximum595. Since595<2*405, odd part405
would require F=405; none occurs. P(2) values1012,1072,1976 have odd parts
253,67,247. Therefore every possible H71 parent form misses the required
405*2^k, proving equality s=71 impossible for every k>=1. The conclusion
s>=72 follows because s is an integer.

Normal and optimized runs agree; both retain explicit runtime checks. They
verify the original census SHA256 and6880-line count, but do not regenerate
its isomorphism representatives or rerun the full census.

## Scope clarification

The closing hypothetical balanced state (h,71,71,h) would repeat the scalar
ray by retained-root R. Its retained root is0-ready and each newborn would
be71-ready; even that hypothetical would not itself be unit-ready renewal.
The original prose says 'scalar ray' and does not claim otherwise. The proved
obstruction is stronger than the elementary unit-ready failure, while leaving
larger-defect and multi-birth possibilities open.

Reproduce from this directory:

```
python -B birth_budget_audit.py
python -O -B birth_budget_audit.py
```
