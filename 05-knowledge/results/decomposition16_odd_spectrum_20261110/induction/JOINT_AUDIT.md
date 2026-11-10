# Independent review of the named deletion-cube reduction

**Verdict: PASS for JOINT_CORES.md J1–J3.** Review by the odd-induction
agent, independent of that note's author. No spectrum-completeness claim is
supported or made. The proof review covered the exact mathematical text
pinned in `audit_joint.json`; later edits require a pin refresh.

## Proof review

- **J1:** Every unmarked SCC is unchanged by every named deletion, so its
  count factors from every coordinate even when marked SCCs split. Hence
  the product divides the coordinate gcd. At gcd1 those components are
  singletons. Removing precisely those vertices preserves the entire cube.
  The integer inequality `5^(core_order-k)<=N^3` gives exactly the stated
  bound `core_order<=k+floor(log_5(N^3))`. It includes marked singleton
  components and `N=1` correctly. The Busch input is cited explicitly.
- **J2:** Distinguished names are fixed individually. Ordered SCCs make the
  core and all `k+1` unmarked transitive gaps unique, including adjacent
  marked singleton components. Thus each fixed core contributes
  `z^b/(1-z)^(k+1)`; no automorphism divides this coefficient further and
  no permutation of the component sequence is allowed. The resulting
  numerator is finite by J1, and the eventual degree is at most `r`.
- **J3:** Separate components give four counts `(Rac,Rbc,Rad,Rbd)`, so the
  determinant is zero without the primitive hypothesis. Its contrapositive
  locates both roots in one SCC. At gcd1 every other component is a
  singleton, so the two-gap formula follows. The converse is properly
  rejected by the strong `(15,3,5,1)` example.
- **Information boundary:** The deleted outside vertices can affect future
  nonuniform attachments, even though they do not affect the present cube.
  The note correctly distinguishes the deletion cube from a joint
  attachment observer and does not claim reconstruction or native
  operation compatibility.

## Independent controls

`audit_joint.py` imports only this agent's independent subset-DP and Floyd
reachability helpers from `check.py`; it does not import the producer's
`joint_check.py`. The producer uses literal path permutations and separate
DFS reachability, giving distinct counting and component algorithms.

The auditor checks all1,096 labelled tournaments of orders3..5 with three
fixed named roots, all8,768 deletion-cube entries after kernel reduction,
all584 marked classes, and all109 corresponding rational-series identities
through order5. Canonicalization fixes the three named roots and permutes
only the unmarked vertices. It independently recounts the two-root strong
zero-determinant example. Normal and optimized receipts are byte-identical.

    python3 -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/induction/audit_joint.py
    python3 -O -B 05-knowledge/results/decomposition16_odd_spectrum_20261110/induction/audit_joint.py

These finite controls do not replace the all-r proofs. The whole root code
was also reviewed: its complete labelled universe, fixed-name quotient,
empty-path convention, exact integer bound, and root order in the four
coordinates match the note.

## An explicit failure without primitivity

The note's warning can be strengthened to an exact all-r example. Fix
`r>=1` named roots and require every deletion coordinate to equal3. Every
H3 tournament is a transitive prefix, one directed triangle, and a
transitive suffix. No named root can belong to the triangle, because its
deletion would reduce the count to1. Thus all r roots are singleton SCCs,
and the skeleton consists of those r individually named singletons and
one unmarked triangle. There are `(r+1)!` ordered skeletons and `r+2`
unmarked singleton gaps. Therefore

    G_(constant cube 3)(z) = (r+1)! z^(r+3)/(1-z)^(r+2).

Its eventual coefficient degree is `r+1`, exceeding the primitive bound
`r`. For one root this is `2z^4/(1-z)^3`, and for two roots it is
`6z^5/(1-z)^4`. This is a proof of the necessity of a hypothesis controlling
unmarked nontrivial components, not merely an example of a gcd exceeding1.
The producer incorporated this refinement into the reviewed JOINT_CORES.md.
