# Session138 later checkpoint: exact mixed-profile sign chambers

This packet continues the actual Q4 mixed-profile hostile from the separately
frozen `zero_rigidity` checkpoint. It does not edit that checkpoint.

`NO_ZERO_THEOREM.md` proves, for the three explicitly labelled star-completed
modules M_867,M_995,M_996 and every a,b>=1,

* 867/995: positive response when min(a,b)=1, negative when a,b>=2;
* 996: positive response in all these contexts, and in fact a sharp lower
  bound135618H(A)H(B) for arbitrary nonempty side tournaments A,B;
* thin positivity for867/995 also allows the other side to be arbitrary;
* the negative quadrant extends to sides whose colored path-cover
  polynomials are `(m³+2m)²` times nonzero nonnegative polynomials, including
  ordinal sums containing two C3 components and arbitrary other components.

The mechanism is an exact tridiagonal run-word response, a one-crossing
comparison, and an elementary binomial-profile likelihood-ratio lemma.
No all-parameter statement is inferred from a bounded scan. General strong-
quotient zero rigidity is still OPEN.

## Files and audit

* `mixed_profile_inputs.json`: full literal twelve-vertex original and
  completed modules, frozen differences, and source receipt hash.
* `audit.py`: standalone stdlib implementation; validates tournament typing
  and the same deletion card, independently counts all induced Hamiltonian
  paths backwards, then counts unordered path covers by anchored partitions.
  It verifies all input differences, tridiagonal sign structure, explicit
  base log-concavity and negative base values.
* `audit.out`, `audit_optimized.out`: identical exact receipts.

Run `python3 -B audit.py` and `python3 -B -O audit.py`. No source imports or
network access are needed. The finite range checks of profile ratios are
hostile controls only; the all-order lemma is proved in the manuscript.

The independent Session138 reconciliation proof audit accepted the complete
stated scope, including the broader negative chamber and boundary sharpenings.
Its separate audit packet is `independent_zero_balance_later`; this packet's
curated bytes are fixed by `public_manifest.json` and `SHA256SUMS`.
