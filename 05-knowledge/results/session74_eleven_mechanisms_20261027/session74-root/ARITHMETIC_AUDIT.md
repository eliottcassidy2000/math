# Independent root audit: weighted bracket and marked Paley transfer

Status: PASS of the explicit elementary claims, with bounded independent checks.
The auditor used independently authored scripts, not imports of the producer checker.

## Proof audit

Read the arithmetic packet in full. For q>=2, a bracket contributes exactly
8m-ceil((2m+1)^2/q) integers when 1<=m<=2q-2. The concave nonemptiness
polynomial is positive at both endpoints and negative and decreasing thereafter.
The two omitted formal endpoints each contribute -1. Restoring them reduces
the sum to a polynomial plus two complete periods of the ceiling tax, proving
the all-q formula, including even q. No primality is needed for this step.

For odd prime p, the square pushforward has multiplicities1 at zero and
1+chi(a) elsewhere, proving R_p=p(p-1)/2+sum a chi(a). For p^e, unit lifts
have all p high digits because 2x is invertible; nonunit roots repeat lower
roots p times and multiply their residues by p^2. Thus
R_e=U_e+p^3 R_(e-2). The zero count p^floor(e/2) is retained. Summing the
recurrence gives the stated geometric coefficient for the weighted character
sum. This is an exact elementary derivation, not a guess from the sample data.

At341, direct root enumeration independently recovers all four conductor
strata, the75-element unit image, and its equality with the subgroup generated
by4 and81. The two axes have total ceiling tax6138, so removing them produces
the exact count error2*6138/341=36. Root multiplicity and least-residue weights
cannot be replaced by unweighted support.

## Independent computation

- check_bracket_residue.py: all q2..200 against direct interval counts;
  literal integer membership q2..10; every odd prime power<=2000 with base<=101;
  cutoff endpoints, q11 weighted/unweighted distinction, and q121 false-zero
  count76998 versus true77008.
- check_marked_paley.py: all odd primes below100; verifies square moments and
  the marked-origin translation carry; independent full341 conductor census;
  exact generated-subgroup equality; explicit eleven-root moment vector.
- Both scripts passed normally and with Python -O; paired receipts are identical.
- The arithmetic producer separately read and accepted
  WEIGHTED_PALEY_INTERFACE.md, independently rederived its translation formula,
  and verified its source-weight convention. Its one wording clarification,
  n>=2 rather than unqualified positive sources, has been applied.

These checks do not assert a prime-count formula or a Collatz/tournament
classification. The proof audit establishes the algebraic formulas; the finite
checks verify implementation and conventions only. No literature-priority claim.
